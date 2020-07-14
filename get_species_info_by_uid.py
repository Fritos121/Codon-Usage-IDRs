# added to prevent error (comment to see it)
from gevent import monkey as curious_george
curious_george.patch_all(thread=False, select=False)

from Bio import Entrez
from time import sleep
import re
import sys

def get_protein_info(uniprot_ids):
    from bioservices import UniProt

    refseq_vers = []
    taxonomy_ids = []
    missing_refseq = []
    missing_taxid = []
    orthos_map = {}
    u = UniProt()

    uniprot_records = list(map(lambda x: x.decode("utf-8"), u.retrieve(uniprot_ids, frmt='txt')))
    # uniprot_records = u.retrieve(uniprot_ids, frmt='txt')

    refseq_pattern = re.compile(r"DR\s+RefSeq.*(XM_.*).")
    taxid_pattern = re.compile(r"OX\s+NCBI_TaxID=(\d+)")

    for i, record in enumerate(uniprot_records):
        # record = "fffffffffffffffffff"
        # print(record)
        # quit()

        refseq_ver = get_match(refseq_pattern, record, missing_refseq, uniprot_ids, i)
        taxonomy_id = get_match(taxid_pattern, record, missing_taxid, uniprot_ids, i)

        refseq_vers.append(refseq_ver)  # refseq version number for coding seq
        taxonomy_ids.append(taxonomy_id)  # tax_id of organism protein belongs to
        orthos_map[uniprot_ids[i]] = [refseq_ver, taxonomy_id]  # map protein info to its uid

    if missing_refseq:
        print('{} Proteins Missing RefSeq Version Numbers: '.format(len(missing_refseq)) + ', '.join(missing_refseq))
    if missing_taxid:
        print('{} Proteins Missing NCBI TaxIDs: '.format(len(missing_taxid)) + ', '.join(missing_taxid))

    return refseq_vers, taxonomy_ids, orthos_map


# keep generic function even though on len() == 1 used in main()
def get_match(pattern, search_str, missing_list, id_list, index):
    match = re.search(pattern, search_str)
    try:
        if len(match.groups()) == 1:
            wanted_item = match.group(1)
        else:
            wanted_item = match.groups()

    except AttributeError:
        missing_list.append(id_list[index])
        # add placeholder value to keep indexing between refseq and uid 1:1, zipped later
        return None

    else:
        # print(wanted_item)
        return wanted_item


# set email in main()
def get_assembly_id(tx_id):

    search_term = 'txid{}'.format(tx_id)  # search term format used by NCBI Genome db
    handle1 = Entrez.esearch(db='genome', term=search_term) # gets internal ID for genome in genome db
    record1 = Entrez.read(handle1)
    # print(record1)
    org_id = record1['IdList'][0]    # what if there are more than 1... not encountered yet so don't worry
    # print(org_id)
    handle2 = Entrez.esummary(db='genome', id=org_id) # gives assembly of genome org
    record2 = Entrez.read(handle2)
    # print(record2)
    assem_id = record2[0]["AssemblyID"]
    # print(assem_id)

    # depending on design of program later, make this get entire report for assembly?

    return assem_id


if __name__ == "__main__":

    # ensure proper arguments are passed. Allows for piping and standalone use of program
    if len(sys.argv) != 3:
        exit("positional arguments: {} <infile> <outfile>".format(sys.argv[0]))

    infile = sys.argv[1]
    outfile = sys.argv[2]

    # if no infile name passed, then assume contents of _ortholog_msa.txt file passed to stdin
    if infile == "-":
        in_fh = sys.stdin
    else:
        in_fh = open(infile)

    # get list of uids (has header line from _ortholog_msa.txt file)
    uids = [x.strip().replace(">", "") for x in in_fh.readlines()]

    # name will be created by program in a directory determined from file's header line, otherwise use given name
    if outfile == "-":
        base_dir = uids.pop(0)
        fn_base = "_ortholog_cds.fasta"
        outfile = base_dir + fn_base

    else:
        uids = uids[1:]

    Entrez.email = 'mlowry2@mymail.vcu.edu'

    # uids = ['A0NAQ1', 'A8DWE3']
    refseq_versions, tax_ids, ortho_mapping = get_protein_info(uids)
    # print(refseq_versions, tax_ids)

    # remove placeholders
    tax_ids = [x for x in tax_ids if x is not None]
    unique_tax_ids = list(set(tax_ids))

    ftp_pattern = re.compile(r"<FtpPath_RefSeq>\S+(/genomes\S+GCF_\S+)<")
    # get partial ftp links to be used by get_assembly_seqs_from_ftp.py
    ftp_partials = []
    missing_ftp = []
    for i, tax_id in enumerate(unique_tax_ids):
        assembly_id = get_assembly_id(tax_id)
        handle = Entrez.esummary(db='assembly', id=assembly_id, report='full')
        string = handle.read().decode("utf-8")
        # print(string)

        # gets ftp (partial) link for all cds in organism
        ftp_link = get_match(ftp_pattern, string, missing_ftp, unique_tax_ids, i)
        ftp_partials.append((ftp_link, tax_id))
        sleep(0.4)  # prevent timeouts from ncbi


    if missing_ftp:
        print('{} Proteins Missing FTP Link for Assembly: '.format(len(missing_ftp)) + ', '.join(missing_ftp))

    # refseq_versions, uids, and tax_ids should have same length
    # get CDS for gene
    for uid, refseq_ver in zip(uids, refseq_versions):
        if refseq_ver is None:
            continue

        else:
            handle = Entrez.efetch(db='nuccore', id=refseq_ver, rettype='fasta', retmode='text')
            cds = ''.join([x.strip() for x in handle.readlines()[1:]])
            ortho_mapping[uid].append(cds)
            handle.close()

    # write cds to file
    # needs to know where
    with open(outfile, 'w') as fh:
        for uid, info in ortho_mapping.items():
            if uid is None or None in info:
                continue
            fh.write(">uid=" + uid + ";refseq_ver=" + info[0] + ";tax_id=" + info[1] + "\n")
            fh.write(info[2] + "\n")

    # pipe out ftp_partials to get_assembly_seqs_from_ftp.py
    print(len(ftp_partials), ftp_partials)
