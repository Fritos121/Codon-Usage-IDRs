# added to prevent error (comment out to see it)
from gevent import monkey as curious_george
curious_george.patch_all(thread=False, select=False)

from Bio import Entrez
from time import sleep
import re
import requests
import argparse


def get_protein_info(uniprot_ids):
    """
    Retrieves EMBL accession numbers and taxonomy ids for list of proteins. Creates a dict to map each protein's
    uid to its EMBL accession number and tax id.
    :param uniprot_ids: List of Uniprot IDs, e.g., "A5HYD5"
    :return: List of EMBL accession numbers, list of taxonmy ids, and dictionary mapping each uid to its info
    """
    from bioservices import UniProt

    embl_accs = []
    taxonomy_ids = []
    missing_embl = []
    missing_taxid = []
    orthos_map = {}
    u = UniProt()

    uniprot_records = list(map(lambda x: x.decode("utf-8"), u.retrieve(uniprot_ids, frmt='txt')))  # WSL cmdline
    # uniprot_records = u.retrieve(uniprot_ids, frmt='txt')  # PyCharm

    embl_pattern = re.compile(r"DR\s+EMBL;.*?;\s+(.*?);")
    taxid_pattern = re.compile(r"OX\s+NCBI_TaxID=(\d+)")

    for i, record in enumerate(uniprot_records):
        # record = "fffffffffffffffffff"    # test placeholder Nones

        embl_acc = get_match(embl_pattern, record, uniprot_ids, missing_embl, i)
        taxonomy_id = get_match(taxid_pattern, record, uniprot_ids, missing_taxid, i)

        embl_accs.append(embl_acc)  # EMBL accession number for coding seq
        taxonomy_ids.append(taxonomy_id)  # tax_id of organism protein belongs to
        orthos_map[uniprot_ids[i]] = [embl_acc, taxonomy_id]  # map protein info to its uid

    if missing_embl:
        print('\n{} Protein(s) Missing EMBL Accession Number: '.format(len(missing_embl)) + ', '.join(missing_embl))

    if missing_taxid:
        print('\n{} Protein(s) Missing NCBI TaxID: '.format(len(missing_taxid)) + ', '.join(missing_taxid))

    return embl_accs, taxonomy_ids, orthos_map


# keep generic function even though only len() == 1 used in main()
def get_match(pattern, search_str, id_list, missing_list, index):
    """
    Find first match using RegEx; Allows for tracking of what IDs can't get the specific pattern match and placeholders
    to maintain indexing.
    :param pattern: compiled RegEx obj
    :param search_str: string to search
    :param id_list: list of identifiers
    :param missing_list: list to add identifier if pattern not found
    :param index: int; to show current spot in list of ids, for error handling only
    :return: Either the wanted pattern found (string), OR appends to list if error and returns placeholder None
    """

    match = re.findall(pattern, search_str)
    match = [x for x in match if x != '-']  # embl_acc are '-' if no Genomic DNA translation available

    if len(match) == 0:
        missing_list.append(id_list[index])
        # add placeholder value to keep indexing between refseq and uid 1:1, zipped later
        return None

    else:
        # grab first embl_acc from list
        return match[0]


def get_assembly_id(tx_id):
    """
    Retrieves NCBI Assembly ID for a given NCBI Taxonomy ID. Email for Entrez must be set in main().
    :param tx_id: string; NCBI tax id
    :return: string; NCBI Assembly ID
    """
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

    return assem_id


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=argparse.FileType('r'), help="Fasta file containing protein uids in header")
    parser.add_argument("outfile", type=argparse.FileType('w'), help="Name of fasta file to write coding seqs to")
    ftp = parser.add_argument_group('ftp', 'Optional arguments Used to create csv file with partial ftp links')
    ftp.add_argument("-f", "--ftp_out", type=argparse.FileType('w'), help="Filename for csv file with partial ftp links")
    ftp.add_argument('-e', '--email', help='Email address required to use during Entrez lookups')

    args = parser.parse_args()
    in_fh = args.infile
    fasta_fh = args.outfile
    csv_fh = args.ftp_out
    email = args.email

    if csv_fh and not email:
        exit('Email address required to use Entrez. Please provide a valid email address.')

    # get list of uids
    uids = [x.strip().replace(">", "") for x in in_fh.readlines() if x.startswith(">")]

    embl_accessions, tax_ids, ortho_mapping = get_protein_info(uids)
    # print(embl_accessions, tax_ids)

    # embl_acc, uids, and tax_ids should have same length
    # get CDS for gene
    embl_base_url = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
    for uid, embl_acc in zip(uids, embl_accessions):
        if embl_acc is None:
            continue

        else:
            r = requests.get(embl_base_url + embl_acc)
            cds = ''.join(r.content.decode('utf-8').split("\n")[1:])
            # account for any http errors... will remove ortho from data for now; write uid info to std error?
            if 'status=' in cds:
                print('\n', uid, "removed from orthologs due to error.")
                del ortho_mapping[uid]
            else:
                ortho_mapping[uid].append(cds)

    # write cds to file
    for uid, info in ortho_mapping.items():
        # do not write any ortholog with incomplete information
        if uid is None or None in info:
            continue

        else:
            fasta_fh.write(">uid=" + uid + ";embl_acc=" + info[0] + ";tax_id=" + info[1] + "\n")
            fasta_fh.write(info[2] + "\n")

    # only run if argument passed in CLI
    if csv_fh:
        Entrez.email = email
        # remove placeholder Nones
        tax_ids = [x for x in tax_ids if x is not None]
        unique_tax_ids = list(set(tax_ids))

        # get partial ftp links to be used by get_assembly_seqs_from_ftp.py
        ftp_pattern = re.compile(r"<FtpPath_RefSeq>\S+(/genomes\S+GCF_\S+)<")
        ftp_partials = []
        missing_ftp = []
        for i, tax_id in enumerate(unique_tax_ids):
            assembly_id = get_assembly_id(tax_id)
            handle = Entrez.esummary(db='assembly', id=assembly_id, report='full')
            string = handle.read().decode("utf-8")  # for WSL
            # string = handle.read()
            # print(string)

            # gets ftp (partial) link for all cds in organism
            ftp_link = get_match(ftp_pattern, string, missing_ftp, unique_tax_ids, i)
            ftp_partials.append((ftp_link, tax_id))
            sleep(0.4)  # prevent timeouts from ncbi

        if missing_ftp:
            print('\n{} Protein(s) Missing FTP Link for Assembly: '.format(len(missing_ftp)) + ', '.join(missing_ftp))

        # write to csv file
        for link, tax_id in ftp_partials:
            csv_fh.write(link + ', ' + tax_id + '\n')

        # print(len(ftp_partials), ftp_partials)


