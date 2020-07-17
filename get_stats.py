# get statistics
from Bio.SeqUtils import GC
from statistics import mean, stdev, median

# for length dist, never used originally
def get_len(seq_list):
    lens = [len(x) for x in seq_list]
    return lens


def parse_fasta(filename):
    seq = ''
    gene_list = []
    gc_list = []
    with open(filename, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                continue
            else:
                line = line.strip()
                seq += line
                gene_list.append(line)
                gc_list.append(GC(line))

    return seq, gene_list, gc_list


def parse_data(filename):
    count = 1
    seq = ''
    seq_list = []
    order_list = []
    disorder_list = []
    want = False
    with open(filename, 'r') as fh:
        for line in fh:
            if count == 3:
                order_AA = line.count('0')
                disorder_AA = line.count('1')
                if order_AA >= 20 and disorder_AA >= 20:
                    want = True
                    scores = list(line)

            if count == 4:
                line = line.strip()
                line = line.upper()
                #does not account for open reading frames... dont want to
                #add chcek for empy line
                #check for reading frame (modulus)
                if want:
                    codons = [line[x:x+3] for x in range(0, len(line), 3)]
                    zipped = zip(scores, codons)
                    order_seq = ''
                    disorder_seq = ''
                    for x in zipped:
                        if x[0] == '0':
                            order_seq += x[1]
                        elif x[0] == '1':
                            disorder_seq += x[1]
                        else:
                            raise ValueError('Invalid prediction score. Need 0 or 1.')

                    order_list.append(order_seq)
                    disorder_list.append(disorder_seq)
                    want = False

                seq += line
                seq_list.append(line)
                count = 0

            count += 1

    return seq, seq_list, order_list, disorder_list

# takes list of protein gene seqs, gets gc
def list_GC(protein_list):
    GC_list = [GC(seq) for seq in protein_list]

    return GC_list


if __name__ == '__main__':
    full_cds, cds_list, cds_gc = parse_fasta('CDS.fasta')
    full_pseudo, pseudo_list, pseudo_gc = parse_fasta('pseudogene.fasta')
    wanted_rna, rna_list, rna_gc = parse_fasta('RNA.fasta')
    full_ncds, ncds_list, ncds_gc = parse_fasta('non_CDS.fasta')
    # data file only has coding seqs...?
    data_cds, data_list, ordered_list, disordered_list = parse_data(
        'Bacteria_Proteobacteria_Gammaproteobacteria_Escherichia.83334.data')

    print('CDS GC:', GC(full_cds))
    print('Number of protein coding bases:', len(full_cds))
    print('Number of protein coding genes:', len(cds_list))
    print('Mean GC and SD of protein coding sequences:', mean(cds_gc), stdev(cds_gc), '\n')

    print('Pseudoprotein GC:', GC(full_pseudo))
    print('Number of pseudoprotein coding bases:', len(full_pseudo))
    print('Number of pseudoprotein coding genes:', len(pseudo_list))
    print('Mean GC and SD of pseudogene coding sequences:', mean(pseudo_gc), stdev(pseudo_gc), '\n')

    print('r/tRNA GC:', GC(wanted_rna))
    print('Number of t/rRNA coding bases:', len(wanted_rna))
    print('Number of t/rRNA coding genes:', len(rna_list))
    print('Mean GC and SD of t/rRNA coding sequences:', mean(rna_gc), stdev(rna_gc), '\n')

    print('non-coding GC:', GC(full_ncds))
    print('Number of non-coding bases:', len(full_ncds))
    print('Number of non-coding segments:', len(ncds_list))
    print('Mean GC and SD of non-coding sequences:', mean(ncds_gc), stdev(ncds_gc), '\n')

    print('data_CDS GC:', GC(data_cds))
    print('Number of protein coding bases:', len(data_cds))
    print('Number of protein coding genes:', len(data_list))
    print('Median GC and SD of protein coding sequences:', median(list_GC(data_list)), stdev(list_GC(data_list)), '\n')
    # print('Total Ordered Residue GC:',

    print(data_list[0])

# #expected to be > len of genome b/c of overlap b/w feature seqs, but need all bases for GC calc
# if len(gb_seq) == ((len(full_cds) + len(wanted_rna) + len(full_ncds) + len(full_pseudo)) - overlap):
# print("Overlap in genes accounts for all extra bases")
# else:
# x = abs(len(gb_seq) - (len(full_cds) + len(wanted_rna) + len(full_ncds) + len(full_pseudo)))
# print(x, "Unaccounted bases in seqs")
