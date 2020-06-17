# create fastas for each class (CDS, pseudogene, RNA, non-CDS)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation


class Record(object):

    def __init__(self, gb_file):

        # only reads 1 record, if 0 or 2+, error
        self.record = SeqIO.read(gb_file, 'genbank')
        # Tuple with record attributes and seqs
        self.genbank = (self.record.id, self.record.name, self.record.seq.strip(), self.record.seq.reverse_complement(),
                   self.record.seq.complement())
        self.features = self.record.features # list constructor didnt work in this style... already a list

        self.gb_seq = self.genbank[2]
        self.rev_c = self.genbank[3]
        self.comp = self.genbank[4]
        # sequence to be replaced with x's where gene of any type is found; used to make non-coding seq
        self.replaced_seq = str(self.genbank[2])

        print('Genome Length: ' + str(len(self.gb_seq)))

        #print(self.features[0:3])

        # check for IUPAC AmbiguousDNA Bases, and gaps (make function)
        chars = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '.', '-']
        present_chars = [x for x in chars if x in self.gb_seq]

        if not present_chars:
            print('No Ambiguity or gaps in sequence \n')

    def get_types(self):
        self.types = {}
        i = 0
        for x in self.features:
            if self.features[i].type in self.types:
                self.types[self.features[i].type] += 1
            else:
                self.types[self.features[i].type] = 1
            i += 1

        print(str(self.types) + '\n')

    def create_fastas(self):

        fh_rna = open("RNA.fasta", 'w')
        fh_dict = {
            ('CDS', False): open("CDS.fasta", 'w'),
            ('CDS', True): open("pseudogene.fasta", 'w'),
            ('tRNA', False): fh_rna,
            ('rRNA', False): fh_rna
        }

        for x in self.features:
            if x.type == 'CDS' or x.type == 'tRNA' or x.type == 'rRNA':
                type = x.type
                pseudo_flag = 'pseudogene' in x.qualifiers
                self.__add_gene(x, fh_dict[(type, pseudo_flag)])
            else:
                continue

        with open('non_CDS.fasta', 'w') as fh:
            start, end = 0, 0
            write = False
            for x in self.replaced_seq:
                if x == 'x':
                    if write:
                        fh.write('>location=' + str(start) + ':' + str(end) + '\n' + self.replaced_seq[start:end] + '\n')
                        write = False
                        start = end + 1
                    else:
                        start += 1

                else:
                    if write:
                        end += 1
                    else:
                        write = True
                        end = start + 1
            else:
                if write:
                    fh.write('>location=' + str(start) + ':' + str(end) + '\n' + self.replaced_seq[start:end] + '\n')

        print("Fasta files created")

    def __add_gene(self, feature, fh):
        fh.write('>type=' + feature.type + ';locus_tag=' + feature.qualifiers['locus_tag'][0] + '\n')
        start = feature.location.start.position
        end = feature.location.end.position
        # build_seq += feature.location.extract(gb_seq)
        feature_loc = FeatureLocation(start, end)
        # need str obj to write into file
        fh.write(str(feature_loc.extract(self.gb_seq)) + '\n')

        self.replaced_seq = self.replaced_seq[:start] + 'x' * len(self.replaced_seq[start:end]) + self.replaced_seq[end:]


def main():
    bacteria = Record('E_coli_O157_H7.gb')
    bacteria.get_types()
    bacteria.create_fastas()


if __name__ == '__main__':
    main()