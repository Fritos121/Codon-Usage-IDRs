# make plots
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import seaborn as sbn
from get_statsV3 import parse_fasta, parse_data, list_GC
from codon_distV15 import codon_iter, get_Org_counts, get_Prot_counts, get_GC, change_counts

#create function to read tt_11

tt_11 = {
			'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
			'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
			'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
			'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
			'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
			'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
			'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
			'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
			'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
			'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
			'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
			'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
			'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
			'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
			'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'O',
			'TGC':'C', 'TGT':'C', 'TGA':'U', 'TGG':'W',
		}

full_cds, cds_list, cds_gc = parse_fasta('CDS.fasta')
full_pseudo, pseudo_list, pseudo_gc = parse_fasta('pseudogene.fasta')
wanted_rna, rna_list, rna_gc = parse_fasta('RNA.fasta')
full_ncds, ncds_list, ncds_gc = parse_fasta('non_CDS.fasta')
data_cds, data_list, ordered_list, disordered_list = parse_data('Bacteria_Proteobacteria_Gammaproteobacteria_Escherichia.83334.data')

# uses original way
disorder_gc = list_GC(disordered_list)
order_gc = list_GC(ordered_list)
data_gc = list_GC(data_list)

# uses codon dist
codon_counts, aa_counts = get_Org_counts(data_list, tt_11)
codon_proportions = change_counts(codon_counts, tt_11)
# in org_counts function, removed cause I only want to check once
# print(codon_counts['ambig_base'])

disordered_counts = get_Prot_counts(disordered_list, tt_11)
disordered_gc = [get_GC(codon_proportions, aa_dict, tt_11) for aa_dict in disordered_counts]
#print(len(disordered_gc), len(disordered_list))
ordered_counts = get_Prot_counts(ordered_list, tt_11)
ordered_gc = [get_GC(codon_proportions, aa_dict, tt_11) for aa_dict in ordered_counts]


for x in [pseudo_gc, rna_gc, ncds_gc, cds_gc, data_gc, order_gc, disorder_gc, ordered_gc, disordered_gc]:
	print(len(x))

plt.grid(True, axis='y')
#plt.hlines(range(30, 70, 5), 0, 10, linewidths=.5, linestyles='dashed')
boxplot = plt.boxplot([pseudo_gc, rna_gc, ncds_gc, cds_gc, data_gc, order_gc, disorder_gc, ordered_gc, disordered_gc], patch_artist=True, whis=[5, 95], showfliers=False)
plt.xticks([1,2,3,4,5,6,7,8,9], ['Pseudoprotein GC', 'r/tRNA GC', 'non-coding GC', 'CDS GC', 'Data_CDS GC', 'Observed Ordered GC', 'Observed Disordered GC', 'Expected Ordered GC', "Expected Disordered GC"], rotation=45)
plt.ylabel('Percent GC')
plt.xlabel('Gene Class')
plt.title('Distribution of Gene GC by Class in E. coli')

colors = ['w', 'lightblue', 'w', 'w', 'w', 'w', 'lightblue', 'w', 'lightblue']
for patch, color in zip(boxplot['boxes'], colors):
	patch.set_facecolor(color)

plt.tight_layout()
#plt.legend()
#plt.show()
plt.savefig('data.png')
plt.close()

bandwidth = .25
cuts = 14
sbn.kdeplot(rna_gc, shade=True, bw=bandwidth, cut=cuts, label='t/rRNA')
sbn.kdeplot(disorder_gc, shade=True, bw=bandwidth, cut=cuts, label='Observed Disordered GC')
sbn.kdeplot(disordered_gc, shade=True, bw=bandwidth, cut=cuts, label='Expected Disordered GC')
plt.ylabel('Density')
plt.xlabel('Percent GC')
plt.title('Density Plot of Gene GC By Gene Class in E. coli')
plt.legend()
plt.show()
#plt.savefig('Density_Plot.png')
#plt.close()

