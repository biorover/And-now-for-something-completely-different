#!/usr/bin/env python

#CompareOverlap: A tool for looking at overlap between gtf files at the nucleotide, exon, and gene levels

import argparse
import numpy

parser = argparse.ArgumentParser(description='CompareOverlap: A tool for looking at overlap between gtf \
                                 files at the nucleotide, exon, and gene levels')

parser.add_argument('--gtf1', help = 'first gtf')
parser.add_argument('--gtf2', help = 'second gtf')
parser.add_argument('--geneNameLog', default = False, help = 'file destination for logging the names of non-overlapping genes')
parser.add_argument('--report', default = 'missing', help = 'what to report to the geneNameLog, "missing" or "fis-fus" (fissions and fussions)')

args = parser.parse_args()

def read_in_gtf(gtf_file_loc):
    """goes through each line of a gtf to add all coordinates overlapped by features to a dictionary keyed by locus.
    Dictionary format is "dict[<locus>] = [{<dictionary_of_genes>},[<list_of_exons>],[<list_of_coords>]]"""
    locus_dict = {}
    gene_names = []
    exon_names = []
    nucleotide_names = []
    for line in open(gtf_file_loc):
        if line[0] != "#" and line.count('\t') > 7:
            fields = line.replace('\r','').replace('\n','').split('\t')
            locus = fields[0]
            if not locus in locus_dict:
                locus_dict[locus]= {}
            for def_field in fields[8].split(';'):
                if "gene_id" in def_field:
                    gene_id = def_field.split()[1]
            coords = range(int(fields[3]),int(fields[4]) + 1)
            if gene_id in locus_dict[locus]:
                locus_dict[locus][gene_id][0].extend(coords[:])
                locus_dict[locus][gene_id][1].append(coords[:])
            else:
                locus_dict[locus][gene_id] = [coords[:],[coords[:]]]
                gene_names.append(gene_id)
            exon_names.append(gene_id + "exon" + str(len(locus_dict[locus][gene_id][1]) - 1))
            for i in coords:
                nucleotide_names.append(locus + str(i))            
    return locus_dict,gene_names,exon_names,nucleotide_names

locdict_one,gene_names_one,exon_names_one,nucleotide_names_one = read_in_gtf(args.gtf1)
locdict_two,gene_names_two,exon_names_two,nucleotide_names_two = read_in_gtf(args.gtf2)

nucleotide_count = [len(set(nucleotide_names_one) & set(nucleotide_names_two)),
                    len(set(nucleotide_names_one) - set(nucleotide_names_two)),
                    len(set(nucleotide_names_two) - set(nucleotide_names_one))]

overlapping_genes_one = []
overlapping_genes_two = []
perfect_genes_one = []
perfect_genes_two = []
overlapping_exons_one = []
overlapping_exons_two = []
perfect_exons_one = []
perfect_exons_two = []
count = 0 # debug

fis_fus_tracker = []

for locus in set(locdict_one.keys()) & set(locdict_two.keys()):
    cor_loc = locdict_two[locus]
    this_loc = locdict_one[locus]
    for gene_name in this_loc:
        gene = this_loc[gene_name]
        for cor_gene_name in cor_loc:
            cor_gene = cor_loc[cor_gene_name]
            if set(gene[0]) & set(cor_gene[0]) != set([]):
                if gene_name in overlapping_genes_one:
                    fis_fus_tracker.append(gene_name+"\t"+cor_gene_name)
                if cor_gene_name in overlapping_genes_two:
                    fis_fus_tracker.append(cor_gene_name+'\t'+gene_name)
                overlapping_genes_one.append(gene_name)
                overlapping_genes_two.append(cor_gene_name)
                if set(gene[0]) ^ set(cor_gene[0]) == set([]):
                    perfect_genes_one.append(gene_name)
                    perfect_genes_two.append(cor_gene_name)
                for exon_num in range(len(gene[1])):
                    exon_set = set(gene[1][exon_num])
                    for cor_exon_num in range(len(cor_gene[1])):
                        if set(gene[1][exon_num]) & set(cor_gene[1][cor_exon_num]) != set([]):
                            overlapping_exons_one.append(gene_name + 'exon' + str(exon_num))
                            overlapping_exons_two.append(cor_gene_name + 'exon' + str(cor_exon_num))
                            if set(gene[1][exon_num]) ^ set(cor_gene[1][cor_exon_num]) == set([]):
                                perfect_exons_one.append(gene_name + 'exon' + str(exon_num))
                                perfect_exons_two.append(cor_gene_name + 'exon' + str(cor_exon_num))

exons_count = [len(set(overlapping_exons_one)) ,len(set(exon_names_one) - set(overlapping_exons_one)),
               len(set(overlapping_exons_two)) ,len(set(exon_names_two) - set(overlapping_exons_two)),]
exons_perfect = [len(set(perfect_exons_one)),len(set(exon_names_one) - set(perfect_exons_one)),
                 len(set(exon_names_two) - set(perfect_exons_two))]
gene_count = [0,0,0,0,0,0]
for gene in gene_names_one:
    if overlapping_genes_one.count(gene) == 1:
        gene_count[0] += 1
    elif overlapping_genes_one.count(gene) > 1:
        gene_count[1] += 1
    else:
        gene_count[2] += 1
for gene in gene_names_two:
    if overlapping_genes_two.count(gene) == 1:
        gene_count[3] += 1
    elif overlapping_genes_two.count(gene) > 1:
        gene_count[4] += 1
    else:
        gene_count[5] += 1

gene_perfect = [len(set(perfect_genes_one)),len(set(gene_names_one) - set(perfect_genes_one)),
                len(set(gene_names_two) - set(perfect_genes_two))]


if args.geneNameLog:
    out = open(args.geneNameLog,'w')
    if "missing" in args.report:
        out.write("gtf1 not gtf2\n" + "\n".join(list( set(gene_names_one) - set(overlapping_genes_one) )) + "\n\ngtf2 not gtf1\n" +
                  "\n".join(list( set(gene_names_two) - set(overlapping_genes_two) )))
    if "fis-fus" in args.report:
        out.write('\n'.join(fis_fus_tracker))
    out.close()
        

print nucleotide_count
print ('sensitivity',100.0 * nucleotide_count[0]/(nucleotide_count[1] + nucleotide_count[0]), 'specificity', 100.0 * nucleotide_count[0]/(nucleotide_count[2] + nucleotide_count[0]))
print exons_count
print exons_perfect
print gene_count
print gene_perfect