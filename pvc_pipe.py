#!/usr/bin/env python

#PVC_pipe, phased variant consensus pipeline

import argparse
import sys
import subprocess

parser = argparse.ArgumentParser(description='PVC_pipe, a phased variant consensus pipeline for converting assemblies into \
                                 pseudo-phased assemblies')

parser.add_argument('-g','--genome',dest='genome', help = 'Reference genome in fasta format')

parser.add_argument('-v','--phased_vcf', dest = 'pvcf', default = None, help = 'Phased vcf file of variants')

parser.add_argument('-p','--phase', dest = 'phase', default = 1, type = int)

parser.add_argument('-s','--sample', dest = 'sample', default = 1, type = int)

args = parser.parse_args()


#sets up genome dictionary
ref_dict = {}
input_scaffold_order = []
working_entry = ""
for line in open(args.genome):
    if line[0] == ">":
        working_entry = line[1:-1]
        input_scaffold_order.append(line[1:-1])
        ref_dict[working_entry] = ""
    elif line[0] == "#":
        continue
    else:
        ref_dict[working_entry] = ref_dict[working_entry] + line[:-1]

pvc_dict = {}


#goes through phased vcf and builds pseudo-haploid sequence
working_entry = ""
for line in open(args.pvcf):
    if line[0] == "#":
        continue
    else:
        fields = line.split('\t')
        if fields[0] != working_entry:
            working_entry = fields[0]
            last_position = 0
        pos = int(fields[1])
        alleles = [fields[3]]
        for alt_allele in fields[4].split(','):
            if alt_allele != ".":
                alleles.append(alt_allele)
        gt_pos = fields[8].split(':').index('GT')
        if '|' in fields[8+args.sample].split(':')[gt_pos]:
            replace_allele = alleles[int(fields[8+args.sample].split(':')[gt_pos].split('|')[args.phase - 1])]
        else:
            replace_allele = alleles[int(fields[8+args.sample].split(':')[gt_pos].split('/')[0])]
        line2add = ref_dict[working_entry][last_position:pos-1] + replace_allele
        last_position = pos + len(alleles[0]) - 1
        if working_entry in pvc_dict:
            pvc_dict[working_entry] = pvc_dict[working_entry] + line2add
        else:
            pvc_dict[working_entry] = line2add


for scaf in input_scaffold_order:
    if scaf in pvc_dict:
        print '>' + scaf + '_pvc\n' + pvc_dict[scaf] + '\n'
    else:
        print '>' + scaf + '_pvc\n' + ref_dict[scaf] + '\n'

