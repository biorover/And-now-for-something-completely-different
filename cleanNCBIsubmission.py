#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description = "takes a contaminant file from \
    NCBI, a genome fasta, and an NCBI formatted tbl file and masks internal \
    contaminants, trims contaminants at the edges of scaffolds/contigs, and \
    adjusts coordinates of genes in the tbl file accordingly")

parser.add_argument('--fasta', help = 'genome fasta file')
parser.add_argument('--tbl', help = 'annotation tbl file')
parser.add_argument('--cont', help = 'NCBI conaminants report')

args = parser.parse_args()

actions = {}
readline = False
for line in open(args.cont):
    if readline:
        locus = line.split()[0]
        seqlen = int(line.split()[1])
        coords = []
        for i in line.split()[2].split('..'):
            coords.append(int(i))
        if not locus in actions:
            actions[locus] = []
        actions[locus].append(coords + [seqlen])

print("contaminated region:")
for locus in actions:
    for action in actions[locus]:
        print([locus] + action)

fasta_out = open(args.fasta.replace('.fsa','') + ".decontaminated.fsa",'w')
tbl_out = open(args.tbl.replace('.tbl','') + ".decontaminated.tbl",'w')

fasta_seqs = {}
fasta_headers = []
for line in open(args.fasta):
    if line[0] == ">":
        fasta_seqs[line[1:-1]] = ""
        fasta_headers.append(line[1:-1])
    else:
        fasta_seqs[fasta_headers[-1]] += line[:-1]


for locus in actions:
    adjust = 0
    for action in actions[locus]:
        if action[0] == 1:
            fasta_seqs[locus] = fasta_seqs[locus][action[1]:]
            adjust += actions[locus][1]
        elif action[1] == action[2]:
            fasta_seqs[locus] = fasta_seqs[locus][:action[0] + 1 - adjust]
        else:
            newseq = ""
            for i in range(len(fasta_seq[locus])):
                if action[0] <= i + adjust <= action[1]:
                    newseq += "N"
                else:
                    newseq += fasta_seq[locus][i]
            fasta_seq[locus] = newseq

for locus in fasta_headers:
    fasta_out.write(">" + locus + "\n" + fasta_seqs[locus] + "\n")


for line in open(args.tbl):
    if ">Feature " in line:
        locus = line.split()[1]
        start_offset = 0
        end_clip = len(fasta_seqs[locus])
        if locus in actions:
            for action in actions[locus]:
                if action[0] == 1:
                    start_offset = action[1]
                elif action[1] == action[2]:
                    end_clip = action[0] - start_offset
        tbl_out.write(line)
        write_features = True
    elif line[0] == "\t" and write_features:
        tbl_out.write(line)
    else:
        fields = line.split('\t')
        start = int(fields[0])
        end = int(fields[1])
        if len(fields) < 2:
            if "gene" in fields[2]:
                if start < start_offset or end > end_clip:
                    write_features = False
                else:
                    write_features = True
            elif "REFERENCE" in fields[2]:
                end = end_clip
                fields[0] = str(start)
                fields[1] = str(end)
                tbl_out.write('\t'.join(fields))
                continue
        if write_features:
            start = start + start_offset
            end = end + start_offset
            fields[0] = str(start)
            fields[1] = str(end)
            if len(fields) == 2:
                fields[1] += "\n"
            tbl_out.write('\t'.join(fields))
