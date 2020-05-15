#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description='takes "interactions" file from \
    SAINT run (classic, not express) and uses it to remove spurious interactions\
    from a maxquant proteinGroups file')
parser.add_argument('--pg',help = 'maxquant proteinGroups file')
parser.add_argument('--int', help = 'SAINT interactions file')
parser.add_argument('--ctl', default = 'ctl', help = 'text designation in \
    control columns. Default = "ctl"')
parser.add_argument('--rep_sep', default = "_", help = 'separator between\
    bait name and replicate number. Set to "none" if no biological replicates. \
    Default = "_"')
args = parser.parse_args()

false_positives = set()
for line in open(args.int):
    fields = line.split('\t')
    if float(fields[5]) < 0.95:
        false_positives.add(fields[1] + "\t" + fields[2])


intensity_cols = []
baits = []
for line in open(args.pg):
    fields = line.split('\t')
    if intensity_cols == []:
        for i in range(len(fields)):
            field = fields[i]
            if (field[:10] == "Intensity "):
                intensity_cols.append(i)
                if not args.ctl in field:
                    intensity_cols.append(i)
                    if args.rep_sep == "none":
                        baits.append(field[field.find(" ") + 1:])
                    else:
                        baits.append(args.rep_sep.join(field[field.find(" ") + \
                            1:].split(args.rep_sep)[:-1]))
    elif len(fields) < len(baits):
        continue
    else:
        if len(fields[1]) > 99:
            fields[1] = fields[1][:99] #these names are truncated in my code to convert maxquant to saint input because otherwise saint errors. This should be removed if names are not truncated in saint results
        for i in range(len(intensity_cols)):
            intensity_index = intensity_cols[i]
            if fields[1] + '\t' + baits[i] in false_positives:
                fields[intensity_index] = "0"
    sys.stdout.write("\t".join(fields))
