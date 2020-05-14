#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='converts maxquant proteinGroups \
    tab-delimited output to files suitable for inputs into SAINTexpress')
parser.add_argument('--pg', help = 'proteinGroups table from maxquant')
parser.add_argument('--ctl', default = 'ctl', help = 'text designation in \
    control columns. Default = "ctl"')
parser.add_argument('--rep_sep', default = "_", help = 'separator between\
    bait name and replicate number. Set to "none" if no biological replicates. \
    Default = "_"')
parser.add_argument('--out_prefix', default = "saint", help = 'prefix for \
    output files')
parser.add_argument('--data_type', default = "intensity", help = '"intensity"\
    or "spectralCount" (default = "intensity")')
parser.add_argument('--saint_version', default = 'express', help = '"express" or "v2" (default = "express")')
args = parser.parse_args()

interaction_file = open(args.out_prefix + "_interactions.tab",'w')
prey_file = open(args.out_prefix + "_prey.tab",'w')
bait_file = open(args.out_prefix + "_bait.tab",'w')

ip_names = []
baits = []
intensity_cols = []
for line in open(args.pg):
    fields = line.split('\t')
    if ip_names == []:
        for i in range(len(fields)):
            field = fields[i]
            if (field[:10] == "Intensity " and args.data_type == 'intensity') \
            or (field[:9] == "Peptides " and args.data_type == 'spectralCount'):
                ip_names.append(field[field.find(" ") + 1:])
                intensity_cols.append(i)
                if args.ctl in field:
                    control_or_treatment = "C"
                else:
                    control_or_treatment = "T"
                if control_or_treatment == "C" or args.rep_sep == "none":
                    baits.append(field[field.find(" ") + 1:])
                else:
                    baits.append(args.rep_sep.join(field[field.find(" ") + 1:].split(args.rep_sep)[:-1]))
                bait_file.write('\t'.join([ip_names[-1],baits[-1],control_or_treatment]) + '\n')
            elif "Sequence length" == field:
                length_index = i
    elif len(fields) < len(ip_names):
        continue
    else:
        for i in range(len(intensity_cols)):
            intensity_index = intensity_cols[i]
            interaction_file.write("\t".join([ip_names[i],baits[i],fields[1],fields[intensity_index]]) + '\n')
        if args.saint_version == "express":
            prey_file.write('\t'.join([fields[1],fields[length_index],fields[6]]) + '\n')
        elif args.saint_version == "v2":
            prey_file.write('\t'.join([fields[1],fields[6]]) + '\n')
        else:
            print('saint_version must be one of "express" or "v2"')
            exit()
