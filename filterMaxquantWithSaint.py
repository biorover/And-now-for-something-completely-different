#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description='takes "interactions" file from \
    SAINT run (classic, not express) and uses it to remove spurious interactions\
    from a maxquant proteinGroups file. Control sample names (as featured in \
    intensity column names) must all have a string not present in bait sample \
    names (specified via the "--ctl" argument, default is "ctl")')
required = parser.add_argument_group('required arguments')
required.add_argument('--pg',help = 'maxquant proteinGroups file', required = True)
required.add_argument('--int', help = 'SAINT interactions file', required = True)
parser.add_argument('--ctl', default = 'ctl', help = 'text designation in \
    control columns. Default = "ctl"')
parser.add_argument('--rep_sep', default = "_", help = 'separator between\
    bait name and replicate number. Set to "none" if no biological replicates. \
    Default = "_"')
parser.add_argument('--mode', default = "bait-specific", help = '"bait-specific" \
    or "all" (default = "bait-specific"): "bait-specific" mode sets bait-prey \
    intensities to zero if SAINT \
    found that specific interaction to be spurious. "all" mode drops row if \
    SAINT found that all prey interactions were spurious.')
parser.add_argument('--subtract', default = False, type = bool, help =
    '"True" or "False" (default = "False"): whether to subtract ')
parser.add_argument('--name_truncate_length', default = 99, type = int,
    help = 'int (default = 99): length after which to truncate protein group names \
    (added because SAINT errors on long protein names so these are usually \
    truncated before running SAINT)')
parser.add_argument('--pro_hits_out', default = None,
                    help = 'File to write results in a format for pro-hits visualization (default = None)')
args = parser.parse_args()

####
#
# Reads each line in SAINT interactions file and adds "bait<tab>prey" to
# false_positives set if probability of true interaction < 0.95
#
####

false_positives = set()
prob_dict = {}

for line in open(args.int):
    fields = line.split('\t')
    if fields[0] != "IP" and fields[1] != "Bait": #in order to skip header
        if float(fields[5]) < 0.95:
            false_positives.add(fields[1] + "\t" + fields[2])
        prob_dict[fields[1] + "\t" + fields[2]] = float(fields[5])

####
#
# Reads through MaxQuant proteinGroups table and either drops false-positive
# bait-prey interaction intensities to zero or drops rows where all interactions
# were deemed false positives, depending on mode. Also performs control
# subtraction if specified
#
####

if args.pro_hits_out:
    phout = open(args.pro_hits_out,'w')
    phout.write('BAIT\tPreyGene\tAvgIntensity\tfauxFDR\n')


intensity_cols = []
baits = []
control_cols = []
for line in open(args.pg):
    fields = line.split('\t')
    fieldscopy = fields[:]
    if intensity_cols == []: #checks if this is the header row (in which case intensity cols not yet filled in)
        for i in range(len(fields)):
            field = fields[i]
            if (field[:10] == "Intensity "):
                if not args.ctl in field:
                    intensity_cols.append(i)
                    if args.rep_sep == "none":
                        baits.append(field[field.find(" ") + 1:])
                    else:
                        baits.append(args.rep_sep.join(field[field.find(" ") + \
                            1:].split(args.rep_sep)[:-1]))
                else:
                    control_cols.append(i)
    elif len(fields) < len(baits): #checks that row is not empty
        continue
    else: #this is a row with data in it!
        if len(fields[1]) > args.name_truncate_length:
            fields[1] = fields[1][:args.name_truncate_length]
        if args.mode == 'bait-specific':
            keep_row = True
            for i in range(len(intensity_cols)):
                intensity_index = intensity_cols[i]
                if baits[i] + "\t" + fields[1] in false_positives:
                    fields[intensity_index] = "0"
        elif args.mode == 'all':
            keep_row = False #sets default to false, then we'll set it to true if any bait is missing from false positive list
            for bait in list(set(baits)):
                if not bait + "\t" + fields[1] in false_positives:
                    keep_row = True
        if args.pro_hits_out:
            bait_intensities = {}
            for i in range(len(intensity_cols)):
                intensity_index = intensity_cols[i]
                bait = baits[i]
                if not bait in bait_intensities:
                    bait_intensities[bait] = []
                bait_intensities[bait].append(float(fieldscopy[intensity_index]))
            for bait in bait_intensities:
                try:
                    score = str(1 - prob_dict[bait + "\t" + fields[1]])
                except:
                    score = "."
                avgint = sum(bait_intensities[bait]) / len(bait_intensities[bait])
                phout.write("\t".join([bait,fields[6],str(avgint),score]) + '\n')
        if not keep_row: #skips subtraction and writing output if dropping row
            continue
        if args.subtract: #subtracts mean control intensity from bait intensities, or sets bait intensities to zero if smaller than control
            ctrl_mean = 0
            for i in control_cols:
                ctrl_mean += float(fields[i]) / len(control_cols)
            for i in intensity_cols:
                if int(fields[i]) > int(ctrl_mean):
                    fields[i] = str(int(fields[i]) - int(ctrl_mean))
                else:
                    fields[i] = "0"

    sys.stdout.write("\t".join(fields))

sys.stdout.close()

if args.pro_hits_out:
    phout.close()
