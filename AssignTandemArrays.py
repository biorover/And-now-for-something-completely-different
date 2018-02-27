#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='analyzes a gff3 (really only looks at the "gene" line) to find \
                                 continuous stretches of related genes (with some sort of identifying \
                                 characters in their names)')

parser.add_argument('-g','--gff',dest='gff',help = 'gff3 file with gene names', type=argparse.FileType('r'))
parser.add_argument('-s','--string', dest='string', help = 'string present in gene family ID but not other genes\' IDs')
parser.add_argument('-i','--interruptions',dest='interruptions', help = 'number of interrupting genes to break\
                    tandem arrays (default 2)',type = int, default = 2)


args = parser.parse_args()

chromedict = {}

for line in args.gff:
    if line[0] != '#':
        fields = line.replace('\n','').replace('\r','').split('\t')
        if len(fields) > 7:
            if fields[2] == 'gene':
                coords = (int(fields[3]),int(fields[4]))
                ID = None
                for attribute in fields[8].split(';'):
                    x = attribute.split('=')
                    if x[0] == "ID":
                        ID = x[1]
                if ID == None:
                    print "crud, there's a line without an ID"
                    exit()
                chrome = fields[0]
                if chrome in chromedict:
                    chromedict[chrome].append([coords,ID])
                else:
                    chromedict[chrome] = [[coords,ID]]

for chrome in chromedict:
    chromedict[chrome].sort()

outlist = []

interruption_counter = 0
working_array = 0
last_in = False

for chrome in chromedict:
    new_chrome = True
    for gene in chromedict[chrome]:
        if args.string in gene[1]:
            if new_chrome and last_in:
                working_array += 1
            interruption_counter = 0
            outlist.append((gene[1],str(working_array)))
            last_in = True
            new_chrome = False
        else:
            interruption_counter += 1
            if interruption_counter == args.interruptions:
                working_array += 1
                new_chrome = False
                last_in = False



for i in outlist:
    print "\t".join(i)




