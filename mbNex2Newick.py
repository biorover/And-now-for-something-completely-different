#!/usr/bin/env python
import sys

nex = open(sys.argv[1]).read().split(';')
translate = nex[1].split('\n')[2:]
tdict = {}
for line in translate:
    fields = line.split()
    #print fields
    taxnum = fields[0]
    taxname = fields[1]
    if taxname[-1] ==",":
        taxname = taxname[:-1]
    tdict[taxnum] = taxname

if len(sys.argv) > 2:
    burnin = int(sys.argv[2])
else:
    burnin = 0

trees = nex[2 + burnin:-2]
newtrees = []
for tree in trees:
    newick = tree.split()[-1]
    newtrees.append(newick)

newtrees = ";\n".join(newtrees) + ';'
for taxnum in tdict:
    newtrees = newtrees.replace(taxnum + ':', tdict[taxnum] + 'dummySeperator:')
newtrees = newtrees.replace('dummySeperator',"")

print newtrees

