#!/usr/bin/env python
#Sean McKenzie, graduate fellow, The Rockefeller University. License can be found in git repository (https://github.com/biorover/And-now-for-something-completely-different)

#Script for taking Mr. Bayes tree output (Nexus files containing multiple newick trees with taxa symbolized with numbers and key
# given in Nexus header) and converting it to a newick file (only newick trees with taxa names).
help_prompt = """
Usage: python mbNex2Newick.py MrBayesOutputTrees [burnin]


MrBayesOutputTrees should be the Nexus file containing a series of trees, usually ending with ".run#.t"

burnin is the number of trees to discard from the beginning (NOT a percentage, you can calculate that yourself
by using 'grep -c "tree gen" fileName'
"""


import sys

if len(sys) < 2:
    print help_prompt
elif sys.argv[1] in ['h','-h','help','-help','--help','--h']:
    print help_prompt

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

