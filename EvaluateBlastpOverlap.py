#!/usr/bin/env python
#Evaluate blastp to examine ratio of overlap between a
#query gene and it's best hit, and vice-versa

import sys
import numpy

target_dict = {}
target_list = []
query_list = []

#debug
line_buffer = []
#/debug

first_query = True
skipline = False
for rawline in open(sys.argv[1]):
    #debug
    line_buffer.append(rawline)
    if len(line_buffer) > 500:
        line_buffer.pop(0)
    #/debug
    line = rawline.replace('\n','').replace('\r','')
    if line == "":
        continue
    elif line[:6] == "Query=":
        if first_query:
            pass
        else:
            if hit_coords:
                hit_coords[1] = last_coord
                qrange[hit_coords[0]:hit_coords[1]] = numpy.zeros(hit_coords[1] - hit_coords[0]) + 1
                thit_coords[1] = last_tcoord
                trange[thit_coords[0]:thit_coords[1]] = numpy.zeros(thit_coords[1] - thit_coords[0]) + 1
            overlap_len = qrange.sum()
            toverlap_len = trange.sum()
            qoverlap = 100.0 * overlap_len / qlen
            toverlap = 100.0 * toverlap_len / tlen
            if first_hit:
                query_list[-1][1] = qoverlap
                if toverlap > target_dict[tname][1]:
                    target_dict[tname][1] = toverlap
                target_list.append([tname,toverlap])
                #debug
                #if tname == "CG10186":
                #    print toverlap_len
                #/debug
        qname = line.split()[1]
        query_list.append([qname,0])
        first_hit = False
        skipline = False
        hit_coords = False
        thit_coords = False
    elif skipline:
        continue
    elif "> " in line and not first_hit:
        first_query = False
        tname = line.split()[1]
        first_hit = True
        if not tname in target_dict:
            target_dict[tname] = [tname,0]
    elif "> " in line and first_hit:
        skipline = True
    elif line[:7] == "Length=":
        if not first_hit:
            qlen = int(line.split('=')[1])
            qrange = numpy.zeros(qlen)
        else:
            tlen = int(line.split('=')[1])
            trange = numpy.zeros(tlen)
    elif "Score =" in line:
        if hit_coords:    
            skipline = True
            hit_coords[1] = last_coord
            qrange[hit_coords[0]:hit_coords[1]] = numpy.zeros(hit_coords[1] - hit_coords[0]) + 1
            thit_coords[1] = last_tcoord
            trange[thit_coords[0]:thit_coords[1]] = numpy.zeros(thit_coords[1] - thit_coords[0]) + 1
        else:
            hit_coords = [0,0]
            thit_coords = [0,0]
    elif line[:7] == "Query  ":
        if hit_coords[0] == 0:
            hit_coords[0] = int(line.split()[1]) - 1
        last_coord = int(line.split()[3])
    elif line[:7] == "Sbjct  ":
        if thit_coords[0] == 0:
            thit_coords[0] = int(line.split()[1]) - 1
        last_tcoord = int(line.split()[3])

percentile_dict = {}
percentiles = [10,20,30,40,50,60,70,80,90,95,98,99,100]
for percentile in percentiles:
    percentile_dict[percentile] = [0,0,0]

for query in query_list:
    for percentile in percentiles:
        if query[1] > percentile:
            percentile_dict[percentile][0] += 1

for target in target_dict.values():
    for percentile in percentiles:
        if target[1] > percentile:
            percentile_dict[percentile][1] += 1

for target in target_list:
    for percentile in percentiles:
        if target[1] > percentile:
            percentile_dict[percentile][2] += 1

print len(query_list)
for percentile in percentiles:
    print str(percentile) + ":\t" + str(percentile_dict[percentile][0]) + '\t' + str(percentile_dict[percentile][1]) + '\t' + str(percentile_dict[percentile][2])


    