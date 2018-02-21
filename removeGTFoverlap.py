#!/usr/bin/env python

#!/usr/bin/env python

#selects longest non-overlapping transcripts from gtf

#assumes sane GTF (only CDS features, all features belonging to same transcript in order,
#defline format: "transcript_id <trancript_id> [; etc.]")

import sys

chunksize = 100000
seqlens = []
lines_dict = {}
coords_dict = {}


for line in open(sys.argv[1]):
    if line[0] != '#' and line.count('\t') > 6:
        fields = line.split('\t')
        transcript_id = fields[-1].split(';')[0].split()[1]
        start, stop = int(fields[3]), int(fields[4])
        locus = fields[0]
        coords = (start,stop)
        start_key = start / chunksize
        stop_key = stop / chunksize
        keys_list = [(locus,start_key)]
        if start_key != stop_key:
            for i in range(stop_key - start_key):
                keys_list.append((locus,keys_list[-1][-1] + 1))
        for key in keys_list:
            if not key in coords_dict:
                coords_dict[key] = []
            for coord in coords:
                if coord / chunksize == key[-1]:
                    coords_dict[key].append(coords)

for line in open(sys.argv[2]):
    if line[0] != '#' and line.count('\t') > 6:
        fields = line.split('\t')
        transcript_id = fields[-1].split(';')[0].split()[1]
        start, stop = int(fields[3]), int(fields[4])
        locus = fields[0]
        start_key = start / chunksize
        stop_key = stop / chunksize
        keys_list = [(locus,start_key)]
        if start_key != stop_key:
            for i in range(stop_key - start_key):
                keys_list.append((locus,keys_list[-1][-1] + 1))
        overlap = False
        for key in keys_list:
            if key in coords_dict:
                for occupied_coord in coords_dict[key]:
                    if occupied_coord[0] <= start <= occupied_coord[1] or \
                    occupied_coord[0] <= stop <= occupied_coord[1] or \
                    start < occupied_coord[0] < stop:
                        overlap = True
                        break
        if not overlap:
            print line[:-1]
