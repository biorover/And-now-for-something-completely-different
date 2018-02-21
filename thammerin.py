#!/usr/bin/env python

#Wraps hmmsearch to give the user a crude ability to conduct an HMMER search of a protein query
#against a nucleotide database.

#One day, one glorious day, the HMMER developers will write a tool to do this correctly. If you're frustrated
#with this script, join with me to bug them to add this feature themselves. I hear rumors that HMMER4 will be out
#soon- would it be too much to hope for that it will contain this function?

import subprocess
import genome
import argparse

parser = argparse.ArgumentParser(description='thammerin\nPronounced Tee hammerin\', so pronounced identically \
                                 to tHMMERn but does NOT infinge on the copyrights of HHMI (HMMER) or \
                                 NCBI (tblastn)\nWraps hmmsearch to give the user a crude ability to conduct \
                                 a HMMER search of a protein query against a nucleotide database.\n\nProgram\
                                 dependencies:python (obviously),hmmer suite v3,genome library from MAGOT')

parser.add_argument('-n','--nucleotide_seqs',dest='target_nucdb',
                    help = 'Nucleotide sequences to be searched (in fasta format)')

parser.add_argument('-p','--protein_hmm', dest = 'hmm_file',
                    help = 'Query protein hmms (created using the hmmbuild program)')

parser.add_argument('-e','--iEvalue_threshold', dest = 'iEval_thresh', default = float('inf'),
                    help = 'threshold "i-Evalue" for inclusion in results')

parser.add_argument('-s','--min_orf_size', dest = 'min_orf_size', default = 10,
                    help = 'threshold "i-Evalue" for inclusion in results')

parser.add_argument('-f','--hmmsearch_filepath',dest = 'hmmsearch_filepath', default = 'hmmsearch',
                    help = 'optional path to hmmsearch program in not in your PATH variable')

args=parser.parse_args()

#sets up default program paths, overwritten by any program paths passes with -f or --program_filepaths
hmmsearch = args.hmmsearch_filepath


target_nucdb = genome.Genome(args.target_nucdb)

#Gets ORFs from the genome and hmmers them

for seq_id in target_nucdb.genome_sequence:
    frameonef = genome.Sequence(target_nucdb.genome_sequence[seq_id])
    frameoner = genome.Sequence(target_nucdb.genome_sequence[seq_id]).reverse_compliment()
    frames = [frameonef,genome.Sequence(frameonef[1:]), genome.Sequence(frameonef[2:]),
              frameoner, genome.Sequence(frameoner[1:]),genome.Sequence(frameoner[2:])]
    fasta_list = []
    for frame_num in (0,1,2,3,4,5):
        frame = frames[frame_num]
        if frame_num < 3:
            frame_offset = frame_num
        else:
            frame_offset = len(frame) - frame_num + 3
        orfs = frame.translate().split('*')
        last_orf_end = 0
        for orf in orfs:
            orf_start = last_orf_end
            last_orf_end = last_orf_end + 1 + len(orf)
            if len(orf) > args.min_orf_size:
                fasta_list.append('>' + seq_id + "_frameOffset-" + str(frame_offset) + "_orfStart-" + str(orf_start) +
                                  '\n' + orf)
    out = open('tmp_thammerin_frames.fa','w')
    out.write('\n'.join(fasta_list))
    out.close()
    subprocess.call(hmmsearch + " --domtblout tmp_thammerin_hmmresults.tab " + args.hmm_file + " tmp_thammerin_frames.fa >> tmp_thammer_log.log", shell=True)
    hmmresults = open('tmp_thammerin_hmmresults.tab')
    for line in hmmresults:
        if line[0] != "#":
            fields = line.split()
            if float(fields[12]) < args.iEval_thresh:
                tname = fields[0]
                frame_offset = int(tname.split('_frameOffset-')[1].split('_orfStart-')[0])
                orf_start = int(tname.split('_orfStart-')[1])
                if frame_offset < 3:
                    strand = "+"
                    start = orf_start * 3 + frame_offset + int(fields[17]) * 3 - 2
                    stop = orf_start * 3 + frame_offset + int(fields[18]) * 3
                else:
                    strand = "-"
                    start = frame_offset - 3 * orf_start - 3 * int(fields[17]) + 4
                    stop = frame_offset - 3 * orf_start - 3 * int(fields[18]) + 2
                print "\t".join([fields[3],tname.split('_frameOffset-')[0],fields[21],'.','.',strand,
                                 fields[15],fields[16],str(start),str(stop),fields[12],fields[13]])




