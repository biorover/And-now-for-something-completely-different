#!/usr/bin/env python
import pandas as pd
import seaborn as sns
import sys
import numpy as np
import copy
import matplotlib.pyplot as plt
import os
import copy
import ete3
ncbi = ete3.NCBITaxa()
import argparse

parser = argparse.ArgumentParser(description='Blobology-inspired taxon annation tool')

parser.add_argument('--out_prefix',help = 'file name prefix for output files', default = 'BITAT')
parser.add_argument('--hit_table',help = 'tab-delimited diamond/blast output')
parser.add_argument('--coverage_file',default = None, help = 'tab-delimited file with coverage for each contig')
#
##Functionality to be added later: build depth table from bam instead of coverage file
#parser.add_argument('--bam_file', default = None help = 'bam file or aligned reads (if coverage file not alread available)')
#
parser.add_argument('--genome_fasta',help = 'fasta file of the genome assembly', default = None)
parser.add_argument('--qid2taxid', help = 'tab-delimited file from the NCBI taxonomy ftp site converting protein query IDs to ')
parser.add_argument('--custom_hit_hierarchy',default = None, help = 'linean taxanomic hierarchy of hits from custom database: \
                    "kingdom,phylum,class,order,family,genus"')
parser.add_argument('--custom_hit_id', default = None, help = 'tag in gene id that indicates the hit is from a custom database')
parser.add_argument('--filtered_hit_table', help = 'comma-delimited diamond/blast output that has already undergone filtering to keep only one hit per window', default = None)
parser.add_argument('--filtered_annotated_hit_table', help = 'comma-delimited diamond/blast out that has already undergone filtering and taxon annotation',default = None)
parser.add_argument('--plot_max_depth', default = None, type = int, help = 'max depth for plot (integer; default = {max depth in data})')
parser.add_argument('--plot_taxon_level',default = 'phyl', help = 'taxonomic at which plot points are colored. Options are \
                    "king","phyl","class", "order", "family", and "genus" (default = "phyl")')
parser.add_argument('--min_depth',default = 0,help = 'minimum sequencing depth for a contig to be included on plot (default = 0)',type = int)

args = parser.parse_args()

def filter_hit_table(hit_table,out_prefix):
    blast_df_raw = pd.read_csv(hit_table,sep="\t",header=None)
    blast_df_raw.columns = ['contig','hit','pid','yeah','hum','dontknow',
        'qstart','qend','sstart','ssend','evalue','score']
    sys.stderr.write('raw blast table shape: ' + str(blast_df_raw.shape) + '\n')

    blast_df = copy.deepcopy(blast_df_raw)

    for contig in set(blast_df['contig']):
        contig_df = blast_df[blast_df['contig'] == contig]
        max_val = max(contig_df['qstart'])
        for window in range(0,max_val,2000):
            window_df = contig_df[(contig_df['qstart'] > window) & \
                (contig_df['qstart'] < window + 2000)]
            chop_list = list(window_df.index)
            if len(chop_list) > 1:
                chop_list.remove(window_df['score'].idxmax())
                blast_df = blast_df.drop(chop_list)
    blast_df.to_csv(out_prefix + '.filtered_hits.csv')
    return blast_df

def annotate_hit_table(blast_df, out_prefix, qid2taxid, custom_hit_hierarchy = None, custom_hit_id = None):
    sys.stderr.write('filtered blast table shape: ' + str(blast_df.shape) + '\n')
    
    seen_lineages = {}
    seen_ranks = {}
    seen_names = {}

    qid2taxid_file = open(qid2taxid)
    qid2taxid_dict = {}
    for line in qid2taxid_file:
        fields = line.split()
        if len(fields) == 3:
            if fields[1] in list(blast_df['hit']):
                qid2taxid_dict[fields[1]] = int(fields[2])
        elif len(fields) == 2:
            if fields[0] in list(blast_df['hit']):
                qid2taxid_dict[fields[0]] = int(fields[1])

    if custom_hit_hierarchy:
        custom_hits_are_from = custom_hit_hierarchy.split(',')

    if not custom_hit_id or not custom_hit_hierarchy:
        custom_hit_id = "NoWayThisIsInAnyFastaHeader$#@($*&^%!"

    blast_df['king'] = "Unk"
    blast_df['phyl'] = "Unk"
    blast_df['class'] = "Unk"
    blast_df['order'] = "Unk"
    blast_df['family'] = "Unk"
    blast_df['genus'] = "Unk"
    blast_df['is_custom'] = False
    error_counts = 0
    for index,row in blast_df.iterrows():
        if custom_hit_id in row['hit']:
            blast_df.at[index,'king'],blast_df.at[index,'phyl'],blast_df.at[index,'class'], blast_df.at[index,'order'],blast_df.at[index,'family'],blast_df.at[index,'genus'] = custom_hits_are_from
            blast_df.at[index,'is_custom'] = True
        else:
            qid = row['hit']
            if qid in qid2taxid_dict:
                if qid2taxid_dict[qid] in seen_lineages:
                    lineage = seen_lineages[qid2taxid_dict[qid]]
                    ranks = seen_ranks[qid2taxid_dict[qid]]
                    names = seen_names[qid2taxid_dict[qid]]
                else:
                    lineage = ncbi.get_lineage(qid2taxid_dict[qid])
                    seen_lineages[qid2taxid_dict[qid]] = lineage
                    ranks = ncbi.get_rank(lineage)
                    seen_ranks[qid2taxid_dict[qid]] = ranks
                    names = ncbi.get_taxid_translator(lineage)
                    seen_names[qid2taxid_dict[qid]] = names
                superkingdom = False
                kingdom = False
                phylum = False
                taxclass = False
                order = False
                family = False
                genus = False
                for i in lineage:
                    if ranks[i] == "superkingdom":
                        superkingdom = names[i]
                    elif ranks[i] == 'kingdom':
                        kingdom = names[i]
                    elif ranks[i] == 'phylum':
                        phylum = names[i]
                    elif ranks[i] == 'class':
                        taxclass = names[i]
                    elif ranks[i] == 'order':
                        order = names[i]
                    elif ranks[i] == 'family':
                        family = names[i]
                    elif ranks[i] == 'genus':
                        genus = names[i]
                if kingdom:
                    blast_df.at[index,'king'] = kingdom
                elif superkingdom:
                    blast_df.at[index,'king'] = superkingdom
                else:
                    blast_df.at[index,'king'] = "kingdom_incertae_sedis"
                if phylum:
                    blast_df.at[index,'phyl'] = phylum
                else:
                    blast_df.at[index,'phyl'] = blast_df.at[index,'king'] + "_pIS"
                if taxclass:
                    blast_df.at[index,'class'] = taxclass
                else:
                    blast_df.at[index,'class'] = blast_df.at[index,'phyl'] + "_cIS"
                if order:
                    blast_df.at[index,'order'] = order
                else:
                    blast_df.at[index,'order'] = blast_df.at[index,'class'] + "_oIS"
                if family:
                    blast_df.at[index,'family'] = family
                else:
                    blast_df.at[index,'family'] = blast_df.at[index,'order'] + "_fIS"
                if genus:
                    blast_df.at[index,'genus'] = genus
                else:
                    blast_df.at[index,'genus'] = blast_df.at[index,'family'] + "_gIS"
            else:
                error_counts += 1

    blast_df.to_csv(out_prefix + '.filtered_taxon_annotated_hits.csv')
    return blast_df

def build_tig_table(genome_fasta,out_prefix,blast_df):
    tigs = []
    tigindi = []
    for line in open(genome_fasta):
        if line[0] == ">":
            tigs.append([line[1:-1],""])
            tigindi.append(line[1:-1].split()[0])
        else:
            tigs[-1][1] += line[:-1]

    rdf = pd.DataFrame(tigs,tigindi)
    rdf.rename(columns={0:"name",1:"sequence"},inplace=True)
    rdf['GC'] = 0.0
    rdf['tiglen'] = 0
    for indi, entries in rdf.iterrows():
        index = indi
        rdf.at[index,'tiglen'] = len(rdf['sequence'][index])
        rdf.at[index,'GC'] = (1.0 * rdf['sequence'][index].count('C') + rdf['sequence'][index].count('G') ) / len(rdf['sequence'][index])

    rdf['king'] = "Unk"
    rdf['phyl'] = "Unk"
    rdf['class'] = "Unk"
    rdf['order'] = "Unk"
    rdf['family'] = 'Unk'
    rdf['genus'] = 'Unk'
    for contig in list(rdf.index):
        if contig in list(blast_df['contig']):
            contig_df = blast_df[blast_df['contig'] == contig]
            if len(contig_df) > 2:
                if len(contig_df['king'].value_counts()) > 0:
                    rdf.at[contig,'king'] = contig_df['king'].value_counts().index[0]
                    rdf.at[contig,'phyl'] = contig_df[contig_df['king'] == rdf.at[contig,'king']]['phyl'].value_counts().index[0]
                    rdf.at[contig,'class'] = contig_df[contig_df['phyl'] == rdf.at[contig,'phyl']]['class'].value_counts().index[0]
                    rdf.at[contig,'order'] = contig_df[contig_df['class'] == rdf.at[contig,'class']]['order'].value_counts().index[0]
                    rdf.at[contig,'family'] = contig_df[contig_df['order'] == rdf.at[contig,'order']]['family'].value_counts().index[0]
                    rdf.at[contig,'genus'] = contig_df[contig_df['family'] == rdf.at[contig,'family']]['genus'].value_counts().index[0]
            elif len(contig_df) > 0:
                hit_index = contig_df['score'].idxmax()
                rdf.at[contig,'king'] = contig_df.at[hit_index,'king']
                rdf.at[contig,'phyl'] = contig_df.at[hit_index,'phyl']
                rdf.at[contig,'class'] = contig_df.at[hit_index,'class']
                rdf.at[contig,'order'] = contig_df.at[hit_index,'order']
                rdf.at[contig,'family'] = contig_df.at[hit_index,'family']
                rdf.at[contig,'genus'] = contig_df.at[hit_index,'genus']

    rdf.to_csv(out_prefix + '.annotated_tig_table.tab')
    return rdf

def add_depths(tig_table,coverage_file):
    tig_table["Depth"] = 0.0
    for line in open(coverage_file):
        fields = line.split()
        if len(fields) > 1:
            tig_table.at[fields[0],'Depth'] = float(fields[1])
    return tig_table

def plot(tig_table, taxa_level, max_depth,out_prefix,min_depth):
    if max_depth:
        plotdf = tig_table[tig_table["Depth"] < max_depth]
        plotdf = tig_table[tig_table["Depth"] > min_depth]
    else:
        plotdf = tig_table[tig_table["Depth"] > min_depth]
    grid = sns.JointGrid(x='GC', y='Depth', data=plotdf)
    taxalist = list(set(list(plotdf[taxa_level])))
    taxalist.sort()
    g = grid.plot_joint(sns.scatterplot, hue=taxa_level, data=plotdf,markers=["."])
    for taxon in taxalist:
        sns.distplot(plotdf.loc[plotdf[taxa_level] == taxon,"Depth"],kde=False, ax=g.ax_marg_y,vertical=True,bins = range(0,int(max(plotdf["Depth"])),3),axlabel=False)
        sns.distplot(plotdf.loc[plotdf[taxa_level]==taxon, 'GC'], ax=g.ax_marg_x, bins = np.arange(0.1,0.8,0.01),kde=False,axlabel=False)
    g.savefig(out_prefix + ".plot.png")

def main():
    if not args.filtered_hit_table and not args.filtered_annotated_hit_table:
        blast_df = filter_hit_table(args.hit_table,args.out_prefix)
    elif not args.filtered_annotated_hit_table:
        blast_df = pd.read_csv(args.filtered_hit_table)
    if not args.filtered_annotated_hit_table:
        blast_df = annotate_hit_table(blast_df,args.out_prefix,args.qid2taxid,args.custom_hit_hierarchy, args.custom_hit_id)
    else:
        sys.stderr.write("reading filtered annotated hit table\n")
        blast_df = pd.read_csv(args.filtered_annotated_hit_table)
    if args.genome_fasta:
        rdf = build_tig_table(args.genome_fasta,args.out_prefix,blast_df)
    if args.coverage_file:
        rdf = add_depths(rdf,args.coverage_file)
        plot(rdf,args.plot_taxon_level,args.plot_max_depth,args.out_prefix,args.min_depth)

if __name__ == "__main__":
    main()
