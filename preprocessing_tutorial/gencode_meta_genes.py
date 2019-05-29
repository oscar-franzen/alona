#!/usr/bin/env python3
"""
  Merges exon intervals of transcripts in a GENCODE GTF genome annotation file. 
 
  How to use:
  1. download annotations
     wget ftp://ftp.ebi.ac.uk/pub/ \
     databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz
  2. gunzip gencode.vM21.annotation.gtf.gz
  3. python3 gencode_meta_genes.py > gencode.vM21.annotation.meta_genes.gtf
 
  Contact:
  Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import re
from collections import defaultdict

def merge_intervals(intervals):
    """ merge a list of intervals """
    loi = []

    for item in intervals:
        loi.append([ item['anno_start'], item['anno_stop'] ])

    loi.sort()

    res = []

    while (len(loi) > 0):
        if len(loi) == 1:
            res.append(loi[0])
            loi.pop(0)
            continue

        if loi[0][1] >= loi[1][0]:
            tmp = [loi[0][0], max(loi[0][1], loi[1][1])]
            loi[0] = tmp
            loi.pop(1)
            continue

        res.append(loi[0])
        loi.pop(0)

    return res

def main():
    coordinates = defaultdict(list)
    tx_types = defaultdict(defaultdict)

    with open('gencode.vM21.annotation.gtf', 'r') as f:
        for line in f:
            # ignore comments
            if line[0] == '#':
                continue
            # ignore exons with retained introns
            if re.search('retained_intron', line):
                continue

            foo = line.split('\t')

            chrm,anno_type,anno_start,anno_stop,anno_strand = foo[0],foo[2],foo[3],foo[4],foo[6]

            if (int(anno_start) - int(anno_stop)) == 0:
                continue

            if anno_type == 'exon':
                # collect identifiers
                gene_symbol = re.search('gene_name "(\S+)"', line).group(1)
                gene_ens = re.search('gene_id "(\S+)\..+"', line).group(1)

                gene_id = '%s_%s' % (gene_symbol, gene_ens)
                
                tx_type = re.search('transcript_type "(.*?)"', line).group(1)
                tx_types[gene_id][tx_type] = 1
                
                coordinates[gene_id].append({'anno_start' : anno_start,
                                             'anno_stop' : anno_stop,
                                             'chr' : chrm,
                                             'anno_strand' : anno_strand})
    
    for gene_id in coordinates:
        merged = merge_intervals(coordinates[gene_id])

        for item in merged:
            #print("$chr\tGENCODE.v27.custom\texon\t$start\t$stop\t1000\t
            #$strand\t.\tgene_id \"$gene\"; gene_name \"$gene\"; transcript_id
            #\"$gene\"; gene_type \"$gene_types\"; transcript_types \"$ttt\"\n");
            
            # transcript_id is not applicable here, however this parameter seems
            # necessary for the file to be viewable in IGV
            print('%s\tGENCODE.metagenes\texon\t%s\t%s\t1000\t%s\t.\tgene_id "%s"; \
gene_name "%s"; transcript_id "%s"; transcript_types "%s"' %
                (coordinates[gene_id][0]['chr'],
                item[0],
                item[1],
                coordinates[gene_id][0]['anno_strand'],
                gene_id,
                gene_id,
                gene_id,
                ','.join(sorted(tx_types[gene_id].keys()))))

if __name__ == '__main__':
    main()
