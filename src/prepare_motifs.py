#!/usr/bin/env python

import os

# Get ENCODE motif matches
"wget http://compbio.mit.edu/encode-motifs/matches.txt.gz"  # get motifs
"gzip -d matches.txt.gz"  # extract
"tr ' ' \\t < matches.txt > matches.tsv"  # replace spaces with tabs
"""perl -ane 'print "$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[0]\n"' matches.tsv > matches.bed"""  # make bed file
"""bedtools merge -c 4 -o distinct -i matches.bed > matches.merged.bed"""  # Merge overlapping motifs

# Filter out non-interesting chromosomes

# Get window around center of motifs
