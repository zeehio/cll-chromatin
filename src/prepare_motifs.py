#!/usr/bin/env python

import os

cmds = list()

# Get ENCODE motif matches
cmds.append("wget http://compbio.mit.edu/encode-motifs/matches.txt.gz")  # get motifs
cmds.append("gzip -d matches.txt.gz")  # extract
cmds.append("tr ' ' \\t < matches.txt > matches.tsv")  # replace spaces with tabs
cmds.append("""perl -ane 'print "$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[0]\n"' matches.tsv > matches.bed""")  # make bed file
cmds.append("""bedtools merge -c 4 -o distinct -i matches.bed > matches.merged.bed""")  # Merge overlapping motifs

# Filter out non-interesting chromosomes/sites

# Run
for cmd in cmds:
    os.system(cmd)
