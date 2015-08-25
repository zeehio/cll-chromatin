#!/usr/bin/env python

import os

# Get blacklist regions
blacklist = "wgEncodeDacMapabilityConsensusExcludable.bed.gz"

os.system("wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/{0}".format(blacklist))
os.system("gzip -d {0}".format(blacklist))
os.system("mv wgEncodeDacMapabilityConsensusExcludable.bed ../data/wgEncodeDacMapabilityConsensusExcludable.bed")


# Get roadmap CD19 perypheral blood HMM state annotation
roadmap_url = "http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/"
roadmap_15statesHMM = "all.mnemonics.bedFiles.tgz"
roadmap_15statesHMM_CD19 = "E032_15_coreMarks_mnemonics.bed.gz"

os.system("wget {0}{1}".format(roadmap_url, roadmap_15statesHMM))
os.system("tar zxvf {0} {1}".format(roadmap_15statesHMM, roadmap_15statesHMM_CD19))
os.system("gzip -d {0}".format(roadmap_15statesHMM_CD19))
os.system("mv E032_15_coreMarks_mnemonics.bed ../data/E032_15_coreMarks_mnemonics.bed")
