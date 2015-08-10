#! /usr/bin/env python

import pandas as pd
import os
import urllib2
import StringIO
import gzip


blueprint_url = "http://ftp.ebi.ac.uk/pub/databases/blueprint/data/"
blueprint_base = "resources/regions/customRegionDB/hg19/blueprint"

if not os.path.exists(blueprint_base):
    os.makedirs(blueprint_base)
    os.makedirs(os.path.join(blueprint_base, "regions"))

# read in table with all files
df = pd.read_csv("http://ftp.ebi.ac.uk/pub/databases/blueprint/releases/current_release/homo_sapiens/20150128.data.index", sep="\t")

df2 = df[
    (df['LIBRARY_STRATEGY'] == "ChIP-Seq") &
    (df['FILE'].str.contains(".bed.gz")) &
    (df['EXPERIMENT_TYPE'].str.contains("H3|H2|Chromatin Accessibility"))
].drop_duplicates()

# get server url
df2['links'] = df2['FILE'].apply(lambda x: blueprint_url + x.replace("blueprint/data/", ""))

# get local filenames
df2['filename'] = df2['links'].apply(lambda x: os.path.basename(x.replace("bed.gz", "bed")))

# set species
df2['species'] = "Human"

# subset, make index for lola
df3 = df2[["species", "CELL_TYPE", "CELL_LINE", "EXPERIMENT_TYPE", "TREATMENT", 'filename']]
df3.columns = ["species", "cellType", "tissue", "antibody", "treatment", "filename"]
df3 = df3.replace("-", "")
df3.to_csv(os.path.join(blueprint_base, "index.txt"), index=False, sep="\t")


# download bed files
# decompress
# save in lola db
for link in df2['links']:
    bed_file = df2[df2['links'] == link]['filename'].tolist()[0]

    response = urllib2.urlopen(link)
    compressedFile = StringIO.StringIO()
    compressedFile.write(response.read())
    compressedFile.seek(0)

    decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='rb')

    with open(os.path.join(blueprint_base, "regions", bed_file), 'w') as handle:
        handle.write(decompressedFile.read())


# Make collection file
lines = list()
lines.append("collector\tdate\tsource\tdescription\n")
lines.append("André F. Rendeiro\t2015-08-10\tBlueprint Consortia\tDownloaded from %s. (Script in get_blueprint_regions.py)\n" % blueprint_url)

with open(os.path.join(blueprint_base, "collection.txt"), 'w') as handle:
    for line in lines:
        handle.write(line)
