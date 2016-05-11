### Chromatin accessibility maps of chronic lymphocytic leukemia identify subtype-specific epigenome signatures and transcription regulatory networks

André F. Rendeiro<sup>\*</sup>, Christian Schmidl<sup>\*</sup>, Jonathan C. Strefford<sup>\*</sup>, Renata Walewska, Zadie Davis, Matthias Farlik, David Oscier, Christoph Bock
Chromatin accessibility maps of chronic lymphocytic leukemia identify subtype-specific epigenome signatures and transcription regulatory networks.

<sup>\*</sup>Shared first authors

**Paper**: [http://dx.doi.org/10.1038/nmeth.3542](http://dx.doi.org/10.1038/nmeth.3542)

**Website**: [cll-chromatin.computational-epigenetics.org](http://cll-chromatin.computational-epigenetics.org)

This repository contains scripts used in the analysis of the data in the paper.

<br>

#### Manuscript
The manuscript is written in [scholarly markdown](http://scholarlymarkdown.com/), therefore you need [Scholdoc](https://github.com/timtylin/scholdoc) to render the markdown into a pdf, rst, word or html manuscript.

A rendered version is available [here](manuscript/main.pdf).

You can [see the raw manuscript here](manuscript/main.md) along with the [figures](manuscript/figures/).

To render the pdf version of the manuscript, run:
```
make manuscript
```
this requires in general a full latex installation.

<br>

#### Analysis

The [`data`](data/external/) directory contains most of the output of the whole analysis.

Here are a few steps needed to reproduce it (more than I'd want to, I admit):

1. Clone the repository: `git clone git@github.com:epigen/cll-chromatin.git`
2. Install required software for the analysis:`make requirements` or `pip install -r requirements.txt`

If you wish to reproduce the processing of the raw data (access has to be requested through [EGA](https://www.ebi.ac.uk/ega/datasets/)), run these steps:

1. Apply for access to the raw data from [EGA](https://www.ebi.ac.uk/ega/datasets/).
2. Download the data localy.
3. Prepare [Looper](https://github.com/epigen/looper) configuration files similar to [these](metadata/project_config_file.yaml) that fit your local system.
4. Run samples through the pipeline: `make preprocessing` or `looper -c metadata/project_config_file.yaml`
5. Get external files (genome annotations mostly): `make external_files` or use the files in [data/external](data/external/).
6. Run the analysis: `make analysis`

Additionaly, processed (bigWig and narrowPeak files together with a gene expression matrix) are available from [GEO with accession number GSE81274](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81274).

If you wish to reproduce the plots from the analysis you can, in principle:

1. run `python src/analysis.py`

Not all parts of the analysis are possible to run *as is*, though. The TF network interence is based on a R package (PIQ) which is really hard to script runs in a system-independent way.
