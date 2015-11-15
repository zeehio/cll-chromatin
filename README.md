cll-patients
============

## Repository organization

A basic review of CLL literature can be found in `/metadata/literature/`.

Clinical, wetlab and sequenced sample annotations are in `/metadata`.

Code used in the analysis is in `/src`

### Data
Cohort clinical data: [CSV](metadata/clinical_annotation.csv)

See annotation sheet here: [CSV](metadata/sequencing_sample_annotation.csv)

See tracks here: [UCSC](http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&hgt.customText=http://www.biomedical-sequencing.at/bocklab/arendeiro/cll-patients/trackHub_hg19.txt)

### Sample names
Sample names are a concatenation of any of the existing annotation fields: `cellLine, numberCells, technique, ip, patientID, treatment, condition, biologicalReplicate, technicalReplicate, genome`.

### Project structure
As defined in [`pipelines`](https://github.com/afrendeiro/pipelines):

`projectsroot`: /data/groups/lab_bock/shared/projects/

`htmlroot`: /data/groups/lab_bock/public_html/arendeiro/


# Analysis

## Reproduce the analysis

1. Clone the repository: `git clone git@github.com:epigen/cll-patients.git`
2. Install required software for the analysis:`make requirements`
3. Run samples through the pipeline: `make preprocessing`
4. Get external files (genome annotations mostly): `make external_files`
5. Run the analysis: `make analysis`

This will produce intermediate files and the manuscript's panel figures.

# Manuscript
The manuscript is written in [scholarly markdown](http://scholarlymarkdown.com/), therefore you need [Scholdoc](https://github.com/timtylin/scholdoc) to render the markdown into a pdf, rst, word or html manuscript.

You can [see it here](manuscript/main.md) along with the figures [here](manuscript/figures/).

To render the pdf version of the manuscript, run:
```
make manuscript
```
this requires in general a full latex installation.
