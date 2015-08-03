cll-patients
============

# Repository organization

A basic review of CLL literature can be found in `/metadata/literature/`.

Clinical, wetlab and sequenced sample annotations are in `/metadata`.

Code used in the analysis is in `/src`

### Data
See annotation sheet here: [CSV](metadata/sequencing_sample_annotation.csv)

See tracks here: [UCSC](http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&hgt.customText=http://www.biomedical-sequencing.at/bocklab/arendeiro/cll-patients/trackHub_hg19.txt)

### Sample names
Sample names are a concatenation of any of the existing annotation fields: `cellLine, numberCells, technique, ip, patientID, treatment, condition, biologicalReplicate, technicalReplicate, genome`.

### Project structure
As defined in [`pipelines`](https://github.com/afrendeiro/pipelines):

For all my projects:

`projectsroot`: /data/groups/lab_bock/shared/projects/

`htmlroot`: /data/groups/lab_bock/public_html/arendeiro/


# Todo

## Patient description
See [annotation](metadata/patient_clinical_annotation.csv)

## Data production
+ Open-chromatin
    + ATAC-seq
+ Histone mods (possible):
    + H3K27AC
    + H3K4ME1
    + H3K4ME3?
+ RNA?

## Analysis (ideas)
+ Accounting for genetic variability:
    + http://www.nature.com/nmeth/journal/v12/n5/full/nmeth.3326.html
    + http://biorxiv.org/content/early/2015/04/30/018788

##### Aditional ideas using only ATAC-seq data:
+ Call nucleosomes and dyads using NucleoATAC
+ Look for global nucleosome positioning differences between clusters and compared with controls
+ Look for global nucleosome positioning differences within groups of gene classes (e.g. genes involved in cellular proliferation, cell cycle, B-cell biology and other classes).

#### Several ideas for downstream:
+ GRNs
    + Enhancer-TF pairing
    + With evolution build GRN
+ De novo enhancers
    + Comparison with other cell lines
+ Variants
    + Call variants
    + check GWAS enrichments
+ Look at time-resolved data for cases available:
    + Poised-enhancer activation later and vice-versa
    + Clonal evolution within patients
+ Global evolution patterns
    + Do certain type of genes get enriched in time over all/clusters of samples? (test all functional classes, correct p*values)
+ Finding new biomarkers
    + Look at earliest points (if possible before diagnosis) for consistent activation

