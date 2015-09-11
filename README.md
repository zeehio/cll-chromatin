cll-patients
============

# Repository organization

A basic review of CLL literature can be found in `/metadata/literature/`.

Clinical, wetlab and sequenced sample annotations are in `/metadata`.

Code used in the analysis is in `/src`

### Data
Cohort clinical data: [CSV](metadata/patient_clinical_annotation.csv)

See annotation sheet here: [CSV](metadata/sequencing_sample_annotation.csv)

See tracks here: [UCSC](http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&hgt.customText=http://www.biomedical-sequencing.at/bocklab/arendeiro/cll-patients/trackHub_hg19.txt)

### Sample names
Sample names are a concatenation of any of the existing annotation fields: `cellLine, numberCells, technique, ip, patientID, treatment, condition, biologicalReplicate, technicalReplicate, genome`.

### Project structure
As defined in [`pipelines`](https://github.com/afrendeiro/pipelines):

`projectsroot`: /data/groups/lab_bock/shared/projects/

`htmlroot`: /data/groups/lab_bock/public_html/arendeiro/


# Analysis

### To-do list
#### 1
- Characterize u/mCLL:
    - LOLA, GO, KEGG, motifs
    - differential expression between mCLL and uCLL (log2 fold-change + p-value)
    - isolate the iCLLs: repeat above
- New comparisons (after each, do as above comparison):
    - Untreated/Treated
    - Relapse
    - (mutations: p53, del13q, del11, tri12)
- Characterize CLL (Cancer)
    - Get ATAC-seq/DNase data from T-cells, PBMCs, Myoblasts, Dendritic cells, MSCs, GM12878, K562, Raji cells, RJ2.2.5 cells, IMR90, Ewing line SKMNC,
    - Get CLL-specific sites
    - Correlate sites
    - characterize CLL-specific sites
        - lola, go, etc...
        - Build network with these (genes with r >=0.7)

#### 2
- Disease progression
    - Get same patient
    -Use time-series method to cluster patterns depending on 'direction' with time
- Footprinting
    - close loop with PIQ
    - try using Nextera background
    - Test on DNAse (reproduce Saeed et al 2014)
    - Try uCLL/mCLL/iCLL footprinting
    - Filter motifs from not expressed stuff in CLL + other criteria
- Predict:
    - Minimal number of sites to detect a trait
    - ROC curves
    - in-silico drug prescription (based on genes under regulation by these peaks)
- Correlate with expression
    - Get spanish gene expresion data
    - Correlate with differentially open sites

#### 3
- Nucleosome position differences:
    - Call nucleosomes and dyads using NucleoATAC
    - Look for global nucleosome positioning differences between clusters and compared with controls
    - Look for global nucleosome positioning differences within groups of gene classes (e.g. genes involved in cellular proliferation, cell cycle, B-cell biology and other classes).
- Once B-cell population data is here:
    - get CLL-specific sites, characterize again
    - revisit drug comparisons (see if cells become more like progenitors after treatment)
    - use cell-type-specific signatures to characterize CLL-specific
    + *de novo* enhancer detection
        + compare with other Naive/Memory B-cells from patients/donors

#### future
+ Co-regulated cis-regulatory modules
    + correlate (or use other measurement) peaks
    + if unsuccessful, do the same with tiling regions instead of peaks
+ genetic variants/QTLs:
    + http://www.nature.com/nmeth/journal/v12/n5/full/nmeth.3326.html
    + http://biorxiv.org/content/early/2015/04/30/018788
    + enrichment in GWAS variants
