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

For all my projects:

`projectsroot`: /data/groups/lab_bock/shared/projects/

`htmlroot`: /data/groups/lab_bock/public_html/arendeiro/


# Todo

## Data production
+ Open-chromatin
    + ATAC-seq
+ Histone mods (possible):
    + H3K27AC
    + H3K4ME1
    + H3K4ME3?
+ RNA?

## Ideas
+ get differentialy open regions between groups:
    + correlate each peak a variable from patient data, plot p-values for each variable along heatmap
    + test enrichment of each groups of regions with:
        + LOLA
        + nearest gene: KEGG pathways, GO, OMIM
+ GRNs
    + Footprint-based
        + Footprint each patient or groups of patients
        + Establish TF -> gene relation
        + Compare networks
+ Co-regulated cis-regulatory modules
    + correlate (or use other measurement) peaks
    + if unsuccessful, do the same with tiling regions instead of peaks
+ Correlate chromatin openness along genes with gene expression
    + use spanish data
+ genetic variants/QTLs:
    + http://www.nature.com/nmeth/journal/v12/n5/full/nmeth.3326.html
    + http://biorxiv.org/content/early/2015/04/30/018788
    + enrichment in GWAS variants
+ nucleosome position differences:
    + Call nucleosomes and dyads using NucleoATAC
    + Look for global nucleosome positioning differences between clusters and compared with controls
    + Look for global nucleosome positioning differences within groups of gene classes (e.g. genes involved in cellular proliferation, cell cycle, B-cell biology and other classes).
+ *de novo* enhancer detection
    + compare with other Naive/Memory B-cells from patients/donors
+ *Time-lapse* of same-patient data:
    + Treat each sample as independent, test differences, measure enrichments
    + Clonal evolution within patients (if genotyped)
+ Finding new biomarkers
    + Look at earliest points (if possible before diagnosis) for consistent activation
