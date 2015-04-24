cll-patients
============

**Nice paper title**

-----------


# Repository organization

Wetlab and sequenced sample annotations are in `/metadata`.

Code used in the analysis will be in `/src`

### Data
See annotation sheet here: [CSV](metadata/cll-patients.sample_annotation.csv)

See tracks here: [UCSC](http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&hgt.customText=http://www.biomedical-sequencing.at/bocklab/arendeiro/cll-patients/bigWig/trackHub_hg19.curated.txt)

### Sample names
Sample names are a concatenation of the annotation fields: `cellLine, numberCells, technique, ip, patient, treatment, biologicalReplicate, technicalReplicate, genome`. 

Samples where all of the above are the same except for `technicalReplicate` (same sample sequenced twice independently), are merged in a new sample labeled with `technicalReplicate=0`, representing one biological replicate which was sequenced more than once. 

Further, Samples where all of the above are the same except for `technicalReplicate` AND `biologicalReplicate` , are merged in a new sample labeled with `technicalReplicate=0`, `biologicalReplicate=0`, representing a concatenation of all biological replicates of this type.

For downstream analysis I've been using mostly the later case: `technicalReplicate=0`, `biologicalReplicate=0`.

This creates a [***sheet with new samples*** containing these merged ones](https://github.com/ComputationalEpigenetics/cll-patients/blob/master/metadata/cll-patients.replicates.annotation_sheet.csv), and paired control samples, which I use for downstream.

### Project structure
As defined in [`chipseq-pipelines`](https://github.com/afrendeiro/chipseq-pipelines):

For all my projects:

`projectsroot`=/fhgfs/groups/lab_bock/shared/projects/
`htmlroot`=/fhgfs/groups/lab_bock/public_html/arendeiro/

```
projectsroot
|__ cll-patients
    |__ runs
    |__ data
    |   |__ fastq
    |   |__ fastqc
    |   |__ raw
    |   |__ mapped
    |   |__ coverage
    |   |__ peaks
    |   |__ motifs
    |__ results
         |__ plots

htmlroot
|__ cll-patients
    |__ bigwig
```
JSON description [here](https://github.com/ComputationalEpigenetics/cll-patients/blob/master/metadata/projectPaths.json).

I will document the content of these folders in due time.


# Todo

## Patient description
+ B-cell CLL or other chronic leukaemias?
+ Patient history
    + ZAP-70 (T-cell r) tested?
    + Blood cells counts? Yes
    + Remissions annotated?
    + Medication history?
+ Time-linked
+ 94 samples
+ ~50 patients (+ ~50 coming)
+ 2 vials with ~20M cells
+ are we getting also samples from normal donors?

## Data production
+ Open-chromatin
    + ATAC-seq
+ Histone mods:
    + H3K27AC
    + H3K4ME3?
    + ChIP-seq or CM? (probably dependent on CM paper coming out)
    + IgG or Input?

Could we have ATAC-seq done on several hematopoietic lineages as a way of seeing samples in light of development and maybe classifying them according to where in differentiation things went wrong?

Could we have gene expression for (at least some of the samples) to show that chromatin opening does change things, and that they correlate?

#### Timeline
2-3 months for ATAC-seq data production (end-to-end)

## Analysis (ideas)
+ Preprocessing
+ QC metrics
+ Stratification:
    + Clustering
    + Pareto
    + Global *vs* only changing regions
    + Compare with physician stratification/remission/disease outcome (if any available)
    + GO/OMIM/LOLA enrichment
    + Interaction networks within clusters
+ Machine learning on stratified patients
    1. Test several approaches
    2. Get metrics (ROC, AUCs, ...)

#### Decision time:
+ Decide at this point if we proceed with histone mods profiling
    + Use imputation with rest of data
+ Decide on the *spin* of the paper:
    + focus on GRNs, biomarkers or deeper characterization

##### Aditional ideas using only ATAC-seq data:
+ Call nucleosomes and dyads using NucleoATAC
+ Look for global nucleosome positioning differences between clusters and compared with controls
+ Look for global nucleosome positioning differences within groups of gene classes (e.g. genes involved in cellular proliferation, cell cycle, B-cell biology and other classes).

#### Several ideas for downstream:
+ GRNs
    + Motif-search
    + Enhancer-TF pairing
    + With evolution build GRN
+ Show chromatin opening in enhancers though normal differentiation and the different states of chromatin opening during disease.
+ De novo enhancers
    + Comparison with Roadmap + Encode
+ Super-enhancers
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

