---
title: CLL literature review
author:
        - André Rendeiro
subject: 
date: August 02, 2015
company: CeMM Research Centre for Molecular Medicine
keywords: cll, chronic lymphocytic leukaemia, review, literature, genomics
bibliography: Documents/collection.bib
---

## What is this:
This is a short overview of fundamental literature on chronic lymphocytic leukaemia (CLL). While this meant both as a review for me and a *get-up-to-speed-quickly* work for others, this should not encourage you not to read the actual papers.


# Basic review

## B-cell biology
1. In secondary lymphoid organs, mature peripheral B cells (follicular (Fo) B cells) acquire antigen from follicular dendritic cells (FDCs) and in turn presents it to cognate CD4+ TfH cells. During TfH-B cell conjugation, the B cell receives critical signals to undergo class switching and to commence monoclonal expansion.
2. After several rounds of cellular division, the B cells go through somatic hypermutation, generating a diversity of clones in the germinal center.
3. Upon receiving an unidentified stimulus, the maturing B cells migrate from the dark zone to the light zone and start to express their antibody on the cell surface. They are in a state of activated apoptosis and compete for survival signals from follicular dendritic cell and TfH cells. This rescue process is dependent on the antibody affinity to the antigen. T cells are also believed to prevent the generation of autoreactive germinal center B cells.
4. Maturing B cells receive a final differentiation signal to exit the germinal center as plasma cells or memory B cells.

## Why is there CLL?
Mistargeted somatic hypermutation is a likely mechanism in the development of CLL and other B-cell lymphomas.

## Clinical staging
Based primarily on the presence of a low platelet or red cell count.

#### Rai Staging System
Absolute values of lymphocytosis (>15,000/mm3) are considered in all stages.

- Stage 0: no adenopathy, hepatosplenomegaly, anemia, or thrombocytopenia.
- Stage I: lymphadenopathy without hepatosplenomegaly, anemia, or thrombocytopenia.
- Stage II: either hepatomegaly or splenomegaly with or without lymphadenopathy.
- Stage III: anemia (hemoglobin <11 g/dL) with or without lymphadenopathy, hepatomegaly, or splenomegaly.
- Stage IV: characterized by absolute lymphocytosis and thrombocytopenia (<100,000/mm3) with or without lymphadenopathy, hepatomegaly, splenomegaly, or anemia.

#### Binet Classification

- Clinical stage A: no anemia or thrombocytopenia and fewer than three areas of lymphoid involvement (Rai stages 0, I, and II).
- Clinical stage B: no anemia or thrombocytopenia with three or more areas of lymphoid involvement (Rai stages I and II).
- Clinical stage C: anemia and/or thrombocytopenia regardless of the number of areas of lymphoid enlargement (Rai stages III and IV).


## Molecular classification
Compiled from various publications (see )

- CD5+ B-cells
- IGHV unmutated CLL (u-CLL) more like naïve B-cells
	+ Likely comes from pre-germinal center B cells
	+ Enriched in high CD38 and ZAP70 (T cell receptor) expression
	+ Enriched in NOTCH1 mutations
	+ Poor prognosis
- IGHV mutated CLL (m-CLL) more like memory B-cells
	+ Likely comes from post-germinal center, mature B cells
	+ More patients left untreated
- there are obvoulsly intermediate states as well (i-CLL)

> **personal note** although the intermediate group has obviously features of both, it seems more similar to the memory B cell (m-CLL) group


## Classical markers
- IGHV mutational status (<98% homology is considered mutated)
- cytogenetics
- CD38 expression
- ZAP-70 expression (predictive of IGHV mutation status)

These are nowadays considered **not fundamental** for the initiation or continuation of leukemic process, but are still **descriptive** of it.


# Recent Papers

## [Puente et al, Nature, 2011](http://doi.org/doi:10.1038/nature10113)
> Puente, X. S., Pinyol, M., Quesada, V., Conde, L., Ordóñez, G. R., Villamor, N., … Campo, E. (2011). **Whole-genome sequencing identifies recurrent mutations in chronic lymphocytic leukaemia**. Nature, 475(7354), 101–105. doi:10.1038/nature10113

### Design
- WGS of 4 CLLs:
	- 2 IGVH mutated/ 2 unmutated
	- no TP53 mutations
	- Patients were already in advanced stages - no treatments so far. Choice was made to avoid cases with obvious mutations and focus on discovery, so these 4 cases are probably not be the most representative of the disease.
- Validation cohort:
	- 363 patients
	- samples (> 70% CLL cells) + matched controls (< 5% CLL cells)

### Results
- ~1000 mutations per patient (< 1/Mb)
- A>C/T>G mutations are more frequent in CLL-M. This is likely due to polymerase η (Eta) action in immunoglobulin genes. Pol η causes mutations that occur in clusters with bias to certain bases (generally more A>C/T>G).

> Side note: somatic hypermutation is not rare and genes such as BCL6, MYC and PIM1, are recurrently mutated by somatic hypermutation in different lymphomas [Pasqualucci et al, 2001, Nature](http://doi.org/doi:10.1038/35085588).


#### 4 recurrently mutated genes with direct functional implications for the disease:
- NOTCH1
	+ frequency: 12%
	+ Early stop codon mutations
	+ activating mutation
	+ Notch pathway significantly affected (upregulated)
	+ Clinicaly advanced, less survival
- MYD88
	+ frequency: 2.9%
	+ MYD88 mutations are relevant in other lymphoid neoplasias
	+ activating mutation
	+ higher MYD88 phosphorilation, higher NFKb binding
	+ stimulation of IL-1r or TLR receptors in MYD88-mutated cases yielded 5-150 fold secretion of IL-1r antagonist, IL-6, chemokine ligands when compared to unmutated MYD88. same stimulation in cells mutated in MYD88 with inactivating mutation resulted in no response. Production of these cytokines has been implicated in recruitment of T cells and macrophages. Recruitment of these cells might be creating a favourable environment for CLL survival
	+ MYD88-mutated patients had earlier diagnose data and more advanced disease
	+ almost all MYD88-mutated patients were IGHV-mutated
- XPO1
	+ frequency: 2.4%
	+ exportin 1 (nuclear export)
	+ IGHV-unmutated exclusive
	+ co-occurent with NOTCH1 mutations
- KLHL6
	+ frequency: 1.8%
	+ implicated in the formation of germinal centres (B-cell maturation site)
	+ almost exclusive to IGHV-mutated (mutations have pol eta pattern)


#### Other mutations/aberations
- del 13q14
- del 6q14-q22
- del CCND2 (cyclin D2)


## [Braggio et al, Leukemia, 2012](http://doi.org/doi:10.1038/leu.2012.14)
> Braggio, E., Kay, N. E., VanWier, S., Tschumper, R. C., Smoley, S., Eckel-Passow, J. E., … Fonseca, R. (2012). **Longitudinal genome-wide analysis of patients with chronic lymphocytic leukemia reveals complex evolution of clonal architecture at disease progression and at the time of relapse**. Leukemia, 26(7), 1698–1701. doi:10.1038/leu.2012.14

This paper is about structural aberrations during CLL clonal evolution.

### Design
- 22 patients, at least two timepoints
- aCGH
- samples grouped acording to time-points:
	+ TP1: collected >6 months before starting first-line treatment;
	+ TP2: collected at the time of progression, during treatment;
	+ TP3: collected >6 months after initial treatment but before subsequent treatments

### Results
- Only 6 of 22 patients showed copy-number abnormalities differences between time-points.
- These were validated with FISH
- Trisomy 12 and deletion 11q32 were stable over time
- Increment in the number of cases with deletion 13q14.3 and 17p
- Four patients had subclones. In all cases, the dominant clone at TP1 or TP2 was no longer at TP3:
	+ One patient had 3 subclones, all sharing common aberations like trisomy 12, showing that they occurred probably after the emergence of a clone. Clone "A" dominated in TP1, was only ~70% of cells in TP2 and 20% at TP3. At TP3 subclone "B" was 60% of all cells and subclone "C" 20%.
	+ Another patient had two subclones, one of which had 70% dominance at TP2, but only 10% at TP3.

CLL progression can occur in either a linear or a branching manner, with multiple genetic subclones evolving either in succession or in parallel


## [Kulis et al, Nature Genetics, 2012](http://doi.org/doi:10.1038/ng.2443)
> Kulis, M., Heath, S., Bibikova, M., Queirós, A. C., Navarro, A., Clot, G., … Martín-Subero, J. I. (2012). **Epigenomic analysis detects widespread gene-body DNA hypomethylation in chronic lymphocytic leukemia**. Nature Genetics, 44(11), 1236–1242. doi:10.1038/ng.2443

### Design
- WGBS of two patients (one u-CLL, one m-CLL)
- naïve B-cells, non-class-switched memory B-cells (ncsMBC) and class-switched memory B-cells (csMBC) + low coverage WGS from a single donor
- 139 CLLs in 450k arrays + expression arrays
- 450k array on various B cell populations from age matched healty donors
	+ 19 whole B-cell (CD19)
	+ three each of CD5+ NBC, NBC, ncsMBCs and csMBCs with > 95% purity
- exome sequecing of the 139 CLLs

### Results
- NBCs vs csMBCs: 1,076,208 DMRs in 19125 genes
- 96% of DMRs were hypomethylated
- little difference between ncsMBCs and csMBCs
- most DMRs in gene bodies and intergenic space (outside CGIs)
- PCA distinguishes u-CLL and m-CLL
- poor overal correlation between DNAme and expression
- ~3,265 DMRs between u-CLL and m-CLL detected in the arrays, > 2.5 million in WGBS
- ~50% of the DMRs were not distinguishable between uCLL and NBCs or mCLL and scMBCs
- differences between uCLL and mCLL in array data:
	+ 3,243 hypermethylated and 29,743 hypomethylated CpGs in U-CLL compared with NBCs and CD5+ NBCs
	+ 246 hypermethylated and 4,606 hypomethylated CpGs in M-CLL compared with MBCs
- higher differences in the WGBS set:
	+ 1,838,346 DMRs between U-CLL and NBCs (1,779,168 hypomethylated and 59,178 hypermethylated)
	+ 1,254,527 DMRs between M-CLL and csMBCs (982,683 hypomethylated and 271,844 hypermethylated)
- in general hypermethylation was at 5' regions of genes and CGIs, whereas hypomethylation was mainly at gene bodies
- hypomethylated regions are enriched in signaling processes such as B-cell receptor–, NF-κB– and calcium-activated pathways; T-cell costimulation; and cytokine-cytokine receptor interactions
- CpGs with significant correlation with expression (4%) were enriched in gene bodies
- in a fraction, correlation between gene expression and DNA methylation was simultaneously negative in 5′ regions and positive in gene bodies
- CpGs hypomethylated in U-CLL or M-CLL were enriched for enhancer regions in B cells and IRF and OCT1 binding sites
- enhancers were also enriched in the uCLL/NBC comparison
- gene body CpGs with negative correlation with expression were enriched in enhancers, whereas positively correlated ones had only slight enrichment	
- tried to loook at the relationship between DNAme and splicing aberations - only two genes showed DNA methylation correlated with exon skipping, and lack of methylation correlated with exon inclusion
- Eight recurrently mutated genes in CLL showed aberrant hypomethylation (ASXL1, BRAF, CDH23, EGR2, FAM117A, PHC2, POT1 and SF3B1) - not more than expected by chance - exome data only used here
- Further clustering of samples based on NBC- and MBC-CpGs showed three clusters, the third one with a majority of mCLLs enriched in fewer IGHV mutations. These CLLs could come from "antigen- experienced, germinal center–independent B cell with low levels of somatic hypermutation".
- The three groups showed significant clinical features between them


## [Queirós et al, Leukemia, 2014](http://doi.org/doi:10.1038/leu.2014.252)
> Queirós, a C., Villamor, N., Clot, G., Martinez-Trillos, A., Kulis, M., Navarro, A., … Martín-Subero, J. I. (2014). **A B-cell epigenetic signature defines three biological subgroups of chronic lymphocytic leukemia with clinical impact**. Leukemia, (August 2014), 598–605. doi:10.1038/leu.2014.252

### Design
- 

### Results
- 


## [Baliakas et al, Leukemia, 2015](http://doi.org/doi:10.1038/leu.2014.196)
Baliakas, P., Hadzidimitriou, A., Sutton, L., Rossi, D., Minga, E., Villamor, N., … Rosenquist, R. (2014). **Recurrent mutations refine prognosis in chronic lymphocytic leukemia**. Leukemia, (April), 1–8. doi:10.1038/leu.2014.196

### Design
- 

### Results
- 


## [Puente et al, Nature, 2015](http://doi.org/doi:10.1038/nature14666)
> Puente, X. S., Beà, S., Valdés-Mas, R., Villamor, N., Gutiérrez-Abril, J., Martín-Subero, J. I., … Campo, E. (2015). **Non-coding recurrent mutations in chronic lymphocytic leukaemia**. Nature. doi:10.1038/nature14666

### Design
- 

### Results
- 



------------
Made with [Scholarly Markdown](http://scholarlymarkdown.com/), and [Scholdoc](http://scholdoc.scholarlymarkdown.com/)
