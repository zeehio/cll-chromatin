---
title: The chromatin accessibility landscape of chronic lymphocytic leukemia
author:
        - André F. Rendeiro^1^^,^^\*^
        - Christian Schmidl^1^^,^^\*^
        - Renata Walewska^2^
        - Zadie Davis^2^
        - David Oscier^2^
        - Jonathan Strefford^3^
        - Christoph Bock^1^^,^^4^^,^^5^^,^^6^
        - ^1^CeMM Research Center for Molecular Medicine of the Austrian Academy of Sciences, Vienna, Austria
        - ^2^Department of Molecular Pathology, Royal Bournemouth Hospital, Bournemouth, UK
        - ^3^Cancer Genomics, Cancer Sciences, University of Southampton, Southampton, United Kingdom
        - ^4^Department of Laboratory Medicine, Medical University of Vienna, Vienna, Austria
        - ^5^Max Planck Institute for Informatics, Saarbrücken, Germany
        - ^6^Correspondence should be addressed to C.B. (<cbock@cemm.oeaw.ac.at>)
        - ^\*^equal contribution
tags: [chronic lymphocytic leukaemia, chromatin, ATAC-seq]
bibliography: /home/afr/Documents/library.bib
---

## Abstract:
<!--
250 words or fewer
-->


## Introduction


## Methods
### Sample acquisition
All cases were diagnosed according to the revised guidelines of the International Workshop Chronic Lymphocytic Leukemia/National Cancer Institute (IWCLL/NCI).

IGHV sequencing, CD38 expression  and screening for CNA's were performed in all cases and ZAP70 expression in the majority. CNA's (del13q, del11q, TP53 loss and trisomy 12) were detected by FISH in most cases or MLPA (MRC Holland SALSA MLPA P037 CLL-1 probemix) in a minority. Most cases also had karyotypic analysis. 

### ATAC-seq
ATAC-seq. Open chromatin mapping was performed using the ATAC-seq method as described [@Buenrostro2013] with minor adaptations for K562 cells. In each experiment, 10 5 cells were washed once in 50 $mu$l PBS, resuspended in 50 $mu$l ATAC-seq lysis buffer (10 mM Tris-HCl, pH 7.4, 10 mM NaCl, 3 mM MgCl 2 and 0.01% IGEPAL CA-630) and centrifuged for 10 min at 4 $^{\circ}$C. Upon centrifugation, the pellet was washed briefly in 50 $mu$l MgCl 2 buffer (10 mM Tris, pH 8.0, and 5 mM MgCl 2 ) before incubating in the transposase reaction mix (12.5 $mu$l 2x TD buffer, 2 $mu$l transposase (Illumina) and 10.5 $mu$l nuclease-free water) for 30 min at 37 $^{\circ}$C. After DNA purification with the MinElute kit, 1 $mu$l of the eluted DNA was used in a qPCR reaction to estimate the optimum number of amplification cycles. Library amplification was followed by SPRI size selection to exclude fragments larger than 1,200 bp. DNA concentration was measured with a Qubit fluorometer (Life Technologies).
Sequencing was performed by the Biomedical Sequencing Facility at CeMM using the Illumina HiSeq 2000/2500 platform (see Supplementary Table {#Table1} for details). Library preparation was performed using custom Nextera primers as described previously [@Buenrostro2013].

### ATAC-seq preprocessing
Reads were trimmed using Skewer [@Jiang2014] and the Nextera adapters as input sequences. Trimmed reads were aligned to the GRCh37/hg19 assembly of the human genome using Bowtie2 [@Langmead2012] with the "--very-sensitive" parameter and the aligner’s default. Duplicate reads were removed using sambamba *markdup*, and only properly paired reads with mapping quality over 30 aligning to the nuclear genome were kept. All downstream analyses were performed on these set of filtered reads.

Genome browser tracks were created with the *genomeCoverageBed* command in BEDTools [@Quinlan2010] and normalized such that each value represents the read count per base pair per thousand reads. Finally, the UCSC Genome Browser's *bedGraphToBigWig* tool was used to produce a bigWig file. The track with percentile signal across the cohort was created by measuring the ATAC-seq read coverage at every basepair with BEDTools *coverage* and normalizing it between samples by divinding each position by the total number of mapped and filtered reads over 10 million reads. Then, the mean as well as the 5, 25, 75 and 95 percentiles of signal across the whole cohort were calculated with Numpy, made into bedgraph and subsequently bigwig format using UCSC's *bedGraphToBigWig*.

Peak-calling was performed with MACS2 [@Zhang2008] using the "--nomodel" and "--extsize 147" flags and arguments, and peaks overlaping the list of blacklisted features of the GRCh37/hg19 assembly as defined by the ENCODE project [@Hoffman2013a] were discarded.
<!--
Reads aligning to the plus strand were offset by +4 bp, and reads aligning to the minus strand were offset by -5 bp as described [@Buenrostro2013]. 
-->

### Analysis of chromatin accessibility variation
The CLL cohort ATAC-seq region set was created by merging all peaks across all samples using BEDTools [@Quinlan2010] *merge* command. To produce Figure [#Figure1]b, we counted the number of unique sites after merging peaks for each sample in an iterative fashion by randomizing the sample order 1000 times and computing 95% confidence intervals across all iterations.

We quantified the acessibility of each region over each sample using Pysam [@citeulike:8274107] by counting the number of reads from the filtered bam file overlaping each region. To normalize read counts across samples, we performed quantile normalization using the *normalize.quantiles* function from the preprocessCore package in R.

For each genomic region we calculated how much support it has over the cohort as the fraction of samples with a called peak in the region over all samples; five measures of ATAC-seq signal variation across the cohort: mean signal, standard deviation, dispersion (variance over mean) and the squared quoeficient of variation (the square of the standard deviation over the mean). Additionaly, we annotated each region with the identity of and distance to the nearest TSS of Ensembl gene models, and the type of genomic region overlaping by using BEDTools *intersect* with features (introns and exons) from the same gene annotation and considering promoters to be 2.5 kb upstream of TSSs and intergenic space the remaining space. To annotate each region with the chromatin state of CD19+ cells, we took the 15 state genome segmentation from the Roadmap Epigenomics Project [@Ernst2015] \(ID E032).

Gene-centric variation presented in Figure [#Figure2]a and b, represents the values of closest region the within 1 kb for promoters, and the mean signall of all distal regions away from the promoter (>1 kb). For Figure [#Figure2]a we ploted the intensity range for these regions, defined as the absolute difference between the 5th and 95th percentiles, whereas in Figure [#Figure2]b the intensity values of the respective regions are used directly to produce the violin plots.

### IGHV mutation status analysis
A Random Forest classifier from the scikit-learn [@Pedregosa2011] implementation (*sklearn.ensemble.RandomForestClassifier*) was trained with each sample's IGHV mutation status as label if known and the matrix of 112168 ATAC-seq CLL cohort region set per sample as input. All CLL samples were used for class prediction using leave-one-out cross validation, and a ROC curve was plotted using the scikit-learn's implementation. Important regions to distinguish samples with different IGHV mutation status were selected by averaging the feature importances of the Random Forest classifier over all iterations and selecting all features with non-zero importance.

Region set enrichment was performed using LOLA [@Sheffield2015] and its core databases: TF binding sites from ENCODE, tissue clustered DNase hypersensitive sites [@Sheffield2013], the Codex database, UCSC feature tables, the Cistrome database and data from the BLUEPRINT project.

Gene ontology and pathway enrichment of the IGHV-specific regions was performed using the R seq2pathway [@Wang2015] and TF motif enrichment analysis was performed with MEME's *ame* tool, using MEME's database of human TF motifs (HOCOMOCOv9).

The sample stratification in Figure [#Figure3]g was performed by correlating the ATAC-seq signal of each sample in the above-defined IGHV status regions of importance in a pairwise fashion and plotting a dendrogram of the correlation with Scipy's hierarchical clustering implementation. With the same correlation values, was performed with R's implementation.

### Gene regulatory network inference from TF footprints
Footprinting was performed with PIQ [@Sherwood2014] using a set of 366 Human transcription factor motifs from the Jaspar database [@Mathelier2014].

As described previously [@Qu2015] we chose to retain TFs with at least 500 high purity binding sites (> 0.7) that overlap any ATAC-seq regions of the CLL cohort.

Our scheme for assigning transcription factor binding sites to genes was as follows: ATAC-seq peaks located in the body of the transcription unit, together with the 2.5 kb regions upstream of the TSS, were assigned to all genes - the case of overlapping transcriptional units (including non-coding transcripts). Intergenic peaks were assigned to the gene whose TSS was closest to the peak. For this we used the GRCh37/hg19 Ensembl gene annotation, and we considered non-protein coding genes in the same manner as protein coding.

<!--
If we need to justify this: 
    - Most ATAC-seq peaks are within gene bodies:
    "transcriptionally active promoters asymmetrically interact with the gene body more than with local upstream sequence, suggesting an activating role for intronic enhancers" [Mifsud, B. et al. Mapping long-range promoter contacts in human cells with high-resolution capture Hi-C. Nat. Genet. 47, 598–606 (2015)]
    - Similar approaches have been used:
    [@Gonzalez2015] McLean, C.Y. et al. GREAT improves functional interpretation of cis-regulatory regions. Nat. Biotechnol. 28, 495–501 (2010).
-->

To infer gene regulatory networks, we calculated an interaction score in a similar way as previously done [@Qu2015]: the interaction score between a transcription factor $t$ and a gene $g$ ($S_{t,g}$) is given by the sum of all TF binding sites (of length $n$) from TF $t$ assigned to gene $g$:
~~~math
S_{t,g} = \sum_{i=0}^{n}2 * (P_{i} - 0.5) * 10 ^{-(\frac{d_{i, g}}{1000000})}
~~~
where $P$ is the PIQ purity score and $d_{i, g}$ is the distance of a particular TF binding site $i$ to gene $g$. This establishes a unidirectional (TF to gene), weigthed (based on the interaction score) relationship, which forms the edges of a graph.

We infered gene regulatory networks for each sample individually, for all samples, and for groups of samples depending on its IGHV mutation stats if known. In the case of the later sample groups, this was performed by using a concatenation of all of the group's bam files as input to PIQ. In order to compare the infered IGHV-unmutated and IGHV-mutated networks, we the logarithm of base two of the difference between the degree of each node between the two networks.

To produce Figure [#Figure4]b and Supplementary Figures [#FigureS10] and [#FigureS11], we considered only TF-TF interactions with score higher than 1.

The entirety of the code used in the is available at [github.com/epigen/cll-patients](https://github.com/epigen/cll-patients).

## Results
### The chromatin accessibility landscape of CLL
Figure [#Figure1]a

Figure [#Figure1]b

Figure [#Figure1]c

Figure [#Figure1]d

### Inter-sample variation in the CLL chromatin accessibility landscape
Figure [#Figure2]a

Figure [#Figure2]b

Figure [#Figure2]c


### Supervised machine learning on clinical features
Figure [#Figure3]a

Figure [#Figure3]b

Figure [#Figure3]c

Figure [#Figure3]d

Figure [#Figure3]e

Figure [#Figure3]f

Figure [#Figure3]g

Figure [#Figure3]h


### Gene regulatory network inference
Figure [#Figure4]a

Figure [#Figure4]b

Figure [#Figure4]c

Figure [#Figure4]d

Figure [#Figure4]e


## Discussion


## Acknowledgements
We thank the Biomedical Sequencing Facility at CeMM for assistance with next-generation sequencing and all members of the Bock lab for their help and advice. This work was performed in the context of the BLUEPRINT project (European Union’s Seventh Framework Programme grant agreement no. 282510) and funded in part by the ERA-NET project CINOCA (FWF grant agreement no. I 1626-B22). C.S. was supported by a Feodor Lynen Fellowship of the Alexander von Humboldt Foundation. C.B. was supported by a New Frontiers Group award of the Austrian Academy of Sciences.

## Authorship Contributions
R.W. and Z.D. followed the patients and isolated lymphocytes from periferal blood, D.O. and J.S. contributed the samples, C.S. performed the experiments; A.F.R. analyzed the data; A.F.R., C.S., D.O., J.S. and C.B. planned the study; C.B. supervised the research; all authors contributed to the writing of the manuscript.

## Disclosure of Conflicts of Interest
The authors declare no conflicts of interest.

## Tables
## Figures

### Figure: {#Figure1}
![](figures/Figure1.pdf){height=100%}
Caption: *The chromatin accessibility landscape of chronic lymphocytic leukemia.*
**a)** Schematic representation of the experimental procedures (top) leading to a map of chromatin accessibility on a CLL cohort.
**b)** Cumulative frequency of unique sites of open chromatin detected across CLL samples in the entire cohort. Lines indicate 95% confidence interval.
**c)** Representative ilustration of ATAC-seq signal in individual CLL samples and across the whole cohort (percentiles indicated in color legend) in the ZAP70 promoter and a regulatory element in its second intron. Note the variance in openness across samples in the promoter contrasting with the constitutively open regulatory element in the intron.
**d)** Enrichment of chromatin states across all detected regions in the cohort over a randomized background of same size. Chromatin state overlap was measured across all human tissues (Ernst and Kellis, 2015 Nature Biotechnology). (THIS IS GOING TO BE CHANGED INTO SOMETHING ELSE)

### Figure: {#Figure2}
![](figures/Figure2.pdf){height=100%}
Caption: *Inter-sample variation in the CLL chromatin accessibility landscape.*
**a)** Variation of accessibilty in promoter and distal regulatory element for all genes with assigned accessible regions. Intensity for genes with more than one distal regulatory element was averaged across elements. Intensity range is defined as the absolute difference between the 95% and 5% percentile of ATAC-seq signal in each element across the cohort.
**b)** Violin plots of chromatin accessibility measurements (openness) at the promoter and enhancers (as defined above) of genes relevant to B cell biology.
**c)** ATAC-seq signal across the CLL cohort in examplary loci of genes in panel **b)**.


### Figure: {#Figure3}
![](figures/Figure3.pdf){height=100%}
Caption: *Clinical data-driven analysis of chromatin accessibility provides insights into CLL diversity and biology.*
**a)** Schematic representation of the machine learning aproach based on clinical data (IGHV mutation status in this case) combined with chromatin accessibility data to achieve sample/patient stratification and relevant feature extraction.
**b)** Receiving operating characteristic (ROC) curve of binary sample classification using a random forest classifier trained on 112168 regions and the IGHV mutation status of the samples. The ROC curve was produced after leave-one-out cross validation and the area under the curve is displayed (0.97).
**c)** Heatmap of samples and regions important to distinguish between the status of IGHV mutation in CLL samples during the classification exercise. Rows and columns were clustered hierarchicaly. Note the two clusters of regions each enriched in IGHV unmutated or mutated samples and therefore termed (uCLL or mCLL, respectively).
**d)** Enrichment of the clusters of regions defined in **c)** in the chromatin states of CD19+ B cells over all CLL open chromatin regions across the cohort.
**e)** Enrichment of cluster 1 and 2 (as defined in **c)**) in predefined functional genomic region sets from other cell types or biological conditions using LOLA.
**f)** Enrichment in cellular pathways from the Reactome and KEGG databases of cluster 1 and 2 (as defined in **c)**)-associated genes.
**g)** Hierarchical clustering of CLL samples based on the sample-wise correlation of chromatin accessibility values of all regions defined in **c)**. Samples are colored with the status of IGHV mutation and with the cluster in which they are placed in the dendrogram.
**h)** Principle component analysis of CLL samples with the chromatin accessibility values of all regions defined in **c)**. Sample colors are as defined in **g)**.

### Figure: {#Figure4}
![](figures/Figure4.pdf){height=100%}
Caption: *Gene regulatory network inference from transcription factor footprints in chromatin accessibility data.*
**a)** Schematic representation of transcription factor (TF) footprinting, TF-gene interaction assignment and gene regulatory network inference from such interactions (TO BE DONE).
**b)** Gene regulatory network infered from TF footprinting of all CLL samples combined. Node size reflects the absolute degree (total number of connections) of the node and color the ratio of out-to in-degree (outgoing or ingoing number of connections, respectively) of each node. Only TFs are visualized.
**c)** Degree of all genes with interactions with TF in the network
**d)** Fold change of the number of connections for each node between networks infered for samples with IGHV unmutated or mutated independently. Nodes are ranked by the fold change between unmutated and mutated IGHV and only transcription factors are shown.
**e)** Genomic loci of two of the most differentially regulated genes from **d)** with chromatin accessibility signal for all samples and in quantiles across cohort.

## Supplemental data
### Figures

#### Figure: {#FigureS1}
Cohort statistics

#### Figure: {#FigureS2}
Quality & saturation

#### Figure: {#FigureS3}
Feature of the CLL cohort region set: support, qv2, etc...

#### Figure: {#FigureS4}
Genomic location and chromatin-state enrichment of the CLL cohort region set.

#### Figure: {#FigureS5}
Feature of IGHV mutation status regions: support, qv2, etc...

#### Figure: {#FigureS6}
Genomic location and chromatin-state enrichment of IGHV mutation status regions.

#### Figure: {#FigureS7}
PIQ stats

#### Figure: {#FigureS8}
CLL network with all interactions

#### Figure: {#FigureS9}
IGHV-unmutated CLL network

#### Figure: {#FigureS10}
IGHV-mutated CLL network

#### Figure: {#FigureS11}
Fold-change of degree for all nodes in the two networks, gene ontology of these.


### Tables
#### Table: {#Table1}
ATAC-seq statistics

#### Table: {#Table2}
Saturated consensus calling of CLL regulatory elements (BED file)

### Website
Track hub with consensus (6 tracks: percentile 5, 25, 50, 75, 95, mean) & all the individual tracks

## References
