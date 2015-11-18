---
title: Large-scale chromatin accessibility profiling uncovers heterogeneity of molecular phenotypes and gene regulatory networks of CLL
author:
        - André F. Rendeiro^1,\*^,
        - Christian Schmidl^1,\*^,
        - Renata Walewska^2^,
        - Zadie Davis^2^,
        - Matthias Farlik^1^,
        - Jonathan Strefford^3^,
        - David Oscier^2^,
        - Christoph Bock^1,4,5^
        - ^1^CeMM Research Center for Molecular Medicine of the Austrian Academy of Sciences, Vienna, Austria
        - ^2^Department of Molecular Pathology, Royal Bournemouth Hospital, Bournemouth, UK
        - ^3^Cancer Genomics, Cancer Sciences, University of Southampton, Southampton, UK
        - ^4^Department of Laboratory Medicine, Medical University of Vienna, Vienna, Austria
        - ^5^Max Planck Institute for Informatics, Saarbrücken, Germany
        - ^\*^equal contribution
tags: [chronic lymphocytic leukaemia, chromatin, ATAC-seq, gene regulatory networks]
bibliography: /home/afr/Documents/library.bib
---

#### Key points:
- First chromatin acessibility profiling in a cancer cohort identifies epigenetic heterogeneity of relevance to the disease
- Gene regulatory network inference unveils alternate regulatory interactions between disease subtypes

## Abstract:
<!--
250 words or fewer
-->
Recent advances in fast, low-input and cost-effective epigenome profiling technologies allow more comprehensive enquiry into the diversity of chromatin regulation in large cohorts of patients. Using ATAC-seq on a large cohort of Chronic Lymphocytic Leukemia (CLL) cases, we've characterized the chromatin landscape of CLL, finding substantial epigenetic heterogeneity.
With a supervised machine learning approach guided on clinical annotations, we've uncovered regulatory regions associated with the IGHV mutation status – the largest molecular determinant of CLL classification. Usage of these regions reflects the likely cell-of-origin of each group during the differentiation of normal B cells. Enrichment analysis of these regions reveals differentially used pathways that may be driving disease progression and contribute to the distinction between aggressive and indolent cases. This classification approach allowed us to find additional sample structure, which we used to stratify patients, finding groups with distinct molecular phenotypes, which highlights the heterogeneity of CLL at the chromatin level.
To gain deeper insights into the molecular wiring of CLL, we've inferred gene regulatory networks (GRN) from transcription factor footprints, and describe the first GRN of this kind for a cancer. Comparing the inferred interaction networks from the previously identified molecular subgroups, we uncover novel differentially used network modules which include previously described differential interactions.
This work describes the heterogeneity of the CLL regulatory landscape, finding groups of molecular phenotypes with distinct usage of regulatory regions and signaling pathways, as well as specific regulatory interactions of potential therapeutical relevance.

## Introduction
"Chronic lymphocytic leukemia (CLL) is a B-cell neoplasm with two major clinical and molecular subtypes related to the mutational status of IGHV, which encodes the immunoglobulin heavy chain variable region [@Zenz2010]. CLLs bearing a high level of IGHV somatic hypermutation (mCLLs) have a favorable clinical prognosis, whereas CLLs with a low or absent IGHV mutational load (uCLLs) have a worse outcome [@Zenz2010]. Although the precise identification of the normal cell counterpart of these CLL subtypes is still controversial [@Seifert2012], IGHV mutational status is thought to reflect the differentiation stage of the cells from which these CLL subtypes derive."

After extensive enquiry into the genetic basis of CLL in the form of somatic mutations [@Puente2011; @Quesada2011; @Puente2015; @Landau2015], the number of recurently mutated genes in CLL is relatively low compared with other neoplasias. Furthermore, genetic factors contributing to the origin and progression are even of smaller number, being recently estimated at 44 [@Landau2015].

Epigenetic alterations are ubiquitous in cancer and continue to accumulate during the progress of disease, while in the meanwhile the mechanisms that promote changes in the tumor epigenome are not yet clear.

[epigenetics gives window into past, present and future]

Previous studies of epigenetics in the context of CLL, perhaps for technical reasons, have focused on DNA methylation [@Kulis2012; @Pei2012; @Oakes2014; @Queiros2014] (reviewed in @Cahill2013) in which methylation outside the traditional CpG island promoters context seems crucial in regulating gene expression in cancer but also during B cell maturation.
While, DNA methylation has been shown capable of uncovering CLL groups with clinical significance [@Kulis2012; @Queiros2014], it is still a costly assay if completeness of coverage genome-wide is paramount and therefore not applyiable to large-scale clinical settings.

Recent technological breakthroughs enable chromatin profiling in large cohorts [@Buenrostro2013] or in rare cell populations by requiring less input material [@Lara-Astiaso2014b; @Schmidl2015] or with faster sample-to-result times [@Schmidl2015] and bring great promise for large scale chromatin profiling in clinical settings.

Since CLL shows substantial heterogeneity of molecular and clinical phenotypes (reviewed in by Gui and Wu [@Gui2015]), we decided to explore the heterogeneity of its chromatin regulation. We performed the assay of transposase-acessible chromatin followed by next-generation sequencing (ATAC-seq)[@Buenrostro2013] in primary lymphocytes of 86 CLL cases, establishing a comprehensive map of epigenetic regulation which reveals the epigenetic variation in regulatory elements across the cohort. With a supervised machine learning approach that incorporates clinical annotations, we uncovered regulatory regions associated with the IGHV mutation status, treatment response and relapse of CLL cases. Chromatin patterns at the most predictive regions reflects the differentiation status of normal B cells, the likely cell-of-origin of each group and differentially active pathways that may distinguish between aggressive and indolent cases. By detecting transcription factor (TF) footprints in the ATAC-seq data, we infer gene regulatory networks (GRNs) of TF-gene interations to gain deeper insights into the molecular wiring of CLL and describe the first GRNs of its kind for any cancer. GRN comparison between the IGHV molecular subgroups uncovers differentially used TFs with known differential regulation and novel interactions that may shed light on the molecular basis of the groups of clinical phenotypes.

## Materials and methods
### Sample acquisition
All cases were diagnosed according to the revised guidelines of the International Workshop Chronic Lymphocytic Leukemia/National Cancer Institute (IWCLL/NCI).

IGHV sequencing, CD38 expression  and screening for CNA's were performed in all cases and ZAP70 expression in the majority. CNA's (del13q, del11q, TP53 loss and trisomy 12) were detected by FISH in most cases or MLPA (MRC Holland SALSA MLPA P037 CLL-1 probemix) in a minority. Most cases also had karyotypic analysis.

### ATAC-seq
Open chromatin mapping was performed using the ATAC-seq method as described [@Buenrostro2013] with minor adaptations for K562 cells. In each experiment, 10 5 cells were washed once in 50 $mu$l PBS, resuspended in 50 $mu$l ATAC-seq lysis buffer (10 mM Tris-HCl, pH 7.4, 10 mM NaCl, 3 mM MgCl 2 and 0.01% IGEPAL CA-630) and centrifuged for 10 min at 4 $^{\circ}$C. Upon centrifugation, the pellet was washed briefly in 50 $mu$l MgCl 2 buffer (10 mM Tris, pH 8.0, and 5 mM MgCl 2 ) before incubating in the transposase reaction mix (12.5 $mu$l 2x TD buffer, 2 $mu$l transposase (Illumina) and 10.5 $mu$l nuclease-free water) for 30 min at 37 $^{\circ}$C. After DNA purification with the MinElute kit, 1 $mu$l of the eluted DNA was used in a qPCR reaction to estimate the optimum number of amplification cycles. Library amplification was followed by SPRI size selection to exclude fragments larger than 1,200 bp. DNA concentration was measured with a Qubit fluorometer (Life Technologies).
Sequencing was performed by the Biomedical Sequencing Facility at CeMM using the Illumina HiSeq 2000/2500 platform (see Supplementary Table {#Table1} for details). Library preparation was performed using custom Nextera primers as described previously [@Buenrostro2013].

### ATAC-seq data preprocessing
Reads were trimmed using Skewer [@Jiang2014] and the Nextera adapters as input sequences. Trimmed reads were aligned to the GRCh37/hg19 assembly of the human genome using Bowtie2 [@Langmead2012] with the "--very-sensitive" parameter and the aligner’s default. Duplicate reads were removed using sambamba *markdup*, and only properly paired reads with mapping quality over 30 aligning to the nuclear genome were kept. All downstream analyses were performed on these set of filtered reads.

Genome browser tracks were created with the *genomeCoverageBed* command in BEDTools [@Quinlan2010] and normalized such that each value represents the read count per base pair per thousand filtered reads. Finally, the UCSC Genome Browser's *bedGraphToBigWig* tool was used to produce a bigWig file. The tracks with percentile signal across the cohort was created by measuring the ATAC-seq read coverage at every basepair with BEDTools *coverage* and normalizing it between samples by divinding each position by the total number of filtered reads over 10 million. Then, the mean as well as the 5, 25, 75 and 95 percentiles of signal across the whole cohort were calculated with Numpy, made into bedgraph and subsequently bigwig format using UCSC's *bedGraphToBigWig*.

Peak-calling was performed with MACS2 [@Zhang2008] using the "--nomodel" and "--extsize 147" flags and arguments, and peaks overlaping the list of blacklisted features of the GRCh37/hg19 assembly as defined by the ENCODE project [@Hoffman2013a] were discarded.
<!--
Reads aligning to the plus strand were offset by +4 bp, and reads aligning to the minus strand were offset by -5 bp as described [@Buenrostro2013].
-->

### Analysis of chromatin accessibility variation
The CLL cohort ATAC-seq region set was created by merging all peaks across all samples using BEDTools [@Quinlan2010] *merge* command. To produce Figure [#Figure1]b, we counted the number of unique sites after merging peaks for each sample in an iterative fashion by randomizing the sample order 1000 times and computing 95% confidence intervals across all iterations.

We quantified the acessibility of each region over each sample using Pysam [@citeulike:8274107] by counting the number of reads from the filtered bam file overlaping each region. To normalize read counts across samples, we performed quantile normalization using the *normalize.quantiles* function from the preprocessCore package in *R*.

For each genomic region we calculated how much support it has over the cohort as the fraction of samples with a called peak in the region over all samples; five measures of ATAC-seq signal variation across the cohort: mean signal, standard deviation, dispersion (variance over mean) and the squared quoeficient of variation (the square of the standard deviation over the mean). Additionaly, we annotated each region with the identity of and distance to the nearest TSS of Ensembl gene models, and the type of genomic region overlaping by using BEDTools *intersect* with features (introns and exons) from the same gene annotation and considering promoters to be 2.5 kb upstream of TSSs and intergenic space the remaining space. To annotate each region with the chromatin state of CD19+ cells, we took the 15 state genome segmentation from the Roadmap Epigenomics Project [@Ernst2015] \(ID E032).

Gene-centric variation presented in Figure [#Figure2]a and b, represents the values of closest region the within 1 kb for promoters, and the mean signall of all distal regions away from the promoter (>1 kb). For Figure [#Figure2]a we ploted the intensity range for these regions, defined as the absolute difference between the 5th and 95th percentiles, whereas in Figure [#Figure2]b the intensity values of the respective regions are used directly to produce the violin plots.

### Sample classification on IGHV mutation status and region characterization
A Random Forest classifier from the scikit-learn [@Pedregosa2011] implementation (*sklearn.ensemble.RandomForestClassifier*) was trained with each sample's IGHV mutation status as label if known and the matrix of 112168 ATAC-seq CLL cohort region set per sample as input. All CLL samples were used for class prediction using leave-one-out cross validation, and a ROC curve was plotted using the scikit-learn's implementation. Important regions to distinguish samples with different IGHV mutation status were selected by averaging the feature importances of the Random Forest classifier over all iterations and selecting all features with non-zero importance.

Region set enrichment was performed using LOLA [@Sheffield2015] and its core databases: TF binding sites from ENCODE, tissue clustered DNase hypersensitive sites [@Sheffield2013], the Codex database, UCSC feature tables, the Cistrome database and data from the BLUEPRINT project.

Gene ontology and pathway enrichment of the IGHV-specific regions was performed using the R seq2pathway [@Wang2015] and TF motif enrichment analysis was performed with MEME's *ame* tool, using MEME's database of human TF motifs (HOCOMOCOv9).

The sample stratification in Figure [#Figure3]g was performed by correlating the ATAC-seq signal of each sample in the above-defined IGHV status regions of importance in a pairwise fashion and plotting a dendrogram of the correlation with Scipy's hierarchical clustering implementation. With the same correlation values, was performed with R's implementation.

### Gene regulatory network inference from TF footprints
Footprinting was performed with PIQ [@Sherwood2014] using a set of 366 Human transcription factor motifs from the Jaspar database [@Mathelier2014].

As described previously [@Qu2015] we chose to retain TFs with at least 500 high purity binding sites (> 0.7) that overlap any ATAC-seq regions of the CLL cohort.

Our scheme for assigning transcription factor binding sites to genes was as follows: TFBSs located in the body of a transcription unit or the 2.5 kb region upstream its transcription start site (TSS), were assigned to the respective overlaping gene - in the case of overlapping transcriptional units we assigned the TFBS to all overlaping genes. Intergenic TFBSs were assigned to the gene whose TSS was closest to the peak. For this association, we used the GRCh37/hg19 Ensembl gene annotation and we considered non-protein coding genes in the same manner as the protein coding ones.

<!--
If we need to justify this:
    - Most ATAC-seq peaks are within gene bodies:
    "transcriptionally active promoters asymmetrically interact with the gene body more than with local upstream sequence, suggesting an activating role for intronic enhancers" [Mifsud, B. et al. Mapping long-range promoter contacts in human cells with high-resolution capture Hi-C. Nat. Genet. 47, 598–606 (2015)]
    - Similar approaches have been used:
    [@Gonzalez2015] McLean, C.Y. et al. GREAT improves functional interpretation of cis-regulatory regions. Nat. Biotechnol. 28, 495–501 (2010).
-->

To infer gene regulatory networks, we calculated an interaction score in a similar way as previously done [@Qu2015]: the interaction score between a transcription factor $t$ and a gene $g$ ($S_{t,g}$) is given by the sum of all TF binding sites (of length $n$) from TF $t$ assigned to gene $g$:
~~~math #interaction-score
S_{t,g} = \sum_{i=0}^{n}2 * (P_{i} - 0.5) * 10 ^{-(d_{i, g})}
~~~
where $P$ is the PIQ purity score and $d_{i, g}$ is the distance of a particular TF binding site $i$ to gene $g$ in megabases. This establishes a unidirectional (TF to gene), weigthed (based on the interaction score) relationship, which forms the edges of a graph.

We infered gene regulatory networks for each sample individually, for all samples, and for groups of samples depending on its IGHV mutation stats if known. In the case of the later sample groups, this was performed by using a concatenation of all of the group's bam files as input to PIQ. In order to compare the infered IGHV-unmutated and IGHV-mutated networks, we the logarithm of base two of the difference between the degree of each node between the two networks.

To produce Figure [#Figure4]b and Supplementary Figures [#FigureS10] and [#FigureS11], we considered only TF-TF interactions with score higher than 1.

The entirety of the source code used in the analysis is available at [github.com/epigen/cll-patients](https://github.com/epigen/cll-patients).

## Results
### The chromatin accessibility landscape of CLL
To establish the chromatin acessibility landscape of CLL and investigate its variation across disease subtypes and particular cases, we chose to use the assay of transposase-acessible chromatin followed by next-generation sequencing (ATAC-seq)[@Buenrostro2013] to profile primary lymphocytes of CLL patients and generate patient-specific chromatin acessibility profiles of CLL cases (Figure [#Figure1]a), primarily due to the relative experimental ease but importantly, its high information content. 

We sequenced an average of 11.8 million read pairs per sample - a total of 2.1 billion read across 86 samples (Table [#TableS1] and Supplementary Figure [#FigureS1]), producing sequencing libraries of high quality (Supplementary Figure [#FigureS2]). These were sequenced with enough depth to discover most accessibility regions in each sample (Supplementary Figure [#FigureS3]). With 86 CLL samples across the whole clinical range of CLL (see Supplementary Figure [#FigureS1]), we manage to cumulatively capture all unique sites of open chromatin to saturation (Figure [#Figure1]b).

To establish a defined region set that characterizes all regions of chromatin acessibility in CLL, we created a set of unique regions of chromatin acessibility based on the individual ATAC-seq peaks of all samples, consisting of 112168 regions which are mainly positioned in gene promoters and introns as well as intergenic space (Figure [#Figure1]c, and Supplementary Figure [#FigureS1]).

These regions are located mostly in intergenic space, gene introns and promoters (Figure [#Figure1]d) as has been described due to the nature of the ATAC-seq assay [@Buenrostro2013]. We decided to investigate the chromatin state of this set of regions regarding a whole CD19+ B cell population. Consistent with its genomic location, most peaks are located in known enhancers and transcription start sites (TSSs) of active genes in B cells. We do find however a considerable amount of regions overlaping quiescent chromatin (without any enriched histone modification in particular) and transcriptional domains of B cells.

### Sample heterogeneity in the CLL chromatin accessibility landscape
Since phenotypes are the product of gene expression and this is dependent on gene regulation  - to which chromatin accessibility contributes to [@Natarajan2012] [@Marstrand2014], we sought to characterize the  of variation in the chromatin accessibility of regulatory elements associated with genes. To this end, for each CLL region, we quantified the dynamic range of chromatin accessibility by measuring the difference between the extremes of the distribution (5th to 95th percentiles).

The range of chromatin accessibility values in promoter elements tends to be more narrow than that of distal regulatory elements of genes, although a considerable number of these also show high dynamic range. Most distal regulatory elements of genes tend to exhibit similar variation across the CLL cohort (Figure [#Figure2]a).

To further describe the chromatin heterogeneity of CLL, we investigated the variability of regulatory elements of genes relevant to B cell biology and CLL pathogenesis across the cohort (Figure [#Figure2]b and Supplementary Figure {#FigureS7}A).

Overall, B cell surface marker genes, genes important for BCR signalling (CD79A/B), genes important for neoplastic proliferation (MYCN, ) or genes known to be recursively mutated in CLL (NOTCH1, SF3BP1, XPO1, CDKN1B) [@Puente2015] are mostly in the lower part of distribution of chromatin accessibility dynamic range (Supplementary Figure [#FigureS7]B).

Promoters of genes known to initiate CLL pathogenesis and other important BCR signaling genes were found preferentially in the lower quantile of chromatin range among all elements (being therefore less variable). Such examples are the MYC and NRAS oncogenes (although with high chromatin acessibility) or the LYN kinase, the initiating step in BCR signaling after receptor stimulation and expressed in both naive and memory B cells (Figure [#Figure2]C).

Nonetheless, some genes important to drive B-cell development into more terminally differenciated states (*e.g.* memory B cells) such as the transcription factor BCL6 have a larger dynamic range in their promoters (Figure [#Figure2]C).

Chromatin acessibility in distal regulatory elements is in general more variable across the cohort and following the global trend, the elements of elements of genes relevant to B cell biology and CLL pathogenesis are also more distributed over the range of dynamic ranges across regulatory elements (Supplementary Figure [#FigureS7]B), although like their promoters, these occupy mainly the first percentile of the distribution as well.


### Alternative chromatin use of CLL IGHV subtypes uncovered by supervised machine learning
We next sought to explore differential chromatin acessibility between CLL cases in the hope that these are and reveal principles of sample structure. We pursued an approach using supervised machine learning classification guided by clinical annotation of CLL cases and using the set of all 112168 CLL regions. This approach has the advantage that focusing on differences brought by considering clinical features are more likely revealing of clinically important features.

Since mutation of IGHV locus is the most determinant factor of CLL molecular variation and clinical outcome [@Zenz2010], we decided to focus on discovering chromatin accessibility features associated with this feature. We trained a Random Forest classifier on the IGHV mutation status of the samples and using leave-on-out cross validation, predicted every sample's IGHV mutation status, including if this had not been known or measured before. The classifier's prediction accuracy was extremely high (AUC = 0.97 - Figure [#Figure3]B) even though we did not use any particular subset of regions, but the whole set of of CLL accessible chromatin regions across the whole cohort. Another advantage of this approach was that this learning process allowed us to easily extract features which are important to distinguish these two molecular subgroups without any hard cutoff. We found 1504 regions contributing to this distintion (Figure [#Figure3]b), which were clearly grouped in two upon hierarchical clustering of the chromatin accessibility signal: one cluster of regions displaying stronger chromatin accessibility in the majority of IHGV-mutated cases (cluster 1) and conversely, another with stronger in most IHGV-unmutated cases. This duality is expected since IGHV mutation status was the trait guiding the learning process leading to the discovery of these regions. Importantly, there were no strong deviations from expectation regarding the chromosomal location of these regions (Supplementary Figure S3), revealing that preferential structural chromosomal aberations between samples of different IGHV status could be providing major differences in the location of these differentially open regions.

To characterise these regions, we first investigated what is their broad genomic context when compared with the whole set of CLL chromatin accessibility regions. The most contrasting genomic feature between the two were gene promoters, where the uCLL regions showed to be enriched whereas the mCLL was depleted, but in its turn more enriched in introns than the uCLL regions (Figure [#Figure3]D). A similar pattern was observed regarding the chromatin states of CD19+ B cells: the uCLL regions was enriched in regions near active transcription start sites (TSS), while the mCLL regions were more enriched in distal and intronic enhancer regions of CD19+ cells. While there was agreement in both clusters of regions in their depletion of repressive chromatin domains, perhaps the most striking finding is the divergence of uCLL and mCLL regions on bivalent TSSs and enhancers of CD19+ cells. These regions are likely enhancers of progenitors, that are poised for use by B cells during their maturation into plasma or memory B cells. Nonetheless, one cannot exclude that since a whole CD19+ population of cells is likely a mixture of B cells in various developmental stages, chromatin domains annotated as bivalent in this case might originate from the mixture of in which a particular elemnent can be used differentially and acquire both H3K4me1 and H3K27me3 histone marks. In the case of these being truly poised enhancers of a pre/pro B cells, it would be interesting to assess when was accessbility in these regions acquired by CLL cases and what is their contribution to the development of the neoplasia, particularly in the group of samples most enriched in these, the uCLL cases.

We therefore turned to a more functional characterization of these clusters of chromatin regions by measuring the overlap with collections of publicly available region sets of functional relevance including genomic and regulatory elements. Regions with higher chromatin accessibility in the uCLL were shown to be enriched in CLL-specific enhancer regions and regions of CD38-negative naive B cells with the enhancer mark H3K4me1, reflecting the likely naive B cell origin of these CLL cells. These regions were also enriched in transcriptional domains (determined by ChIP-seq of H3K36me3) of naive B cells and B cell-derived cell lines such as the pre CSR BL-2 line hinting that uCLL might be using new, possibly exclusive enhancers from the transcription sites of B cells. The enrichment of either DNase hypersensitivity regions or transcriptional domains of other hematopoietic cell types likely reflects the more undifferentiated status of uCLL in the B cell lineage compared to mCLL. Regions more accessible in mCLL, on the other hand, are domainated by enhancer regions of lymphocyte-derived cell lines (SU-DHL-5, JVM-2, GM12878, KARPAS-422) and TF binding sites (Figure [#Figure3]E). One such example is BATF - known to regulate class-switch recombination (CSR) [@Ise2011], but binding sites of BCL6 in germinal centers (GC), as well as EBF1, BCL3 and 11A - both  associated with B-cell malignancies - are also significantly enriched in the mCLL regions.

Being cellular signaling pathways a very important and their activation often a determinant of cell proliferation and clinical outcome in cancer, we performed pathway enrichment of genes associated with IGHV-specific chromatin acessibility regions (Figure [#Figure3]F). We found striking differences between the two region clusters in terms of their enriched molecular pathways: the mCLL region cluster was enriched in pathways relevant for lymphocytes, but perhaps not particularly prominent in CLL pathogenesis (CTLA4 inhibitory signaling, FCERI pathway and FC) with the exception of BCR signaling; in the uCLL-specific regions we find NOTCH pathway signaling proeminently enriched, as well as FGF signaling, together with several pathways that promote cytoskeletal activity and cell movement.

The sample structure created by hierarchical clustering of all samples in all of the 1504 IGHV-specific regions (Figure [#Figure3]C) is mostly as expected, with cases with the same IGHV mutation status being found together. Nonetheless, a gradient of signal across samples (in inverse directions dependent on the cluster of regions) can be detected, suggesting that the sample structure is not hermetic but heterogenious, and that further molecular groups of the disease may exist. To further structure CLL cases in the previously discovered IGHV-specific regions (Figure [#Figure3]C), we performed hierarchical clustering on the correlation of chromatin acessibility signal (Figure [#Figure3]G) in these regions, as well as principal component analysis (PCA) on the chromatin acessibility signal (Figure [#Figure3]H). The sample structure reveals three major groups, corresponding to uCLL and mCLL and a third intermediate group (iCLL) stemming mostly from the mCLL group, but having a mixture of IGHV mutation status. Interestingly, another intermediate group of samples stemming from the uCLL group is observed, with most samples bearing a mutated IGHV locus, despite its very small size (n = 3).


### Gene regulatory network inference from TF footprints
Due to the tagmentation-based nature of the ATAC-seq assay, chromatin locations bound by proteins can hinder the action of the transposase enzyme, reducing the acessibility signal in a locally. This protection therefore leaves an impression of the binding of TFs on the chromatin called footprints which can be detected through specialized analysis of the accessibility signal in specific positions containing TF motifs [@Buenrostro2013]. We applyied a TF footprint detection algorithm to a set of 366 human transcription factors and interogated the binding sites of these TFs that overlap the set of all CLL accessible sites for the ocupation of the respective TF (Figure [#Figure4]a). This allowed us to establish a relationship between the binding TF and genes in proximity of its binding sites through an interaction score (more in the Methods section), which we use to reconstruct gene regulatory networks (GRN).

We first started with an attempt to reconstruct a GRN representing all TF-gene interations mediated by TFs in CLL, by infering the binding of all detectable TFs from the acessibility profiles of all CLL samples together (Figure [#Figure4]b). In this general CLL GRN, TFs are the most well-connected nodes (Figure [#Figure4]c), among which some general regulators (SP1/2, CTCF), but also families of prominent regulators of proliferation (EGR, KLF) and cell cycle (E2F), and some factors of relevance for B cell maturation (SPI1, PAX5, BCL6).

We reasoned that the same approach could be used to infer GRNs of specific CLL molecular subtypes, by considering regulatory interactions predicted from TF footprints in groups of samples. To gain deeper insights into the molecular wiring of CLL groups with different IGHV mutation status, we therefore infered GRNs specific to each group independently (Supplementary Figure [#FigureS4]a/b). Although the overal structure of the networks are similar (Supplementary Figure [#FigureS4]c), by comparison of the node connectivity between the IGHV-specific networks, we find a range of TFs which interaction degree is changing (Figure [#Figure4]d), among these, prominently a group of Homeobox genes of the HOXA cluster. These and other differentially regulated genes have noticeable variation in chromatin acessibility (Figure [#Figure4]e) consistent with the differential accessibility and likely, differential TF binding in their regulatory elements.


## Discussion

###### fig1
By performing ATAC-seq on a large set of primary CLL samples, we have for the first time characterized the chromatin acessibility landscape of CLL across the whole range of clinical subtypes of the disease.
With this approach, we could reach saturation of chromatin acessibility sites in the genome, in a relatively cost-effective manner, highlighting that chromatin acessibility profiling in a large cohort of cancer samples directly from primary patient samples is feasable, and likely to become more prevalent.

The existence of a considerable amount of regions outside the expected context of promoters and known enhancers of whole B cell populations might indicate that these are either novel, CLL-specific enhancers or enhancers of more specific cell populations throughout B cell differentiation. Such a question would be answered with more complete epigenome profiling of B cell populations through their differentiation course.


###### fig2
With this project, we aimed to investigate the variation in chromatin regulation
Most genes have considerable chromatin variation

The finding that genes relevant in the context of CLL (such as recurently mutated genes[@Puente2015]) display lower variability at the chromatin level may suggest that a likely requirement to achieve the molecular phenotype leading to the manifestation of CLL is the consistent .


###### fig3
The identificatin of IGHV-specific regions highlighted the genomic and epigenetic differences between the two major molecular subtypes of CLL: uCLL and mCLL. The functional investigation into these differentially used regions revealed that these likely reflect the potential cell-of-origin of the CLL cells: the mCLL group was enriched in enhancers of B cell lines and in TF binding sites of known factors required for diferentiation into memory B cells (BATF, BCL6), while the uCLL group was enriched in regions used by naïve B cells and other hematopoietic cells, likely reflecting the more undiferentiated state of these cells. A common for both groups is the enrichemnt in regions used specifically for transcription (with the histone mark H3K36me3) in other cell types. We speculate that these are likely novel, CLL-excluusive enhancer regions, an hypotesis that would be clarified by histone mark profiling in the different CLL groups and normal B cell populations.

The enriched signling pathways between the two CLL groups may lead to new discoveries on the molecular phenotypes of these disease groups: inhibitory signaling by CTLA4 is a known pathway in T cells which leads to overal decrease in T cell population expansion [@Kong2014]. However, when CTLA4 was experimentaly downregulated in CLL it showed increased cell survival through the upregulation of molecules involved in B-cell proliferation or survival signaling [@Mittal2013]. The fact that this pathway is exclusively enriched in the mCLL group could contribute to the generally less agressive outcome of mCLL compared with uCLL. Signaling mediated through the Fc region of immunoglobulin E (FCERI) includes prominent BCR signaling molecules such as the proeminint BCR signaling kinases LYN, SYK and BTK that may lead to NFKB activation. This contribution of Ca+ mediated FCERI signaling in CLL seems largely unknown and could therefore be a potential new area of molecular enquiry in the mCLL group. Also in mCLL, FC gamma receptor (FCGR) activation seems interesting: this molecule binds to the Fc portion of immunoglobulin G (IgG) and it is particularly interesting due to its relation with CD20 - a common clinical target in lymmphoid malignancies: high FCGR expression reduces Rituximab efficiency by causing CD20 (an anti-CD20 antibody prescribed in CD20+ CLL) internalization [@Lim2011], and FCGR polymorphisms have been shown to alter the clinical efficacy of rituximab [@Zhuang2015]. Enrichment of the BCR signaling pathway in mCLL was surprising, particularly in light of its higher activity in uCLL. [ Nonetheless, this may indicate that potentially some antagonists of this pathway might be need to be activated in mCLL. ] Antigen-independent cell-autonomous BCR stimulation is also known to broadly occur in CLL [@Minden2012].

The sample substructure we found when measuring chromatin acessibility in IGHV-specific regions is similar to that previously found in DNA methylation studies of CLL [@Kulis2012; @Queiros2014], with a third group of CLL cases (iCLL) which seems to stem primarily from the mCLL group, having a mixture of the chromatin acessibility signature of the uCLL and mCLL groups in these regions. The existence of a fourth group, of iCLL more close to the uCLL one would be of interest and if expanded upon the analysis of more samples, would certain warant further molecular investigation of these cases.

Supervised machine-learning approaches based on clinical annotations seems therefore valuable in uncovering sample substructure and characterization of regions contributing to the distinshion between sample grouups and therefore should be applyied to new groups of cases with features of clinical relevance (*e.g.* different therapies, overal survival).

###### fig4
Although reverse engineering of regulatory networks in human cells has been a goal for a while (inclusively in B cells [@Basso2005; @Lefebvre2010]) large-scale detection of regulatory interactions in a patient-specific manner has never been achieved. This work presents thus the first application of GRN inference through TF footprinting in a cancer setting cancer.

While the network infered from all CLL samples might reveal general cellular activity or pan-cancer signals, some B cell-specific TFs are also shown to have high regulatory activity. Such examples include important factors for B cell maturation such as BCL6, PAX5 and SPI1, but other factors known to be markers of CLL progression such as E2F4 [@Falt2005].

Among the most differentially regulated TFs between uCLL and mCLL, we have detected a few known cases that have IGHV mutation status-dependent regulation: 
several genes in the HOXA cluster are among the most differentialy regulated genes, although interestingly, in opposing directions. Since this locus has been shown to be hypermethylated in mCLL [@Strathdee2006] we hypotesise that DNA methylation at its locus in mCLL and prevents the binding of TFs either directly or by recruitment of chromatin modifiers which repress the locus' chromatin. 

Inference of GRNs through TF foootprinting seems therefore capable of retrieving known regulatory interactions in disease subgroups and due to the genome-wide nature nature of ATAC-seq data, is likely to reveal many new interactions of potential clinical relevance. This approach should be of particular use in cancer, where where the gene regulatory program of cells is often changed, but also in percision medicine in the context of several diseases since it is possible to infer GRNs with relative confidence albeit with low sensitivity from a single ATAC-seq sample.

###### end
This work systematically dissects the regulatory patterns and inter-patient heterogeneity of CLL based on chromatin data for a large patient cohort – finding groups of molecular phenotypes with distinct usage of regulatory regions and signalling pathways, and specific regulatory interactions of potential therapeutic relevance.


## Acknowledgements
We thank the Biomedical Sequencing Facility at CeMM for assistance with next-generation sequencing and all members of the Bock lab for their help and advice. This work was performed in the context of the BLUEPRINT project (European Union’s Seventh Framework Programme grant agreement no. 282510) and funded in part by the ERA-NET project CINOCA (FWF grant agreement no. I 1626-B22). C.S. was supported by a Feodor Lynen Fellowship of the Alexander von Humbold
t Foundation. C.B. was supported by a New Frontiers Group award of the Austrian Academy of Sciences.


## Authorship
**Contribution:**
R.W. and Z.D. followed the patients and isolated lymphocytes from periferal blood, D.O. and J.S. contributed the samples, C.S. performed the experiments; A.F.R. analyzed the data; A.F.R., C.S., D.O., J.S. and C.B. planned the study; C.B. supervised the research; all authors contributed to the writing of the manuscript.

**Conflict-of-interest disclosure:**
The authors declare no competing financial interests.

**Correspondence:**
Christoph Bock, CeMM Research Center for Molecular Medicine of the Austrian Academy of Sciences, Vienna, Austria; e-mail: <cbock@cemm.oeaw.ac.at>.


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
Sequencing statistics

#### Figure: {#FigureS2}
Library quality

#### Figure: {#FigureS3}
Cohort statistics

#### Figure: {#FigureS4}
Sample-wise sequencing saturation

#### Figure: {#FigureS5}
Genomic location and chromatin-state enrichment of the CLL cohort region set.

#### Figure: {#FigureS6}
Features of the CLL cohort region set: support, qv2, etc...

#### Figure: {#FigureS7}
Heterogenity in chromatin accessibility in B-cell and CLL related genes.

#### Figure: {#FigureS8}
Feature of IGHV mutation status regions: support, qv2, chromosomal locations.

#### Figure: {#FigureS9}
Genomic location and chromatin-state enrichment of IGHV mutation status regions.

#### Figure: {#FigureS10}
PIQ stats

#### Figure: {#FigureS11}
CLL network with all interactions

#### Figure: {#FigureS12}
IGHV-unmutated CLL network

#### Figure: {#FigureS13}
IGHV-mutated CLL network

#### Figure: {#FigureS14}
Fold-change of degree for all nodes in the two networks, gene ontology of these.


### Tables
#### Table: {#Table1}
ATAC-seq statistics

#### Table: {#Table2}
Saturated consensus calling of CLL regulatory elements (BED file)

### Website
Track hub with consensus (6 tracks: percentile 5, 25, 50, 75, 95, mean) & all the individual tracks

## References
