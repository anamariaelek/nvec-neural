Analysis files:

`01_RNAseq.Rmd` analysis of RNA-seq data\
`02_ATACseq.Rmd` analysis of ATAC-seq data\
`03_Combined.Rmd` joint analysis of TF-gene regulatory networks

Directories:

`RNASEQ_QUANTIFICATION/` data and results used for analysis in `01_RNAseq.Rmd`\
`ATACSEQ/nucleosome_free_regions/` data and results used for analysis in `02_ATACseq.Rmd`\
`plots` plots produced in joint analysis of TF-gene regulatory networks\
`results` result of analysis in `03_Combined.Rmd`\
`published` published data used for Pou4 target comparisons in `03_Combined.Rmd`

Analysis workflow:

Genes were clustered in 7 groups based on their expression fold change in different transgenic lines (Fig1A). Similarly, accessible regions (i.e. ATAC peaks) were also clustered, resulting in the 7 groups with accessibility patterns across transgenic lines equivalent to those of gene expression (SFig2). Finally, gene accessibility scores were calculated as weighted sum of the accessibility of peaks assigned to nearby genes, and genes were also clustered in the same 7 clusters based on gene accessibility scores across transgenic lines (Fig1B).

Set of PWMs for archetype motifs was scored across all peaks, and motif hits were determined as those positions in peaks where PWM score was higher than 95th quantile of genomic distribution of motif scores (i.e. motif PWM scores in random genomic positions). Motif enrichment was then calculated for groups of peaks defined earlier (Fig1C).

Footprint scores were calculated for all accessible regions i.e. ATAC peaks (using TOBIAS ATACorrect and FootprintScores). Motif scores and footprint scores were combined (using TOBIAS BINDetect) to determine motif binding in ATAC peaks.

We filtered all genes by expression (log fold change \> 0) and targets by footprint score (bound sites from TOBIAS BINDetect) to construct reporter line-specific Gene Regulatory Networks (GRNs).
