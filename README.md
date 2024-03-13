### Analysis files:

`01_RNAseq.Rmd` analysis of RNA-seq data\
`02_ATACseq.Rmd` analysis of ATAC-seq data\
`03_Combined.Rmd` joint analysis of TF-gene regulatory networks

### Directories:

`RNASEQ_QUANTIFICATION/` data and results used for analysis in `01_RNAseq.Rmd`\
`ATACSEQ/nucleosome_free_regions/` data and results used for analysis in `02_ATACseq.Rmd`\
`plots` plots produced in joint analysis of TF-gene regulatory networks\
`results` result of analysis in `03_Combined.Rmd`\
`published` published data used for Pou4 target comparisons in `03_Combined.Rmd`

### Analysis workflow:

Genes were clustered in 7 groups based on their expression fold change in different transgenic lines ([Fig2A](figures/Fig2.pdf)). Similarly, accessible regions (i.e. ATAC peaks) were also clustered, resulting in the 7 groups with accessibility patterns across transgenic lines equivalent to those of gene expression ([SFig2](figures/SFig2.pdf)). Finally, gene accessibility scores were calculated as weighted sum of the accessibility of peaks assigned to nearby genes, and genes were also clustered in the same 7 clusters based on gene accessibility scores across transgenic lines ([Fig2B](figures/Fig2.pdf)).

Set of PWMs for archetype motifs was scored across all peaks, and motif hits were determined as those positions in peaks where PWM score was higher than 95th quantile of genomic distribution of motif scores (i.e. motif PWM scores in random genomic positions). Motif enrichment was then calculated for groups of peaks defined earlier ([Fig2C](figures/Fig2.pdf)).

Footprint scores were calculated for all accessible regions i.e. ATAC peaks (using `TOBIAS ATACorrect` and `TOBIAS FootprintScores`). Motif scores and footprint scores were combined (using `TOBIAS BINDetect`) to determine motif binding in ATAC peaks ([Fig2D](figures/Fig2.pdf)).

We filtered all genes by expression and targets by footprint score (bound sites from `TOBIAS BINDetect`) ([Fig2D](figures/Fig2.pdf)) to construct reporter line-specific Gene Regulatory Networks (GRNs, [Fig4](figures/Fig4.pdf)). For Pou4, we compared target genes predicted by our computational GRN inference with published [ChIP-seq](https://elifesciences.org/articles/74336) and [mutant](https://doi.org/10.1016/j.celrep.2020.03.031) data ([Fig5](figures/Fig5.pdf)).
