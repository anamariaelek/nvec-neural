TF-target gene binding networks files:  

`network.tsv.gz` Networks with all genes (TPM>100), before filtering by expression.
`network_expression.tsv.gz` Networks after filtering by expression (fold change (pos vs neg) > 1)  
`network_expression.tsv.gz` Same but in binary format
`network_expression.xslx` Same as `network_expression.tsv.gz`, with separate sheets for each cell line; includes all TF-target genes and TF-TF networks  
`network_Elav.tsv.gz`  Same as `network_expression.tsv.gz` but subset only for Elav1 reporter line  
`network_Fox.tsv.gz`  Same as `network_expression.tsv.gz` but subset only for FoxQ2d reporter line  
`network_Ncol.tsv.gz` Same as `network_expression.tsv.gz`but subset only for Ncol3 reporter line  

Network columns:  

|   column  | description |
|-----------| ------------|
| sample | Sample name |
| gene | TF gene ID |
| name | Known TF name or orthology-based name |
| pfam | TF PFAM annotation
| expression_tpm | TF expression TPM
| expression_lfc | TF expression log2 fold change positive vs negative sample
| target_gene | Target gene ID
| target_name | Target gene known name or orthology-based name
| target_pfam | Target gene PFAM annotation
| target_expression_tpm | Target gene expression TPM
| target_expression_lfc | Target gene expression log2 fold change positive vs negative sample
| motif | Motif archetype assigned to TF
| motif_score | Motif PWM score
| footprint_score | Footprint score
| expression_correlation | Correlation of TF-target gene expression
| target_TF | Is target gene a TF?


