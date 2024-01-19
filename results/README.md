TF-target gene binding networks files:  

`network.tsv.gz` All networks in one table  
`network.xslx` All networks in one file, separate sheets for each cell line, all TF-target genes and TF-TF networks  
`network_Elav.tsv.gz`  Network for elav1 line  
`network_Fox.tsv.gz`  Network for foxQ2d line  
`network_Ncol.tsv.gz`  Network for ncol3 line  

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


