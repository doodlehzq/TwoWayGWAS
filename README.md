# TwoWayGWAS
R code for two-way GWAS to study genomic interaction between rice and Xanthomonas oryzae   

A small sample data are provided to illustrate how two-way GWAS was conducted, including
  (1) genotypes for N1 rice accessions
  (2) genotype for N2 Xoo species
  (3) phenotype (leaf lesion length) between each rice accession and Xoo species
 
For each pair of variants (a rice one and a Xoo one), a linear model including N1-1 rice covariants and N2-1 Xoo coraviates is used to test their interaction:

pheno ~ rice_geno + xoo_geno + rice_geno:xoo_geno + rice_covariates + xoo_corariates

