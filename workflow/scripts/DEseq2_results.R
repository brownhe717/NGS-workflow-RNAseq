# setup ------------------------------------------------------------------------
## load necessary packages
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))

# import dds object ------------------------------------------------------------
## import dds2 object
dds <- readRDS(snakemake@input[["dds"]])

# get DEseq results for all contrasts ------------------------------------------
contrast <- c("condition", snakemake@params[["contrast"]])
results <- results(dds, contrast = contrast, alpha = snakemake@params[["padj_cutoff"]]) %>% 
  as.data.frame() %>% 
  arrange(padj) %>% 
  rownames_to_column(var = "gene_id")

# add additional information to results table ----------------------------------
# #read in GTF file
gtf <- rtracklayer::import(snakemake@input[["annotation"]]) %>% 
    mcols() %>%
    as.data.frame() %>%
    filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name)

# add gene symbol to results
out_table <- dplyr::left_join(results, gtf, by = "gene_id")

## add column indicating if gene is differentially expressed with padj < 0.05 and FC > 2
out_table <- out_table %>%
  dplyr::mutate(is_diff = (padj < snakemake@params[["padj_cutoff"]] & (abs(log2FoldChange) > snakemake@params[["FC_cutoff"]]))) %>%
  replace_na(list(is_diff = FALSE))

# write output file ------------------------------------------------------------
write_tsv(out_table, snakemake@output[[1]])
