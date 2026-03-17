# setup ========================================================================
library(GenomicRanges)
library(tidyverse)

# function definitions =========================================================
# define function to compute RPKM
rpkm <- function(count_table, widths) {
  require(tidyverse)
  
  row_names <- tibble(id = rownames(count_table))
  
  
  # calculate rpkm
  kb_widths <- widths / 10^3
  
  column_rpkm <- function(x, widths = kb_widths) {
    cr <- (x / widths) / (sum(x) / 10^6 )
    return(cr)
  }
  
  all_rpkm <- count_table |>
    map(column_rpkm) |>
    as_tibble() |>
    bind_cols(row_names) |>
    column_to_rownames(var = "id")
  
  return(all_rpkm)
  
}

# read in count data ===========================================================
counts <- read_tsv(snakemake@input[[1]], comment = "#") |> 
  select(-c(2:5)) |> 
  column_to_rownames(var = "Geneid")
colnames(counts) <- gsub(".bam", "", basename(colnames(counts)))
# perform RPKM normalization ===================================================
rpkm_table <- rpkm(count_table = select(counts, -c("Length")), widths = counts$Length)

# add gene names to RPKM table =================================================
gene_annotation <- rtracklayer::import(snakemake@input[["annotation"]]) |> 
  mcols() |>
  as.data.frame() |>
  filter(type == "gene") |>
  dplyr::select(gene_id, gene_name)

rpkm_table <- rpkm_table |> 
  rownames_to_column(var = "gene_id") |> 
  left_join(gene_annotation, by = "gene_id") |> 
  select(gene_id, gene_name, everything())

# write output table ===========================================================
rpkm_table |> 
  write_tsv(snakemake@output[[1]])
