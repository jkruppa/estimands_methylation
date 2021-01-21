## ------------------------------------------------------------
## get all packages
library(magrittr)
packages <- c("tidyverse", "magrittr", "readxl", "stringr",
              "plyr", "foreach", "haven", "ArrayExpress",
              "data.table", "lumi")
load_packages <- function(packages){
  flag <- lapply(packages, require, character.only = TRUE) %>%
    unlist %>%
    all
  message("All packages loaded: ", flag)
}
load_packages(packages)
## some function overwriting issues
select <- function(...) dplyr::select(...)
rename <- function(...) dplyr::rename(...)
## we remove missing and set Inf values to 7
prepare_values <- function(mat){
 tmp_mat <- na.omit(mat) %>%
    beta2m
 tmp_mat[is.infinite(tmp_mat)] <- 7 * sign(tmp_mat[is.infinite(tmp_mat)])
 return(tmp_mat)
}

## ------------------------------------------------------------
## E_GEOD_55763
in_dir <- file.path("/path/to/array_express/files")
proc_dir <- file.path("/path/to/store/processed/array_express/files")

E_GEOD_55763_sdrf <- read_delim(file.path(in_dir, "E_GEOD_55763", "E-GEOD-55763.sdrf.txt"), "\t")
E_GEOD_55763_pheno_tbl <- tibble(sample_id = str_replace(E_GEOD_55763_sdrf$`Comment [Sample_title]`,
                                                         "Peripheral blood, ", ""),
                                 sample_description = E_GEOD_55763_sdrf$`Comment [Sample_description]`)
write_rds(E_GEOD_55763_pheno_tbl, file.path(proc_dir, "E_GEOD_55763_pheno.rds"))

## !!! CAUTION this will take time... !!! 
## E_GEOD_55763_mat <- fread(file.path(in_dir, "E_GEOD_55763", "GSE55763_normalized_betas.txt"))
## wanted_col_lgl <- str_detect(names(E_GEOD_55763_mat), "Pval", negate = TRUE)
## E_GEOD_55763_reduced_mat <- E_GEOD_55763_mat[, ..wanted_col_lgl]
## write_rds(E_GEOD_55763_reduced_mat, file.path(proc_dir, "E_GEOD_55763_beta.rds"))
## !!! CAUTION this will take time... !!! 

E_GEOD_55763_tbl <- read_rds(file.path(proc_dir, "E_GEOD_55763_beta.rds")) %>%
  as_tibble()

## generate beta matrix
E_GEOD_55763_mat <- E_GEOD_55763_tbl %>%
  select(-ID_REF) %>%
  as.matrix
row.names(E_GEOD_55763_mat) <- E_GEOD_55763_tbl$ID_REF

## order by pheno_tbl
E_GEOD_55763_ord_mat <- E_GEOD_55763_mat[, E_GEOD_55763_pheno_tbl$sample_id]

## write beta and m_values matrix
write_rds(E_GEOD_55763_ord_mat, file.path(proc_dir, "E_GEOD_55763_beta_values_mat.rds"))
write_rds(prepare_values(E_GEOD_55763_ord_mat), file.path(proc_dir, "E_GEOD_55763_m_values_mat.rds"))

## ------------------------------------------------------------

E_GEOD_55763_beta <- read_rds(file.path(proc_dir, "E_GEOD_55763_beta_values_mat.rds"))
E_GEOD_55763_m <- read_rds(file.path(proc_dir, "E_GEOD_55763_m_values_mat.rds"))

## statistics
dim(E_GEOD_55763_m)

## m values
E_GEOD_55763_m_vec <- E_GEOD_55763_m %>%
  as.vector 
E_GEOD_55763_m_vec %>%
  summary
E_GEOD_55763_m_vec %>%
  sd(na.rm = TRUE)


## overall samples
E_GEOD_55763_pheno_tbl <- read_rds(file.path(proc_dir, "E_GEOD_55763_pheno.rds"))
## study samples
E_GEOD_55763_study_pheno_tbl <- E_GEOD_55763_pheno_tbl %>%
  filter(sample_description == "Population study sample.")
## technical samples
E_GEOD_55763_technical_pheno_tbl <- E_GEOD_55763_pheno_tbl %>%
  filter(sample_description != "Population study sample.") %>%
  mutate(group = as.numeric(str_replace(sample_description, ".*group\\s(\\d),.*", "\\1")),
         sample = as.numeric(str_replace(sample_description, ".*sample\\s(\\d+).*", "\\1"))) %>%
  select(-sample_description)

## ------------------------------------------------------------
## by J.Kruppa on Wednesday, August 28, 2019 (17:05)
## E_GEOD_68379  

E_GEOD_68379_sdrf <- read_delim(file.path(in_dir, "E_GEOD_68379", "E-GEOD-68379.sdrf.txt"), "\t")
E_GEOD_68379_pheno_tbl <- tibble(sample_id = E_GEOD_68379_sdrf$'Assay Name',
                                 primary_site = E_GEOD_68379_sdrf$'FactorValue [primary site]',
                                 primary_histology = E_GEOD_68379_sdrf$'FactorValue [primary histology]',
                                 cell_line = E_GEOD_68379_sdrf$'FactorValue [cell line]') %>%
  distinct(cell_line, .keep_all = TRUE) %>%
  mutate(cell_line = toupper(cell_line))
write_rds(E_GEOD_68379_pheno_tbl, file.path(proc_dir, "E_GEOD_68379_pheno.rds"))

## !!! CAUTION this will take time... !!! 
## this will take time... 
## E_GEOD_68379_mat <- fread(file.path(in_dir, "E_GEOD_68379", "GSE68379_Matrix.processed.txt"))
## wanted_col_lgl <- str_detect(names(E_GEOD_68379_mat), "PVal", negate = TRUE)
## E_GEOD_68379_reduced_mat <- E_GEOD_68379_mat[, ..wanted_col_lgl]
## names(E_GEOD_68379_reduced_mat) <- str_replace(names(E_GEOD_68379_reduced_mat), "_AVG.Beta", "")
## write_rds(E_GEOD_68379_reduced_mat, file.path(proc_dir, "E_GEOD_68379_beta.rds"))
## !!! CAUTION this will take time... !!! 

E_GEOD_68379_tbl <- read_rds(file.path(proc_dir, "E_GEOD_68379_beta.rds")) %>%
  as_tibble()

## generate beta matrix
E_GEOD_68379_mat <- E_GEOD_68379_tbl %>%
  select(-V1, -Row.names) %>%
  as.matrix
row.names(E_GEOD_68379_mat) <- E_GEOD_68379_tbl$Row.names
colnames(E_GEOD_68379_mat) <- toupper(colnames(E_GEOD_68379_mat))

## order by pheno_tbl
E_GEOD_68379_ord_mat <- E_GEOD_68379_mat[, E_GEOD_68379_pheno_tbl$cell_line]

## write beta and m_values matrix
write_rds(E_GEOD_68379_ord_mat, file.path(proc_dir, "E_GEOD_68379_beta_values_mat.rds"))
write_rds(prepare_values(E_GEOD_68379_ord_mat), file.path(proc_dir, "E_GEOD_68379_m_values_mat.rds"))

## ------------------------------------------------------------
## by J.Kruppa on Monday, September 30, 2019 (13:04)
## !!! START HERE !!!

E_GEOD_68379_beta <- read_rds(file.path(proc_dir, "E_GEOD_68379_beta_values_mat.rds"))
E_GEOD_68379_m <- read_rds(file.path(proc_dir, "E_GEOD_68379_m_values_mat.rds"))

## overall samples
E_GEOD_68379_pheno_tbl <- read_rds(file.path(proc_dir, "E_GEOD_68379_pheno.rds"))

## statistics
dim(E_GEOD_68379_m)

## m values
E_GEOD_68379_m_vec <- E_GEOD_68379_m %>%
  as.vector 
E_GEOD_68379_m_vec %>%
  summary
E_GEOD_68379_m_vec %>%
  sd(na.rm = TRUE)

## ------------------------------------------------------------
