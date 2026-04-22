#rm(list=ls())
#library(targets)
#library(tidyverse)
#library(mixOmics)
#library(ggplot2)
#datafile = tar_read(Final_metabs_long)
#metab_mapping = tar_read(Mapping_file)
#set.seed(11042012) 
#id_var   <- "sidno"
#block_size = 50

  
PCA_function <- function(datafile) {
  set.seed(11042012) 
  
  #---------------------Data prep-----------
  

  
  dat_full <- datafile |>
    dplyr::filter(exam==1) |>
    dplyr::left_join(
      read.table("/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_C18Neg_SampleInfo_20250329.txt", header=TRUE, sep="\t") |>
        dplyr::filter(exam==1) |>
        dplyr::rename(sidno = subject_id),
      by = dplyr::join_by(TOM_ID, sidno, exam)) |>
    dplyr::mutate(Column_number = as.factor(Column_number))


   metab_names <- metab_mapping$Metabolite
  
   dat <- dat_full |>
     dplyr::select(dplyr::all_of(id_var), dplyr::any_of(metab_names))
   #---------------------Greedy-algorithm-----------
  
  
  target_ids <- nrow(dat) / 100* 75
  required <- c(id_var)
  preds    <- setdiff(names(dat), required)
  
  #First drop obvious missingprobs
  na_mat  <- is.na(dat[, preds, drop = FALSE])
  miss_rate <- as.data.frame(t(colMeans(na_mat))) 
  names(miss_rate) <- "missing"
  miss_rate$Var <- row.names(miss_rate)
  range(miss_rate$missing)
  
  #none (sob)
  rm(na_mat)
  rm(miss_rate)
  
  #Run algorithm
  
  id <- dat[[id_var]]
  id_index <- as.integer(factor(id))
  
  ok_rows <- complete.cases(dat[, required, drop = FALSE])
  na_mat  <- is.na(dat[, preds, drop = FALSE])
  
  count_valid_ids <- function(ok_rows, id_index) {
    sum(rowsum(as.integer(ok_rows), group = id_index, reorder = FALSE) > 0)
  }
  
  ok_ids_n <- count_valid_ids(ok_rows, id_index)
  if (ok_ids_n < target_ids) stop("Not enough IDs even with required vars.")
  
  chosen <- character(0)
  
  #block_size <- 20   # tune this: 100–500 is usually good
  
  repeat {
    if (ok_ids_n < target_ids) break
    if (length(chosen) == length(preds)) break
    
    remaining <- setdiff(seq_along(preds), match(chosen, preds, nomatch = 0L))
    if (!length(remaining)) break
    
    blocks <- split(remaining, ceiling(seq_along(remaining) / block_size))
    
    loss_list <- lapply(blocks, function(idx) {
      # rows that would survive for each candidate in this block
      ok_block <- (!na_mat[, idx, drop = FALSE]) * ok_rows
      
      # number of surviving rows per ID, for each candidate
      id_counts <- rowsum(ok_block, group = id_index, reorder = FALSE)
      
      # number of IDs with >=1 surviving row
      valid_ids_per_var <- colSums(id_counts > 0)
      
      # loss relative to current number of valid IDs
      ok_ids_n - valid_ids_per_var
    })
    
    losses <- unlist(loss_list, use.names = FALSE)
    
    j_best <- remaining[which.min(losses)]
    
    ok_rows <- ok_rows & !na_mat[, j_best]
    chosen  <- c(chosen, preds[j_best])
    ok_ids_n <- count_valid_ids(ok_rows, id_index)
    
    cat("chosen:", length(chosen), " ids:", ok_ids_n, "\n")
  }
  
  rm(required)
  rm(preds)
  rm(ok_rows)
  rm(ok_ids_n)
  rm(j_best)
  rm(losses)
  
  #--------------Make data  
  
  
  dat_PCA <- dat |>
    dplyr::select(dplyr::all_of(id_var), 
                  dplyr::any_of(chosen)) |>
    tidyr::drop_na()
  

  
  #--------------Check unstructured data  
  
  dat_PCA_preds <- dat_PCA |>
    dplyr::select(dplyr::any_of(chosen)) 
  
  
  pca <- mixOmics::pca(dat_PCA_preds, ncomp = 2, scale = TRUE)
  
plotIndiv(pca, ind.names = TRUE,
            legend = FALSE, 
            title = 'MESA Metabs, PCA comp 1 - 2')
  

          
PCA_res <- cbind(dat_PCA, pca$variates$X) |>
  dplyr::left_join(dat_full |>
                     dplyr::select(-dplyr::any_of(chosen)), dplyr::join_by(sidno)) |>
  dplyr::left_join(foreign::read.dta("/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1FinalLabel02092016.dta"), by = dplyr::join_by(idno))

plotIndiv(pca, ind.names = TRUE, group = as.numeric(PCA_res$site1c),
          legend = TRUE, 
          title = 'MESA Metabs, PCA comp 1 - 2')

plot(PCA_res$PC1, PCA_res$Sample_order)


PCA_res  <- PCA_res |>
  dplyr::mutate(PC1_binary = dplyr::case_when(PC1 >= 40 ~ 1,
                                              PC1 < -40 ~0,
                                              TRUE ~ NA_real_))

write.csv(PCA_res |>
            dplyr::filter(PC1_binary == 1) |>
            dplyr::select(idno, sidno, TOM_ID), "/media/Analyses/DEFINE-T2D_Clustering-project/script-building/data-problem/MESA-IDs_PC1.csv", row.names = FALSE)


  
  loadings_matrix <- as.data.frame(pca$loadings$X) |>
    tibble::rownames_to_column("Metabolite") |>
    dplyr::left_join(tar_read(Metabolite_CV_and_missingness), by = dplyr::join_by(Metabolite))
  
  ggplot(loadings_matrix, aes(x = PC1, fill = Assay)) + 
    geom_histogram()
  
  
  
  #--------------Check unstructured C8 data
  
  C8 <- loadings_matrix |>
    dplyr::filter(Assay == "C8") |>
    dplyr::pull(Metabolite)
  
  dat_PCA_preds_C8 <- dat_PCA_preds |>
    dplyr::select(dplyr::any_of(C18)) 
  
  
  pca_C8 <- mixOmics::pca(dat_PCA_preds_C8, ncomp = 2, scale = TRUE)
  
  plotIndiv(pca_C8, ind.names = TRUE,
            legend = FALSE, 
            title = 'MESA C8 Metabs, PCA comp 1 - 2')
  
  #--------------Check unstructured C18 data
  
  C18 <- loadings_matrix |>
    dplyr::filter(Assay == "C18") |>
    dplyr::pull(Metabolite)
  
  dat_PCA_preds_C18 <- dat_PCA_preds |>
    dplyr::select(dplyr::any_of(C18)) 
  
  
  pca_C18 <- mixOmics::pca(dat_PCA_preds_C18, ncomp = 2, scale = TRUE)
  
  plotIndiv(pca_C18, ind.names = TRUE,
            legend = FALSE, 
            title = 'MESA C18 Metabs, PCA comp 1 - 2')
  
  #--------------Check unstructured Amide data
  
  Amide <- loadings_matrix |>
    dplyr::filter(Assay == "Amide") |>
    dplyr::pull(Metabolite)
  
  dat_PCA_preds_Amide <- dat_PCA_preds |>
    dplyr::select(dplyr::any_of(Amide)) 
  
  
  pca_Amide <- mixOmics::pca(dat_PCA_preds_Amide, ncomp = 2, scale = TRUE)
  
  plotIndiv(pca_Amide, ind.names = TRUE,
            legend = FALSE, 
            title = 'MESA Amide Metabs, PCA comp 1 - 2')
  
  #--------------Check unstructured HILIC data
  
  HILIC <- loadings_matrix |>
    dplyr::filter(Assay == "HILIC") |>
    dplyr::pull(Metabolite)
  
  dat_PCA_preds_HILIC <- dat_PCA_preds |>
    dplyr::select(dplyr::any_of(HILIC)) 
  
  
  pca_HILIC <- mixOmics::pca(dat_PCA_preds_HILIC, ncomp = 2, scale = TRUE)
  
  plotIndiv(pca_HILIC, ind.names = TRUE,
            legend = FALSE, 
            title = 'MESA HILIC Metabs, PCA comp 1 - 2')
  
  
  #--------------Check unstructured Not C18 data
  

  dat_PCA_preds_notC18 <- dat_PCA_preds |>
    dplyr::select(-dplyr::any_of(C18)) 
  
  
  pca_notC18 <- mixOmics::pca(dat_PCA_preds_notC18, ncomp = 2, scale = TRUE)
  
  plotIndiv(pca_notC18, ind.names = TRUE,
            legend = FALSE, 
            title = 'MESA Metabs without C18, PCA comp 1 - 2')
  
  
}
  