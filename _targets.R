library('targets')
library('tarchetypes')
library(crew)
tar_option_set(controller = crew_controller_local(workers = 10)) 
tar_option_set(seed = 11042012)
Sys.setenv(VROOM_CONNECTION_SIZE = as.character(10 * 1024 * 1024)) #For the large metab datasets

####This runs a workflow for cleaning the MESA TOPMed multi-omics project metabolite data, formatting it into the following format: wide for metabolites & long for exams, and providing basic quality control (QC) metrics



tar_option_set(packages = c("dplyr", "tidyr", "tibble", "readr", "data.table", "bit64", "foreign", "quarto", 
                            "rlang", "purrr", "rcompanion", "mixOmics"))

tar_source("/media/Analyses/MESA-Metabolomics-PCA/R")

list(
 
  #---------------------------------------------------------------------------------------#
  #--------------------------------1. Build metabolites-----------------------------------#
  #---------------------------------------------------------------------------------------#
  
  #Abundance tables
  tar_target(path_amide, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/25_0107_TOPMed_MESA_Amide-neg_rev031325.csv", format = "file"),
  tar_target(path_C8, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_C8-pos_checksums_rev031325.csv", format = "file"),
  tar_target(path_C18, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_C18-neg_checksums_rev031325.csv", format = "file"),
  tar_target(path_HILIC,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_HILIC-pos_checksums_rev031325.csv", format = "file"),
  
  #Info tables
  tar_target(path_amide_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_AmideNeg_SampleInfo_20250329.txt", format = "file"),
  tar_target(path_C8_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_C8Pos_SampleInfo_20250329.txt", format = "file"),
  tar_target(path_C18_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_C18Neg_SampleInfo_20250329.txt", format = "file"),
  tar_target(path_HILIC_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_HILIC-Pos_SampleInfo_20250329.txt", format = "file"),
  
  #Bridging file
  tar_target(path_bridge,"/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESA-SHARE_IDList_Labeled.csv", format = "file"),
  
  #Run function
  tar_target(build_metabs_out,   
             build_metabs_function(
               path_amide = path_amide,
               path_C8 = path_C8,
               path_C18 = path_C18,
               path_HILIC = path_HILIC,
               path_amide_info = path_amide_info,
               path_C8_info = path_C8_info,
               path_C18_info = path_C18_info,
               path_HILIC_info = path_HILIC_info,
               path_bridge = path_bridge,
               QC_label = "QC-pooled_ref")
  ),
  
  #Save outputs
  #Initial metab info after cleaning
  tar_target(build_metabs_QC, build_metabs_out$QC_info_out),
  #Cleaned long metabs file (with duplicates)
  tar_target(Metabs_long, build_metabs_out$metabs_long),
  #QC info in mapping file
  tar_target(build_mapping_QC, build_metabs_out$QC_info_mapping_out),
  #Mapping file
  tar_target(Mapping_file, build_metabs_out$mapping_final),
  #CV and missingness info
  tar_target(Metabolite_CV_and_missingness, build_metabs_out$CV_info),
  #File with duplicates info
  tar_target(Duplicates_flag_file, build_metabs_out$Duplicates_file),
  #Final metabs, cleaned and without duplicates
  tar_target(Final_metabs_long, build_metabs_out$metabolite_file_nodupes),
  tar_target(Final_metabs_long_info,build_metabs_out$metabolite_file_nodupes_dims),
  
 

  
  #----------------------------------------------------------------------------#
  #-------------------------2. PCA                     ------------------------#
  #----------------------------------------------------------------------------#

  
  #------------Files-------
  
 # tar_target(Metab_PCA, 
 # PCA_function(datafile = Final_metabs_long,
 # metab_mapping = Mapping_file,
 # id_var   <- "sidno",
 # block_size = 50)
 # ),
  

  #---------------------------------
  #---4B. Run individual sPLS-DA----
  #---------------------------------
  
  #------------Metabs-------
  # T2D in one year vs no T2D in 10 years
  
  
  
  

  #----------------------------------------------------------------------------#
  #--------------------------------Quarto file---------------------------------#
  #----------------------------------------------------------------------------#
  
  tarchetypes::tar_quarto(
   build_proteins_quarto,
   path = "/media/Analyses/MESA-Metabolomics-PCA/MESA-Metabolomics-PCA.qmd",
   quiet = FALSE
  )
  
)
  