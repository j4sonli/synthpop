suppressPackageStartupMessages({
  library(R.utils)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tidycensus)
  library(bnlearn)
  library(foreach)
  library(doSNOW)
})

options(dplyr.summarise.inform=F)

rm(list=ls())
# graphics.off()
# cat("\014")
set.seed(1033)

synthpop_start_time <- Sys.time()
  
cat("Synthetic population generation for the United States.", sep="\n")


##### HANDLE ARGUMENTS FROM COMMAND LINE CALL #####

args <- commandArgs(trailingOnly=T, asValues=T)

if(length(args) == 0 || (nchar(args$puma) != 7 || grepl("\\D", args$puma))){
  stop("You must supply a numeric 7-digit STATE+PUMA code via the -puma= option.", call. = F)
}
pumas <- args$puma

use_prev_marg <- F
if("use_prev_marg" %in% names(args)){
  use_prev_marg <- args$use_prev_marg
}

marg_dir <- paste0("synthpop_data/acs_marginals/", pumas, "/")
if("marg_dir" %in% names(args)){
  marg_dir <- args$marg_dir
}

oversample_by <- 3
if("samp_rate" %in% names(args)){
  oversample_by <- as.numeric(args$samp_rate)/100-1
}

output_suffix <- ""
if("suffix" %in% names(args)){
  output_suffix <- paste0("-", args$suffix)
}

n_cores <- 4
if("n_cores" %in% names(args)){
  n_cores <- args$n_cores
}


##### RUN SETTINGS #####

Rscript_path <- "/Program Files/R/R-4.1.0/bin/Rscript.exe"

julia_path <- "/Program Files/Julia-1.6.1/bin/julia.exe"

gurobi_home <- '"C:\\Program Files\\gurobi910\\win64"'

gurobi_license <- '"C:\\Program Files\\gurobi910\\gurobi.lic"'

my_census_api_key <- "ac6cb3e106c860e52384fe71cf0407a13c25b96c"


prepare_data_file <- "prepare_data.R"
process_pums_file <- "process_pums.R"
process_marg_file <- "process_marg.R"
train_bns_file <- "train_bns.R"
sample_hhs_file <- "sample_hhs.R"
sample_indvs_file <- "sample_indvs.R"
select_synpop_file <- "select_synpop.jl"
assign_spatial_file <- "assign_spatial.jl"

hh_pool_file <- paste0("synthpop_output2/hh_pool_", pumas, output_suffix, ".csv")
indv_pool_file <- paste0("synthpop_output3/indv_pool_", pumas, output_suffix, ".csv")
syn_hhs_file <- paste0("synthpop_output/syn_hhs_", pumas, output_suffix, ".csv")
syn_indvs_file <- paste0("synthpop_output/syn_indvs_", pumas, output_suffix, ".csv")
syn_hhs_spatial_file <- paste0("synthpop_output/syn_hhs_spatial_", pumas, output_suffix, ".csv")
syn_indvs_spatial_file <- paste0("synthpop_output/syn_indvs_spatial_", pumas, output_suffix, ".csv")

gurobi_logfile <- "synthpop_gurobi_log.txt"

puma_tract_equiv_file <- "synthpop_data/2010_Census_Tract_to_2010_PUMA.csv"

tract_block_equiv_file <- "synthpop_data/tract_block_equiv_2010.Rds"


##### CREATE OUTPUT DIRECTORIES #####

for(output_file in c(hh_pool_file, indv_pool_file, syn_hhs_file, syn_indvs_file, syn_hhs_spatial_file, syn_indvs_spatial_file)){
  if(!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = T)
  }
}


##### DEFINE SYSTEM COMMANDS #####

select_synpop_cmd <- paste0('"', julia_path,'" "', select_synpop_file, '" ', 
                            pumas, " ",
                            hh_pool_file, " ",
                            indv_pool_file, " ",
                            syn_hhs_file, " ",
                            syn_indvs_file, " ",
                            puma_tract_equiv_file, " ",
                            marg_dir, " ",
                            gurobi_logfile, " ",
                            gurobi_home, " ",
                            gurobi_license)

assign_spatial_cmd <- paste0('"', julia_path,'" "', assign_spatial_file, '" ',
                             pumas, " ",
                             syn_hhs_file, " ",
                             syn_hhs_spatial_file, " ",
                             syn_indvs_file, " ",
                             syn_indvs_spatial_file, " ",
                             puma_tract_equiv_file, " ",
                             marg_dir, " ",
                             gurobi_logfile, " ",
                             gurobi_home, " ",
                             gurobi_license)


##### RUN PROCEDURE #####

source(prepare_data_file)

source(train_bns_file)

source(sample_hhs_file)

source(sample_indvs_file)

select_synpop_return <- system(select_synpop_cmd)

assign_spatial_return <- system(assign_spatial_cmd)

#########################

synthpop_end_time <- Sys.time()

cat(paste0("Synthetic population generation process for PUMA ", pumas, " completed."), sep="\n")
cat(paste0("Total time: ", round(synthpop_end_time-synthpop_start_time, 2), 
           " ", units(synthpop_end_time-synthpop_start_time)), sep="\n")

