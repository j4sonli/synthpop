library(dplyr)
library(tidyr)
library(bnlearn)
library(foreach)
library(doSNOW)
library(stringr)

old <- Sys.time()

options(dplyr.summarise.inform=F)

rm(list=ls())
graphics.off()
cat("\014")
set.seed(1033)


start_time <- Sys.time()


##### RUN SETTINGS #####

pumas_file <- "PUMAs.csv"

my_census_api_key <- "ac6cb3e106c860e52384fe71cf0407a13c25b96c"

process_pums_file <- "process_pums-20210910-gq.R"

process_marg_file <- "process_marg-20210909.R"
use_prev_marg <- T
marg_dir <- "./Optim_Data/ACS_Marg/cambridge/"

puma_tract_equiv_file <- "./Optim_Data/2010_Census_Tract_to_2010_PUMA.csv"

tract_block_equiv_file <- "./data/tract_block_equivs/COL2010_TAB2010_ST_25_v2/COL2010_TAB2010_ST_25_v2.txt"

output_suffix <- "20210915-gq"
sampled_hhs_file <- paste0("./output/sampled_hhs-", output_suffix, ".csv")
sampled_hh_indvs_file <- paste0("./output/sampled_indvs-", output_suffix, ".csv")
hh_pool_file <- paste0("./Optim_Data/hh_pool-", output_suffix, ".csv")
indv_pool_file <- paste0("./Optim_Data/indv_pool-", output_suffix, ".csv")

oversample_by <- 3

##### PREPARE DATA #####

pumas <- sort(read.csv(pumas_file, colClasses=c("character"))$pumas)
n_pumas <- length(pumas)

puma_tract_equiv <- read.csv(puma_tract_equiv_file)
puma_tract_join <- puma_tract_equiv %>%
  mutate(puma = paste0(str_pad(as.character(STATEFP), 2, pad="0"),
                       str_pad(as.character(PUMA5CE), 5, pad="0"))) %>%
  filter(puma %in% pumas) %>%
  mutate(tract = paste0(str_pad(as.character(STATEFP), 2, pad="0"),
                        str_pad(as.character(COUNTYFP), 3, pad="0"),
                        str_pad(as.character(TRACTCE), 6, pad="0")))

tracts <- puma_tract_join$tract

tract_block_equiv <- read.csv(tract_block_equiv_file) %>%
  select(TAB_STATE_2010, TAB_COUNTY_2010, TRACT_2010, BLK_2010)
tract_block_join <- tract_block_equiv %>%
  mutate(tract = paste0(str_pad(as.character(TAB_STATE_2010), 2, pad="0"),
                        str_pad(as.character(TAB_COUNTY_2010), 3, pad="0"),
                        str_pad(as.character(TRACT_2010), 6, pad="0"))) %>%
  filter(tract %in% tracts) %>%
  group_by(tract) %>%
  summarise(blkgp_ids = list(sort(unique(substr(BLK_2010, 1, 1))))) %>%
  unnest(blkgp_ids) %>% 
  mutate(blkgp = paste0(tract, blkgp_ids)) %>%
  data.frame()

blkgps <- tract_block_join$blkgp
  

trs <- function(w){
  w.int <- floor(w) 
  rems <- w - w.int 
  topup <- sample(length(w), size = round(sum(rems)), prob = rems)
  w.int[topup] <- w.int[topup] + 1
  return(w.int)
}

## get the microdata
cat("Preparing PUMS microdata...")
source(process_pums_file)

## get marginal data
cat("Preparing marginal data...")
source(process_marg_file)

t1 <- Sys.time()

##### HH BNs: LEARN STRUCTURES AND PARAMETERS #####

n_types <- 4

cat("Creating household BNs...")

## prepare training dataset
hh.train <- mapply(function(x, type) {
  x <- x %>%
    select(-c(HHID, hh_type, puma)) %>%
    slice(rep(1:n(), hhwt)) %>%
    select(-hhwt)
  if(type == "SM"){
    x <- select(x, -c(hhsiz, nchil))
  } else if(type == "NF"){
    x <- select(x, -c(nchil))
  }
  return(x)
}, pums_hh.type.puma[1:4], rep(colnames(pums_hh.type.puma)[1:4], each=n_pumas)) %>%
  matrix(nrow=n_types, ncol=n_pumas, byrow=T)

bn_score <- "aic"

## learn structure of each PUMA BN
hh_bns.puma <- lapply(hh.train, function(data) hc(data, score=bn_score)) %>%
  matrix(nrow=n_types, ncol=n_pumas)

## find arc strength
hh_arc_strengths.puma <- mapply(function(bn,data) arc.strength(bn, data, criterion=bn_score), t(hh_bns.puma), t(hh.train), SIMPLIFY=F)

## plot each PUMA BN
struc_names <- c("Married Couples", "No Spouse", "Single-member", "Nonfamily, not living alone", "Group Quarters")[1:n_types]
par(mfrow=c(n_types,n_pumas), mar=c(2,3,2,1))
par(mfrow=c(2,2), mar=c(2,3,2,1))
invisible(mapply(function(bn,data,arc_s,n,struc_name){
    strength.plot(bn,arc_s,main=paste0(struc_name, " (AIC: ", round(AIC(bn, data), digits=2), ") (n=", n, ")"), shape="ellipse")
  }, t(hh_bns.puma), t(hh.train), hh_arc_strengths.puma, sapply(t(hh.train), nrow), rep(struc_names, each=n_pumas)))

## fit parameters of each BN
hh_ordinals <- c("hhinc", "nwork", "ncars", "nchil", "hhsiz")
hh_bns.puma.fit <- mapply(function(bn, data){
    bn.fit(bn, data, method="mle", replace.unidentifiable=T, ordinal=hh_ordinals)
  }, hh_bns.puma, hh.train)


##### HH_INDV BNs: LEARN STRUCTURES AND PARAMETERS #####

cat("Creating individual BNs...")

## prepare training dataset
hh_indv.train <- mapply(function(x, type) {
  x <- x %>%
    select(-c(HHID, INDVID, hhwt, hh_type, puma)) %>%
    slice(rep(1:n(), wt)) %>%
    select(-wt)
  if(type == "SM"){
    x <- select(x, -c(hhsiz, i_inc, nchil))
  } else if(type == "NF"){
    x <- select(x, -c(nchil))
  }
  return(x)
}, pums_hh_indv.type.puma[1:4], rep(colnames(pums_hh.type.puma)[1:4], each=n_pumas)) %>%
  matrix(nrow=n_types, ncol=n_pumas, byrow=T)

## prohibit edges between HH-level variables, which are useless
hh_varss <- lapply(hh.train, names)
hh_indv_varss <- lapply(hh_indv.train, names)
hh_indv_prohibits <- mapply(function(hh_vars, hh_indv_vars) {
  left_join(data.frame(hh_indv_vars), data.frame(hh_vars), by=character()) %>%
    rename(from=hh_indv_vars, to=hh_vars)
}, hh_varss, hh_indv_varss, SIMPLIFY=F)

## learn structure of each PUMA BN
hh_indv_bns.puma <- mapply(function(data, prohib) {
  hc(data, score=bn_score, blacklist=prohib)
}, hh_indv.train, hh_indv_prohibits, SIMPLIFY=F) %>%
  matrix(nrow=n_types, ncol=n_pumas)

## find arc strengths
hh_indv_arc_strengths.puma <- mapply(function(bn,data) arc.strength(bn, data, criterion=bn_score), t(hh_indv_bns.puma), t(hh_indv.train), SIMPLIFY=F)

## plot each PUMA BN
par(mfrow=c(n_types,n_pumas), mar=c(2,3,2,1))
invisible(mapply(function(bn,data,arc_s,n,struc_name){
  strength.plot(bn,arc_s,main=paste0(struc_name, " (AIC: ", round(AIC(bn, data), digits=2), ") (n=", n, ")"), shape="ellipse")
}, t(hh_indv_bns.puma), t(hh_indv.train), hh_indv_arc_strengths.puma, sapply(t(hh_indv.train), nrow), rep(struc_names, each=n_pumas)))

## fit parameters of each BN
hh_indv_ordinals <- c("hhinc", "nwork", "ncars", "nchil", "hhsiz", "i_age", "i_inc")
hh_indv_bns.puma.fit <- mapply(function(bn, data){
  bn.fit(bn, data, method="mle", replace.unidentifiable=T, ordinal=hh_indv_ordinals)
}, hh_indv_bns.puma, hh_indv.train)

## build the GQ BN
gq_data <- pums_hh_indv.type.puma[[5]] %>%
  select(-c(HHID, INDVID, hhwt, hh_type, puma)) %>%
  slice(rep(1:n(), wt)) %>%
  select(-c(wt, tenur, hhinc, nwork, ncars, nchil, hhsiz))
gq_bn <- hc(gq_data, score=bn_score, blacklist=hh_indv_prohibits[[1]] %>% 
              filter(!(from %in% c("tenur", "nwork", "ncars", "nchil", "hhsiz", "hhinc"))) %>%
              filter(!(to %in% c("tenur", "nwork", "ncars", "nchil", "hhsiz", "hhinc"))))
gq_bn_arc_strengths <- arc.strength(gq_bn, gq_data, criterion=bn_score)
strength.plot(gq_bn,gq_bn_arc_strengths, main=paste0("Group Quarters (AIC: ", round(AIC(gq_bn, gq_data), digits=2), ")"), shape="ellipse")
gq_bn.fit <- bn.fit(gq_bn, gq_data, method="mle", replace.unidentifiable=T, ordinal=hh_indv_ordinals)

t2 <- Sys.time()


##### BN: SAMPLE HHs AND INDVs #####

## get HH struc counts according to marginal, and multiply by oversampling amount
hh_strucs_by_struc <- puma_hhtype %>% 
  select(-puma_id) %>% 
  `*`(1 + oversample_by) %>%
  ceiling() %>%
  as.list()

## sample HHs, separately for each PUMA
cat("Sampling ", sum(pop_by_puma$puma_pop), " total households")

cl <- makeCluster(4)
registerDoSNOW(cl)
sampled_hhs.strucs.puma <- 
  foreach(hhstruc_i=1:n_types, hhstruc_pops=hh_strucs_by_struc, .packages=c("dplyr","bnlearn")) %:%
  foreach(puma=factor(pumas, levels=levels(pums_hh_indv$puma)), puma_pop=hhstruc_pops) %dopar% {
    puma_struc_sample <- hh_bns.puma.fit[[hhstruc_i]] %>% 
      rbn(n = puma_pop) %>%
      mutate(hhtype = names(pums_hh_indv.type)[hhstruc_i],
             puma = puma)
    ## fill in missing hhsiz and nchil
    if(hhstruc_i == 3){
      puma_struc_sample <- puma_struc_sample %>%
        mutate(hhsiz = ordered(rep("1", puma_pop), levels=levels(hh.train[[1]]$hhsiz)),
               nchil = ordered(rep("0", puma_pop), levels=levels(hh.train[[1]]$nchil)))
    } else if(hhstruc_i == 4){
      puma_struc_sample <- puma_struc_sample %>%
        mutate(nchil = ordered(rep("0", puma_pop), levels=levels(hh.train[[1]]$nchil)))
    }
    puma_struc_sample
  } %>%
  unlist(recursive=F) %>%
  matrix(nrow=n_pumas)
stopCluster(cl)

## add HHIDs to each PUMA's samples
sampled_hhs.strucs.puma <- 
  lapply(1:n_pumas, function(puma_i){
    n_hh_struc_samples <- sapply(hh_strucs_by_struc, function(type) type[puma_i])
    lapply(1:n_types, function(hhstruc_i) {
      sampled_hhs.strucs.puma[[puma_i,hhstruc_i]] %>% 
        mutate(HHID=seq_len(nrow(.)) + sum(n_hh_struc_samples[1:hhstruc_i]) - n_hh_struc_samples[hhstruc_i]) %>%
        select(HHID, everything())
    })
  }) %>%
  unlist(recursive=F) %>%
  matrix(nrow=n_pumas, byrow=T)

sampled_hhs.puma <- base::apply(sampled_hhs.strucs.puma, 1, bind_rows)

## aggregate by household variables
sampled_hhs.strucs.agg.puma <- lapply(sampled_hhs.strucs.puma, function(sampled_hhs.struc){
    sampled_hhs.struc %>% 
      group_by_at(vars(-c(HHID))) %>%
      summarise(ct=n(), hhids=list(HHID)) %>%
      data.frame()
  }) %>%
  matrix(nrow=n_pumas)

## sample individuals to go in sampled HHs for each PUMA
cat("Sampling individuals")

hhsiz6p_dist <- pums_raw$NP %>% 
  table() %>%
  .[min(which(names(.) >= 6)):length(.)]

cl <- makeCluster(4)
registerDoSNOW(cl)
pb <- txtProgressBar(max=sum(sapply(sampled_hhs.strucs.agg.puma, nrow)), style=3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))
sampled_hh_indvs.puma <- 
  foreach(puma_i=1:n_pumas, .packages=c("dplyr","bnlearn"), .options.snow=opts) %:%
  foreach(hhstruc_i=1:n_types, .combine=bind_rows) %:%
  foreach(row=seq_len(nrow(sampled_hhs.strucs.agg.puma[[puma_i,hhstruc_i]])), combine=rbind) %dopar% {
    hh <- sampled_hhs.strucs.agg.puma[[puma_i,hhstruc_i]][row,]
    evid <- list(tenur=hh$tenur, hhinc=hh$hhinc, nwork=hh$nwork, ncars=hh$ncars)
    if(hhstruc_i != 3){
      evid$hhsiz <- hh$hhsiz
    }
    if(hhstruc_i != 3 && hhstruc_i != 4){
      evid$nchil <- hh$nchil
    }
    # if household size is "6+", resample household size from tail end distribution of PUMS
    true_hhsiz <- as.integer(hh$hhsiz)
    if(true_hhsiz == 6){
      true_hhsiz <- sample(as.integer(names(hhsiz6p_dist)), 1, prob=hhsiz6p_dist)
    }
    hh_indvs <- mutilated(hh_indv_bns.puma.fit[[hhstruc_i]], evidence=evid) %>% 
      rbn(n = hh$ct * true_hhsiz) %>%
      mutate(HHID = rep(unlist(hh$hhids), each=true_hhsiz)) %>%
      select(HHID, everything())
    if(hhstruc_i == 3 || hhstruc_i == 4){
      hh_indvs$nchil <- ordered(rep(0, nrow(hh_indvs)), levels=levels(hh.train[[1]]$nchil))
      if(hhstruc_i == 3){
        hh_indvs$hhsiz <- ordered(rep(1, nrow(hh_indvs)), levels=levels(hh.train[[1]]$hhsiz))
        hh_indvs$i_inc <- ordered(recode(hh_indvs$hhinc,
                                         "101-150K"=">100K",
                                         ">150K"=">100K"), levels=levels(hh_indv.train[[1]]$i_inc))
        }
      hh_indvs <- select(hh_indvs, 1:5, nchil, everything())
    }
    hh_indvs %>%
      mutate(hhtype = names(pums_hh_indv.type)[hhstruc_i],
             puma = pumas[puma_i])
  } %>%
  matrix(nrow=n_pumas)
close(pb)
stopCluster(cl)

## fix IDs
sampled_hh_indvs.puma <- lapply(sampled_hh_indvs.puma, function(sample) {
  sample %>%
  arrange(HHID) %>%
  mutate(INDVID = as.integer(seq_len(nrow(.)))) %>%
    # mutate(HHID = rep(1:length(rle(.$HHID)$lengths), rle(.$HHID)$lengths), 
    #        INDVID = as.integer(seq_len(nrow(.)))) %>%
    select(HHID, INDVID, everything())
})

## sample GQ individuals
n_gq <- 16660*(1+oversample_by)
gq_indvs <- rbn(gq_bn.fit, n=n_gq) %>%
  mutate(tenur = "GQ",
         hhinc = "GQ",
         nwork = "GQ",
         hhtype = "GQ",
         hhsiz = ordered(rep(1, n()), levels=levels(sampled_hhs.puma[[1]]$hhsiz)))

write.csv(bind_rows(sampled_hhs.puma, gq_indvs), file=sampled_hhs_file, row.names=F)
write.csv(bind_rows(sampled_hh_indvs.puma, gq_indvs), file=sampled_hh_indvs_file, row.names=F)

t3 <- Sys.time()


##### OPTIMIZATION #####
## introduce proxy variables matching optimization marginals
hh_pool <- bind_rows(sampled_hhs.puma, gq_indvs) %>%
  mutate(tenur_hhinc_prox = paste0(recode(tenur,
                                          "Owned with mortgage or loan"="Owned",
                                          "Owned free and clear"="Owned"), 
                                   ".",
                                   recode(hhinc, 
                                          "<=0K"="<25K",
                                          "1-25K"="<25K")),
         tenur_hhsiz_prox = paste0(recode(tenur,
                                          "Owned with mortgage or loan"="Owned",
                                          "Owned free and clear"="Owned"), 
                                   ".",
                                   recode(hhsiz, 
                                          "5"="5+", 
                                          "6+"="5+")),
         tenur_prox = recode(tenur,
                             "Owned with mortgage or loan"="Owned",
                             "Owned free and clear"="Owned"))
indv_pool <- bind_rows(sampled_hh_indvs.puma, gq_indvs)

write.csv(hh_pool, file=hh_pool_file, row.names=F)
write.csv(indv_pool, file=indv_pool_file, row.names=F)


end_time <- Sys.time()