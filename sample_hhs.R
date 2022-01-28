set.seed(1033)

start_time <- Sys.time()


##### BN: SAMPLE HHs #####
cat("Sampling household pool (", sep='')

## get HH struc counts according to marginal, and multiply by oversampling amount
hh_strucs_by_struc <- puma_hhtype %>% 
  select(-c(geoid, GQ)) %>% 
  `*`(1 + oversample_by) %>%
  ceiling()

puma_pop <- tract_hhtype %>%
  select(-c(geoid, GQ)) %>% 
  sum()

cat(sum(hh_strucs_by_struc), " households).....")

## compute tail end distribution of household sizes >= 6
hhsiz6p_dist <- pums %>% 
  group_by(HHID) %>%
  summarise(hhwt = first(hhwt),
            hhsiz = first(hhsiz)) %>%
  slice(rep(1:n(), hhwt)) %>% 
  pull(hhsiz) %>%
  table() %>%
  .[min(which(as.integer(names(.)) >= 6)):length(.)]

## compute tail end distribution of number of workers >= 3
## (approximated by sum of employed individuals per HHID)
nwork3p_dist <- pums_hh_indv %>% 
  filter(tenur!="GQ") %>%
  group_by(HHID) %>% 
  summarise(hhwt = first(hhwt),
            emplyd = sum(emply %in% c("Civilian employed", "Armed forces employed"))) %>% 
  slice(rep(1:n(), hhwt)) %>%
  pull(emplyd) %>%
  table() %>%
  .[min(which(as.integer(names(.)) >= 3)):length(.)]

## compute tail end distribution of number of children >= 3
nchil3p_dist <- pums %>%
  group_by(HHID) %>%
  summarise(hhwt = first(hhwt),
            nchil = first(nchil)) %>%
  slice(rep(1:n(), hhwt)) %>%
  pull(nchil) %>%
  table() %>%
  .[min(which(as.integer(names(.)) >= 3)):length(.)]

## sample HHs, separately for each PUMA
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
sampled_hhs <-
  foreach(hhstruc_i=1:n_types, .packages=c("dplyr", "bnlearn"), .combine=bind_rows, .multicombine=T) %dopar% {
    struc_pop <- as.numeric(hh_strucs_by_struc[hhstruc_i])
    struc_sample <- data.frame()
    while(nrow(struc_sample) < struc_pop){
      struc_sample_new <- rbn(hh_bns.puma.fit[[hhstruc_i]], n = struc_pop - nrow(struc_sample)) %>%
        filter(tenur != "GQ")
      to_keep <- rep(T, nrow(struc_sample_new))
      #TODO: batch this process
      true_hhsizs <- if(hhstruc_i != 3) as.integer(substr(struc_sample_new$hhsiz, 1, 1)) else rep(1, nrow(struc_sample_new))
      true_nworks <- as.integer(substr(struc_sample_new$nwork, 1, 1))
      true_nchils <- if(hhstruc_i %in% c(1,2)) as.integer(substr(struc_sample_new$nchil, 1, 1)) else rep(0, nrow(struc_sample_new))
      for(row in seq_len(nrow(struc_sample_new))){
        hh = struc_sample_new[row,]
        evid <- list(tenur=hh$tenur, hhinc=hh$hhinc, nwork=hh$nwork, ncars=hh$ncars)
        if(hhstruc_i != 3){
          evid$hhsiz <- hh$hhsiz
        }
        if(hhstruc_i %in% c(1,2)){
          evid$nchil <- hh$nchil
        }
        ## if household size is "6+", resample household size from tail end distribution of PUMS
        if(true_hhsizs[row] == 6){
          true_hhsizs[row] <- sample(as.integer(names(hhsiz6p_dist)), 1, prob=hhsiz6p_dist)
        }
        true_hhsiz <- true_hhsizs[row]
        ## if nwork or nchil is "3+", resample them from tail end distributions of PUMS
        if(true_nworks[row] == 3){
          true_nworks[row] <- if(length(nwork3p_dist) == 1) as.integer(names(nwork3p_dist)) else sample(as.integer(names(nwork3p_dist)), 1, prob=nwork3p_dist)
        }
        true_nwork <- true_nworks[row]
        if(true_nchils[row] == 3){
          true_nchils[row] <- if(length(nchil3p_dist) == 1) as.integer(names(nchil3p_dist)) else sample(as.integer(names(nchil3p_dist)), 1, prob=nchil3p_dist)
        }
        true_nchil <- true_nchils[row]
        true_nonemployed_nonchild <- true_hhsiz - true_nwork - true_nchil
        ## reject households where:
        ##  nchil + nwork > hhsiz (has child workers);
        ##  nwork > 0 but impossible to sample a worker;
        ##  nchil > 0 but impossible to sample a child;
        ##  none of the above > 0 but impossible to sample a non-employed adult
        p_sample_worker <- cpquery(hh_indv_bns.puma.fit[[hhstruc_i]], method="lw", (emply=="Armed forces employed" | emply=="Civilian employed"), evid)
        p_sample_child <- cpquery(hh_indv_bns.puma.fit[[hhstruc_i]], method="lw", (i_age=="<18"), evid)
        p_sample_nonemployed_nonchild <- cpquery(hh_indv_bns.puma.fit[[hhstruc_i]], method="lw", ((emply == "Not in labor force" | emply == "Unemployed") & i_age != "<18"), evid)
        if((true_nchil + true_nwork > true_hhsiz) ||
           (true_nwork > 0 && (is.na(p_sample_worker) || p_sample_worker == 0)) ||
           (true_nchil > 0 && (is.na(p_sample_child) || p_sample_child == 0)) ||
           (true_nonemployed_nonchild > 0 && (is.na(p_sample_nonemployed_nonchild) || p_sample_nonemployed_nonchild == 0))){
          to_keep[row] <- F
        }
      }
      struc_sample_new <- struc_sample_new %>%
        mutate(true_hhsiz = true_hhsizs,
               true_nwork = true_nworks,
               true_nchil = true_nchils) %>%
        filter(to_keep)
      struc_sample <- rbind(struc_sample, struc_sample_new)
    }
    struc_sample <- struc_sample %>%
      mutate(hhtype = names(pums_hh_indv.type)[hhstruc_i],
             puma = pumas)
    ## fill in missing hhsiz and nchil
    if(hhstruc_i == 3){
      struc_sample <- struc_sample %>%
        mutate(hhsiz = ordered(rep("1", struc_pop), levels=levels(hh.train[[1]]$hhsiz)),
               nchil = ordered(rep("0", struc_pop), levels=levels(hh.train[[1]]$nchil)))
    } else if(hhstruc_i == 4){
      struc_sample <- struc_sample %>%
        mutate(nchil = ordered(rep("0", struc_pop), levels=levels(hh.train[[1]]$nchil)))
    }
    struc_sample
  }
stopCluster(cl)

## add HHIDs to each PUMA's samples
sampled_hhs <- sampled_hhs %>%
  mutate(HHID=seq_len(nrow(.))) %>%
  select(HHID, everything())

## sample GQ individuals
n_gq <- as.integer(nrow(gq_data) * (1 + oversample_by))+1
gq_indvs <- rbn(gq_bn.fit, n=n_gq) %>%
  mutate(HHID = seq(nrow(sampled_hhs) + 1, nrow(sampled_hhs) + n_gq),
         puma = pumas,
         tenur = "GQ",
         hhinc = "GQ",
         nwork = "GQ",
         hhtype = "GQ",
         hhsiz = "GQ")

t1 <- Sys.time()
cat(paste0(round(as.numeric(t1-start_time), 2), " ", units(t1-start_time)), sep="\n")


##### PREPARE FOR OPTIMIZATION #####
cat("Preparing for optimization.....")

## introduce proxy variables matching optimization marginals
hh_pool <- bind_rows(sampled_hhs, gq_indvs) %>%
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
                                   hhsiz),
         tenur_prox = recode(tenur,
                             "Owned with mortgage or loan"="Owned",
                             "Owned free and clear"="Owned"))

write.csv(hh_pool, file=hh_pool_file, row.names=F)

end_time <- Sys.time()
cat(paste0(round(as.numeric(end_time-t1), 2), " ", units(end_time-t1)), sep="\n")

cat("Household pool completed.", sep="\n")

