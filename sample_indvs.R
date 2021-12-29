set.seed(1033)

start_time <- Sys.time()

##### BN: SAMPLE INDVs #####
cat("Sampling individuals.....", sep="\n")

## split household pool by household type
hh_pool.type <- split(hh_pool, hh_pool$hhtype) %>%
  .[c("MC", "NS", "SM", "NF")]

## sample individuals to go in sampled HHs for each PUMA
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max=nrow(filter(hh_pool, hhtype != "GQ")), style=3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))
sampled_hh_indvs <- 
  foreach(hhstruc_i=1:n_types, .packages=c("dplyr", "bnlearn"), .options.snow=opts, .combine=bind_rows, .multicombine=T) %:%
  foreach(row=seq_len(nrow(hh_pool.type[[hhstruc_i]])), .combine=bind_rows, .multicombine=T) %dopar% {
    hh <- hh_pool.type[[hhstruc_i]][row,]
    evid <- list(tenur=hh$tenur, hhinc=hh$hhinc, nwork=hh$nwork, ncars=hh$ncars)
    if(hhstruc_i != 3){
      evid$hhsiz <- hh$hhsiz
    }
    if(hhstruc_i %in% c(1,2)){
      evid$nchil <- hh$nchil
    }
    ## mutilate BN
    hh_indv_bn <- mutilated(hh_indv_bns.puma.fit[[hhstruc_i]], evidence=evid)
    hh_indvs <- data.frame()
    ## sample workers; and if MC or NS type, sample children (rejecting any children that are employed)
    ## use previously sampled "true" hhsiz, nwork, and nchil
    true_hhsiz <- as.numeric(hh$true_hhsiz)
    true_nwork <- as.numeric(hh$true_nwork)
    curr_workers <- 0 
    true_nchil <- as.numeric(hh$true_nchil)
    curr_children <- 0
    true_nonemployed_nonchild <- true_hhsiz - true_nwork - true_nchil
    curr_nonemployed_nonchild <- 0
    while(nrow(hh_indvs) < true_hhsiz){
      candidate <- rbn(hh_indv_bn, n = 1)
      if(curr_workers < true_nwork){
        if(candidate$emply == "Armed forces employed" || candidate$emply == "Civilian employed"){
          hh_indvs <- rbind(hh_indvs, candidate)
          curr_workers <- curr_workers + 1
        }
      } else if(curr_children < true_nchil){
        if(candidate$i_age == "<18"){
          hh_indvs <- rbind(hh_indvs, candidate)
          curr_children <- curr_children + 1
        }
      } else if(curr_nonemployed_nonchild < true_nonemployed_nonchild){
        if((candidate$emply == "Not in labor force" || candidate$emply == "Unemployed") && candidate$i_age != "<18"){
          hh_indvs <- rbind(hh_indvs, candidate)
          curr_nonemployed_nonchild <- curr_nonemployed_nonchild + 1
        }
      }
    }
    hh_indvs <- hh_indvs %>%
      mutate(HHID = hh$HHID) %>%
      select(HHID, everything())
    ## handle variables for SM and NF households
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
             puma = pumas)
  }
close(pb)
stopCluster(cl)

## fix IDs
indv_pool <- sampled_hh_indvs %>%
  arrange(HHID) %>%
  mutate(INDVID = as.integer(seq_len(nrow(.)))) %>%
  select(HHID, INDVID, everything())

t1 <- Sys.time()
cat(paste0(".....", round(as.numeric(t1-start_time), 2), " ", units(t1-start_time)), sep="\n")


##### PREPARE FOR OPTIMIZATION #####
cat("Preparing for optimization.....")

indv_pool <- bind_rows(indv_pool, gq_indvs) %>%
  mutate(i_sex_i_age = paste0(i_sex,".",i_age),
         emply_prox = recode(emply,
                             "Civilian employed"="Employed",
                             "Armed forces employed"="Employed"))

write.csv(indv_pool, file=indv_pool_file, row.names=F)

end_time <- Sys.time()
cat(paste0(round(as.numeric(end_time-t1), 2), " ", units(end_time-t1)), sep="\n")

cat("Individual pool completed.", sep="\n")

