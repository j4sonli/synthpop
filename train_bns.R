set.seed(1033)

start_time <- Sys.time()


##### HH BNs: LEARN STRUCTURES AND PARAMETERS #####
cat("Learning household BNs.....")

n_types <- 4

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

## suppress numerous bnlearn warnings
default_warn_option <- getOption("warn")
options(warn = -1)

## learn structure of each PUMA BN
hh_bns.puma <- lapply(hh.train, function(data) hc(data, score=bn_score)) %>%
  matrix(nrow=n_types, ncol=n_pumas)

## find arc strength
hh_arc_strengths.puma <- mapply(function(bn,data) arc.strength(bn, data, criterion=bn_score), t(hh_bns.puma), t(hh.train), SIMPLIFY=F)

# ## plot each PUMA BN
# struc_names <- c("Married Couples", "No Spouse", "Single-member", "Nonfamily, not living alone", "Group Quarters")[1:n_types]
# par(mfrow=c(2,2), mar=c(2,3,2,1))
# invisible(mapply(function(bn,data,arc_s,n,struc_name){
#     strength.plot(bn,arc_s,main=paste0(struc_name, " (AIC: ", round(AIC(bn, data), digits=2), ") (n=", n, ")"), shape="ellipse")
#   }, t(hh_bns.puma), t(hh.train), hh_arc_strengths.puma, sapply(t(hh.train), nrow), rep(struc_names, each=n_pumas)))

## fit parameters of each BN
hh_ordinals <- c("hhinc", "nwork", "ncars", "nchil", "hhsiz")
hh_bns.puma.fit <- mapply(function(bn, data){
    bn.fit(bn, data, method="mle", replace.unidentifiable=T, ordinal=hh_ordinals)
  }, hh_bns.puma, hh.train)

options(warn = default_warn_option)

t1 <- Sys.time()
cat(paste0(round(as.numeric(t1-start_time), 2), " ", units(t1-start_time)), sep="\n")


##### HH_INDV BNs: LEARN STRUCTURES AND PARAMETERS #####
cat("Learning individual BNs.....")

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

## suppress numerous bnlearn warnings
default_warn_option <- getOption("warn")
options(warn = -1)

## learn structure of each PUMA BN
hh_indv_bns.puma <- mapply(function(data, prohib) {
  hc(data, score=bn_score, blacklist=prohib)
}, hh_indv.train, hh_indv_prohibits, SIMPLIFY=F) %>%
  matrix(nrow=n_types, ncol=n_pumas)

## find arc strengths
hh_indv_arc_strengths.puma <- mapply(function(bn,data) arc.strength(bn, data, criterion=bn_score), t(hh_indv_bns.puma), t(hh_indv.train), SIMPLIFY=F)

# ## plot each PUMA BN
# par(mfrow=c(n_types,n_pumas), mar=c(2,3,2,1))
# invisible(mapply(function(bn,data,arc_s,n,struc_name){
#   strength.plot(bn,arc_s,main=paste0(struc_name, " (AIC: ", round(AIC(bn, data), digits=2), ") (n=", n, ")"), shape="ellipse")
# }, t(hh_indv_bns.puma), t(hh_indv.train), hh_indv_arc_strengths.puma, sapply(t(hh_indv.train), nrow), rep(struc_names, each=n_pumas)))

## fit parameters of each BN
hh_indv_ordinals <- c("hhinc", "nwork", "ncars", "nchil", "hhsiz", "i_age", "i_inc")
hh_indv_bns.puma.fit <- mapply(function(bn, data){
  bn.fit(bn, data, method="mle", replace.unidentifiable=T, ordinal=hh_indv_ordinals)
}, hh_indv_bns.puma, hh_indv.train)

## train the GQ BN; prohibit edges between HH-level variables, which are useless
gq_data <- pums_hh_indv.type.puma[[5]] %>%
  select(-c(HHID, INDVID, hhwt, hh_type, puma)) %>%
  slice(rep(1:n(), wt)) %>%
  select(-c(wt, tenur, hhinc, nwork, ncars, nchil, hhsiz))
gq_bn <- hc(gq_data, score=bn_score, blacklist=hh_indv_prohibits[[1]] %>% 
              filter(!(from %in% c("tenur", "nwork", "ncars", "nchil", "hhsiz", "hhinc"))) %>%
              filter(!(to %in% c("tenur", "nwork", "ncars", "nchil", "hhsiz", "hhinc"))))
gq_bn_arc_strengths <- arc.strength(gq_bn, gq_data, criterion=bn_score)
# strength.plot(gq_bn,gq_bn_arc_strengths, main=paste0("Group Quarters (AIC: ", round(AIC(gq_bn, gq_data), digits=2), ")"), shape="ellipse")
gq_bn.fit <- bn.fit(gq_bn, gq_data, method="mle", replace.unidentifiable=T, ordinal=hh_indv_ordinals)

options(warn = default_warn_option)

end_time <- Sys.time()
cat(paste0(round(as.numeric(end_time-t1), 2), " ", units(end_time-t1)), sep="\n")
