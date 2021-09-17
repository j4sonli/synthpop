library(tidycensus)

## read and filter PUMS data
census_api_key(my_census_api_key, overwrite=T, install=T)
readRenviron("~/.Renviron")

pums_raw <- get_pums(
  variables = c("ST", "NP", "HINCP", "PINCP", "AGEP", "RAC1P",
                "ESR", "SEX", "SCHL", "TEN", "HISP", "VEH", "NOC", "HHT",
                "WIF", "INDP"),
  state = "multiple",
  puma = split(pumas, substr(pumas, 1, 2)) %>%
    lapply(function(st_pumas) paste(st_pumas, collapse=",7950000US") %>% substring(3)) %>% 
    unlist(),
  survey = "acs5",
  year = 2019,
  key = CENSUS_API_KEY,
  show_call = T
) %>% 
  mutate(puma = paste0(paste0(str_pad(ST, 2, pad="0"), 
                              str_pad(PUMA, 5, pad="0")))) %>%
  # filter(TEN != "b") %>%  ## filter out GQs
  select(-c(ST, PUMA, SERIALNO, SPORDER)) %>%
  data.frame()

## compute HHID sequence from household sizes
hhids <- integer(nrow(pums_raw))
row_i <- 1
curr_hhid <- 1
while(row_i < nrow(pums_raw)+1){
  hhids[row_i:(row_i+as.integer(pums_raw$NP[row_i])-1)] <- rep.int(curr_hhid, pums_raw$NP[row_i])
  row_i <- row_i + as.integer(pums_raw$NP[row_i])
  curr_hhid <- curr_hhid + 1
}

pums <- pums_raw %>%
  rename(hhwt=WGTP,
         wt=PWGTP,
         hhsiz=NP,
         hhinc=HINCP,
         i_inc=PINCP,
         i_age=AGEP,
         race_=RAC1P,
         emply=ESR,
         i_sex=SEX,
         edu_q=SCHL,
         tenur=TEN,
         i_eth=HISP,
         ncars=VEH,
         nchil=NOC,
         hh_type=HHT,
         nwork=WIF,
         indus=INDP) %>%
  mutate(HHID=hhids,
         INDVID=seq_len(nrow(.)))

## recode missing data
pums$tenur[pums$tenur == "b"] <- "GQ"
pums$hhinc[pums$hhinc == -60000] <- NA  ## GQ
pums$nwork[pums$nwork == "b"] <- NA  ## GQ, non-family households
pums$ncars[pums$ncars == "b"] <- 0  ## GQ
pums$nchil[pums$nchil == -1] <- 0  ## GQ
pums$hh_type[pums$hh_type == "b"] <- "GQ"
pums$edu_q[pums$edu_q == "bb"] <- "< 3 yrs. old"
pums$emply[pums$emply == "b"] <- "< 16 yrs. old"
pums$indus[pums$indus == 9920] <- "< 16 yrs. old / NILF"  ## < 16 yrs. old / NILF who last worked > 5 yrs. ago or never worked
pums$i_inc[pums$i_inc == -19999] <- 0  ## < 15 yrs. old

pums$hhwt[pums$hhwt == 0] <- 1 ## GQ

## recode all variables
pums_hh_indv <- pums %>% 
  mutate(puma = factor(puma),
         tenur = recode(tenur, 
                        "1"="Owned with mortgage or loan",
                        "2"="Owned free and clear",
                        "3"="Rented",
                        "4"="Rented") %>%
           factor(),
         hhinc = cut(hhinc, breaks=c(-Inf, 0, 25000, 50000, 75000, 100000, 150000, Inf), 
                     labels=c("<=0K", "1-25K", "26-50K", "51-75K", "76-100K", "101-150K", ">150K")) %>%
           ordered(),
         nwork = cut(as.numeric(nwork), breaks=c(-Inf, 0, 1, 2, Inf), labels=c("0", "1", "2", "3+")) %>%
           ordered(),
         ncars = cut(as.numeric(ncars), breaks=c(-Inf, 0, 1, Inf), labels=c("0", "1", "2+")) %>%
           ordered(),
         nchil = cut(nchil, breaks=c(-Inf, 0, 1, Inf), labels=c("0", "1", "2+")) %>%
           ordered(),
         hhsiz = cut(hhsiz, breaks=c(-Inf, 1, 2, 3, 4, 5, Inf), labels=c("1", "2", "3", "4", "5", "6+")) %>%
           ordered(),
         i_age = cut(i_age, breaks=c(-Inf, 17, 35, 50, 70, Inf), 
                     labels=c("<18", "18-35", "36-50", "51-70", ">70")) %>%
           ordered(),
         i_sex = recode(i_sex, 
                        "1"="Male", 
                        "2"="Female") %>%
           factor(),
         race_ = recode(race_,
                        "1"="White",
                        "2"="Black",
                        "3"="American Indian, Alaskan, or Pacific Islander",
                        "4"="American Indian, Alaskan, or Pacific Islander",
                        "5"="American Indian, Alaskan, or Pacific Islander",
                        "6"="Asian",
                        "7"="American Indian, Alaskan, or Pacific Islander",
                        "8"="Other",
                        "9"="Two or more races") %>%
           factor(),
         i_eth = cut(as.integer(i_eth), breaks=c(-Inf, 1, Inf), labels=c("No Hispanic origin", "Hispanic origin")) %>%
           factor(),
         i_inc = cut(i_inc, breaks=c(-Inf, 0, 25000, 50000, 75000, 100000, Inf), 
                     labels=c("<=0K", "1-25K", "26-50K", "51-75K", "76-100K", ">100K")) %>%
           ordered(),
         emply = recode(emply, 
                        "1"="Civilian employed",
                        "2"="Civilian employed",
                        "3"="Unemployed",
                        "4"="Armed forces employed",
                        "5"="Armed forces employed",
                        "6"="Not in labor force") %>%
           factor(),
         indus = ifelse(indus == "< 16 yrs. old / NILF", "< 16 yrs. old / NILF", substr(indus, 1, 1)) %>%
           factor(),
         edu_q = recode(edu_q, 
                        "2"="1", ## no schooling or nursery/preschool
                        "4"="3", "5"="3", "6"="3", "7"="3", ## incomplete grade school
                        "9"="8", "10"="8", ## completed grade school
                        "12"="11", "13"="11", "14"="11", "15"="11", ## completed middle school
                        "17"="16", ## high school graduate or GED
                        "19"="18", ## some college, no degree
                        ) %>% 
           recode("1"="No schooling or nursery/preschool",
                  "3"="Some grade school",
                  "8"="Grade school",
                  "11"="Middle school",
                  "16"="High school or GED",
                  "18"="Some college",
                  "20"="Associate's",
                  "21"="Bachelor's",
                  "22"="Master's",
                  "23"="Professional degree beyond bachelor's",
                  "24"="Doctorate") %>%
           factor(),
         hh_type = recode(hh_type,
                          "1"="MC",
                          "2"="NS",
                          "3"="NS",
                          "4"="SM",
                          "5"="NF",
                          "6"="SM",
                          "7"="NF")) %>%
  select(HHID, INDVID,
         puma, tenur, hhinc, nwork, ncars, nchil, hhsiz,
         i_age, i_sex, race_, i_eth, i_inc, emply, indus, edu_q, 
         hh_type, hhwt, wt)

## split into household types
n_types <- 5
pums_hh_indv.type <- list(MC=filter(pums_hh_indv, hh_type == "MC"),
                          NS=filter(pums_hh_indv, hh_type == "NS"),
                          SM=filter(pums_hh_indv, hh_type == "SM"),
                          NF=filter(pums_hh_indv, hh_type == "NF"),
                          GQ=filter(pums_hh_indv, hh_type == "GQ"))

## for SM, NF, and GQ types, set nwork equal to the number of nonzero incomes in the household
pums_hh_indv.type$SM$nwork <- as.integer(pums_hh_indv.type$SM$i_inc != "<=0K") %>%
  ordered(levels = levels(pums_hh_indv.type$MC$nwork))
pums_hh_indv.type$NF$nwork <- pums_hh_indv.type$NF %>% 
  group_by(HHID) %>% 
  summarise(nonzero=sum(i_inc != "<=0K")) %>% 
  pull(nonzero) %>%
  rep(times=pums_hh_indv.type$NF %>% 
            count(HHID) %>% 
            pull(n)) %>%
  cut(breaks=c(-Inf, 0, 1, 2, Inf), labels=c("0", "1", "2", "3+")) %>%
  ordered(levels = levels(pums_hh_indv.type$MC$nwork))
# if("GQ" %in% names(pums_hh_indv.type)){
#   pums_hh_indv.type$GQ$nwork <- pums_hh_indv.type$GQ %>% 
#     group_by(HHID) %>% 
#     summarise(nonzero=sum(i_inc != "<=0K")) %>% 
#     pull(nonzero) %>%
#     rep(times=pums_hh_indv.type$GQ %>% 
#           count(HHID) %>% 
#           pull(n)) %>%
#     cut(breaks=c(-Inf, 0, 1, 2, Inf), labels=c("0", "1", "2", "3+")) %>%
#     ordered(levels = levels(pums_hh_indv.type$MC$nwork))
#   
#   ## for GQ type, set hhinc equal to i_inc since GQ hhsiz is 1
#   pums_hh_indv.type$GQ$hhinc <- pums_hh_indv.type$GQ$i_inc %>%
#     recode(">100K"="101-150K") %>%
#     ordered(levels=levels(pums_hh_indv.type$MC$hhinc))
# }

pums_hh_indv.type$GQ$nwork <- "GQ"
pums_hh_indv.type$GQ$hhinc <- "GQ"

## split by PUMA
pums_hh_indv.type.puma <- lapply(pums_hh_indv.type, function(dataset) {
  dataset %>%
    split(dataset$puma)
}) %>% 
  unlist(recursive=F) %>%
  matrix(nrow=n_pumas, dimnames=list(pumas, names(pums_hh_indv.type)))

## create the household-level data
pums_hh.type.puma <- lapply(pums_hh_indv.type.puma, function(dataset) {
  dataset %>%
    select(HHID, puma, tenur, hhinc, nwork, ncars, nchil, hhsiz, hh_type, hhwt) %>%
    group_by(HHID) %>%
    slice(1) %>%
    data.frame()
}) %>% 
  matrix(nrow=n_pumas, dimnames=list(pumas, names(pums_hh_indv.type)))
  