set.seed(1033)

start_time <- Sys.time()


##### PREPARE DATA #####
cat(paste0("Preparing PUMS and ACS data for PUMA ", pumas, "....."), sep="\n")

n_pumas <- 1

puma_tract_equiv <- read.csv(puma_tract_equiv_file)
puma_tract_join <- puma_tract_equiv %>%
  mutate(puma = paste0(str_pad(as.character(STATEFP), 2, pad="0"),
                       str_pad(as.character(PUMA5CE), 5, pad="0"))) %>%
  filter(puma %in% pumas) %>%
  mutate(tract = paste0(str_pad(as.character(STATEFP), 2, pad="0"),
                        str_pad(as.character(COUNTYFP), 3, pad="0"),
                        str_pad(as.character(TRACTCE), 6, pad="0")))

tracts <- puma_tract_join$tract

tract_block_join <- read.csv(tract_block_equiv_file) %>%
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
  if(length(w)==1 || all(w == floor(w))){
    return(w)
  }
  w.int <- floor(w) 
  rems <- w - w.int 
  topup <- sample(length(w), size = round(sum(rems)), prob = rems)
  w.int[topup] <- w.int[topup] + 1
  return(w.int)
}

## get the microdata
source(process_pums_file)

## get marginal data
source(process_marg_file)

end_time <- Sys.time()
cat(paste0(".....", round(as.numeric(end_time-start_time), 2), " ", units(end_time-start_time)), sep="\n")
