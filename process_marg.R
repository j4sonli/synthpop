if(!use_prev_marg) {
  
  fetch_from_acs <- function(table_id, geo, geoids) {
    get_acs(geography=geo, state=substr(geoids[1], 1, 2), table=table_id, output="wide") %>%
      filter(GEOID %in% geoids) %>%
      rename(geoid=GEOID) %>%
      arrange(geoid) %>%
      data.frame()
  }

  if(!dir.exists(marg_dir)) {
    cat(paste0("Creating new directory at ", marg_dir), sep="\n")
    dir.create(marg_dir, recursive = T)
  }
  
  ########## UTILITY ##########
  
  ##### tract x gq_pop #####
  tract_gq_pop <- fetch_from_acs("B26001", "tract", tracts) %>%
    mutate("gq_pop" = B26001_001E) %>%
    select(-c(2:4))
  write.csv(tract_gq_pop, file=paste0(marg_dir, "tract_gq_pop.csv"), row.names=F)
  
  blkgp_pop <- fetch_from_acs("B01003", "block group", blkgps) %>%
    mutate(pop = B01003_001E) %>%
    select(-c(2:4))
  
  
  ########## CALIBRATION ##########
  ## add GQ category to all household-level calibration and validation marginals
  
  ##### tract x tenur x hhinc #####
  tract_tenur_hhinc <- fetch_from_acs("B25118", "tract", tracts) %>%
    mutate("Owned.<25K" = B25118_003E + B25118_004E + B25118_005E + B25118_006E + B25118_007E,
           "Owned.26-50K" = B25118_008E + B25118_009E,
           "Owned.51-75K" = B25118_010E,
           "Owned.76-100K" = B25118_011E,
           "Owned.101-150K" = B25118_012E,
           "Owned.>150K" = B25118_013E,
           "Rented.<25K" = B25118_015E + B25118_016E + B25118_017E+ B25118_018E + B25118_019E,  
           "Rented.26-50K" = B25118_020E + B25118_021E,  
           "Rented.51-75K" = B25118_022E,  
           "Rented.76-100K" = B25118_023E,  
           "Rented.101-150K" = B25118_024E,
           "Rented.>150K" = B25118_025E) %>%
    select(-c(2:52)) %>%
    mutate("GQ.GQ" = tract_gq_pop$gq_pop)
  write.csv(tract_tenur_hhinc, file=paste0(marg_dir, "tract_tenur_hhinc.csv"), row.names=F)
  
  ##### blkgp x tenur x hhsiz #####
  blkgp_tenur_hhsiz <- fetch_from_acs("B25009", "block group", blkgps) %>% 
    mutate("Owned.1" = B25009_003E,
           "Owned.2" = B25009_004E,
           "Owned.3" = B25009_005E,
           "Owned.4" = B25009_006E,
           "Owned.5" = B25009_007E,
           "Owned.6+" = B25009_008E + B25009_009E,
           "Rented.1" = B25009_011E,
           "Rented.2" = B25009_012E,
           "Rented.3" = B25009_013E,
           "Rented.4" = B25009_014E,
           "Rented.5" = B25009_015E,
           "Rented.6+" = B25009_016E + B25009_017E) %>%
    select(-c(2:36)) %>%
    mutate(tract = substr(geoid, 1, 11)) %>% 
    left_join(rename(tract_gq_pop, tract=geoid), by='tract') %>%
    group_by(tract) %>% 
    mutate(gq_pop = trs(gq_pop/n())) %>%
    ungroup() %>%
    select(-tract) %>%
    rename(GQ.GQ = gq_pop) %>%
    data.frame() %>%
    rename("Owned.6+"="Owned.6.", "Rented.6+"="Rented.6.")
  write.csv(blkgp_tenur_hhsiz, file=paste0(marg_dir, "blkgp_tenur_hhsiz.csv"), row.names=F)
  
  ##### tract/puma x hhtype #####
  tract_hhtype <- fetch_from_acs("B11001", "tract", tracts) %>%
    mutate(MC = B11001_003E,
           NS = B11001_005E + B11001_006E,
           SM = B11001_008E,
           NF = B11001_009E) %>%
    mutate("GQ" = tract_gq_pop$gq_pop) %>%
    select(-c(2:20))
  write.csv(tract_hhtype, file=paste0(marg_dir, "tract_hhtype.csv"), row.names=F)
  puma_hhtype <-  data.frame(geoid = pumas,
                             MC = sum(tract_hhtype$MC),
                             NS = sum(tract_hhtype$NS),
                             SM = sum(tract_hhtype$SM),
                             NF = sum(tract_hhtype$NF),
                             GQ = sum(tract_hhtype$GQ))
  write.csv(puma_hhtype, file=paste0(marg_dir, "puma_hhtype.csv"), row.names=F)
  
  ##### tract x nwork #####
  tract_nwork <- fetch_from_acs("B08203", "tract", tracts) %>%
    mutate("0" = B08203_007E,
           "1" = B08203_013E,
           "2" = B08203_019E,
           "3+" = B08203_025E) %>%
    select(-c(2:62)) %>%
    mutate("GQ" = tract_gq_pop$gq_pop)
  write.csv(tract_nwork, file=paste0(marg_dir, "tract_nwork.csv"), row.names=F)
  
  ##### tract x i_sex x i_age #####
  tract_i_sex_i_age.raw <- fetch_from_acs("B01001", "tract", tracts)
  tract_i_sex_i_age <- tract_i_sex_i_age.raw %>%
    mutate("Male.<18" = B01001_003E + B01001_004E + B01001_005E + B01001_006E,
           "Male.18-35" = B01001_007E + B01001_008E + B01001_009E + B01001_010E + B01001_011E + B01001_012E,
           "Male.36-50" = B01001_013E + B01001_014E + B01001_015E,
           "Male.51-70" = B01001_016E + B01001_017E + B01001_018E + B01001_019E + B01001_020E + B01001_021E,
           "Male.>70" = B01001_022E + B01001_023E + B01001_024E + B01001_025E,
           "Female.<18" = B01001_027E + B01001_028E + B01001_029E + B01001_030E,
           "Female.18-35" = B01001_031E + B01001_032E + B01001_033E + B01001_034E + B01001_035E + B01001_036E,
           "Female.36-50" = B01001_037E + B01001_038E + B01001_039E,
           "Female.51-70" = B01001_040E + B01001_041E + B01001_042E + B01001_043E + B01001_044E +B01001_045E,
           "Female.>70" = B01001_046E + B01001_047E + B01001_048E + B01001_049E) %>%
    select(-c(2:100))
  write.csv(tract_i_sex_i_age, file=paste0(marg_dir, "tract_i_sex_i_age.csv"), row.names=F)
  
  ## approximate split into block groups, for estimation of blkgp_emply additions
  blkgp_i_sex_i_age <- 
  
  ##### blkgp x emply #####
  blkgp_emply <- fetch_from_acs("B23025", "block group", blkgps) %>%
    mutate("Employed" = B23025_004E + B23025_006E,
           "Unemployed" = B23025_005E,
           "Not in labor force" = B23025_007E) %>%
    select(-c(2:16))
  # need to add category for individuals less than 16 years old
  # the count of which we determine by (actual.blkgp_pop - blkgp_emply.blkgp_pop)
  # i.e., completely independent of tract_i_sex_i_age
  blkgp_emply <- blkgp_emply %>%
    mutate(`< 16 yrs. old` = blkgp_pop$pop - (Employed + Unemployed + `Not in labor force`))
  write.csv(blkgp_emply, file=paste0(marg_dir, "blkgp_emply.csv"), row.names=F)
  
  
  ########## VALIDATION ##########
  
  ##### tract x nwork x ncars #####
  tract_nwork_ncars <- fetch_from_acs("B08203", "tract", tracts) %>%
    mutate("0.0" = B08203_008E,
           "0.1" = B08203_009E,
           "0.2+" = B08203_010E + B08203_011E + B08203_012E,
           "1.0" = B08203_014E,
           "1.1" = B08203_015E,
           "1.2+" = B08203_016E + B08203_017E + B08203_018E,
           "2.0" = B08203_020E,
           "2.1" = B08203_021E,
           "2.2+" = B08203_022E + B08203_023E + B08203_024E,
           "3+.0" = B08203_026E,
           "3+.1" = B08203_027E,
           "3+.2+" = B08203_028E + B08203_029E + B08203_030E) %>%
    select(-c(2:62)) %>%
    mutate("GQ.GQ" = tract_gq_pop$gq_pop)
  write.csv(tract_nwork_ncars, file=paste0(marg_dir, "tract_nwork_ncars.csv"), row.names=F)

  ##### blkgp x i_eth x race_ #####
  blkgp_i_eth_race_ <- fetch_from_acs("B03002", "block group", blkgps) %>%
    mutate("No Hispanic origin.White" = B03002_003E,
           "No Hispanic origin.Black" = B03002_004E,
           "No Hispanic origin.American Indian, Alaskan, or Pacific Islander" = B03002_005E + B03002_007E,
           "No Hispanic origin.Asian" = B03002_006E,
           "No Hispanic origin.Other" = B03002_008E,
           "No Hispanic origin.Two or more races" = B03002_009E,
           "Hispanic origin.White" = B03002_013E,
           "Hispanic origin.Black" = B03002_014E,
           "Hispanic origin.American Indian, Alaskan, or Pacific Islander" = B03002_015E + B03002_017E,
           "Hispanic origin.Asian" = B03002_016E,
           "Hispanic origin.Other" = B03002_018E,
           "Hispanic origin.Two or more races" = B03002_019E) %>%
    select(-c(2:44))
  write.csv(blkgp_i_eth_race_, file=paste0(marg_dir, "blkgp_i_eth_race_.csv"), row.names=F)
  
} else {
  tract_gq_pop <- read.csv(paste0(marg_dir, "tract_gq_pop.csv"), check.names=F)
  
  tract_tenur_hhinc <- read.csv(paste0(marg_dir, "tract_tenur_hhinc.csv"), check.names=F)
  blkgp_tenur_hhsiz <- read.csv(paste0(marg_dir, "blkgp_tenur_hhsiz.csv"), check.names=F)
  tract_hhtype <- read.csv(paste0(marg_dir, "tract_hhtype.csv"), check.names=F)
  puma_hhtype <- read.csv(paste0(marg_dir, "puma_hhtype.csv"), check.names=F)
  tract_nwork <- read.csv(paste0(marg_dir, "tract_nwork.csv"), check.names=F)
  
  tract_nwork_ncars <- read.csv(paste0(marg_dir, "tract_nwork_ncars.csv"), check.names=F)
  tract_i_sex_i_age <- read.csv(paste0(marg_dir, "tract_i_sex_i_age.csv"), check.names=F)
  blkgp_i_eth_race_ <- read.csv(paste0(marg_dir, "blkgp_i_eth_race_.csv"), check.names=F)
  blkgp_emply <- read.csv(paste0(marg_dir, "blkgp_emply.csv"), check.names=F)
}
