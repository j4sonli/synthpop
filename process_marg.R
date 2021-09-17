library(acs)

api.key.install(my_census_api_key)

geo.sets.puma <- split(pumas, substr(pumas, 1, 2)) %>% 
  lapply(function(st_pumas) substr(st_pumas, 3, 7)) %>% 
  mapply(function(st, puma_codes) 
    geo.make(state=as.numeric(st), 
             puma=puma_codes), 
    names(.), .)
geo.sets.tract <- split(tracts, substr(tracts, 1, 2)) %>%
  lapply(function(st_tracts) substr(st_tracts, 3, 11)) %>% 
  mapply(function(st, cty_tract) 
    geo.make(state=as.numeric(st), 
             county=as.numeric(substr(cty_tract, 1, 3)), 
             tract=substr(cty_tract, 4, 9)), 
    names(.), .)
geo.sets.blkgp <- split(blkgps, substr(blkgps, 1, 2)) %>%
  lapply(function(st_blkgps) substr(st_blkgps, 3, 12)) %>% 
  mapply(function(st, cty_tract_blkgp) 
    geo.make(state=as.numeric(st), 
             county=as.numeric(substr(cty_tract_blkgp, 1, 3)), 
             tract=substr(cty_tract_blkgp, 4, 9),
             block.group=substr(cty_tract_blkgp, 10, 10)), 
    names(.), .)

fetch_from_acs <- function(table_id, geo.sets, geoids) {
  lapply(geo.sets, function(geo.set) 
    acs.fetch(2019, span=5, geo.set, table_id, col.names="pretty")@estimate %>%
      data.frame()
  ) %>%
    bind_rows() %>%
    `rownames<-`(NULL) %>%
    mutate(geoid=geoids)
}

if(!use_prev_marg) {

##### tract x tenur x hhinc #####
tract_tenur_hhinc <- fetch_from_acs("B25118", geo.sets.tract, tracts) %>%
  mutate("Owned.<25K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied..Less.than..5.000 +
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...5.000.to..9.999 +
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...10.000.to..14.999 + 
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...15.000.to..19.999 +
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...20.000.to..24.999,
         "Owned.26-50K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...25.000.to..34.999 +
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...35.000.to..49.999,
         "Owned.51-75K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...50.000.to..74.999,
         "Owned.76-100K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...75.000.to..99.999,
         "Owned.101-150K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...100.000.to..149.999,
         "Owned.>150K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Owner.occupied...150.000.or.more,
         "Rented.<25K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied..Less.than..5.000 +    
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...5.000.to..9.999 +    
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...10.000.to..14.999 +  
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...15.000.to..19.999 +  
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...20.000.to..24.999,  
         "Rented.26-50K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...25.000.to..34.999 +  
           B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...35.000.to..49.999,  
         "Rented.51-75K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...50.000.to..74.999,  
         "Rented.76-100K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...75.000.to..99.999,  
         "Rented.101-150K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...100.000.to..149.999,
         "Rented.>150K"=B25118..Tenure.by.Household.Income.in.the.Past.12.Months..in.2015.Inflation.Adjusted.Dollars...Renter.occupied...150.000.or.more) %>%
  select(-c(1:25))
write.csv(tract_tenur_hhinc, file=paste0(marg_dir, "tract_tenur_hhinc.csv"), row.names=F)

##### blkgp x tenur x hhsiz #####
blkgp_tenur_hhsiz <- fetch_from_acs("B25009", geo.sets.blkgp, blkgps) %>% 
  mutate("Owned.1"=TENURE.BY.HOUSEHOLD.SIZE..Owner.occupied..1.person.household,
         "Owned.2"=TENURE.BY.HOUSEHOLD.SIZE..Owner.occupied..2.person.household,
         "Owned.3"=TENURE.BY.HOUSEHOLD.SIZE..Owner.occupied..3.person.household,
         "Owned.4"=TENURE.BY.HOUSEHOLD.SIZE..Owner.occupied..4.person.household,
         "Owned.5+"=TENURE.BY.HOUSEHOLD.SIZE..Owner.occupied..5.person.household + 
           TENURE.BY.HOUSEHOLD.SIZE..Owner.occupied..6.person.household +
           TENURE.BY.HOUSEHOLD.SIZE..Owner.occupied..7.or.more.person.household,
         "Rented.1"=TENURE.BY.HOUSEHOLD.SIZE..Renter.occupied..1.person.household,
         "Rented.2"=TENURE.BY.HOUSEHOLD.SIZE..Renter.occupied..2.person.household,
         "Rented.3"=TENURE.BY.HOUSEHOLD.SIZE..Renter.occupied..3.person.household,
         "Rented.4"=TENURE.BY.HOUSEHOLD.SIZE..Renter.occupied..4.person.household,
         "Rented.5+"=TENURE.BY.HOUSEHOLD.SIZE..Renter.occupied..5.person.household + 
           TENURE.BY.HOUSEHOLD.SIZE..Renter.occupied..6.person.household +
           TENURE.BY.HOUSEHOLD.SIZE..Renter.occupied..7.or.more.person.household) %>%
  select(-c(1:17))
write.csv(blkgp_tenur_hhsiz, file=paste0(marg_dir, "blkgp_tenur_hhsiz.csv"), row.names=F)

##### blkgp x tenur x hhsiz #####
blkgp_tenur <- fetch_from_acs("B25009", geo.sets.blkgp, blkgps) %>% 
  mutate("Owned"=TENURE.BY.HOUSEHOLD.SIZE..Owner.occupied.,
         "Rented"=TENURE.BY.HOUSEHOLD.SIZE..Renter.occupied.) %>%
  select(-c(1:17))
write.csv(blkgp_tenur, file=paste0(marg_dir, "blkgp_tenur.csv"), row.names=F)

##### puma x hhtype #####
hhtype_marg <- fetch_from_acs("B11001", geo.sets.puma, pumas)
puma_hhtype <- hhtype_marg %>%
  mutate(MC=Household.Type..including.Living.Alone...Family.households..Married.couple.family,
         NS=Household.Type..including.Living.Alone...Family.households..Other.family..Male.householder..no.wife.present +
           Household.Type..including.Living.Alone...Family.households..Other.family..Female.householder..no.husband.present,
         SM=Household.Type..including.Living.Alone...Nonfamily.households..Householder.living.alone,
         NF=Household.Type..including.Living.Alone...Nonfamily.households..Householder.not.living.alone) %>%
  select(-c(1:90)) %>%
  rename(puma_id=geoid)
write.csv(hhtype_marg.puma, file=paste0(marg_dir, "hhtype_marg.csv"), row.names=F)
write.csv(puma_hhtype, file=paste0(marg_dir, "puma_hhtype.csv"), row.names=F)

##### tract x hhtype #####
tract_hhtype <- fetch_from_acs("B11001", geo.sets.tract, tracts) %>%
  mutate(MC=Household.Type..including.Living.Alone...Family.households..Married.couple.family,
         NS=Household.Type..including.Living.Alone...Family.households..Other.family..Male.householder..no.wife.present +
           Household.Type..including.Living.Alone...Family.households..Other.family..Female.householder..no.husband.present,
         SM=Household.Type..including.Living.Alone...Nonfamily.households..Householder.living.alone,
         NF=Household.Type..including.Living.Alone...Nonfamily.households..Householder.not.living.alone) %>%
  select(-c(1:90))
write.csv(tract_hhtype, file=paste0(marg_dir, "tract_hhtype.csv"), row.names=F)

##### tract x nwork #####
tract_nwork <- fetch_from_acs("B08203", geo.sets.tract, tracts) %>%
  mutate("0"=NUMBER.OF.WORKERS.IN.HOUSEHOLD.BY.VEHICLES.AVAILABLE..No.workers.,
         "1"=NUMBER.OF.WORKERS.IN.HOUSEHOLD.BY.VEHICLES.AVAILABLE..1.worker.,
         "2"=NUMBER.OF.WORKERS.IN.HOUSEHOLD.BY.VEHICLES.AVAILABLE..2.workers.,
         "3+"=NUMBER.OF.WORKERS.IN.HOUSEHOLD.BY.VEHICLES.AVAILABLE..3.or.more.workers.) %>%
  select(-c(1:30))
write.csv(tract_nwork, file=paste0(marg_dir, "tract_nwork.csv"), row.names=F)

} else {
  tract_tenur_hhinc <- read.csv(paste0(marg_dir, "tract_tenur_hhinc.csv"), check.names=F)
  blkgp_tenur_hhsiz <- read.csv(paste0(marg_dir, "blkgp_tenur_hhsiz.csv"), check.names=F)
  hhtype_marg <- read.csv(paste0(marg_dir, "hhtype_marg.csv"), check.names=F)
  puma_hhtype <- read.csv(paste0(marg_dir, "puma_hhtype.csv"), check.names=F)
  tract_hhtype <- read.csv(paste0(marg_dir, "tract_hhtype.csv"), check.names=F)
  tract_nwork <- read.csv(paste0(marg_dir, "tract_nwork.csv"), check.names=F)
}

##### households by PUMA #####
pop_by_puma <- hhtype_marg %>%
  mutate(puma_pop = Household.Type..including.Living.Alone...Total.) %>%
  select(-c(1:90)) %>%
  rename(puma=geoid)

