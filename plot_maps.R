library(R.utils)
library(dplyr)
library(tidyr)
library(stringr)
library(tidycensus)
library(ggplot2)
library(gridExtra)
library(tidycensus)
library(tmap)
library(sf)
library(viridis)

options(dplyr.summarise.inform=F)

plot_study_areas <- function(map_pumas) {
  our_msas <- c("Los Angeles-Long Beach-Santa Ana, CA",
                "Phoenix-Mesa-Glendale, AZ",
                "Houston-Sugar Land-Baytown, TX",
                "Jacksonville, FL")
  us_msas <- st_read("map_data/tl_2010_us_cbsa10/tl_2010_us_cbsa10.shp") %>%
    filter(NAME10 %in% our_msas)
  
  puma_sfs <- rbind(
    st_read("map_data/tl_2019_04_puma10/tl_2019_04_puma10.shp"),
    st_read("map_data/tl_2019_06_puma10/tl_2019_06_puma10.shp"),
    st_read("map_data/tl_2019_12_puma10/tl_2019_12_puma10.shp"),
    st_read("map_data/tl_2019_25_puma10/tl_2019_25_puma10.shp"),
    st_read("map_data/tl_2019_48_puma10/tl_2019_48_puma10.shp"),
    st_read("map_data/tl_2019_56_puma10/tl_2019_56_puma10.shp")
  ) %>%
    # filter out islands in LA MSA for better visualization
    filter(GEOID10 != "0603768")
  
  study_area_pumas <- st_intersection(puma_sfs, us_msas) %>%
    select(STATEFP10, PUMACE10, GEOID10, NAME10, geometry)
  
  # add back the PUMAs in MAPC, MA
  mapc_pumas <- sort(c(2500503,2500504,2500505,2500506,2500507,2500508,2500703,2500704,
                       2501000,2501400,2502400,2502800,2503301,2503302,2503303,2503304,
                       2503305,2503306,2503400,2503500,2503601,2503602,2503603,2503900,
                       2504903))
  study_area_pumas <- rbind(study_area_pumas, puma_sfs %>% 
                              filter(GEOID10 %in% as.character(mapc_pumas)) %>%
                              select(STATEFP10, PUMACE10, GEOID10, geometry) %>%
                              mutate(NAME10 = "Metro Boston, MA *"))
  # add back the PUMAs in WY
  wy_pumas <- c(5600100, 5600200, 5600300, 5600400, 5600500)
  study_area_pumas <- rbind(study_area_pumas, puma_sfs %>% 
                              filter(GEOID10 %in% as.character(wy_pumas)) %>%
                              select(STATEFP10, PUMACE10, GEOID10, geometry) %>%
                              mutate(NAME10 = "State of WY"))
  
  # highlight the low, median, and high densities
  study_area_pumas$density <- NA
  density_vals <- c("Lowest", "Median", "Highest")
  density_counter <- 1
  for(map_puma in map_pumas[1:18]) {
    if(density_counter > 3) {
      density_counter <- 1
    }
    study_area_pumas[study_area_pumas$GEOID10 == map_puma, "density"] = density_vals[density_counter]
    density_counter <- density_counter + 1
  } 
  study_area_pumas$density <- factor(study_area_pumas$density, levels=rev(density_vals))
  study_area_pumas <- study_area_pumas %>%
    mutate(NAME10 = factor(NAME10, levels=c("Los Angeles-Long Beach-Santa Ana, CA",
                                            "Houston-Sugar Land-Baytown, TX",
                                            "Phoenix-Mesa-Glendale, AZ",
                                            "Metro Boston, MA *",
                                            "Jacksonville, FL",
                                            "State of WY")))
  
  dev.new(width=7.2, height=9, noRStudioGD = T)
  
  tm_shape(study_area_pumas) + 
    tm_polygons(col="density",
                palette = "RdYlBu",
                title="PUMA Density",
                showNA=F,
                colorNA="#f0f0f0",
                legend.is.portrait=F,
                border.alpha=0.5) +
    tm_compass(position=c("right", "top")) +
    tm_scale_bar(breaks=c(0, 100, 200),
                 text.size=0.8) +
    tm_facets(by="NAME10",
              ncol=2) +
    tm_layout(main.title="Study areas",
              main.title.fontface="bold",
              main.title.size=1.4,
              main.title.position="center",
              panel.label.size=1.3,
              legend.outside=T,
              legend.outside.position="bottom",
              legend.outside.size=0.15,
              legend.title.size=1.175,
              legend.title.fontface="bold",
              legend.text.size=1,
              outer.margins=c(0, 0.02, 0, 0.02),
              inner.margins=c(0.06, 0.06, 0.06, 0.06))  
}

error_heatmap <- function(synpop, marginal, dim1, dim2, plot_title, limits=NULL) {
  pop <- synpop %>% 
    mutate(tenur_prox = recode(tenur,
                               "Owned with mortgage or loan"="Owned",
                               "Owned free and clear"="Owned"),
           hhinc_prox = recode(hhinc, 
                               "<=0K"="<25K",
                               "1-25K"="<25K"),
           hhsiz_prox = hhsiz,
           emply_prox = {if(!all(is.na(emply))) recode(emply,
                                                       "Civilian employed"="Employed",
                                                       "Armed forces employed"="Employed") else NULL}) %>%
    mutate(tenur_hhinc_prox = paste0(tenur_prox, ".", hhinc_prox),
           tenur_hhsiz_prox = paste0(tenur_prox, ".", hhsiz_prox),
           nwork_ncars = paste0(nwork, ".", ncars),
           i_sex_i_age = paste0({if("i_sex" %in% names(.)) i_sex else NULL}, 
                                ".", 
                                {if("i_age" %in% names(.)) i_age else NULL})) %>%
    count(!!dim1, !!dim2) %>%
    rename(synpop_ct=n)
  
  marg <- marginal %>%
    group_by(!!dim1, !!dim2) %>%
    summarise(marg_ct = sum(marg_ct)) %>%
    data.frame()
  marg[[1]] = as.character(marg[[1]])
  marg[[2]] = as.character(marg[[2]])
  
  heatmap_data <- left_join(pop, marg, by=c(quo_name(dim1), quo_name(dim2))) %>%
    mutate(percent_error = (synpop_ct - marg_ct)/(1+marg_ct) * 100,
           error = synpop_ct - marg_ct)
  
  p <- ggplot(heatmap_data, aes(!!dim1, !!dim2, fill=percent_error)) + 
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red", limits=limits) + 
    geom_tile() +
    ggtitle(plot_title)
  
  return(list(p=p, heatmap_data=heatmap_data))
}

rmse <- function(errors) {
  return(sqrt(sum(errors^2)/length(errors)))
}

get_rmses <- function(heatmap_data, spatial_level) {
  return(heatmap_data %>%
           group_by(geoid) %>%
           mutate(share_error=(synpop_ct/sum(synpop_ct)-marg_ct/sum(marg_ct))*100) %>%
           ungroup() %>%
           pivot_wider(id_cols=geoid, ., names_from=names(.)[1], values_from=share_error) %>%
           # replace NAs that arise because of 0 count in synpop/marg (mostly GQ)
           replace(is.na(.), 0) %>%
           mutate(rmse = apply(., 1, function(x) rmse(as.numeric(x[2:length(x)])))) %>%
           select(geoid, rmse))
}

process_marginals <- function(marg_dir) {
  source("process_marg.R")
  tract_tenur_hhinc.l <<- gather(tract_tenur_hhinc, tenur_hhinc_prox, marg_ct, -geoid) %>%
    rename(tract=geoid) %>%
    mutate(tract = str_pad(tract, 11, "left", "0")) %>%
    mutate(tenur_prox = unlist(strsplit(tenur_hhinc_prox, "\\."))[c(T, F)],
           hhinc_prox = unlist(strsplit(tenur_hhinc_prox, "\\."))[c(F, T)])
  blkgp_tenur_hhsiz.l <<- gather(blkgp_tenur_hhsiz, tenur_hhsiz_prox, marg_ct, -geoid) %>%
    rename(blkgp=geoid) %>%
    mutate(blkgp = str_pad(blkgp, 12, "left", "0")) %>%
    mutate(tract = substr(blkgp, 1, 11),
           tenur_prox = unlist(strsplit(tenur_hhsiz_prox, "\\."))[c(T, F)],
           hhsiz_prox = unlist(strsplit(tenur_hhsiz_prox, "\\."))[c(F, T)])
  tract_nwork.l <<- gather(tract_nwork, nwork, marg_ct, -geoid) %>%
    rename(tract=geoid) %>%
    mutate(tract = str_pad(tract, 11, "left", "0"))
  tract_hhtype.l <<- gather(tract_hhtype, hhtype, marg_ct, -geoid) %>%
    rename(tract=geoid) %>%
    mutate(tract = str_pad(tract, 11, "left", "0"))
  tract_nwork_ncars.l <<- gather(tract_nwork_ncars, nwork_ncars, marg_ct, -geoid) %>%
    rename(tract=geoid) %>%
    mutate(tract = str_pad(tract, 11, "left", "0")) %>%
    mutate(nwork = unlist(strsplit(nwork_ncars, "\\."))[c(T, F)],
           ncars = unlist(strsplit(nwork_ncars, "\\."))[c(F, T)])
  tract_i_sex_i_age.l <<- gather(tract_i_sex_i_age, i_sex_i_age, marg_ct, -geoid) %>%
    rename(tract=geoid) %>%
    mutate(tract = str_pad(tract, 11, "left", "0")) %>%
    mutate(i_sex = unlist(strsplit(i_sex_i_age, "\\."))[c(T, F)],
           i_age = unlist(strsplit(i_sex_i_age, "\\."))[c(F, T)])
  blkgp_emply.l <<- gather(blkgp_emply, emply, marg_ct, -geoid) %>%
    rename(blkgp=geoid) %>%
    mutate(blkgp = str_pad(blkgp, 12, "left", "0")) %>%
    mutate(tract = substr(blkgp, 1, 11),
           emply_prox = recode(emply,
                               "Civilian employed"="Employed",
                               "Armed forces employed"="Employed"))
  tract_i_inc.l <<- tract_i_inc %>%
    mutate(`<=25K`=`age<15`+`<=25K`) %>%
    select(-`age<15`) %>%
    gather(i_inc, marg_ct, -geoid) %>%
    rename(tract=geoid) %>%
    mutate(tract = str_pad(tract, 11, "left", "0"),
           i_inc = recode(i_inc, "age<15"="<=25K")) %>%
    rename(i_inc_prox = i_inc)
  
  return(NULL)
}

plot_error_heatmaps <- function(syn_hhs_spatial_file, syn_indvs_spatial_file, marg_dir) {
  syn_hhs <- read.csv(syn_hhs_spatial_file) %>%
    mutate(geoid = str_pad(geoid, 12, "left", "0")) %>%
    mutate(tract = substr(geoid, 1, 11),
           blkgp = substr(geoid, 1, 12))
  syn_indvs <- read.csv(syn_indvs_spatial_file) %>%
    mutate(geoid = str_pad(geoid, 12, "left", "0")) %>%
    mutate(tract = substr(geoid, 1, 11),
           blkgp = substr(geoid, 1, 12)) %>%
    mutate(i_inc_prox = recode(i_inc,
                               "<=0K"="<=25K",
                               "1-25K"="<=25K",
                               "76-100K"=">75K",
                               ">100K"=">75K"))
  process_marginals(marg_dir)
  heatmap1 <- error_heatmap(syn_hhs, tract_tenur_hhinc.l,
                            quo(tenur_hhinc_prox), quo(tract),
                            "Tract X Tenure X HH Income Errors")
  heatmap2 <- error_heatmap(syn_hhs, blkgp_tenur_hhsiz.l,
                            quo(tenur_hhsiz_prox), quo(blkgp),
                            "Block Group X Tenure X HH Size Errors")
  heatmap3 <- error_heatmap(syn_hhs, tract_nwork.l,
                            quo(nwork), quo(tract),
                            "Tract X Num. Workers Errors")
  heatmap4 <- error_heatmap(syn_hhs, tract_hhtype.l,
                            quo(hhtype), quo(tract),
                            "Tract X HH Type Errors")
  heatmap5 <- error_heatmap(syn_indvs, tract_i_sex_i_age.l,
                            quo(i_sex_i_age), quo(tract),
                            "Tract X Sex X Age Errors")
  heatmap6 <- error_heatmap(syn_indvs, blkgp_emply.l,
                            quo(emply_prox), quo(blkgp),
                            "Block Group X Employment Errors")
  heatmap8 <- error_heatmap(syn_hhs %>% filter(hhtype != "GQ"), tract_nwork_ncars.l,
                            quo(ncars), quo(tract),
                            "Tract X Num. Cars Errors")
  heatmap9 <- error_heatmap(syn_indvs, tract_i_inc.l,
                            quo(i_inc_prox), quo(tract),
                            "Tract X Indv. Income Errors")
  gl <- list(heatmap1$p, heatmap2$p, heatmap3$p, heatmap4$p, heatmap5$p,
             heatmap6$p, heatmap8$p, heatmap9$p)
  
  grid.arrange(
    grobs = gl,
    nrow = 3,
    layout_matrix = rbind(c(1, 1, 1, 2, 2, 2),
                          c(3, 3, 4, 4, 7, 7),
                          c(5, 5, 6, 6, 8, 8))
  )
}

get_error_map_data <- function(map_pumas, variable, spatial_level, synpop_suffix) {
  sf_geos.puma <- data.frame()
  for(pumas in map_pumas){
    print(pumas)
    syn_hhs_file <- paste0("synthpop_output/",pumas,"/syn_hhs_spatial_",pumas,"-",synpop_suffix,".csv")
    syn_indvs_file <- paste0("synthpop_output/",pumas,"/syn_indvs_spatial_",pumas,"-",synpop_suffix,".csv")
    syn_hhs <- read.csv(syn_hhs_file) %>%
      mutate(geoid = str_pad(geoid, 12, "left", "0")) %>%
      mutate(tract = substr(geoid, 1, 11),
             blkgp = substr(geoid, 1, 12))
    syn_indvs <- read.csv(syn_indvs_file) %>%
      mutate(geoid = str_pad(geoid, 12, "left", "0")) %>%
      mutate(tract = substr(geoid, 1, 11),
             blkgp = substr(geoid, 1, 12)) %>%
      mutate(i_inc_prox = recode(i_inc,
                                 "<=0K"="<=25K",
                                 "1-25K"="<=25K",
                                 "76-100K"=">75K",
                                 ">100K"=">75K"))
    marg_dir <<- paste0("synthpop_data/acs_marginals/", pumas, "/")
    process_marginals(marg_dir)
    marg_name <- if(variable!="ncars") paste0(spatial_level,"_",variable,".l") else "tract_nwork_ncars.l"
    sym_variable <- if(variable %in% c("tenur_hhinc", "tenur_hhsiz", "emply", "i_inc")) paste0(variable,"_prox") else variable
    heatmap_data <- error_heatmap(if(variable %in% c("tenur_hhinc", "tenur_hhsiz", "nwork", "hhtype", "ncars")) syn_hhs else syn_indvs, 
                                  eval(parse(text=marg_name)),
                                  sym(sym_variable), sym(spatial_level),
                                  "")$heatmap_data %>%
      rename(geoid=!!sym(spatial_level))
    
    ############
    sf_geos <- readRDS(file=paste0("map_data/geos/",pumas,"_",spatial_level,"_geos.Rds"))
    sf_geos <- left_join(sf_geos, get_rmses(heatmap_data, spatial_level), by=c("geoid"))
    sf_geos$puma <- pumas
    
    # total census spatial unit sizes
    hh_pops <- (if(spatial_level=="tract") tract_hhtype.l else blkgp_tenur_hhsiz.l) %>% 
      rename(geoid=!!sym(spatial_level), hh_pop=marg_ct) %>% 
      group_by(geoid) %>%
      summarise(hh_pop = sum(hh_pop)) %>%
      ungroup()
    indv_pops <- (if(spatial_level=="tract") tract_i_sex_i_age.l else blkgp_emply.l) %>%
      rename(geoid=!!sym(spatial_level), indv_pop=marg_ct) %>% 
      group_by(geoid) %>%
      summarise(indv_pop = sum(indv_pop)) %>%
      ungroup()
    
    sf_geos <- left_join(left_join(sf_geos, hh_pops, by=c("geoid")), indv_pops, by=c("geoid"))
    
    sf_geos.puma <- rbind(sf_geos.puma, sf_geos)
  }
  return(sf_geos.puma)
}

plot_error_map <- function(main_title, map_pumas, densities, variable, spatial_level, synpop_suffix, 
                           gray_hh_pop_under, gray_indv_pop_under, NA_label) {
  puma_labels <- c(paste0("PUMA ", substr(map_pumas[1:15], 3, 7), " (", round(densities[1:15]), " people/sq.mi.)"),
                   paste0("State of WY (", round(densities[16]), " people/sq.mi.)"))
  sf_geos.puma <- get_error_map_data(map_pumas,
                                     variable=variable,
                                     spatial_level=spatial_level,
                                     synpop_suffix=synpop_suffix)
  # convert to sf and apply HH and individual size filters
  sf_geos.puma.sf <- st_sf(sf_geos.puma) %>%
    filter(!st_is_empty(geometry)) %>%
    mutate(rmse = ifelse(indv_pop >= gray_indv_pop_under & hh_pop >= gray_hh_pop_under, rmse, NA))
  
  # combine WY PUMAs together and factor PUMA column to ensure plot ordering
  sf_geos.puma.sf$puma <- factor(recode(sf_geos.puma.sf$puma,
                                        "5600100"="State of WY",
                                        "5600200"="State of WY",
                                        "5600300"="State of WY",
                                        "5600400"="State of WY",
                                        "5600500"="State of WY"), 
                                 levels=c(map_pumas[1:15], "State of WY"))
  
  dev.new(width=8, height=8, noRStudioGD = T)
  
  spatial_name <- if(spatial_level=="tract") "Tract" else "Block Group"
  
  tm_shape(sf_geos.puma.sf) + 
    tm_polygons(col="rmse",
                style="cont",
                title=paste0("RMSE % (Census ",spatial_name,")"),
                textNA=NA_label,
                colorNA="#e6e6e6",
                legend.is.portrait=F) +
    tm_compass(position=c("right", "top")) +
    tm_scale_bar(text.size=0.8) +
    tm_facets(by="puma",
              ncol=3) +
    tm_ylab("                        Jacksonville, FL      Boston, MA        Phoenix, AZ       Houston, TX      Los Angeles, CA", 
            size=1.2, 
            space=1) +
    tm_layout(main.title=main_title,
              main.title.fontface="bold",
              main.title.size=1.4,
              main.title.position="center",
              panel.labels=puma_labels,
              panel.label.size=1.3,
              legend.outside=T,
              legend.outside.position="bottom",
              legend.outside.size=0.2,
              legend.title.size=1.175,
              legend.title.fontface="bold",
              legend.text.size=1,
              outer.margins=c(0, 0.02, 0, 0.02),
              inner.margins=c(0.06, 0.06, 0.06, 0.06)) 
}

plot_tract_hhtype_diffs_map <- function(main_title, map_pumas, densities) {
  puma_labels <- c(paste0("PUMA ", substr(map_pumas[1:15], 3, 7), " (", round(densities[1:15]), " people/sq.mi.)"),
                   paste0("State of WY (", round(densities[16]), " people/sq.mi.)"))
  
  sf_geos.puma <- data.frame()
  for(pumas in map_pumas){
    print(pumas)
    marg_dir <<- paste0("synthpop_data/acs_marginals/", pumas, "/")
    process_marginals(marg_dir)
    
    # using the column names synpop_ct and marg_ct just nominally, so that
    # we can reuse the rmse functions
    distr_by_tract <- tract_hhtype %>% 
      select(-geoid) %>% 
      as.matrix() %>% 
      prop.table(margin=1) %>%  
      `*`(100) %>%
      data.frame() %>% 
      cbind(tract_hhtype$geoid, .) %>% 
      rename(geoid=`tract_hhtype$geoid`) %>%
      mutate(geoid = str_pad(as.character(geoid), 11, "left", "0")) %>%
      pivot_longer(cols=c("MC", "NS", "SM", "NF", "GQ"), names_to="hhtype", values_to="synpop_ct") %>%
      data.frame() %>%
      select(hhtype, geoid, everything()) %>%
      mutate(synpop_ct = replace_na(synpop_ct, 0),
             marg_ct = rep(tract_hhtype %>% 
                             select(-geoid) %>% 
                             colSums() %>% 
                             prop.table() %>%
                             `*`(100), length(unique(tract_hhtype$geoid))))
    
    sf_geos <- readRDS(file=paste0("map_data/geos/",pumas,"_tract_geos.Rds"))
    sf_geos <- left_join(sf_geos, get_rmses(distr_by_tract, "tract"), by=c("geoid"))
    sf_geos$puma <- pumas
    
    sf_geos.puma <- rbind(sf_geos.puma, sf_geos)
  }
  
  # convert to sf
  sf_geos.puma.sf <- st_sf(sf_geos.puma) %>%
    filter(!st_is_empty(geometry))
  
  # combine WY PUMAs together and factor PUMA column to ensure plot ordering
  sf_geos.puma.sf$puma <- factor(recode(sf_geos.puma.sf$puma,
                                        "5600100"="State of WY",
                                        "5600200"="State of WY",
                                        "5600300"="State of WY",
                                        "5600400"="State of WY",
                                        "5600500"="State of WY"), 
                                 levels=c(map_pumas[1:15], "State of WY"))
  
  dev.new(width=8, height=8, noRStudioGD = T)
  
  spatial_name <- "Tract"
  
  tm_shape(sf_geos.puma.sf) + 
    tm_polygons(col="rmse",
                style="cont",
                palette="BuPu",
                title=paste0("RMSE % (Census ",spatial_name,")"),
                colorNA="#e6e6e6",
                legend.is.portrait=F) +
    tm_compass(position=c("right", "top")) +
    tm_scale_bar(text.size=0.8) +
    tm_facets(by="puma",
              ncol=3) +
    tm_ylab("                        Jacksonville, FL      Boston, MA        Phoenix, AZ       Houston, TX      Los Angeles, CA", 
            size=1.2, 
            space=1) +
    tm_layout(main.title=main_title,
              main.title.fontface="bold",
              main.title.size=1.25,
              main.title.position="center",
              panel.labels=puma_labels,
              panel.label.size=1.3,
              legend.outside=T,
              legend.outside.position="bottom",
              legend.outside.size=0.2,
              legend.title.size=1.175,
              legend.title.fontface="bold",
              legend.text.size=1,
              outer.margins=c(0, 0.02, 0, 0.02),
              inner.margins=c(0.06, 0.06, 0.06, 0.06)) 
}

plot_pie <- function(main_title, variable, spatial_level, puma) {
  marg_dir <<- paste0("synthpop_data/acs_marginals/", puma, "/")
  process_marginals(marg_dir)
  marg_name <- if(variable!="ncars") paste0(spatial_level,"_",variable,".l") else "tract_nwork_ncars.l"
  sym_variable <- if(variable %in% c("tenur_hhinc", "tenur_hhsiz", "emply", "i_inc")) paste0(variable,"_prox") else variable
  
  pie_data <- eval(parse(text=marg_name)) %>% 
    mutate(var = factor(!!sym(sym_variable), levels=unique(!!sym(sym_variable)))) %>%
    group_by(var) %>% 
    summarise(marg_ct = sum(marg_ct))
  pie_data.vec <- setNames(pie_data$marg_ct, pie_data$var)
  
  # save the figure
  png(file=paste0(variable,"_",puma,".png"), width=250, height=250)
  pie(pie_data.vec,
      col=viridis(length(pie_data.vec)),
      main=main_title)
  dev.off()
}

############################### PLOT SPECS #####################################

use_prev_marg <- T
download_geometries <- F

map_pumas <- c("0603701", "0603716", "0603733",
               "4805000", "4804901", "4804614",  
               "0400134", "0400105", "0400122",
               "2501400", "2503603", "2503302",
               "1208900", "1203105", "1203103",
               "5600200", "5600100", "5600300", "5600400", "5600500")

densities <- c(383.7, 9049.6, 42788.1,
               66.1, 3457.2, 11736.8,
               112.0, 3947.0, 8075.6,
               885.9, 5817.3, 28710.4,
               160.5, 1989.8, 3285.1,
               # last entry is density for entire state of WY
               5.97)

############################## STUDY AREA MAPS #################################

plot_study_areas(map_pumas)

###################### TRACT LEVEL HHTYPE DIFFERENCES ##########################

plot_tract_hhtype_diffs_map("Tract-level deviation of household type composition from PUMA mean", 
                            map_pumas, densities)

################################# PIE CHARTS ###################################

var_names <- c("tenur_hhinc", "tenur_hhsiz", "nwork", "hhtype", "i_sex_i_age", "emply", "ncars", "i_inc")
spatial_levels <- c("tract", "blkgp", "tract", "tract", "tract", "blkgp", "tract", "tract")
pie_titles <- c("Tenure by Household Income", "Tenure by Household Size", "Number of Workers", "Household Type", 
                "Sex by Age", "Employment Status", "Number of Cars", "Income")
for(v in seq_along(var_names)) {
  for(puma in map_pumas) {
    plot_pie(pie_titles[v], 
             variable=var_names[v],
             spatial_level=spatial_levels[v],
             puma=puma)
  }
}

################################## ERROR MAPS ##################################

synpop_suffix <- "20220127"

##### DOWNLOAD GEOMETRIES #####
# get and save tract- and block group-level geometries
# do this once if you do not have geometries downloaded; then, set download_geometries <- FALSE
if(download_geometries) {
  if(!dir.exists("map_data/geos")) {
    dir.create("map_data/geos", recursive=T)
  }
  for(pumas in map_pumas){
    marg_dir <- paste0("synthpop_data/acs_marginals/", pumas, "/")
    my_census_api_key <- "ac6cb3e106c860e52384fe71cf0407a13c25b96c"
    puma_tract_equiv_file <- "synthpop_data/2010_Census_Tract_to_2010_PUMA.csv"
    tract_block_equiv_file <- "synthpop_data/tract_block_equiv_2010.Rds"
    process_pums_file <- "process_pums.R"
    process_marg_file <- "process_marg.R"
    source("prepare_data.R")
    
    sf_geos.tract <- get_acs(geography="tract", state=substr(tracts[1], 1, 2), table="B11001", output="wide", geometry = TRUE) %>%
      filter(GEOID %in% tracts) %>%
      select(GEOID, NAME, geometry) %>%
      rename(geoid=GEOID) %>%
      arrange(geoid) %>%
      data.frame()
    
    saveRDS(sf_geos.tract, file=paste0("map_data/geos/",pumas,"_tract_geos.Rds"))
    
    sf_geos.blkgp <- get_acs(geography="block group", state=substr(tracts[1], 1, 2), table="B25009", output="wide", geometry = TRUE) %>%
      filter(GEOID %in% blkgps) %>%
      select(GEOID, NAME, geometry) %>%
      rename(geoid=GEOID) %>%
      arrange(geoid) %>%
      data.frame()
    
    saveRDS(sf_geos.blkgp, file=paste0("map_data/geos/",pumas,"_blkgp_geos.Rds"))
  }
}

##### ERROR HEATMAPS #####

hm <- plot_error_heatmaps(syn_hhs_spatial_file="synthpop_output/2503302/syn_hhs_spatial_2503302-example.csv", 
                          syn_indvs_spatial_file="synthpop_output/2503302/syn_indvs_spatial_2503302-example.csv", 
                          marg_dir="synthpop_data/acs_marginals/2503302/")

##### CALIBRATION MARGINALS #####

plot_error_map("Spatial errors in tenure by houshold income distribution",
               map_pumas, densities,
               variable="tenur_hhinc", 
               spatial_level="tract", 
               synpop_suffix=synpop_suffix, 
               gray_hh_pop_under=30, 
               gray_indv_pop_under=50, 
               NA_label="too few households")

plot_error_map("Spatial errors in tenure by houshold size distribution",
               map_pumas, densities,
               variable="tenur_hhsiz", 
               spatial_level="blkgp", 
               synpop_suffix=synpop_suffix, 
               gray_hh_pop_under=30, 
               gray_indv_pop_under=50, 
               NA_label="too few households")

plot_error_map("Spatial errors in number of workers distribution",
               map_pumas, densities,
               variable="nwork", 
               spatial_level="tract", 
               synpop_suffix=synpop_suffix, 
               gray_hh_pop_under=30, 
               gray_indv_pop_under=50, 
               NA_label="too few households")

plot_error_map("Spatial errors in household type distribution",
               map_pumas, densities,
               variable="hhtype", 
               spatial_level="tract", 
               synpop_suffix=synpop_suffix, 
               gray_hh_pop_under=30, 
               gray_indv_pop_under=50, 
               NA_label="too few households")

plot_error_map("Spatial errors in sex by age distribution",
               map_pumas, densities,
               variable="i_sex_i_age", 
               spatial_level="tract", 
               synpop_suffix=synpop_suffix, 
               gray_hh_pop_under=30, 
               gray_indv_pop_under=50, 
               NA_label="too few individuals")

plot_error_map("Spatial errors in employment distribution",
               map_pumas, densities,
               variable="emply", 
               spatial_level="blkgp", 
               synpop_suffix=synpop_suffix, 
               gray_hh_pop_under=30, 
               gray_indv_pop_under=50, 
               NA_label="too few individuals")

##### VALIDATION MARGINALS #####

plot_error_map("Spatial errors in number of cars distribution",
               map_pumas, densities,
               variable="ncars", 
               spatial_level="tract", 
               synpop_suffix=synpop_suffix, 
               gray_hh_pop_under=30, 
               gray_indv_pop_under=50, 
               NA_label="too few households")

plot_error_map("Spatial errors in income distribution",
               map_pumas, densities,
               variable="i_inc", 
               spatial_level="tract", 
               synpop_suffix=synpop_suffix, 
               gray_hh_pop_under=30, 
               gray_indv_pop_under=50, 
               NA_label="too few individuals")



