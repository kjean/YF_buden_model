rm(list=ls(all=TRUE))

model = "FOI" # model = "FOI" or "R0"

library(maptools)
library(sp) 
library(shapefiles)
library(fields)
library(Hmisc)
library(EnvStats)


computer<-3
if (computer==0) { homedir<-"Y:/Kevin/"
} else if (computer==3) { homedir<-"C:/Users/Kevin JEAN/Desktop/Recherche/Yellow_Fever/Kevin_YF_DIDE_share_drive/GitDir/"  
}

run_or_load = "run"
save_work_env=T



shpdir = commondir = currdir = homedir

if(computer==5) currdir=homedir

setwd(currdir)




if(run_or_load == "run"){
 
  ####### ####### ####### ####### ####### 
  ####### load shapefiles
  
  
  shp0 = readShapePoly(paste(shpdir, "Africa_adm0.shp",sep=""))
  shp1 = readShapePoly(paste(shpdir, "Africa_adm1.shp",sep=""))
  
  shp1$adm0_adm1<-paste(shp1$ISO, shp1$ID_1, sep="_")
  shp1 = shp1[order(shp1$adm0_adm1),]
  
  
  c34 = c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ERI", "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GNQ", "KEN", "LBR", "MLI", "MRT", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZMB")
  country34 = c("Angola","Burundi","Benin","Burkina Faso","Central African Republic","Cote d'Ivoire","Cameroon","Democratic Republic of the Congo","Republic of Congo","Eritrea","Ethiopia","Gabon","Ghana","Guinea","Gambia","Guinea-Bissau","Equatorial Guinea","Kenya","Liberia","Mali","Mauritania","Niger","Nigeria","Rwanda","Sudan","Senegal","Sierra Leone","Somalia", "South Sudan", "Chad","Togo","Tanzania","Uganda", "Zambia") 
  
  
  ####### ####### ####### ####### ####### 
  ####### load serological data
  load("../serological_data.Rdata")
  
  
  
  ###########
  source(paste0(model, "_functions/", model, "_functions_launch_data.R"))
  
  ####### ####### ####### ####### ####### 
  ####### load environmental data
  env_table = paste(homedir,"environment/Africa_adm1_dat.csv", sep="")
  envdat = launch_env_dat(env_table)
  names(envdat)
  dim(envdat$dat)
  dat = envdat$dat
  
  
  #### fit glm
  objet_glm = fit_glm(dat =envdat$dat, depi = envdat$depi)
  names(objet_glm)
  
  beta0 = objet_glm[[1]]
  x = objet_glm[[2]]
  y = objet_glm[[3]]
  
  
  ####### ####### ####### ####### ####### 
  ####### population and vaccination data
  source(paste0(model, "_functions/", model, "_functions_format_pop_data.R"))
  
  pop_path = paste(homedir, "population/", sep="")
  
  pop1 = import_pop_data_for_FOI(path =pop_path , c_country=c34, pop_source="LS2014", adm="adm1");dim(pop1)
  pop1 = pop1[pop1$adm0_adm1 %in% envdat$dat$adm0_adm1,];dim(pop1)
  pop.full = pop1
  
  # latest vaccination data csv file
  vc.dir = paste(commondir, "vaccination_scenarios/", sep="")
  latest_vac2d_csv = "vaccination_coverage_year_cut_2017_best_estimate_skew0.csv"
  vc1 = read.csv(paste0(vc.dir,latest_vac2d_csv ), h=T, stringsAsFactors = F); dim(vc1)
  names(vc1)
  names(vc1)[names(vc1)=="adm1"]= "adm0_adm1"
  vc1 = vc1[,-which(names(vc1)=="country")]
  vc1 = vc1[vc1$adm0_adm1 %in% envdat$dat$adm0_adm1,];dim(vc1)
  all(vc1$adm0_adm1 == pop1$adm0_adm1)
  vc1 = repair_vc_data(vc1)
  
  
  vc.full_100 = vc1 # vaccine coverage with vaccine efficacy
  vc3d = transform_into_vc3d(vc1)
  
  
  ## life expectation data
  life_exp_path =  paste0(homedir, "population/life expectancy_by_country_year_1950-2100_kj_2015.csv")
  life0 = import_life_exp(life_exp_path)
  
  
  ## create pop30.agg, pop30.agg and pop.vc.moments
  list_moment = create_pop30.agg_vc30.agg(pop1, vc1)
  pop30.agg = list_moment$pop30.agg
  vc30.agg = list_moment$vc30.agg
  pop.vc.moments = list_moment$pop.vc.moments
  vc30.agg_100 = vc30.agg # vaccine coverage with vaccine efficacy
  
  
  #### creat pop vc data for YFSD years
  pop1.afro.years = pop1[pop1$year>=2005 & pop1$year<=2013,]
  vc1.afro.years = vc1[vc1$year>=2005 & vc1$year<=2013,]
  pop1.afro.years[is.na(pop1.afro.years) | is.na(vc1.afro.years)] = 0
  vc1.afro.years[is.na(pop1.afro.years) | is.na(vc1.afro.years)] = 0
  
  
  #### creat pop vc data for serosurveys provinces
  
  
  
  #### creat pop.agg vc.agg data for serosurveys provinces
  list.agg.ag = create_pop_and_vc.agg.ag(pop1, vc1, n.serosurveys, adm1s)
  pop.agg = list.agg.ag$pop.agg
  vc.agg = list.agg.ag$vc.agg
  vc.agg.ag = list.agg.ag$vc.agg.ag
  
  
  
  
  #########################################################################
  ####### Create status0
  source(paste0(model, "_functions/", model, "_functions_create_status0.R"))
  list_status0 = create_status0(x, beta0, n.serosurveys)
  status0 = list_status0$status0
  sdj = list_status0$sdj
  sd.prior = list_status0$sd.prior
  mu.prior= list_status0$mu.prior
  varsin.nc = list_status0$varsin.nc
  transform = list_status0$transform
  pars.ini = list_status0$pars.ini
  ii = list_status0$ii
  
  
  
  #########################################################
  #### save workspace
  
  if(save_work_env){
    save(file="workspace_FOI_model_mcmc.Rdata", list=setdiff(ls(), c("homedir", "currdir", "scenar_type", "commondir", "run_or_load",
                                                                     "computer")))
  }
  
} else if (run_or_load == "load"){ # end if run_or_load == "run"  
  load("workspace_FOI_model_mcmc.Rdata")
}

##### need to reset the paths that have changed if loaded






#########################################################################
####### run MCMC
source(paste0(model, "_functions/", model, "_main_functions_for_MCMC.R"))
amat = matrix(0:100,nrow=nrow(pop30.agg),ncol=101,byrow=T)

## tweak MCMC
status = status0
set.seed(101)
results_tweak = tweak.mcmc.FOI.model(chain.length=5,  #500
                                     bb.max=3,  #20
                                     omit=1,  #1
                                     sdj=sdj,
                                     mu.prior=mu.prior, sd.prior=sd.prior, 
                                     transform=transform,varsin.nc=varsin.nc)
names(results_tweak)

mcmc.params = results_tweak$mcmc.params
sd_mat = results_tweak$sd_mat
accept = results_tweak$accept
status = results_tweak$status
sdj = results_tweak$sdj



## run MCMC
index_code = as.numeric(commandArgs(TRUE))
if(identical(index_code, numeric(0))) index_code=1

t= as.numeric(Sys.time())
seed= (t - floor(t)) * 1e8 
set.seed(seed)

results_mcmc_foi = run.mcmc.FOI.model(chain.length = 3, #400
                                      omit = 1, #100
                                      bb.tot=2, #250
                                      nb_codes= 2, #25
                                      index_code = index_code,
                                      status = status,
                                      sdj=sdj,
                                      mu.prior=mu.prior, sd.prior=sd.prior, 
                                      transform=transform,varsin.nc=varsin.nc,
                                      print_predictions = T, print_p.detect.nb_inf = T, print_burden = F)


