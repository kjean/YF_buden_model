rm(list=ls(all=TRUE))

model = "R0" # model = "FOI" or "R0"

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
  avec=c(0:100)
  
  
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
  all_res_pop_3d = get_pop_data_3d(path =pop_path, c_country=c34, pop_source="LS2014", adm="adm1", dat=dat)
  
  pop1 = all_res_pop_3d$pop1
  pop3d = all_res_pop_3d$pop3d
  P_tot_2d = all_res_pop_3d$P_tot_2d
  p_prop3d = all_res_pop_3d$p_prop3d
  dn1 = dimnames(pop3d)[[1]]
  dn2 = as.numeric(dimnames(pop3d)[[2]])
  dn3 = dimnames(pop3d)[[3]]
  dn1_survey = sero.studies
  #save(p_prop3d, file="p_prop3d.Rdata")
  #save(P_tot_2d, file="P_tot_2d.Rdata")
  
  ## life expectation data
  life_exp_path =  paste0(homedir, "population/life expectancy_by_country_year_1950-2100_kj_2015.csv")
  life0 = import_life_exp(life_exp_path)
  
  ### vaccination data
  vc.dir = paste(commondir, "vaccination_scenarios/", sep="")
  latest_vac2d_csv = "vaccination_coverage_year_cut_2017_best_estimate_skew0.csv"
  
  vc2d = read.csv(paste0(vc.dir,latest_vac2d_csv),stringsAsFactors = F, h=T); dim(vc2d)
  names(vc2d)[names(vc2d)=="country"]= "adm0"
  names(vc2d)[names(vc2d)=="adm1"]= "adm0_adm1"
  vc2d = repair_vc_data(vc2d, dat=dat)
  
  vc3d = transform_into_vc3d(vc2d,  adm="adm1")
  
  ## calculate inc_v3d
  inc_v3d = calc_incidence_vac(vc3d)
  dim(inc_v3d)
  dimnames(inc_v3d)
  
  ### calculate t0_vac_africa
  t0_vac_africa = calc_t0_vac_africa(vc3d)
  
  ### 
  #vc2d = vc2d[,-which(names(vc2d)=="adm0")]
  pop.moments.whole = calc.pop.moments.whole(p_prop3d, t0_vac_africa)
  
  
  #### aggregate pop and vc at serosurvey sites
  list_pop.agg_vc.agg =create_pop.agg_vc.agg(pop1=pop1, vc2d=vc2d)
  pop.agg3d = list_pop.agg_vc.agg$pop.agg3d
  vc.agg3d = list_pop.agg_vc.agg$vc.agg3d
  
  
  inc_v3d.agg = calc_inc_v3d.agg(vc.agg3d);dim(inc_v3d.agg)
  pop.moments.agg = calc.pop.moments.agg(pop.agg3d); dim(pop.moments.agg)
  
  
  ##########
  # launch MCMC function to create R0 lookup table
  source(paste0(model, "_functions/", model, "_main_functions_for_MCMC.R"))
  
  
  ### needs to adapt fun.recurrence.seroprev.whole
  ###### create the R0 look.up table
  if(!file.exists("R0_lookup_table.Rdata")){
    ptm = proc.time()
    create_R0.lookup()
    proc.time()-ptm
  } else {
    load("R0_lookup_table.Rdata")
  }
  
  
  
  ###### create status R0
  foi_const_surv = c(0,1e-6,0,0,0,0,rep(0,25))
  source(paste0(model, "_functions/", model, "_functions_create_status0.R"))
  list_status0 = create_status0(x, beta0, n.serosurveys)
  status0 = list_status0$status0
  sdj = list_status0$sdj
  sd.prior = list_status0$sd.prior
  varsin.nc = list_status0$varsin.nc
  transform = list_status0$transform
  pars.ini = list_status0$pars.ini
  ii = list_status0$ii
  
  list.pop.at.survey = create.pop.at.survey(pop.agg3d)
  p_at_survey_3d = list.pop.at.survey$p_at_survey_3d
  P_tot_survey_2d = list.pop.at.survey$P_tot_survey_2d
  
  #########################################################
  #### save workspace
  
  if(save_work_env){
    save(file="workspace_R0_model_mcmc.Rdata", list=setdiff(ls(), c("homedir", "currdir", "scenar_type", "commondir", "run_or_load",
                                                                    "computer")))
  }

} else if (run_or_load == "load"){ # end if run_or_load == "run"  
  load("workspace_R0_model_mcmc.Rdata")
}




#########################################################################
####### tweak MCMC
source(paste0(model, "_functions/", model, "_main_functions_for_MCMC.R"))

status = status0
if(!file.exists("sdj_after_tweaking_nb_run=12500.csv")){
  results_tweak = tweak.mcmc.R0.model(chain.length=500,  #500
                                      bb.max=25,   #25
                                      omit=1,   #1
                                      sdj=sdj,
                                      sd.prior=sd.prior, 
                                      transform=transform,varsin.nc=varsin.nc)

} else if(file.exists("sdj_after_tweaking_nb_run=12500.csv")){ # if sdj has already been calculated thus run 2000 chains for burnin
  sdj = as.numeric(read.csv(paste0(currdir,"sdj_after_tweaking_nb_run=12500.csv" ), h=T)[,1])
  results_tweak = tweak.mcmc.R0.model(chain.length=2,  #2000
                                      bb.max=1,   #1
                                      omit=1,   #1
                                      sdj=sdj,
                                      sd.prior=sd.prior, 
                                      transform=transform,varsin.nc=varsin.nc)
}

mcmc.params = results_tweak$mcmc.params
sd_mat = results_tweak$sd_mat
accept = results_tweak$accept
status = results_tweak$status
sdj = results_tweak$sdj


##################################
## run MCMC
index_code = as.numeric(commandArgs(TRUE))
if(identical(index_code, numeric(0))) index_code=1

t= as.numeric(Sys.time())
seed= (t - floor(t)) * 1e8 
set.seed(seed)

results_mcmc_r0 = run.mcmc.r0.model(chain.length = 2, #400
                                      omit = 1, #100
                                      bb.tot=1, #250
                                      nb_codes= 1, #25
                                      index_code = index_code,
                                      status=status,
                                      sdj=sdj,
                                      sd.prior=sd.prior, 
                                      transform=transform,varsin.nc=varsin.nc,
                                      print_predictions = T, print_p.detect.nb_inf = F, print_burden = T)