
import_pop_data = function(path, c_country, pop_source, adm){
 
  pop_year_age = NULL          
  if(pop_source == "LS2014")   {
    for(adm0 in c_country){
      if(!grepl("NFS", homedir)) print(adm0)
      # tmp = read.csv(paste0("W:/Data/processed_Afpop/pop_size/from_LandScan2014/", adm, "/pop_year_age_", adm, "_", adm0, "_LandScan2014.csv"),
      #                h=T)
      tmp = read.csv(paste0(path, "/from_LandScan2014/", adm, "/pop_year_age_", adm, "_", adm0, "_LandScan2014.csv"),
                     h=T)
      pop_year_age = rbind(pop_year_age, tmp)
    }
  } 
  pop_year_age[,adm] = paste(pop_year_age[,"adm0"], pop_year_age[,adm], sep="_")
  colnames(pop_year_age)[colnames(pop_year_age)==adm] = paste("adm0", adm, sep = "_")
  return(pop_year_age)
}



###########
add_1940_1950 = function(pop2d){
  if(min(pop2d$year) != 1950)  stop("pop2d should start in 1950")
  
  for(y in 1949:1940) {
    pop1.early = pop2d[pop2d$year==y+1,]
    pop1.early$year = y
    pop2d = rbind(pop1.early,pop2d)
  }
  
  return(pop2d)
}

###########
transform_into_pop3d = function(pop2d, adm){
 
  adm0_adm = paste0("adm0_", adm)
  
  pop3d=rep(NA, nrow(pop2d)*(ncol(pop2d)-3))
  dim(pop3d)=c(length(table(pop2d[,adm0_adm])), length(table(pop2d$year)), ncol(pop2d)-3)
  # 1st dim = adm0_adm1
  # 2nd dim = year
  # 3rd dim = age
  
  dn1 = names(table(pop2d[,adm0_adm]))
  dn2 = as.numeric(names(table(pop2d$year))) # 111 annees
  dn3 = names(pop2d)[4:length(names(pop2d))] # 101 classes d'age
  
  for(a in 1:length(dn3)) { # dn3 = nb de classes d age =101
    for(y in min(dn2):max(dn2)) { 
      mm = match(pop2d[,adm0_adm][pop2d$year==y],dn1) 
      pop3d[mm,y-min(dn2)+1,a]=pop2d[pop2d$year==y,a+3]
    }
  }
  
  dimnames(pop3d)[[1]] <- dn1
  dimnames(pop3d)[[2]] <- dn2
  dimnames(pop3d)[[3]] <- dn3
  
  return(pop3d)
  
}



#################
get_P_tot_2d= function(pop3d){
  
  dn1 = dimnames(pop3d)[[1]]
  dn2 = dimnames(pop3d)[[2]]
  
  P_tot_2d = rep(NA, length(dn1)*length(dn2) )
  dim(P_tot_2d)=c(length(dn1),length(dn2) )
  for (adm in 1:length(dn1)){
    for (year in  1:length(dn2) ){
      P_tot_2d[adm,year] = sum(pop3d[adm,year,], na.rm=T)
    }
  } 
  rownames(P_tot_2d)=dn1
  colnames(P_tot_2d)=dn2
  
  return(P_tot_2d)
}


################
get_p_prop3d = function(pop3d, P_tot_2d){
  
  dn1 = dimnames(pop3d)[[1]]
  dn2 = dimnames(pop3d)[[2]]
  dn3 = dimnames(pop3d)[[3]]
  
  p_prop3d=rep(NA, dim(pop3d)[1]*dim(pop3d)[2]*dim(pop3d)[3])
  dim(p_prop3d)=dim(pop3d)
  for(adm in 1:length(dn1)){
    for (year in 1:length(dn2)){
      for (a in 1:length(dn3)){
        p_prop3d[adm,year,a] = pop3d[adm, year, a]/P_tot_2d[adm, year]
      }
    }
  }
  dimnames(p_prop3d)[[1]] <- dn1
  dimnames(p_prop3d)[[2]] <- dn2
  dimnames(p_prop3d)[[3]] <- dn3
  
  return(p_prop3d)
}




################
# overall functions

get_pop_data_3d = function(path, c_country, pop_source, adm, dat=dat){
  if (!pop_source %in% c("LS2011", "LS2014", "WP2010")) {
    stop("pop source has to be LS2011, LS2014, WP2010")
  }
  if (!adm %in% c("adm1", "adm2")) {
    stop("adm source has to be adm1 or adm2")}
  adm0_adm = paste0("adm0_", adm)
  
  pop1 = import_pop_data(path=path, c_country, pop_source, adm)
 # pop1[is.na(pop1)]=0 # to avoid the part with na before 1995 that put the mess in inc_v3d - pas sur que ce soit utile, le pb est aussi dans les données vaccinales
  
  if(!grepl("NFS", homedir)) print("add 1940 to 1950")
  pop1 = add_1940_1950(pop1)
  
  # restrict to lines in dat
  pop1 = pop1[pop1[,adm0_adm] %in% dat[,adm0_adm],]
  
  if(!grepl("NFS", homedir)) print("create pop3d")
  pop3d = transform_into_pop3d(pop2d=pop1, adm)
  
  if(!grepl("NFS", homedir)) print("create P_tot_2d")
  P_tot_2d = get_P_tot_2d(pop3d=pop3d)
  
  if(!grepl("NFS", homedir)) print("create p_prop3d")
  p_prop3d = get_p_prop3d(pop3d=pop3d, P_tot_2d=P_tot_2d)
  
  res = list(pop1= pop1, pop3d= pop3d, P_tot_2d=P_tot_2d, p_prop3d=p_prop3d)
  
  return(res)
}





#################################
# vc data
repair_vc_data = function(vc2d, dat){# before 1995, we have NA values for those aged >75. We affect the value of the latest age group
  for (a in 2:length(dn3)){
    vc2d[,a] = ifelse(is.na(vc2d[,a]), vc2d[,a-1], vc2d[,a])
  }
  # restrict to lines in dat
  vc2d = vc2d[vc2d[,"adm0_adm1"] %in% dat[,"adm0_adm1"],]
  return(vc2d)
}


transform_into_vc3d = function(vc2d, adm){
  
  # just transform adm1 in adm0_adm1 /adm2 in adm0_adm2
  # temp =  paste0(vc2d[,"adm0"], "_", vc2d[,adm] )
  # names(vc2d)[names(vc2d) == adm ] =  paste0("adm0_", adm)
  # vc2d[, paste0("adm0_", adm)] = temp
  
  adm0_adm = paste0("adm0_", adm)
  
  vc3d=rep(NA, nrow(vc2d)*(ncol(vc2d)-3))
  dim(vc3d)=c(length(table(vc2d[,adm0_adm])), length(table(vc2d$year)), ncol(vc2d)-3)
  # 1st dim = adm0_adm1
  # 2nd dim = year
  # 3rd dim = age
  
  dn1 = names(table(vc2d[,adm0_adm]))
  dn2 = as.numeric(names(table(vc2d$year))) # 111 annees
  dn3 = names(vc2d)[4:length(names(vc2d))] # 101 classes d'age
  
  for(a in 1:length(dn3)) { # dn3 = nb de classes d age =101
    for(y in min(dn2):max(dn2)) { 
      mm = match(vc2d[,adm0_adm][vc2d$year==y],dn1) 
      vc3d[mm,y-min(dn2)+1,a] = vc2d[vc2d$year==y,a+3]
    }
  }
  
  dimnames(vc3d)[[1]] <- dn1
  dimnames(vc3d)[[2]] <- dn2
  dimnames(vc3d)[[3]] <- dn3
  
  return(vc3d)
  
}



#################################
# life expectation data
import_life_exp = function(life_data_path){
  life0 = read.csv(life_data_path,stringsAsFactors=FALSE)
  life0 = life0[life0$country %in% country34,]
  life0$country = c34[match(life0$country,country34)] # converting from country to adm0.
  
  life0[,-(1:2)] = life0[,-(1:2)]-matrix(0:100,nrow=nrow(life0),ncol=101,byrow=TRUE)
  life0[,ncol(life0)] = 0
  
  return(life0)
}





#################################
# incidence vaccination 
calc_incidence_vac = function(vc3d){
  inc_v3d=rep(NA, dim(vc3d)[1]*dim(vc3d)[2]*dim(vc3d)[3] )
  dim(inc_v3d)= dim(vc3d)
  dimnames(inc_v3d)=dimnames(vc3d)
  
  for (adm in 1:dim(vc3d)[1] ) { #
    inc_v3d[adm,1,]= 0 # no vaccination in 1940
    
    for(year in 2:(dim(vc3d)[2]-1)){
      inc_v3d[adm,year,1] = vc3d[adm,year+1,1] # for a0, incidence = coverage in the a1 age group the next year
      
      for(age in 2:(dim(vc3d)[3]-1)) { 
        if(!is.na(vc3d[adm,year,age]) & vc3d[adm,year,age]==1){ #if vc is already =1, incidence =0
          inc_v3d[adm,year,age]=0
        } else {                             
          inc_v3d[adm,year,age]=  ( vc3d[adm,year+1,age+1] - vc3d[adm,year,age] ) / (1 - vc3d[adm,year,age] ) # among those suscpetible at y-1
          inc_v3d[adm,year,age]=ifelse(inc_v3d[adm,year,age]<10e-15,0,inc_v3d[adm,year,age])
        }# with rounding, some values become negative
      }
      inc_v3d[adm,year,dim(vc3d)[3]]=0 # incidence vaccination = 0 for age=100
    }
    
    inc_v3d[adm,dim(vc3d)[2],]=0# incidence vaccination = 0 for year=2050
    inc_v3d[adm,dim(vc3d)[2],1] = vc3d[adm,dim(vc3d)[2],1] # except among new borns
  }
  
  return(inc_v3d)
}


###########################
# first year of vaccination
calc_t0_vac_africa = function(vc3d){
  t0_vac_africa = rep(NA,dim(vc3d)[1])
  for (i in 1:dim(vc3d)[1]){
    if(!grepl("NFS", homedir)) print(i)
    year_i = 1940
    sum_vac=0
    while(sum_vac == 0 & year_i<=2050) {
      sum_vac = sum( vc3d[i, year_i-1940+1,], na.rm=T)
      year_i = year_i + 1
    }
    t0_vac_africa[i]=year_i-2
  }
  names(t0_vac_africa)=dimnames(vc3d)[[1]]
  return(t0_vac_africa)
}



#############
## calculate pop.moments.whole =  moments of the Taylor expansion
calc.pop.moments.whole = function(p_prop3d, t0_vac_africa){

  avec=c(0:100)# a matrix of age with the same dimension that pop.agg for 1 fixed year
  
  pop.moments.whole = rep(NA, length(dn1)*6)
  #pop.tot.whole= rowSums(pop1[,-c(1:3)], na.rm=T)
  dim(pop.moments.whole)=c(length(dn1),6)
  
  # Before, pop.moment was calculated for each year, then I selected the year of vaccine intro for Taylor expansion
  # that's not relevant, especially for provinces where t0_vac is 2050
  # Now, I calculate the mean pop.moments for between 1940 and t0_vac if t0_vac <2013
  # else I calculate the mean pop.moment between 1984 and 2013
  for (adm in 1:length(dn1)){
    if (t0_vac_africa[adm]<2013){
      vec = which(dn2==1940):which(dn2==t0_vac_africa[adm])
      pop.moments.whole_tmp = rep(NA, length(vec)*6)
      dim(pop.moments.whole_tmp) = c(length(vec), 6)
      for(year in vec){
        pop.moments.whole_tmp[year-min(vec)+1, 1] = sum(p_prop3d[adm, year,], na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 2] = sum(p_prop3d[adm, year,]*avec, na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 3] = sum(p_prop3d[adm, year,]*avec*avec, na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 4] = sum(p_prop3d[adm, year,]*avec*avec*avec, na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 5] = sum(p_prop3d[adm, year,]*avec*avec*avec*avec, na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 6] = sum(p_prop3d[adm, year,]*avec*avec*avec*avec*avec, na.rm=T)
      }
      
    } else if (t0_vac_africa[adm]>=2013){
      pop.moments.whole_tmp = rep(NA, 30*6)
      dim(pop.moments.whole_tmp) = c(30, 6)
      vec = which(dn2==1984):which(dn2==2013)
      for(year in vec){
        pop.moments.whole_tmp[year-min(vec)+1, 1] = sum(p_prop3d[adm, year,], na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 2] = sum(p_prop3d[adm, year,]*avec, na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 3] = sum(p_prop3d[adm, year,]*avec*avec, na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 4] = sum(p_prop3d[adm, year,]*avec*avec*avec, na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 5] = sum(p_prop3d[adm, year,]*avec*avec*avec*avec, na.rm=T)
        pop.moments.whole_tmp[year-min(vec)+1, 6] = sum(p_prop3d[adm, year,]*avec*avec*avec*avec*avec, na.rm=T)
      }
    }
    pop.moments.whole[adm,] = apply(pop.moments.whole_tmp,2,mean)
  }
  return(pop.moments.whole)
}






##############################################
##### function aggregating pop and vc at serosurvey sites

#### first function creating pop.agg and vc.agg 

create_pop.agg_vc.agg = function(pop1, vc2d){
  
  pop.agg = vc.agg = NULL
  for(year in 1940:2100) { # if needed remplace 2050 by max(unlist(study.years)) 
    for(i in 1:n.serosurveys) {
      #if (year <= max(unlist(study.years[[i]])) ){
      pop.tmp = pop1[pop1$year == year & pop1$adm0_adm1 %in% adm1s[[i]],]
      vc.tmp = vc2d[vc2d$year == year & vc2d$adm0_adm1 %in% adm1s[[i]],]
      if (nrow(pop.tmp)==1) {
        pop.tmp$adm0 = sero.studies[i]
        vc.tmp$adm0 = sero.studies[i]
      } else
        if(nrow(pop.tmp)>1) { # need to do some summing/averaging...
          pop.tmp[is.na(vc.tmp)] = NA
          vc.tmp = colSums(pop.tmp[,-(1:3)]*vc.tmp[,-(1:3)],na.rm=T)/colSums(pop.tmp[,-(1:3)],na.rm=T) 
          dim(vc.tmp) = c(1,length(vc.tmp))
          colnames(vc.tmp) = paste("a",0:100,sep="")
          vc.tmp = data.frame(adm0 = sero.studies[i], adm0_adm1 = NA, year = year, vc.tmp)
          
          pop.tmp = colSums(pop.tmp[,-(1:3)])
          dim(pop.tmp) = c(1,length(pop.tmp))
          colnames(pop.tmp) = paste("a",0:100,sep="")
          pop.tmp = data.frame(adm0 = sero.studies[i], adm0_adm1 = NA, year = year,pop.tmp)
        }
      pop.agg = rbind(pop.agg,pop.tmp)
      vc.agg = rbind(vc.agg,vc.tmp)
      #}
    }
  }
  rm(pop.tmp, vc.tmp)
  
  pop.agg[,-(1:3)][is.na(pop.agg[,-(1:3)])] = 0
  vc.agg[,-(1:3)][is.na(vc.agg[,-(1:3)])] = 0
  
  
  ## pass in 3d
  pop.agg3d=rep(NA, nrow(pop.agg)*(ncol(pop.agg)-3))
  dim(pop.agg3d)=c(length(sero.studies), length(table(pop.agg$year)), ncol(pop.agg)-3)
  vc.agg3d=rep(NA, nrow(vc.agg)*(ncol(vc.agg)-3))
  dim(vc.agg3d)=c(length(table(vc.agg$adm0)), length(table(vc.agg$year)), ncol(vc.agg)-3)
  
  dn1_survey = sero.studies # 31 surveys
  dn2 = as.numeric(names(table(pop.agg$year))) # 111 annees
  dn3 = names(pop1)[4:length(names(pop.agg))] # 101 classes d'age
  
  for(a in 1:length(dn3)) { # dn3 = nb de classes d age =101
    for(y in min(dn2):max(dn2)) { 
      mm = match(pop.agg$adm0[pop.agg$year==y],dn1_survey) 
      pop.agg3d[mm,y-min(dn2)+1,a]=pop.agg[pop.agg$year==y,a+3]
      vc.agg3d[mm,y-min(dn2)+1,a]=vc.agg[vc.agg$year==y,a+3]
    }
  }
  dim(pop.agg3d)
  dimnames(pop.agg3d)[[1]] = dimnames(vc.agg3d)[[1]] <- dn1_survey
  dimnames(pop.agg3d)[[2]] = dimnames(vc.agg3d)[[2]]  <- dn2
  dimnames(pop.agg3d)[[3]] = dimnames(vc.agg3d)[[3]] <- dn3
  
  
  return(list(pop.agg3d = pop.agg3d, vc.agg3d=vc.agg3d))
  
}

#### incidence of vaccination aggregated by serosite
calc_inc_v3d.agg = function(vc.agg3d){

  inc_v3d.agg=rep(NA, dim(vc.agg3d)[1]*dim(vc.agg3d)[2]*dim(vc.agg3d)[3] )
  dim(inc_v3d.agg)= dim(vc.agg3d)
  dimnames(inc_v3d.agg)=dimnames(vc.agg3d)
  
  for (adm in 1:dim(vc.agg3d)[1] ) { #
    inc_v3d.agg[adm,1,]=0 # no vaccination in 1940
    
    for(y in 2:(dim(vc.agg3d)[2]-1)){
      inc_v3d.agg[adm,y,1] = vc.agg3d[adm,y+1,1] # for a0, incidence = coverage in the a0 age group
      
      for(age in 2:(dim(vc.agg3d)[3]-1)) {
        if( !is.na(vc.agg3d[adm,y,age]) & vc.agg3d[adm,y,age]==1){
          inc_v3d.agg[adm,y,age]=0
        } else { 
          inc_v3d.agg[adm,y,age]= ( vc.agg3d[adm,y+1,age+1]-vc.agg3d[adm,y,age] ) / ( 1 - vc.agg3d[adm,y,age] ) # incidence among those susceptibles at year y-1
          inc_v3d.agg[adm,y,age]=ifelse(inc_v3d.agg[adm,y,age]<10e-15,0,inc_v3d.agg[adm,y,age]) # with rounding, some values become negative
        }
      }
      inc_v3d.agg[adm,y,dim(vc.agg3d)[3]]=0 # incidence vaccination = 0 for age=100
    }
    inc_v3d.agg[adm,dim(vc.agg3d)[2],]=0# incidence vaccination = 0 for year=2050
    inc_v3d.agg[adm,dim(vc.agg3d)[2],1]= vc.agg3d[adm,dim(vc.agg3d)[2],1] # except amon newborns
  }
  return(inc_v3d.agg)
}

####### moment of the Taylor expansion

# Before, pop.moment was calculated for each year, then I selected the year of vaccine intro for Taylor expansion
# that's not relevant, especially for provinces where t0_vac is 2050
# Now, I calculate the mean pop.moments for between 1940 and t0_vac if t0_vac <2013
# else I calculate the mean pop.moment between 1984 and 2013

calc.pop.moments.agg = function(pop.agg3d){
  pop.moments.agg = rep(NA, n.serosurveys*6)
  #pop.tot.whole= rowSums(pop1[,-c(1:3)], na.rm=T)
  dim(pop.moments.agg)=c(n.serosurveys,6)
  avec=c(0:100)
  for (index_survey in 1:n.serosurveys){
    
    if (t0_vac[index_survey]<study.years[[index_survey]]){
      vec = which(dn2==1940):which(dn2==t0_vac[index_survey])
      pop.moments.agg_tmp = rep(NA, length(vec)*6)
      dim(pop.moments.agg_tmp) = c(length(vec), 6)
      for(year in vec){
        p_prop_agg = pop.agg3d[index_survey, year,]/(sum(pop.agg3d[index_survey, year,], na.rm=T))
        pop.moments.agg_tmp[year, 1] = sum( p_prop_agg, na.rm=T)
        pop.moments.agg_tmp[year, 2] = sum(p_prop_agg*avec, na.rm=T)
        pop.moments.agg_tmp[year, 3] = sum(p_prop_agg*avec*avec, na.rm=T)
        pop.moments.agg_tmp[year, 4] = sum(p_prop_agg*avec*avec*avec, na.rm=T)
        pop.moments.agg_tmp[year, 5] = sum(p_prop_agg*avec*avec*avec*avec, na.rm=T)
        pop.moments.agg_tmp[year, 6] = sum(p_prop_agg*avec*avec*avec*avec*avec, na.rm=T)
      }
      
    } else if (t0_vac[index_survey]>=2013){
      vec = which(dn2==1984):which(dn2==2013)
      pop.moments.agg_tmp = rep(NA, length(vec)*6)
      dim(pop.moments.agg_tmp) = c(length(vec), 6)
      for(year in vec){
        #print(year)
        p_prop_agg = pop.agg3d[index_survey, year,]/(sum(pop.agg3d[index_survey, year,], na.rm=T))
        pop.moments.agg_tmp[year-min(vec)+1, 1] = sum( p_prop_agg, na.rm=T)
        pop.moments.agg_tmp[year-min(vec)+1, 2] = sum(p_prop_agg*avec, na.rm=T)
        pop.moments.agg_tmp[year-min(vec)+1, 3] = sum(p_prop_agg*avec*avec, na.rm=T)
        pop.moments.agg_tmp[year-min(vec)+1, 4] = sum(p_prop_agg*avec*avec*avec, na.rm=T)
        pop.moments.agg_tmp[year-min(vec)+1, 5] = sum(p_prop_agg*avec*avec*avec*avec, na.rm=T)
        pop.moments.agg_tmp[year-min(vec)+1, 6] = sum(p_prop_agg*avec*avec*avec*avec*avec, na.rm=T)
      }
    }
    pop.moments.agg[index_survey,] = apply(pop.moments.agg_tmp,2,mean)
  }  
  return(pop.moments.agg)
}



###### some pop tables/ matrix of pop at survey sites

create.pop.at.survey = function(pop.agg3d){
  # Create p_at_survey and P_tot, population structure and number, aggregated at survey provinces
  p_at_survey=rep(NA, dim(pop.agg3d)[1]*dim(pop.agg3d)[2]*dim(pop.agg3d)[3])
  dim(p_at_survey)=dim(pop.agg3d)
  dimnames(p_at_survey) =dimnames(pop.agg3d)
  dim(p_at_survey)
  P_tot_survey=  matrix(rep(NA, dim(pop.agg3d)[1]*dim(pop.agg3d)[2]), nrow=length(dn1_survey))
  rownames(P_tot_survey) = dn1_survey
  colnames(P_tot_survey) = dn2
  for (i in 1:length(dn1_survey)){
    for (Y in  1:length(dn2) ){
      P_tot_survey[i,Y] = sum(pop.agg3d[i,Y,], na.rm=T)
      p_at_survey[i,Y,] = pop.agg3d[i,Y,]/P_tot_survey[i,Y]
    }
  }  
  p_at_survey_3d = p_at_survey
  P_tot_survey_2d = P_tot_survey

  return(list(p_at_survey_3d=p_at_survey_3d, P_tot_survey_2d=P_tot_survey_2d))
}
###### R0 lookup table

#########################################################################################
#########################################################################################
#################  Creation of Look-Up table for R0 and Ninf30 ##########################
#########################################################################################
#########################################################################################
create_R0.lookup = function(){
  dn_R0 = c(seq(1, 1.099 , by=0.001), seq(1.1, 1.198, b=0.002), seq(1.2, 1.495, by=0.005), seq(1.5, 1.99, by=0.01), 
            seq(2, 2.9, by= 0.05 ), seq(3,3.9, by = 0.1), seq(4,10, by=0.2), c(11:20), seq(25, 50, by=5))
  length(dn_R0)
  R0_lookup = rep(NA, length(dn1)*length(dn_R0))
  dim(R0_lookup) = c(length(dn1),length(dn_R0))
  colnames(R0_lookup) = dn_R0
  rownames(R0_lookup) = dat$adm0_adm1
  adm2 = 1:length(dn1)
  
  ptm = proc.time()
  for (r in 1:length(dn_R0)){
    print(paste(r,"/", length(dn_R0), ";   R0=",dn_R0[r]))
    R0_rep = rep(dn_R0[r], length(dn1))
    tmp = fun.recurrence.seroprev.whole (adm=adm2, R0=R0_rep, t0_vac_adm=t0_vac_africa)$Ninf_t_province # this is for each year from 1940 to 2050
    R0_lookup[,r] = rowSums(tmp[, which(dn2==1984):which(dn2==2013)])
  }
  save(R0_lookup, file="R0_lookup_table.Rdata")

}


