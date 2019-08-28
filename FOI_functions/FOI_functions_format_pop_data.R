
import_pop_data = function(path, c_country, pop_source, adm){
 
  pop_year_age = NULL          
  if(pop_source == "LS2014")   {
    for(adm0 in c_country){
      print(adm0)
      # tmp = read.csv(paste0("W:/Data/processed_Afpop/pop_size/from_LandScan2014/", adm, "/pop_year_age_", adm, "_", adm0, "_LandScan2014.csv"),
      #                h=T)
      tmp = read.csv(paste0(path, "/from_LandScan2014/", adm, "/pop_year_age_", adm, "_", adm0, "_LandScan2014.csv"),
                     h=T)
      pop_year_age = rbind(pop_year_age, tmp)
    }
  } 
  
  return(pop_year_age)
  
}



import_pop_data_for_FOI = function(path, c_country, pop_source, adm){
  pop1 = import_pop_data(path, c_country, pop_source, adm)
  pop1 = add_1940_1950(pop1);dim(pop1)
  pop1$adm0_adm1<-paste(pop1$adm0,pop1$adm1, sep="_")
  pop1 = pop1[,c(1,2, ncol(pop1), 3:(ncol(pop1)-1))]
  pop1 = pop1[,-c(2)] # retrieve adm1, information is in  adm0_adm1
  pop1 = pop1[,-c(1)]
  pop1 = pop1[order(pop1$adm0_adm1),]
  
  return(pop1)
}


#############################
add_1940_1950 = function(pop2d){
  if(min(pop2d$year) != 1950)  stop("pop2d should start in 1950")
  
  for(y in 1949:1940) {
    pop1.early = pop2d[pop2d$year==y+1,]
    pop1.early$year = y
    pop2d = rbind(pop1.early,pop2d)
  }
  
  return(pop2d)
}


transform_into_pop3d = function(pop2d, adm){
 
  # just transform adm1 in adm0_adm1 /adm2 in adm0_adm2
  temp =  paste0(pop2d[,"adm0"], "_", pop2d[,adm] )
  names(pop2d)[names(pop2d) == adm ] =  paste0("adm0_", adm)
  pop2d[, paste0("adm0_", adm)] = temp

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

get_pop_data_3d = function(path, c_country, pop_source, adm){
  if (!pop_source %in% c("LS2011", "LS2014", "WP2010")) {
    stop("pop source has to be LS2011, LS2014, WP2010")
  }
  if (!adm %in% c("adm1", "adm2")) {
    stop("adm source has to be adm1 or adm2")}
  
  
  pop2d = import_pop_data(path, c_country, pop_source, adm) # check here if ... is sufficient or if need to tell path=path
  
  print("add 1940-1950")
  pop2d = add_1940_1950(pop2d)
  
  print("build pop3d, P_tot_2d, p_prop3d")
  pop3d = transform_into_pop3d(pop2d=pop2d, adm)
  P_tot_2d = get_P_tot_2d(pop3d=pop3d)
  p_prop3d = get_p_prop3d(pop3d=pop3d, P_tot_2d=P_tot_2d)
  
  res = list(pop3d= pop3d, P_tot_2d=P_tot_2d, p_prop3d=p_prop3d)
  
  return(res)
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
# vc data
repair_vc_data = function(vc2d){ # before 1995, we have NA values for those aged >75. We affect the value of the latest age group
  for (a in 3:ncol(vc2d)){
    vc2d[,a] = ifelse(is.na(vc2d[,a]), vc2d[,a-1], vc2d[,a])
  }
  return(vc2d)
}



transform_into_vc3d = function(vc2d, adm=NA){
  
  # just transform adm1 in adm0_adm1 /adm2 in adm0_adm2
  # temp =  paste0(vc2d[,"adm0"], "_", vc2d[,adm] )
  # names(vc2d)[names(vc2d) == adm ] =  paste0("adm0_", adm)
  # vc2d[, paste0("adm0_", adm)] = temp
  if(is.na(adm)) adm="adm1"
  adm0_adm = paste0("adm0_", adm)
  
  vc3d=rep(NA, nrow(vc2d)*(ncol(vc2d)-2))
  dim(vc3d)=c(length(table(vc2d[,adm0_adm])), length(table(vc2d$year)), ncol(vc2d)-2)
  # 1st dim = adm0_adm1
  # 2nd dim = year
  # 3rd dim = age
  
  dn1 = names(table(vc2d[,adm0_adm]))
  dn2 = as.numeric(names(table(vc2d$year))) # 111 annees
  dn3 = names(vc2d)[3:length(names(vc2d))] # 101 classes d'age
  
  for(a in 1:length(dn3)) { # dn3 = nb de classes d age =101
    for(y in min(dn2):max(dn2)) { 
      mm = match(vc2d[,adm0_adm][vc2d$year==y],dn1) 
      vc3d[mm,y-min(dn2)+1,a]=vc2d[vc2d$year==y,a+2]
    }
  }
  
  dimnames(vc3d)[[1]] <- dn1
  dimnames(vc3d)[[2]] <- dn2
  dimnames(vc3d)[[3]] <- dn3
  
  return(vc3d)
  
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
    print(i)
    year_i = 1940
    sum_vac=0
    while(sum_vac == 0 & year_i<=2050) {
      sum_vac = sum( vc3d[i, year_i-1940+1,], na.rm=T)
      year_i = year_i + 1
    }
    t0_vac_africa[i]=year_i-1
  }
  names(t0_vac_africa)=dimnames(vc3d)[[1]]
  return(t0_vac_africa)
}




############################
# create pop30.agg
create_pop30.agg_vc30.agg = function(pop1, vc1){
    pop30 = pop1[pop1$year>1983 & pop1$year<2014,] # need to change the years here
    pop30 = pop30[order(pop30$adm0_adm1),]
    
    vc30 = vc1[vc1$year>1983 & vc1$year<2014,]
    
    pop30[is.na(pop30)] = 0
    pop30[is.na(pop30) | is.na(vc30)] = 0 
    vc30[is.na(vc30)] = 0
    
    pop30.agg = aggregate(pop30[,3],by=list(adm0_adm1=pop30$adm0_adm1),sum) # a0 is the 4th column, 
    # pop by age aggregated over 30 years
    vc30.agg = aggregate(vc30[,3]*pop30[,3],by=list(adm0_adm1=vc30$adm0_adm1),sum) # number of vaccinated people by age, aggregated obver 30y
    vc30.agg$x = vc30.agg$x/pop30.agg$x
    
    
    names(pop30.agg)[2] = names(vc30.agg)[2] = names(pop30)[3]
    
    for(i in 3:ncol(pop30)) { # then follow to a1 and others
      pop30.agg = cbind(pop30.agg, aggregate(pop30[,i],by=list(adm0_adm1=pop30$adm0_adm1),sum)$x)
      
      vc30.agg = cbind(vc30.agg, aggregate(vc30[,i]*pop30[,i],by=list(adm0_adm1=vc30$adm0_adm1),sum)$x/aggregate(pop30[,i],by=list(adm0_adm1=pop30$adm0_adm1),sum)$x)
      names(pop30.agg)[ncol(pop30.agg)] = names(vc30.agg)[ncol(vc30.agg)] = names(pop30)[i]
    }
    vc30.agg[is.na(vc30.agg)] = 0
    
    amat = matrix(0:100,nrow=nrow(pop30.agg),ncol=101,byrow=T)
    pop.vc.moments = data.frame(adm0_adm1 = pop30.agg$adm0_adm1, 
                                m0 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1]),na.rm=T),
                                m1 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat,na.rm=T), 
                                m2 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat*amat,na.rm=T),
                                m3 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat*amat*amat,na.rm=T),
                                m4 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat*amat*amat*amat,na.rm=T),
                                m5 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat*amat*amat*amat*amat,na.rm=T)) 
    
    return(list(pop30.agg=pop30.agg, vc30.agg=vc30.agg, pop.vc.moments=pop.vc.moments))
}



##### 
create_pop_and_vc.agg.ag = function(pop1, vc1, n.serosurveys, adm1s){
  
  pop_sero = pop1[pop1$adm0_adm1 %in% unlist(adm1s),] 
  vc1_sero = vc1[vc1$adm0_adm1 %in% unlist(adm1s),]
  
  for(i in 1:n.serosurveys) {
    pop_sero = pop_sero[!(pop_sero$adm0_adm1 %in% adm1s[[i]] & !pop_sero$year %in% study.years[[i]]),]
    vc1_sero = vc1_sero[!(vc1_sero$adm0_adm1 %in% adm1s[[i]] & !vc1_sero$year %in% study.years[[i]]),]
  }
  
  
  pop.agg = vc.agg = NULL
  for(i in 1:n.serosurveys) {
    pop.tmp = pop_sero[pop_sero$year %in% study.years[[i]] & pop_sero$adm0_adm1 %in% adm1s[[i]],]
    vc.tmp = vc1_sero[vc1_sero$year %in% study.years[[i]] & vc1_sero$adm0_adm1 %in% adm1s[[i]],]
    if (nrow(pop.tmp)==1) {
      pop.tmp$adm0 = vc.tmp$adm0 = sero.studies[i]
      pop.tmp = pop.tmp[,c(ncol(pop.tmp), 1:(ncol(pop.tmp)-1))] # just need to reput adm0 in the first place
      vc.tmp = vc.tmp[,c(ncol(vc.tmp), 1:(ncol(vc.tmp)-1))]
    } else if(nrow(pop.tmp)>1) { # need to do some summing/averaging...
      pop.tmp[is.na(vc.tmp)] = NA
      vc.tmp = colSums(pop.tmp[,-(1:2)]*vc.tmp[,-(1:2)],na.rm=T)/colSums(pop.tmp[,-(1:2)],na.rm=T) 
      dim(vc.tmp) = c(1,length(vc.tmp))
      colnames(vc.tmp) = paste("a",0:100,sep="")
      vc.tmp = data.frame(adm0 = sero.studies[i], adm0_adm1 = NA, year = mean(study.years[[i]]),vc.tmp)
      
      pop.tmp = colSums(pop.tmp[,-(1:2)])/length(study.years[[i]])
      dim(pop.tmp) = c(1,length(pop.tmp))
      colnames(pop.tmp) = paste("a",0:100,sep="")
      pop.tmp = data.frame(adm0 = sero.studies[i], adm0_adm1 = NA, year = mean(study.years[[i]]),pop.tmp)
      
    }
    pop.agg = rbind(pop.agg,pop.tmp)
    vc.agg = rbind(vc.agg,vc.tmp)
  }
  pop.agg[,-(1:3)][is.na(pop.agg[,-(1:3)])] = 0
  vc.agg[,-(1:3)][is.na(vc.agg[,-(1:3)])] = 0
  
  vc.agg.ag = NULL
  for(i in 1:n.serosurveys) { # averages vaccin coverage over the age groups
    ag = findInterval(0:100,age.min[[i]])
    vc.agg.ag[[i]] = aggregate(c(as.matrix(vc.agg[i,-(1:3)]*pop.agg[i,-(1:3)])),by=list(ag=ag),sum)$x/aggregate(c(as.matrix(pop.agg[i,-(1:3)])),by=list(ag=ag),sum)$x 
  }
  
  
  return(list(pop.agg=pop.agg, vc.agg=vc.agg, vc.agg.ag=vc.agg.ag))
  
}