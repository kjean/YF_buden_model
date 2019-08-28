
##### launch environmental data

launch_env_dat = function(env_table, delete_surv_AGO=T){
  
  dat.full = read.csv(env_table, stringsAsFactors=F)
  
  depvar = "cas.or.out" # which presence/abence outcome to consider. Alternatives : "cases" # "outbreaks" # 
  depi = match(depvar,names(dat.full)) # column number for the chosen outcome
  
  # excluding LC1,3,5,15 as they only ever cover <5% - doesn't make sense to use these...
  ex.i = match(c("LC1","LC3","LC5","LC15"),names(dat.full))
  dat.full = dat.full[,-ex.i]
  
  # adding a categorical "dominant land cover" variable:
  LC.i = grep("LC",names(dat.full))
  LC.dom = names(dat.full)[LC.i][apply(dat.full[,LC.i],1,which.max)]
  dat.full = cbind(dat.full, LC.dom) # 37 potential covariates/levels
  
  
  ###############
  dat.full$surv.qual.adm0[is.na(dat.full$surv.qual.adm0)] = 0
  dat.full$surv.qual.adm1[is.na(dat.full$surv.qual.adm1)] = 0
  
  if(delete_surv_AGO==T){# setting the surv.qual.adm0 for AGO to 0 (very few cases reported): 
    dat.full$surv.qual.adm0[dat.full$adm0=="AGO"] = 0 
  }
  
  #adding log(surv.qual.adm0)
  dat.full = cbind(dat.full, log.surv.qual.adm0 = log(dat.full$surv.qual.adm0)) # log.surv better fits data than surv
  dat.full$log.surv.qual.adm0[is.infinite(dat.full$log.surv.qual.adm0)] = 0  # A way not to apply beta(log.surveillance) to coutries uotside of YFSD
  
  adm05 = dat.full$adm0
  adm05[dat.full$surv.qual.adm0>0] = "AFR" # a same categorical variable for all countries within the YFSD
  dat.full = cbind(dat.full, adm05=adm05) # for countries outside YFSD, 1 categorical variable per country
  
  dat = dat.full[dat.full$adm0 %in% c34,]
  dat$adm05 = as.factor(as.character(dat$adm05))
  #table(dat$adm05)
  
  #summary(dat)
  #str(dat)
  v1 = apply(dat,2,var)
  for(i in 8:(ncol(dat))) {
    if(!is.factor(dat[,i]) & !is.character(dat[,i])) {
      dat[,i] = dat[,i]/sqrt(v1[i])
      dat.full[,i] = dat.full[,i]/sqrt(v1[i])  
    } 
    # Explanation:
    # If we fit the model on dat and if we want to project estimates on dat.full, variabble from dat.full need to be expressed on the same scale than those from dat , thus we normalize dat.full relatively to dat
    
  }

  # dat = dat[,names(dat)!="adm0"] ligne enlevee au 04/09/15, a remettre si pb
  dat = dat[,names(dat)!="surv.qual.adm0"]
  depi = match(depvar,names(dat))
  
  
  dim(dat)
  class(dat$adm0_adm1)
  
  # I do that because pop is ordered that way and we need to match both pop and dat
  dat = dat[ order(dat$adm0_adm1), ]
  
  return(list(depi = depi, dat = dat))
}


######## 
## fit GLM 

################################################################################
################################################################################

# Preparing the models to fit

################################################################################
################################################################################

# fm.best = "cas.or.out ~ lon+logpop+EVI.max+surv.qual.adm1"
# replacing the surveillance quality with a country factor:
# the first 10 models are the best 10 models according to bic,
# but with surv.qual.adm1 replaced by adm0,
# models 11 to 20 are the best 10 models that include adm0,
# but not surv.qual.adm1 or surv.qual.adm0

fit_glm = function(dat, depi){
  fm.list = c(
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+logpop+LC14+EVI.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+logpop+EVI.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+logpop+MIR.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+lat+logpop+LC14+EVI.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+logpop+LC7",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+lat+logpop+MIR.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+logpop+LC14+MIR.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lat+logpop+EVI.mean+MIR.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lat+logpop+EVI.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+lat+logpop+LC16+MIR.mean",
    
    "cas.or.out ~log.surv.qual.adm0+adm05+ lat+logpop+LC14+EVI.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lat+logpop+MIR.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+logpop+LC4+LC14+EVI.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lat+logpop+LC16+MIR.mean",
    "cas.or.out ~log.surv.qual.adm0+adm05+ lon+logpop+LC7+LC14+EVI.mean"
  )
  
  model=1
  model = as.numeric(model)
  fm.best = fm.list[model]
  
  bm = glm(as.formula(fm.best), data=dat, family=binomial(link="cloglog")) #bm = best model
  
  # setting up the evaluation of the likelihood:
  beta = coefficients(bm)
  vl=NULL
  for(i in 1:ncol(dat)) if(length(grep(names(dat)[i], names(beta)))>0) vl = c(vl,i) # select the variables used in the GLM
  x = cbind(Intercept=1,dat[,vl])
  j.expand = !sapply(1:ncol(x), function(i) is.numeric(x[,i]) & !is.factor(x[,i]))
  x.num = NULL
  for(j in 1:ncol(x)) {   # create indicative variables for adm05
    if(j.expand[j]==T) {
      tab = table(1:nrow(x),x[,j])
      colnames(tab) = paste(names(x)[j],colnames(tab),sep="")
      x.num = cbind(x.num,tab[,-1])
    } else {
      x.num = cbind(x.num,x[,j])
      colnames(x.num)[ncol(x.num)] = names(x)[j]
    }
  }
  x = x.num # covariate matrix
  rm(x.num)
  y = dat[,depi]  # dependant variable vector
  
  
  beta0 = coefficients(bm)
  names(beta0)[names(beta0)=="(Intercept)"] = "Intercept"
  
  mm = match(names(beta0),colnames(x))
  x = x[,mm]
  
  beta0 = beta0[match(colnames(x),names(beta0))]
  nn = names(beta0)

 return(list(beta0=beta0, x=x, y=y))
   
}

