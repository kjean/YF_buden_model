##################################333
# Nouvelles fonctions YF Kevin





################################################################################
# Calculating lnL of glm regression:
fun.glm.lnL = function(beta, x, y) {
  # beta are the model coefficients,
  # x the independent covariates (needs one column for the Intercept),
  # and y the binary outcome.
  # model predictions pi = p, 1-pi = q
  eta = as.numeric(x %*% beta)
  logq = -exp(eta) # -exp(X*beta) (= log(1-q) )
  logp = log(1-exp(logq)) # 
  #q = exp(logq)
  #p = 1-q
  logl = sum(logp[y==1]) + sum(logq[y==0])
  return (logl)
}


################################################################################
################################################################################

# calculating my own model predictions from the generalised linear regression
# model with family=binomial(link="cloglog"):
# une fois les coefs de GLM fixes, calcule la proba predite de report de YF

fun.calcPred = function(coefs, newdata, type="response",varsin=NA) {
  if(!(type %in% c("link","response"))) stop("fun.calcPred: invalid type of predictions specified.\n")
  if(is.na(varsin[1])) varsin = 1:length(coefs)
  eta = newdata[,varsin] %*% coefs[varsin]
  if(type=="link") {
    preds = eta #X*beta
  } else if(type=="response") {
    preds = 1-exp(-exp(eta)) # q = 1 - exp (- exp(X*beta))
  }
  return(preds)
}

################################################################################
# fun.lnL: for estimating the foi from serosurveys

# First calculate expected seroprevalence in each age group for a given FOI
fun.calcSeroprev.ag = function(foi,age.groups,pop,vc=0) {
  if(length(vc)==1) {
    vc = rep(vc,length(age.groups))
  } else {
    if(length(vc)!=length(age.groups)) {
      stop("fun.calcSeroprev.ag: incompatible length of vaccination coverage vc.\n")
    }
  }
  if(!all(!is.na(vc)) |min(vc)<0 |max(vc)>1){
    stop("fun.calcSeroprev: invalid values for vaccination coverage vc.\n")
  }
  sero1 = c(as.matrix((1-exp(-foi*0:100))*pop[,-(1:3)])) # the first 3 columns are  year, adm0, adm0_adm1
  ag = findInterval(0:100,age.groups)
  sero.ag = aggregate(sero1,by=list(age.group=ag),sum)
  sero.ag$x = sero.ag$x/aggregate(c(as.matrix(pop[,-(1:3)])),by=list(age.group=ag),sum)$x
  # proportion immune due to either vaccination or previous infection:
  sero.out = 1 - (1-sero.ag$x[sero.ag$age.group>0])*(1-vc)
 
  return(sero.out)
}

# then calculate the lnL for a given serosurvey and a given FOI
fun.lnL = function(foi,obs.tot,obs.pos,age.groups,pop,vc=0) {
  if(length(obs.pos)!=length(obs.tot) | length(obs.pos) != length(age.groups)){
    stop("fun.lnL: incompatible dimensions.\n")
  }
    serop = fun.calcSeroprev.ag(foi,age.groups,pop,vc)
  
    # log_serop and log_1_minus_serop lines added on 6 Oct 2015, to be retrieved if problems.
    log_serop = ifelse(serop==0, -50, log(serop)) # done to avoid Infinite lnL
    log_1_minus_serop = ifelse(serop==1, -50, log(1-serop)) # done to avoid Infinite lnL
    lnL = sum(lgamma(obs.tot+1)-lgamma(obs.pos+1)-lgamma(obs.tot-obs.pos+1) + obs.pos*log_serop + (obs.tot-obs.pos)*log_1_minus_serop)
  # pourquoi pas sous la forme simple d'une lnL binomiale? est-ce equivalent?
    return(lnL)
}
# that the same than log of binomial likelihood



################################################################################
# Calculating lnL of of the prior: (to be validated by Tini)
fun.lnL.prior = function(params, sd.prior, transform){
  if(length(sd.prior)!= length(params) | length(transform)!=length(params))
    stop("fun.lnL.prior: invalid sd.prior input.\n")
  
  lnLPrior = rep(NA, 5) # there are 5 different cases
  
  jj = grep("^log.adm05",names(status$params)) # select country parameters
  lnLPrior[1] = -0.5*sum( (params[jj]/sd.prior[jj])^2 ) # that s the first cases, for country parameters (mu.prior=0 for country factors)
  
  lnLPrior[2] = sum( dnorm(params[transform == 0 & grepl("^log.adm05",names(status$params))==F], mean=0, sd=30, log=TRUE))
  # this term is for normaly distributed non-country parameters : normal distrib with high sd (flat distrib)
  
  lnLPrior[3] = sum( dexp(params[transform == 1], rate=0.001,log=TRUE))
  # this term is for parameters distributed over [0,+Inf], ie FOI: quasi flat distrib
  
  lnLPrior[4] = 0 # here, we use an uniform prior which is constant, thus no need to include it, it will cancel we comparing to the previous MCMC step
  #lnLPrior[4] = sum( dexp(params[transform == 2], rate=0.001,log=TRUE))
  # this term is for parameters distributed over [0,1], ie vc.factor : quasi flat distrib
  
  lnLPrior[5] = dlogis(log(params[transform == 3]/(1-params[transform == 3])), location =3.999, scale = 0.128*sqrt(3)*pi, log=T)
  
  return(sum(lnLPrior))
}
  



#################################################################################
# function for transition in the MCMC


fun.mcmc.step = function(status, sdj, mu.prior, sd.prior, transform, index=NA, varsin.nc) {
  # index indexes the parameters  
  #print(index)
  status.new = status
  
  if(is.na(index[1])) {
    stop("in proposal: invalid index")
    } 
  
  ############
  # Proposals calculation
  ############
  if(transform[index]==0) { ## use a normal prior
    status.new$params[index] =  rnorm(1,mean=status$params[index],sd=sdj[index]) # normal distrib for glm covar
  
  } else if (transform[index]==1) { ## use a log normal prior
    status.new$params[index] = rlnorm(1,meanlog=log(status$params[index]),sdlog=sdj[index])  # lognormal distrib for FOI (varying from 0 to Inf)
    
  } else if (transform[index]==2 | transform[index]==3) {## use a truncated log normal prior
    status.new$params[index] = rlnormTrunc(1, meanlog= log(status$params[index]) , sdlog=sdj[index], min=0, max=1) 
    # logistic distrib for vc.factor.CMRs and vaccine efficacy(varying from 0 to 1)
        
  } else { stop("in proposal: invalid transform")    }
  
  # les 4 lignes suivantes sont a faire sauter en cas de pb
  # if (index == 18) { tmp = (status$params[index] + status$params[1]) / status$params[1]
  # tmp.new = rnorm(1,mean=tmp,sd=sdj[index])
  # status.new$params[index] = status$params[1] * (tmp.new - 1)
  #}
  
    
  ######
  # Update Likelyhood calc
  #ii = grep("^log.",names(status$params)) # ii indexes the GLM parameters
  if (index == 1){
    #change vacc coverage
    vac_eff = status.new$params[1]
    mm = grep("foi", names(status.new$params))
    lnL_tmp = NULL
    for (i in 1:length(mm)) {
      if (vc.factor[i]==1 | is.na(vc.factor[i])){ # need to recalculate lnL only for surveys including vaccinted people
        foi = status.new$params[mm[i]]
        vcfac = vc.factor[i]
        if(is.na(vcfac)) vcfac = status.new$params[which(names(status.new$params)==paste("vc.factor.",sero.studies[i],sep=""))] # if vc is unknown, it is fitted with MCMC
        lnL_tmp = c(lnL_tmp, fun.lnL(foi=foi,obs.tot=sero.dat[[i]][,1], obs.pos=sero.dat[[i]][,2], age.groups=age.min[[i]], pop=pop.agg[i,], 
                                      vc=vac_eff*vcfac*vc.agg.ag[[i]]) )
      }
    }
    status.new$lnL[1] = sum(lnL_tmp)
    
  } else if (index %in%  ii) { # if the index parameter is a GLM parameter, update GLM likelyhood (otherwise, not)
    # select GLM parameters
    status.new$lnL[2] = fun.glm.lnL(status.new$params[ii], x, y) # lnL  GLM
    
  } else if (max(ii)<index & index<=(max(ii)+n.serosurveys) ){ # if the index parameter is one of the FOI parameter but not the CMRs one
    i = index - max(ii)
    foi = status.new$params[max(ii)+i] # takes the FOI for the study i
    vac_eff = status.new$params[1]
    vcfac = vc.factor[i]
    if(is.na(vcfac)) vcfac = status.new$params[which(names(status.new$params)==paste("vc.factor.",sero.studies[i],sep=""))] # if vc is unknown, it is fitted with MCMC
    status.new$lnL[i+2] = fun.lnL(foi=foi,obs.tot=sero.dat[[i]][,1], obs.pos=sero.dat[[i]][,2], age.groups=age.min[[i]], pop=pop.agg[i,], 
                                  vc=vac_eff*vcfac*vc.agg.ag[[i]])
   
    
    # theline below has been changed on 5th Nov 2015 in oredr to allow leave-one-out validation
  #} else if (index==length(status$params) ) { # fitting FOI for CMRs or vc.factor.CMRs
  } else if  ( grepl("vc.factor", names(status.new$params)[index]) )  { # fitting FOI for CMRs or vc.factor.CMRs 
    i = which(is.na(vc.factor))  # WARNING: this works only if we have 1 vc.factor  = NA !
    
    foi = status.new$params[max(ii)+which(is.na(vc.factor))]
    #vcfac = status.new$params[which(names(status.new$params)==paste("vc.factor.",sero.studies[i],sep=""))]
    vcfac = status.new$params[index]
    vac_eff = status.new$params[1]
    status.new$lnL[i+2] = fun.lnL(foi=foi,obs.tot=sero.dat[[i]][,1], obs.pos=sero.dat[[i]][,2], age.groups=age.min[[i]], pop=pop.agg[i,], 
                                  vc=vac_eff*vcfac*vc.agg.ag[[i]])
  }
    

    # priors likelyhood
  status.new$lnL[length(status.new$lnL)] = fun.lnL.prior(status.new$params, sd.prior, transform)
  
  
  # accept/reject
  ratio.post = sum(status.new$lnL, na.rm=T) - sum(status$lnL, na.rm=T)
  # 
  if(transform[index]==0) {
    correction = 0 
  } else if (transform[index]==1) {
    correction = log(status.new$params[index]) - log(status$params[index]) # correction for lognormal distribution
  } else if(transform[index]==2 | transform[index]==3) { # correction factor for Trunc lognormal
    tmp1 = - log (status$params[index]) / sdj[index]
    tmp2 = - log (status.new$params[index]) / sdj[index]
    correction = log(status.new$params[index]) + log(pnorm(tmp1)) - log(status$params[index]) - log(pnorm(tmp2)) # correction for truncated lognormal distribution
    if(is.nan(correction)) correction = 0
  } 
  
  if(is.finite(ratio.post)) {
    p.accept = ratio.post + correction 
  } else {p.accept = status.new$lnL}
  
  ## accept/reject step:
  tmp = log(runif(1))
  if(tmp<p.accept) { # accept:
    status.new$accept=1
  } else { # reject:
    status.new = status
    status.new$accept = 0
  }
  
  if( index == 1 & status.new$accept ==1){
    vac_eff = status.new$params[1]
    mm = grep("foi", names(status.new$params))
    for (i in 1:length(mm)) {
      if (vc.factor[i]==1 | is.na(vc.factor[i])){
        foi = status.new$params[mm[i]]
        vcfac = vc.factor[i]
        if(is.na(vcfac)) vcfac = status.new$params[which(names(status.new$params)==paste("vc.factor.",sero.studies[i],sep=""))] # if vc is unknown, it is fitted with MCMC
        status.new$lnL[i+2] = fun.lnL(foi=foi,obs.tot=sero.dat[[i]][,1], obs.pos=sero.dat[[i]][,2], age.groups=age.min[[i]], pop=pop.agg[i,], 
                                      vc=vac_eff*vcfac*vc.agg.ag[[i]])
      }
    } 
  }
  
  return(status.new)
}




####################################
# Functions still to adapt

#############
# A chaque iteration de MCMC, on doit calculer un burden...
# c'est comme ca seulement qu'on arrive a des CI 95%


################################################################################
################################################################################
fun.calcFoi.Africa = function(status,polydeg=5) {
  # pourquoi sous la forme d'un polygone?
  
  mypreds.nc = fun.calcPred(coefs = status$params[ii],newdata=x,type="link",varsin=varsin.nc)
  p.detect.link = fun.calc.pdetect(params= status$params,mypreds.nc=mypreds.nc)
  
  # using pop.vc.moments as a global variable
  if(polydeg>ncol(pop.vc.moments)) error("fun.calcFoi.Africa: invalid value for polydeg.\n")
  # c'est la solution numerique de l equation 7: Tini la resout par developpement de Taylor en prenant les 5 premiers degres..
  
  #  z = log(1-status$mypreds.nc)*exp(-status$params[which(names(status$params)=="p.detect")]) # problem if mypreds.nc==1 on the response scale
  #z = -exp(status$mypreds.nc)*exp(-status$p.detect.link) # taking into account thtat mypreds are now stored on the link scale.
  z = -exp(mypreds.nc)*exp(-p.detect.link)
  
  if(polydeg>0) for(i in 1:polydeg) {
    z = cbind(z,(-1)^(i+1)*pop.vc.moments[,i+1]/factorial(i-1))
  }
  out = sapply(1:nrow(x), function(i) polyroot(z[i,]))
  out[abs(Arg(out))<=1e-10] = Re(out)[abs(Arg(out))<=1e-10]
  out[abs(Arg(out))>1e-10] = NA
  dt = dim(out)
  out = as.numeric(out)
  dim(out) = dt
  out = apply(out,2,min,na.rm=T)
  return(out)
}
################################################################################

fun.calcBurden = function(fois,years,adm1 = NA, vac_eff_arg = 1, seed = NA, total_run, current_run) {
  # using vc.full and pop.full as global variables.
  # using life0 as global variable.
  # total_run and current_run are just here to be sure that we use the same CFR with the same R0 estimates btw different scenarios
  # total_run : nb of posterior samples generated
  # current_run = index of the current sample we're estimating
  
  if(is.na(adm1[1])) adm1 = pop.full$adm0_adm1[pop.full$year==years[1]]
  
  if (!is.na(seed)) {set.seed(seed) 
  } else {set.seed(101) 
  }
  prop.severe_all = rlnorm(total_run,meanlog=-2.222, sdlog=0.427)
  prop.severe_all[prop.severe_all>1]=1
  prop.death.all.cases_all = rlnorm(total_run, meanlog= -3.00, sdlog=0.442)
  prop.death.all.cases_all[prop.death.all.cases_all>1]=1
  
  prop.severe = prop.severe_all[current_run]
  prop.death.all.cases = prop.death.all.cases_all[current_run]
  
  # here we have in ~10% of the cases prop.death.all.cases > prop.severe
  # this is not a technical problem, but may be a theoretical one
  
  cases = severe = deaths = life.years.left = life.years.left_central =NULL
  for(year in years){
    mm = match(adm1,pop.full$adm0_adm1[pop.full$year==year])
    foi.mat = matrix(fois,nrow=length(fois),ncol=101,byrow=F)
    a.mat = matrix(0:100,nrow=length(fois),ncol=101,byrow=T)
    burden.by.age = foi.mat*exp(-foi.mat*a.mat)*(1-vac_eff_arg*vc.full[vc.full$year==year,-(1:2)])*pop.full[pop.full$year==year,-(1:2)]
    
    # calculating the remaining life expectancy for all infecteds - remember to scale this down to deaths later!
    adm0 = substr(pop.full$adm0_adm1[pop.full$year==year][mm],1,3)
    life1 = life0[life0$year==year,]
    mm.adm0 = match(adm0,life1$country)
    life1 = as.matrix(life1[mm.adm0,-(1:2)])
    
    cases_tmp = rowSums(burden.by.age,na.rm=T)
    cases = cbind(cases, cases_tmp)
    
    severe_tmp = prop.severe*cases_tmp
    severe = cbind(severe, severe_tmp)
    
    deaths_tmp = prop.death.all.cases*cases_tmp
    deaths = cbind(deaths, deaths_tmp)
    
    
    YLL = prop.death.all.cases*rowSums(burden.by.age*life1,na.rm=TRUE)
    # times and disability weights are those from LaBeaud et al
    if(prop.death.all.cases>prop.severe)prop.death.all.cases = prop.severe # otherwise i can have negative value for (prop.severe-prop.death.all.cases)
    YLD = 17.8/365.25*0.172 *prop.severe*rowSums(burden.by.age*life1,na.rm=TRUE) + #acute phase among all severe cases
      28/365.25*0.024*(prop.severe-prop.death.all.cases)*rowSums(burden.by.age*life1,na.rm=TRUE)
    
    life.years.left_tmp = YLL + YLD
    life.years.left = cbind(life.years.left, life.years.left_tmp)
    
    # YLL_central = 0.0548*rowSums(burden.by.age*life1,na.rm=TRUE)
    # YLD_central = 17.8/365.25*0.172 *0.118*rowSums(burden.by.age*life1,na.rm=TRUE) + #acute phase among all severe cases
    #   28/365.25*0.024*(0.118-0.0548)*rowSums(burden.by.age*life1,na.rm=TRUE)
    # life.years.left_central_tmp = YLL_central + YLD_central # this is for GAVI central outputs - need a constant fraction of deaths
    # life.years.left_central = cbind(life.years.left_central, life.years.left_central_tmp)
  }
  
  colnames( cases) = colnames(severe ) = colnames(deaths ) = colnames( life.years.left) =  years
  
  
  return(list(cases=cases,severe=severe,deaths=deaths,life.years.left=life.years.left))
}
################################################################################
fun.aggBurden = function(fois, burden.years, vac_eff_arg=1, total_run, current_run) {
  burden.adm1 = fun.calcBurden(fois,burden.years, vac_eff_arg=vac_eff_arg,
                               total_run = total_run, current_run=current_run)
  
  cases = aggregate(burden.adm1$cases,by=list(adm0 = dat$adm0),sum)
  severe = aggregate(burden.adm1$severe,by=list(adm0 = dat$adm0),sum)
  deaths = aggregate(burden.adm1$deaths,by=list(adm0 = dat$adm0),sum)
  lyleft = aggregate(burden.adm1$life.years.left,by=list(adm0 = dat$adm0),sum)
 # lyleft_central = aggregate(burden.adm1$life.years.left_central,by=list(adm0 = dat$adm0),sum)
  #names(cases)[2] = names(severe)[2] =names(deaths)[2] =names(lyleft)[2] = names(lyleft_central)[2] =paste("b.",burden.years[1],sep="")
  
  return(list(cases=cases,severe.cases=severe, deaths=deaths, life.years.left=lyleft))
}
################################################################################
fun.calc.pdetect = function(params,mypreds.nc){#,mypreds,mypreds.nc) {
  # calcul le terme b d apres la partie gauche de l'equation 6
  
  # using adm1s, n.serosurveys, sero.studies, dat , vc30.agg, pop30.agg as global variables
  mm = match(unlist(adm1s),dat$adm0_adm1) # attention, avec la nouvelle variable adm0_adm1, il faut que ce soit bien trie
  
  vec.vc.factor = vec.foi = NULL
  for(i in 1:n.serosurveys) {
    vec.foi = c(vec.foi,rep(params[which(names(params)==paste("foi.",sero.studies[i],sep=""))],length(adm1s[[i]])))
    if(!is.na(vc.factor[i])) {
      vec.vc.factor = c(vec.vc.factor, rep(vc.factor[i],length(adm1s[[i]])))
    } else {
      vec.vc.factor = c(vec.vc.factor, rep(params[which(names(params)==paste("vc.factor.",sero.studies[i],sep=""))],length(adm1s[[i]])))
    }
  }
  vec.foi = vec.foi[match(vc30.agg[mm,1],unlist(adm1s))]
  vec.vc.factor = vec.vc.factor[match(vc30.agg[mm,1],unlist(adm1s))]
  
  vac_eff = params[1] #vaccine efficacy
  popvc = (1-vac_eff*vec.vc.factor*vc30.agg[mm,-(1:2)])*pop30.agg[mm,-(1:2)]
  expfoi = exp(outer(-vec.foi, 0:100))
  ninf.30 = vec.foi*rowSums(expfoi*popvc)
  
  b = mypreds.nc[mm]-log(ninf.30)
  return(mean(b))
}

##########################################################################
# calculating factor b - same but returns the vector un-averaged
fun.calc.pdetect.multi = function(params,mypreds.nc){#,mypreds,mypreds.nc) {
  # calcul le terme b d apres la partie gauche de l'equation 6
  
  # using adm1s, n.serosurveys, sero.studies, dat , vc30.agg, pop30.agg as global variables
  mm = match(unlist(adm1s),dat$adm0_adm1) # attention, avec la nouvelle variable adm0_adm1, il faut que ce soit bien trie
  
  vec.vc.factor = vec.foi = NULL
  for(i in 1:n.serosurveys) {
    vec.foi = c(vec.foi,rep(params[which(names(params)==paste("foi.",sero.studies[i],sep=""))],length(adm1s[[i]])))
    if(!is.na(vc.factor[i])) {
      vec.vc.factor = c(vec.vc.factor, rep(vc.factor[i],length(adm1s[[i]])))
    } else {
      vec.vc.factor = c(vec.vc.factor, rep(params[which(names(params)==paste("vc.factor.",sero.studies[i],sep=""))],length(adm1s[[i]])))
    }
  }
 # vec.foi = vec.foi[match(vc30.agg[mm,1],unlist(adm1s))]  # ENLEVEES AU 08 sept
 # vec.vc.factor = vec.vc.factor[match(vc30.agg[mm,1],unlist(adm1s))]
 vec_eff = params[1]
  popvc = (1-vec_eff*vec.vc.factor*vc30.agg[mm,-(1:2)])*pop30.agg[mm,-(1:2)]
  expfoi = exp(outer(-vec.foi, 0:100))
  ninf.30 = vec.foi*rowSums(expfoi*popvc)
  
  b = mypreds.nc[mm]-log(ninf.30)
  names(b) = unlist(adm1s)
  return(b)
}





######## main functions for MCMC


## tweak MCMC
tweak.mcmc.FOI.model = function(chain.length, bb.max, omit, 
                                sdj, mu.prior, sd.prior, transform,varsin.nc,
                                print_tweak = F){
  
  #status=status.input
  mcmc.params = matrix(NA,nrow=chain.length/omit, ncol=length(status$params)+length(status$lnL))
  lnLserosurvey=NULL
  for (i in 1:n.serosurveys) {lnLserosurvey[i] = paste("lnL serosurvey", i, sep="_")}
  colnames(mcmc.params) = c(names(status$params),"lnL vac_eff","lnL GLM", lnLserosurvey, "lnL Prior")#names(status$lnL))
  
  accept = matrix(0,nrow=bb.max,ncol=length(pars.ini))
  colnames(accept) = names(status$params)
  sd_mat = matrix(0,nrow=bb.max,ncol=length(pars.ini))
  colnames(sd_mat) = names(status$params)
  
  n.to.tweak=n.non.static =rep(NA,bb.max) #  si le burn in marche bien, permet de ne pas lancer les 20 batch (bb.max)
  i.ex = integer(0)
  
  ar.min = 0.2
  ar.max = 0.4
  ar.low = 0.1
  ar.high = 0.6
  
  pars.ini = status$params
  
  if(print_tweak){
    windows()
    windows()
    windows()
    windows()
    windows()
  }
  
  set.seed(1)
  for(bb in 1:bb.max) { # tweak obtained after 13-15 batches 
    if(!grepl("NFS", homedir)) print(paste("tweak bb=", bb)) # means if we are not working on the cluster
    for(i in 1:chain.length) {
      #print(i)
      for(k in 1:length(pars.ini)) { # k index les parametres
       # print(k)
        status = fun.mcmc.step(status, sdj, mu.prior, sd.prior, transform, index=k, varsin.nc)
        accept[bb,k] = accept[bb,k]+status$accept
      } # end for(k in 1:length(par))
      if(i%%omit==0) {
        mcmc.params[i/omit,] = c(status$params,status$lnL) # une fois que j'ai fait toutes les chaines pour 1 batch, j'enregistre les valeurs de param dans la tab mcmc.params
      } # end if(i%%omit==0)
    }  # end for(i in 1:chain.length)
    ## tweaking the AR:
    n.to.tweak[bb] = length(pars.ini) #-length(i.ex)
    n.non.static[bb] = 0
    for(k in 1:length(pars.ini)) {  # k index les parametres
      if(!(k %in% i.ex)) {
        if(accept[bb,k]/chain.length<ar.low) sdj[k] = sdj[k]/1.3
        else if(accept[bb,k]/chain.length<ar.min) sdj[k] = sdj[k]/1.1
        else if(accept[bb,k]/chain.length>ar.high) sdj[k] = sdj[k]*1.3
        else if(accept[bb,k]/chain.length>ar.max) sdj[k] = sdj[k]*1.1
        else n.to.tweak[bb] = n.to.tweak[bb]-1
        #if(t.test(mcmc.params[1:(chain.length/2),k],mcmc.params[(1:(chain.length/2))+chain.length/2,k],alternative="two.sided")$p.value<0.05) n.non.static[bb] = n.non.static[bb]+1
      } # end if(!(k %in% i.ex))
    } # end for(k in 1:length(pars.ini))
    
    
    
    # plotting the tweaking
    if(print_tweak){
      dev.set(2)
      accept1 = accept[,1:(trunc(ncol(accept)/2))]
      barplot(accept1/chain.length,beside=TRUE,las=2)
      lines(c(0,1000),rep(ar.min,2),col=2)
      lines(c(0,1000),rep(ar.max,2),col=2)
      
      dev.set(3)
      accept2=accept[,(trunc(ncol(accept)/2)+1):ncol(accept)]
      barplot(accept2/chain.length,beside=TRUE,las=2)
      lines(c(0,1000),rep(ar.min,2),col=2)
      lines(c(0,1000),rep(ar.max,2),col=2)
      
      dev.set(4)
      par(mfcol=c(5,5),mar=c(2,3,1,1)+0.1,mgp=c(2,1,0),oma=c(0,0,2,0))
      for(i in 1:max(ii) ) {
        plot(mcmc.params[,i],type="l",xlab="",ylab=colnames(mcmc.params)[i])
      }
      
      dev.set(5)
      par(mfcol=c(5,5),mar=c(2,3,1,1)+0.1,mgp=c(2,1,0),oma=c(0,0,2,0))
      for(i in (max(ii)+1):(max(ii)+21) ) {
        plot(mcmc.params[,i],type="l",xlab="",ylab=colnames(mcmc.params)[i])
      }
      
      dev.set(6)
      par(mfcol=c(5,4),mar=c(2,3,1,1)+0.1,mgp=c(2,1,0),oma=c(0,0,2,0))
      for(i in (max(ii)+21):length(pars.ini)) {
        plot(mcmc.params[,i],type="l",xlab="",ylab=colnames(mcmc.params)[i])
      }
    
    }
    
    sd_mat[bb,] = sdj
    #write.csv(mcmc.params,paste("mcmc_bb",bb,".csv",sep=""),row.names=FALSE)
    if(n.to.tweak[bb]==0 & n.non.static[bb]==0) break
  } ## end for(bb in 1:bb.max)
  
  ### print acept + params
  date=format(Sys.time(),"%Y%m%d")
  nb.runs = chain.length*bb.max
  name.dir = paste(currdir, "tweak_MCMC_FOI_nb_runs=", nb.runs, "_", date, sep='')
  
  if(!dir.exists(name.dir)) dir.create(name.dir,  showWarnings = TRUE)
  setwd(name.dir)
  accept=accept/chain.length
  write.csv(accept, file="tweaking_accept.csv")
 
  # write the last burn-in stauts
  # to_write=c(mcmc.params[dim(mcmc.params)[1],], status$accept)
  # names(to_write)[length(to_write)]="accept"
  # write.csv(to_write, file=paste("last_mcmc_params_after_burnin-bbmax=",bb.max,"_nbchains=",chain.length, ".csv", sep="" ))
  # 
  
  return(list(mcmc.params=mcmc.params, sd_mat=sd_mat, accept=accept,
              status = status, sdj = sdj))
}








### function to run MCMC after tweaking
run.mcmc.FOI.model = function(chain.length, omit, bb.tot, nb_codes,index_code,
                              status,
                              sdj, mu.prior, sd.prior, transform,varsin.nc,
                              print_predictions, print_p.detect.nb_inf, print_burden){
  
  nb.runs = chain.length*bb.tot/omit
  batch_by_code=bb.tot/nb_codes
  bb.max= index_code*batch_by_code
  bb.min = bb.max-batch_by_code+1  
  
  #date=gsub(' ', '_', date())
  date=format(Sys.time(),"%Y%m%d")
  name.dir = paste(currdir, "MCMC_FOI_nb_runs=", nb.runs, "_", date, sep='')
  
  if(!dir.exists(name.dir)) dir.create(name.dir,  showWarnings = TRUE)
  setwd(name.dir)
  
  load(paste0(currdir,"p_prop3d.Rdata"))
  load(paste0(currdir,"P_tot_2d.Rdata"))
  
  # keep MCMC params
  mcmc.params = matrix(NA,nrow=chain.length/omit, ncol=length(status$params)+length(status$lnL)) # 1 table for each batch
  lnLserosurvey=NULL
  for (i in 1:n.serosurveys) {lnLserosurvey[i] = paste("lnL serosurvey", i, sep="_")}
  colnames(mcmc.params) = c(names(status$params),"lnL vac_eff","lnL GLM", lnLserosurvey, "lnL Prior")#names(status$lnL))
  rm(lnLserosurvey)
  
  accept = matrix(0,nrow=(bb.max-bb.min+1),ncol=length(status$params))
  colnames(accept) = names(status$params)
  
  # keep FOI/R0
  FOI_runs = matrix(NA, nrow=length(dat$adm0_adm1), ncol=chain.length/omit  )
  rownames(FOI_runs) = dat$adm0_adm1
  
  # keep burden: cases and DALYs
  years = 1984:2013
  burden_cases = matrix(NA, nrow=length(c34)*length(years), ncol=(1+chain.length/omit))
  burden_cases [,1] = rep(years, each=length(c34))
  rownames(burden_cases) =rep(c34, length(years))
  burden_DALY = burden_cases
  
  # keep p.detect.link
  p.detect.link_runs= rep(NA, chain.length/omit)
  
  # keep p.detect.multi and nb_inf_30
  p.detect.multi= matrix(NA, nrow=length(unlist(adm1s)), ncol=chain.length/omit)
  #nb_inf_30 = p.detect.multi
  
  # keep seroprev_predict_survey & sero.fit.glm
  seroprev_predict_survey = matrix(NA, nrow=(chain.length/omit)*n.serosurveys, ncol=length(0:100) )
  seroprev_predict_survey = data.frame(seroprev_predict_survey)
  sero.fit.glm = matrix(NA, nrow=(chain.length/omit)*n.serosurveys, ncol=length(0:100) )
  sero.fit.glm = data.frame(sero.fit.glm)
  
  t= as.numeric(Sys.time())
  seed= (t - floor(t)) * 1e8 
  #print(seed)
  set.seed(seed)
  
  
  ii = grep("^log.",names(status$params))
  for(bb in bb.min:bb.max) {
    if(!grepl("NFS", homedir)) print(paste("bb=", bb)) # means if we are not working on the cluster
    for(i in 1:chain.length) {
      #print(paste("i chain=", i))
      for(k in (1:length(status$params))){
        #print(k)
        status = fun.mcmc.step(status, sdj, mu.prior, sd.prior, transform, index=k, varsin.nc)
        accept[bb-bb.min+1,k] = accept[bb-bb.min+1,k]+status$accept
      }
      if(i%%omit==0) {
        mcmc.params[i/omit,] = c(status$params,status$lnL)
        
        # update vaccine coverage with vaccine efficacy
        vacc_eff = status$params[1]
        vc.full = vc.full_100
        vc.full[,4:length(vc.full)] = vacc_eff*vc.full[,4:length(vc.full_100)] # the 3 first columns are adm0, adm0_adm1, year
        vc30.agg= vc30.agg_100
        vc30.agg[,3:length(vc30.agg)]= vacc_eff*vc30.agg_100[,3:length(vc30.agg_100)] # the 2 first columns are adm0_adm1, year
        
        pop.vc.moments = data.frame(adm0_adm1 = pop30.agg$adm0_adm1, 
                                    m0 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1]),na.rm=T),
                                    m1 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat,na.rm=T), 
                                    m2 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat*amat,na.rm=T),
                                    m3 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat*amat*amat,na.rm=T),
                                    m4 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat*amat*amat*amat,na.rm=T),
                                    m5 = rowSums(pop30.agg[,-1]*(1-vc30.agg[,-1])*amat*amat*amat*amat*amat,na.rm=T)) 
        
        
        foi_vec = fun.calcFoi.Africa(status)
        FOI_runs[,i/omit] = foi_vec
        
        mypreds.nc = fun.calcPred(coefs = status$params[ii],newdata=x,type="link",varsin=varsin.nc)
        p.detect.link_runs[i/omit] = fun.calc.pdetect(params= status$params,mypreds.nc=mypreds.nc)
        
        if(print_p.detect.nb_inf==T){
          p.detect.multi[,i/omit] = fun.calc.pdetect.multi(params= status$params,mypreds.nc=mypreds.nc)
          
        }
        
        if(print_burden==T){
          for(YYYY in years){
            foi_vec = fun.calcFoi.Africa(status)
            burden_cases[burden_cases[,1]==YYYY,i/omit+1] = 0.1 * as.numeric(fun.aggBurden(foi_vec, YYYY, vac_eff_arg=status$params[1])$cases[,2]) #i/omit+1 because the first column is year
            burden_DALY[burden_DALY[,1]==YYYY,i/omit+1] = as.numeric(fun.aggBurden(foi_vec, YYYY)$life.years.left[,2])
            
          }
        }
        
        if(print_predictions==T){
          
          ## d abord les prediction serosurveys
          sero.out= NULL
          vac_eff = status$params[1]
          
          for (index_survey in 1:n.serosurveys){
            vcfact = ifelse(!is.na(vc.factor[index_survey]), vc.factor[index_survey], status$params[grep("vc.fac", names(status0$params))])
            vc= c(as.matrix(vac_eff*vcfact*vc.agg[index_survey,-(1:3)]))
            
            foi = status$params[length(ii)+index_survey+1]
            sero1 = c(as.matrix(1-exp(-foi*0:100)))
            sero.out_tmp = 1 - (1-sero1)*(1-vc)
            sero.out = rbind(sero.out, sero.out_tmp)  
          }
          #sero.out$sero.studies = sero.studies
          seroprev_predict_survey[ (1+((i/omit)-1)*n.serosurveys):(((i/omit))*n.serosurveys) ,] = sero.out
          
          ## puis les prediction modele complet
          sero.fit_glm_withoutherd_adm1 = list()
          for(index_survey in 1:n.serosurveys){
            sero.fit_glm_withoutherd = sero.fit_glm_withoutherd_indiv_tmp= NULL
            for (j in 1:length(adm1s[[index_survey]])){
              tadm=match(adm1s[[index_survey]][j], dat$adm0_adm1)
              foi_pred = FOI_runs[tadm,i/omit]
              
              s_tadm =  c(as.matrix(1-exp(-foi_pred*0:100)))
              vcfact = ifelse(!is.na(vc.factor[index_survey]), vc.factor[index_survey], status$params[grep("vc.fac", names(status0$params))])
              vc= c(as.matrix(vac_eff*vcfact*vc3d[tadm, study.years[[index_survey]]-1939,]))
              
              sero.fit_glm_withoutherd_indiv_tmp = 1 - (1-s_tadm)*(1-vc)
              sero.fit_glm_withoutherd = rbind(sero.fit_glm_withoutherd,  sero.fit_glm_withoutherd_indiv_tmp)
            }
            rownames(sero.fit_glm_withoutherd) = adm1s[[index_survey]]
            sero.fit_glm_withoutherd_adm1[[index_survey]] = sero.fit_glm_withoutherd
          }
          
          sero.fit_glm_withoutherd_agg = NULL
          for (index_survey in 1:n.serosurveys){
            tadm=match(adm1s[[index_survey]], dat$adm0_adm1)
            vec_tmp = sero.fit_glm_withoutherd_adm1[[index_survey]][,]*p_prop3d[tadm, study.years[[index_survey]]-1939,]*
              P_tot_2d[tadm, study.years[[index_survey]]-1939]
            if( length(adm1s[[index_survey]])>1 ) vec_tmp = colSums(vec_tmp)
            
            vec_tmp = vec_tmp/pop.agg[index_survey,-c(1:3)]
            sero.fit_glm_withoutherd_agg = rbind(sero.fit_glm_withoutherd_agg, vec_tmp)
          }
          #sero.fit_glm_withoutherd_agg$sero.studies = sero.studies
          sero.fit.glm[(1+((i/omit)-1)*n.serosurveys):(((i/omit))*n.serosurveys),] = sero.fit_glm_withoutherd_agg
        }# if(print_predictions==T)
      }# if(i%%omit==0)
    }#for(i in 1:chain.length)
    
    write.csv(mcmc.params,paste("mcmc_bb",bb,".csv",sep=""),row.names=FALSE)
    write.csv(FOI_runs,paste("FOI_bb",bb,".csv",sep=""),row.names=TRUE)
    write.csv(p.detect.link_runs,paste("p.detect.link_bb",bb,".csv",sep=""),row.names=FALSE)
    if(print_p.detect.nb_inf==T) write.csv(p.detect.multi,paste("p.detect.multi_bb",bb,".csv",sep=""),row.names=FALSE)
    if(print_burden==T){
      write.csv(burden_cases,paste("burden_cases_bb",bb,".csv",sep=""),row.names=FALSE)
      write.csv(burden_DALY,paste("burden_DALY_bb",bb,".csv",sep=""),row.names=FALSE)
    }
    if(print_predictions==T){
      seroprev_predict_survey_print = data.frame(sero.studies=rep(sero.studies, times = chain.length/omit), seroprev_predict_survey)
      sero.fit.glm_print = data.frame(sero.studies=rep(sero.studies, times = chain.length/omit), sero.fit.glm)
      
      write.csv(seroprev_predict_survey_print,paste("seroprev_predict_survey_bb",bb,".csv",sep=""),row.names=FALSE)
      write.csv(sero.fit.glm_print,paste("sero.fit.glm_bb",bb,".csv",sep=""),row.names=FALSE)
    }
    
    
  } ## end for(bb in 1:bb.max)
  accept=accept/chain.length
  write.csv(accept, file=paste("accept_bb_", bb.max, ".csv", sep="" ))
}
