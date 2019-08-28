##################################
# Functions for herd immunity model

# These functions include new serosurveys, vaccine efficacy (estimated against prior distribution only)
# CAF serosurvey is not accounted for right now
# This new version includes in the recurrence_seroprev function the 4 compartiments i+/_ v+/_





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



###########################
# calculate FOI for the prevac era from the R0 value

foi_survey_prevac = function(R0, pop.moments, polydeg=6){ # returns the foi of the pre_vacc period
  
  if(polydeg>length(pop.moments)) error("fun.calcFoi.Africa: invalid value for polydeg.\n")
  
  lambda = NULL
  if(polydeg>0) for(deg in 1:polydeg) {
    lambda = cbind(lambda,(-1)^(deg+1)*pop.moments[deg]/factorial(deg-1)) # attention, ici, le i fait reference a la survey sur laquelle on fait tourner MCMC
  }
  lambda[,1]=lambda[,1]-1/R0 # we have to put the equation =0, so the term of order 0 (first column) should integrate -1/R0
  
  out = sapply(1:nrow(lambda), function(j) polyroot(lambda[j,]) )
  
  out[abs(Arg(out))<=1e-10] = Re(out)[abs(Arg(out))<=1e-10]
  out[abs(Arg(out))>1e-10] = NA # here we have a resolution problem : in case of survey=5 and R0=5 (for instance), we have no solution in Real
  dt = dim(out)
  out = as.numeric(out)
  dim(out) = dt
  if (polydeg>2){out = apply(out,2,min,na.rm=T)}
  #names(out)= paste("FOI", sero.studies[i],sep="_") #pb now that i is no more used
  return(out)
  # consider to check whether the first FOI is the one of the first survey...
}


######################################
### Calculate the post vaccination fraction of susceptibles and fill the different tables (s_atau3_survey, etc)
# return a list with s(a,t), N_inf (a,t)...
## a specific function for the survey provinces (with aggregated population)

fun.recurrence.seroprev.survey = function(index_survey, R0, foi_const=0, t0_vac, p_at_survey, P_tot_survey,  inc_v3d, pop.moments,
                                          vac_eff_arg=1) { # t0_vac=t0_vac[I], p_at_survey=p_at_survey[I,,]
  # P_tot_survey=P_tot_survey[I,]
  
  # foi_const is a constant FOI for environmental transmissions
  
  # first create the table to keep the results
  s_atau1_survey=rep(NA, length(dn2)*length(dn3))
  dim(s_atau1_survey)= c(length(dn2), length(dn3))
  Stot_at_survey = Ninf_t_survey = Nenv_t_survey = rep(NA, length(dn2))
  n_inf_at_survey = s_atau3_survey = s_atau2_survey=  s_atau1bis_survey= n_env_at_survey = s_atau1_survey 
  
  # create the 4 compartiments accounting for natural infection i (+/-) and vaccination v (+/-)
  # I need these compartiments at the 3 tau times
  i_neg_v_neg_tau1 = i_neg_v_neg_tau2 = i_neg_v_neg_tau3 = matrix( rep(NA, length(dn2)*length(dn3)), nrow=length(dn2)) # i- v-
  i_pos_v_neg_tau1 = i_pos_v_neg_tau2 = i_pos_v_neg_tau3 = matrix( rep(NA, length(dn2)*length(dn3)), nrow=length(dn2)) # i+ v-
  i_neg_v_pos_tau1 = i_neg_v_pos_tau2 = i_neg_v_pos_tau3 = matrix( rep(NA, length(dn2)*length(dn3)), nrow=length(dn2)) # i- v+
  i_pos_v_pos_tau1 = i_pos_v_pos_tau2 = i_pos_v_pos_tau3 = matrix( rep(NA, length(dn2)*length(dn3)), nrow=length(dn2)) # i+ v+
  
  i_neg_v_neg_tau1_bis = i_neg_v_pos_tau1_bis = i_pos_v_neg_tau1_bis = i_pos_v_pos_tau1_bis = matrix( rep(NA, length(dn2)*length(dn3)), nrow=length(dn2)) # tis is for environmental infection
  
  # then calculate the pre_vacc prevalence
  foi = foi_survey_prevac(R0=R0, pop.moments=pop.moments, polydeg=6) # ici au besoin mettre R0 = round(R0, 6)
  
  # then calculate the rest of the prevalence
  
  for (y in  1:length(dn2) ){
    
    if ( dn2[y] < t0_vac ){ #il faudra plus t0_vac=t0_vac[I]  # SI PB CHANGER < pour <=
      # equation (4) at t+tau1
      i_neg_v_neg_tau1[y,] = exp( - (foi+foi_const)*avec )
      i_pos_v_neg_tau1[y,] = 1- exp( - (foi+foi_const)*avec )
      i_neg_v_pos_tau1[y,] = 0 # no vaccination
      i_pos_v_pos_tau1[y,] = 0 # no vaccination
      
      # aging
      i_neg_v_neg_tau2[y,] = i_neg_v_neg_tau3[y,] = exp( - (foi+foi_const)*(avec+1) )
      i_pos_v_neg_tau2[y,] = i_pos_v_neg_tau3[y,] = 1- exp( - (foi+foi_const)*(avec+1) )
      i_neg_v_pos_tau2[y,] = i_neg_v_pos_tau3[y,] = 0
      i_pos_v_pos_tau2[y,] = i_pos_v_pos_tau3[y,] = 0
      
      #record infections
      n_inf_at_survey[y,]= ( i_pos_v_neg_tau2[y,] -  i_pos_v_neg_tau1[y,])*p_at_survey[y,]*P_tot_survey[y]
      Ninf_t_survey[y] = sum(n_inf_at_survey[y,], na.rm=T)
      n_env_at_survey[y,] = (1-exp(-foi_const))*i_neg_v_neg_tau1[y,]*p_at_survey[y,]*P_tot_survey[y] # not sure it's the good formula...
      Nenv_t_survey[y] = sum(n_env_at_survey[y,], na.rm=T) 
      Stot_at_survey[y] = sum(i_neg_v_neg_tau1[y,]*p_at_survey[y,], na.rm=T)
      
    } else if (dn2[y] == t0_vac){
      
      i_neg_v_neg_tau1[y,] = exp( - (foi+foi_const)*avec )
      i_pos_v_neg_tau1[y,] = 1- exp( - (foi+foi_const)*avec )
      i_neg_v_pos_tau1[y,] = 0 # no vaccination
      i_pos_v_pos_tau1[y,] = 0 # no vaccination
      
      # aging
      i_neg_v_neg_tau2[y,] = exp( - (foi+foi_const)*(avec+1) )
      i_pos_v_neg_tau2[y,] = 1- exp( - (foi+foi_const)*(avec+1) )
      i_neg_v_pos_tau2[y,] = 0
      i_pos_v_pos_tau2[y,] = 0
      
      #record infections
      n_inf_at_survey[y,]= ( i_pos_v_neg_tau2[y,] -  i_pos_v_neg_tau1[y,])*p_at_survey[y,]*P_tot_survey[y]
      Ninf_t_survey[y] = sum(n_inf_at_survey[y,], na.rm=T)
      n_env_at_survey[y,] = (1-exp(-foi_const))*i_neg_v_neg_tau1[y,]*p_at_survey[y,]*P_tot_survey[y] # not sure it's the good formula...
      Nenv_t_survey[y] = sum(n_env_at_survey[y,], na.rm=T) 
      Stot_at_survey[y] = sum(i_neg_v_neg_tau1[y,]*p_at_survey[y,], na.rm=T)
      
      #tau3: vaccination arise
      i_neg_v_neg_tau3[y,] = i_neg_v_neg_tau2[y,]*(1-vac_eff_arg*inc_v3d[y,]) # eq (8)
      i_pos_v_neg_tau3[y,] = i_pos_v_neg_tau2[y,]*(1-vac_eff_arg*inc_v3d[y,]) # eq (8)
      i_neg_v_pos_tau3[y,] = i_neg_v_pos_tau2[y,]+i_neg_v_neg_tau2[y,]*(vac_eff_arg*inc_v3d[y,])
      i_pos_v_pos_tau3[y,] = i_pos_v_pos_tau2[y,]+i_pos_v_neg_tau2[y,]*(vac_eff_arg*inc_v3d[y,])
           
    } else if (dn2[y] > t0_vac){ # SI PB CHANGER >= pour >
      
      # children born susceptible and unvaccinated
      i_neg_v_neg_tau1[y,1] = 1
      i_pos_v_neg_tau1[y,1] = 0
      i_neg_v_pos_tau1[y,1] = 0 
      i_pos_v_pos_tau1[y,1] = 0 
      
      for (a in 2:length(dn3)){
        # tau1: aging   eq(4)
        i_neg_v_neg_tau1[y,a] = i_neg_v_neg_tau3[y-1,a-1]
        i_pos_v_neg_tau1[y,a] = i_pos_v_neg_tau3[y-1,a-1]
        i_neg_v_pos_tau1[y,a] = i_neg_v_pos_tau3[y-1,a-1]
        i_pos_v_pos_tau1[y,a] = i_pos_v_pos_tau3[y-1,a-1]
      }
      # tau1_bis: environmental infection
      i_neg_v_neg_tau1_bis[y,] = exp( - (foi_const)*avec )*i_neg_v_neg_tau1[y,]
      i_pos_v_neg_tau1_bis[y,] = i_pos_v_neg_tau1[y,] + (i_neg_v_neg_tau1[y,] - i_neg_v_neg_tau1_bis[y,])
      i_neg_v_pos_tau1_bis[y,] = i_neg_v_pos_tau1[y,] # no vaccination
      i_pos_v_pos_tau1_bis[y,] = i_pos_v_pos_tau1[y,]  # no vaccination
      
      n_env_at_survey[y,] = (1-exp(-foi_const))*i_neg_v_neg_tau1[y,]*p_at_survey[y,]*P_tot_survey[y] # not sure it's the good formula...
      Nenv_t_survey[y] = sum(n_env_at_survey[y,], na.rm=T)  
      
      if(dn2[y] < 1990) { # changed the 24/11/16
        Stot_at_survey[y] = sum(i_neg_v_neg_tau1_bis[y,]*p_at_survey[y,], na.rm=T) # total susceptibles are only those uninfected and unvaccinated
      } else if (dn2[y] >= 1990){ # from 1990, we have ages btw 77 and 100y, what we don't have before
        Stot_at_survey[y] = sum(i_neg_v_neg_tau1_bis[y,!is.na(i_neg_v_neg_tau1_bis[y,])]*p_at_survey[y,!is.na(i_neg_v_neg_tau1_bis[y,])])/
          sum(p_at_survey[y,!is.na(i_neg_v_neg_tau1_bis[y,])])
      }
      Ninf_t_survey[y] = max (0, P_tot_survey[y]*(Stot_at_survey[y] - 1/R0) ) # eq (5) 
      n_inf_at_survey[y,] = Ninf_t_survey[y]*i_neg_v_neg_tau1_bis[y,]*p_at_survey[y,]/Stot_at_survey[y] # eq (6) 
      
      #tau2: new infections arise
      i_neg_v_neg_tau2[y,] = i_neg_v_neg_tau1_bis[y,] - n_inf_at_survey[y,]/(p_at_survey[y,]*P_tot_survey[y] ) # eq (7)  
      i_pos_v_neg_tau2[y,] = i_pos_v_neg_tau1_bis[y,] + n_inf_at_survey[y,]/(p_at_survey[y,]*P_tot_survey[y] )
      i_neg_v_pos_tau2[y,] = i_neg_v_pos_tau1_bis[y,]
      i_pos_v_pos_tau2[y,] = i_pos_v_pos_tau1_bis[y,]
      
      #tau3: vaccination arise
      i_neg_v_neg_tau3[y,] = i_neg_v_neg_tau2[y,]*(1-vac_eff_arg*inc_v3d[y,]) # eq (8)
      i_pos_v_neg_tau3[y,] = i_pos_v_neg_tau2[y,]*(1-vac_eff_arg*inc_v3d[y,]) # eq (8)
      i_neg_v_pos_tau3[y,] = i_neg_v_pos_tau2[y,]+i_neg_v_neg_tau2[y,]*(vac_eff_arg*inc_v3d[y,])
      i_pos_v_pos_tau3[y,] = i_pos_v_pos_tau2[y,]+i_pos_v_neg_tau2[y,]*(vac_eff_arg*inc_v3d[y,])
    }
       
  } # end for (y in  1:length(dn2) )
  liste_resu= list(i_neg_v_neg_tau3=i_neg_v_neg_tau3, i_pos_v_neg_tau3=i_pos_v_neg_tau3, i_neg_v_pos_tau3=i_neg_v_pos_tau3, i_pos_v_pos_tau3=i_pos_v_pos_tau3,
                   Ninf_t_survey=Ninf_t_survey, Stot_at_survey=Stot_at_survey, Nenv_t_survey=Nenv_t_survey)
  #n_env_at_survey=n_env_at_survey, Nenv_t_survey=Nenv_t_survey)
  return(liste_resu)
}


################################################################################
### Calculate the aggregated seropositivity among the age groups of the serosurveys
fun.herd.calcSeroprev.ag = function(index_survey, R0, foi_const = 0, age.groups,vcfac_arg=1, vac_eff_arg=1) { # pourquoi ici je donne en argument s_atau3_survey???
#   if(length(vc)==1) {
#     vc = rep(vc,length(age.groups))
#   } else {
#     if(length(vc)!=length(age.groups)) {
#       stop("fun.calcSeroprev.ag: incompatible length of vaccination coverage vc.\n")
#     }
#   }
#   if(!all(!is.na(vc)) |min(vc)<0 |max(vc)>1)
#     stop("fun.calcSeroprev: invalid values for vaccination coverage vc.\n")
#   
  res = fun.recurrence.seroprev.survey(index_survey=index_survey, R0=R0, foi_const=foi_const, t0_vac=t0_vac[index_survey],
                                       p_at_survey= p_at_survey_3d[index_survey,,],
                                       P_tot_survey= P_tot_survey_2d[index_survey,], inc_v3d=inc_v3d.agg[index_survey,,], 
                                       pop.moments=pop.moments.agg[index_survey,],
                                       vac_eff_arg=vac_eff_arg)
  
  if (vcfac_arg==1){
    res_out = 1 - res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]
  } else if(vcfac_arg==0){
    res_out = 1 - res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]/(res$i_pos_v_neg_tau3[study.years[[index_survey]]-1939,] +res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,])
  } else if(vcfac_arg>0 & vcfac_arg<1 ){# for CMRs
    res_out = vcfac_arg*(1 - res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]) +
                (1-vcfac_arg)*(1 - res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]/(res$i_pos_v_neg_tau3[study.years[[index_survey]]-1939,] +res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]))
  } else  stop("fun.calcSeroprev: impossible value for vcfac\n")
  
  ag = findInterval(0:100,age.groups)  
  preval1 = res_out* p_at_survey_3d[index_survey, study.years[[index_survey]]-1939,]
  preval.ag = aggregate(preval1,by=list(age.group=ag),sum, na.rm=T)
  preval.ag$x = preval.ag$x/aggregate(p_at_survey_3d[index_survey, study.years[[index_survey]]-1939,],by=list(age.group=ag),sum)$x
  
  # here we calculate fractions of seropositive (accounting from vc) from the susceptibles
  #sero.out = 1 - (1-vc)*suscep.ag$x # esct-ce cela ou     sero.out = 1 - suscep.ag$x/(1-vc)
  # la formule change de celle de Tini car moi je travaille a partir de la fraction SUSCPETIBLe de la pop totale.
  
  sero.out = preval.ag$x # actually I should not correct vc by vac_eff cause people are excluded if they DECLARE having been vaccinated, which is not affected by vac_eff 
  
  # here if needed, make a condition for vc=1
  sero.out = ifelse( sero.out< 1e-10, 0, sero.out) # rounding may lead to very low (1e-16) negative values
  
  
  return(sero.out)
}



# then calculate the lnL for a given serosurvey and a given FOI
fun.herd.lnL.survey = function(index_survey, R0, foi_const=0, obs.tot,obs.pos,age.groups,vcfac_arg=1, vac_eff_arg=1) {
  if(length(obs.pos)!=length(obs.tot) | length(obs.pos) != length(age.groups))
    stop("fun.lnL: incompatible dimensions.\n")
  serop = fun.herd.calcSeroprev.ag(index_survey=index_survey, R0=R0, foi_const=foi_const, age.groups=age.groups,vcfac_arg=vcfac_arg, vac_eff_arg=vac_eff_arg)
  
  # log_serop and log_1_minus_serop lines added on 19 Oct 2015, to be retrieved if problems.
  log_serop = ifelse(serop==0, -50, log(serop)) # done to avoid Infinite lnL
  log_1_minus_serop = ifelse(serop==1, -50, log(1-serop)) # done to avoid Infinite lnL
  lnL = sum(lgamma(obs.tot+1)-lgamma(obs.pos+1)-lgamma(obs.tot-obs.pos+1) + obs.pos*log_serop + (obs.tot-obs.pos)*log_1_minus_serop)
  
  #lnL_groupe=lgamma(obs.tot+1)-lgamma(obs.pos+1)-lgamma(obs.tot-obs.pos+1) + obs.pos*log(serop) + (obs.tot-obs.pos)*log(1-serop)
  #lnL_groupe= ifelse( is.nan(lnL_groupe), -Inf, lnL_groupe)
  # I have to do this beacause if obs.pos=0 and serop=0, R can't evaluate -Inf*0 and returns NaN
  
  lnL = sum(lnL)
  
  return(lnL)
}



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
  
  lnLPrior[3] = sum( dexp(params[transform == 1]-1, rate=0.001, log=TRUE))
  # this term is for parameters distributed over [1,+Inf], ie R0:. We substract 1 to be distributed on [0, Inf] and follow a dexp
  
  lnLPrior[4] = 0 # here, we use an uniform prior which is constant, thus no need to include it, it will cancel we comparing to the previous MCMC step
  #lnLPrior[4] = sum( dexp(params[transform == 2], rate=0.001,log=TRUE))
  # this term is for parameters distributed over [0,1], ie vc.factor : quasi flat distrib                  
  
  lnLPrior[5] = dlogis(log(params[transform == 3]/(1-params[transform == 3])), location =3.999, scale = 0.128*sqrt(3)*pi, log=T)
  lnLPrior[5] = ifelse(params[transform == 3]<0.6, (lnLPrior[5]-15), lnLPrior[5]) # this is done to penalize values of vac_eff <0.6
  
  return(lnLPrior)
}





###############################
# fit MCMC - update status

fun.herd.mcmc.step = function(status, sdj, index=NA) {
  # index indexes the parameters  
  
  status.new = status
  
  if(is.na(index[1])) {
    stop("in proposal: invalid index")
  } 
  
  ############
  # Proposals calculation
  ############
  if(transform[index]==0) { ## use a normal prior
    status.new$params[index] =  rnorm(1,mean=status$params[index],sd=sdj[index]) # normal distrib for glm covar
    
  } else if (transform[index]==1){    
    status.new$params[index] = rlnormTrunc(1, meanlog= log(status$params[index]) , sdlog=sdj[index], min=1, max=Inf) 
    
  } else if (transform[index]==2 | transform[index]==3) {## use a truncated log normal prior
    status.new$params[index] = rlnormTrunc(1, meanlog= log(status$params[index]) , sdlog=sdj[index], min=0, max=1) 
    # log-truncated distrib for vc.factor.CMRs and vac_eff(varying from 0 to 1)
    
  } else { stop("in proposal: invalid transform")    }
  
  ############
  # Update Likelyhood calc
  ############
  # ii indexes the GLM parameters
  if( index == 1){# if param = vac_eff
    #     vac_eff = status.new$params[1]
    #     mm = grep("R0", names(status.new$params))
    #     lnL_tmp = NULL
    #     for (iii in 1:length(mm)) {
    #       if (iii != 2){ # just exclude CAF for lnL of vaccine efficacy - look if KEN_zone2 puts the mess
    #         R0 = status.new$params[mm[iii]]
    #         vcfac = vc.factor[iii]
    #         if(is.na(vcfac)) vcfac = status.new$params[which(names(status.new$params)==paste("vc.factor.",sero.studies[iii],sep=""))] # if vc is unknown, it is fitted with MCMC
    #         index_survey = iii
    #         lnL_tmp = c(lnL_tmp, fun.herd.lnL.survey(index_survey=index_survey, R0=R0, foi_const = foi_const_surv[index_survey],
    #                                                  obs.tot=sero.dat[[index_survey]][,1] , 
    #                                                  obs.pos=sero.dat[[index_survey]][,2], age.groups = age.min[[index_survey]], 
    #                                                  vc=(1-vcfac)*vc.agg.ag[[index_survey]], vac_eff_arg = vac_eff))
    #       }
    #     }
    # actually nothing to do: lnL is in the prior_lnL term
    #    status.new$lnL[1] = sum(lnL_tmp)
    
    
  } else if (index>1 & index<=max(ii)) { # if the index parameter is a GLM parameter, update GLM likelyhood (otherwise, not)
    # select GLM parameters
    status.new$lnL[1] = fun.glm.lnL(status.new$params[ii], x, y) # lnL  GLM
    
  } else if (max(ii)<index & index<=(max(ii)+n.serosurveys) ) { # if the index parameter is one of the R0 parameter but not the CMRs one
    vac_eff = status.new$params[1]
    index_survey = index - max(ii)
    vcfac = vc.factor[index_survey]
    if(is.na(vcfac)) vcfac = status.new$params[which(names(status.new$params)==paste("vc.factor.",sero.studies[index_survey],sep=""))] # if vc is unknown, it is fitted with MCMC
    R0=status.new$params[index]
    status.new$lnL[index_survey+1] = fun.herd.lnL.survey(index_survey=index_survey, R0=R0, foi_const = foi_const_surv[index_survey],
                                                         obs.tot=sero.dat[[index_survey]][,1] , 
                                                         obs.pos=sero.dat[[index_survey]][,2], age.groups = age.min[[index_survey]], 
                                                         vcfac_arg=vcfac, vac_eff_arg = vac_eff)
    # here I changed , previously it was  vc=vcfac*vc.agg.ag[[index_survey]])
  } else if (index==length(status$params)) { # fitting R0 for CMRs or vc.factor.CMRs
    vac_eff = status.new$params[1]
    i = which(is.na(vc.factor))  # WARNING: this works only if we have 1 vc.factor  = NA !
    index_survey = i
    R0 = status.new$params[max(ii)+which(is.na(vc.factor))]
    vcfac = status.new$params[length(status$params)]
    status.new$lnL[index_survey+1] = fun.herd.lnL.survey(index_survey=index_survey, R0=R0, foi_const = foi_const_surv[index_survey],
                                                         obs.tot=sero.dat[[index_survey]][,1] , 
                                                         obs.pos=sero.dat[[index_survey]][,2], age.groups = age.min[[index_survey]], 
                                                         vcfac_arg=vcfac, vac_eff_arg = vac_eff)
    # here I changed , previously it was  vc=vcfac*vc.agg.ag[[index_survey]])
  } else { stop("error in the index")    }
  
  # priors likelyhood
  status.new$lnL[length(status.new$lnL)] = sum(fun.lnL.prior(status.new$params, sd.prior, transform))
  
  # these lines are added to penalize low values of vac_eff for CAF only
#   if (index==grep("CAF", (names(status$params))) ) {
#     index_survey = index - max(ii)
#     vac_eff =  status.new$params[1]
#     prior_vac= dlogis(log(vac_eff/(1-vac_eff)), location =3.999, scale = 0.128*sqrt(3)*pi, log=T)
#     status.new$lnL[index_survey+1] = status.new$lnL[index_survey+1] + 5*prior_vac
#   }
  
  # accept/reject
  ratio.post = sum(status.new$lnL, na.rm=T) - sum(status$lnL, na.rm=T)
  
  if(transform[index]==0) {
    correction = 0 
  } else if (transform[index]==1) {
    # code modified the 20 Oct 2015
    tmp1 =  log (status$params[index]) / sdj[index]
    tmp2 =  log (status.new$params[index]) / sdj[index]
    #tmp1 = - log (status$params[index]) / sdj[index]
    #tmp2 = - log (status.new$params[index]) / sdj[index]
    correction = log(status.new$params[index]) + log(pnorm(tmp1)) - log(status$params[index]) - log(pnorm(tmp2)) # correction for truncated lognormal distribution
    if(is.nan(correction)) correction = 0
  } else if(transform[index]==2 | transform[index]==3){
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
  
  #the following is done to update the survey lnL if a new value of vac_eff is accepted
  # Pb : then the term lnl_surv changes when vac_eff changes, so not good... I have to create another status$lnL_bis...
  if( index == 1 & status.new$accept ==1){
    vac_eff = status.new$params[1]
         mm = grep("R0", names(status.new$params))
    #     lnL_tmp = NULL
         for (iii in 1:length(mm)) {
    #       if (iii != 2){ # just exclude CAF for lnL of vaccine efficacy - look if KEN_zone2 puts the mess
            R0 = status.new$params[mm[iii]]
            vcfac = vc.factor[iii]
            if(is.na(vcfac)) vcfac = status.new$params[which(names(status.new$params)==paste("vc.factor.",sero.studies[iii],sep=""))] # if vc is unknown, it is fitted with MCMC
            index_survey = iii
            
            status.new$lnL[index_survey+1] = fun.herd.lnL.survey(index_survey=index_survey, R0=R0, foi_const = foi_const_surv[index_survey],
                                                                 obs.tot=sero.dat[[index_survey]][,1] , 
                                                                 obs.pos=sero.dat[[index_survey]][,2], age.groups = age.min[[index_survey]], 
                                                                 vcfac_arg=vcfac, vac_eff_arg = vac_eff)
    }
  }
  return(status.new)
}






####################################################
####################################################
####################################################
####### Functions for Burden calculation
####################################################
####################################################
####################################################


###########################
# calculate FOI for the prevac era from the R0 value 
# function vectorized for the calculation among all provinces


foi_whole_prevac = function(adm, R0, pop.moments, polydeg=6){ # returns the foi of the pre_vacc period
  
  if(polydeg>ncol(pop.moments)) error("fun.calcFoi.Africa: invalid value for polydeg.\n")
  
  lambda = NULL
  if(polydeg>0) for(deg in 1:polydeg) {
    lambda = cbind(lambda,(-1)^(deg+1)*as.numeric(pop.moments[,deg+1])/factorial(deg-1)) # attention, ici, le i fait reference a la survey sur laquelle on fait tourner MCMC
  }
  lambda[,1]=lambda[,1]-1/R0 # we have to put the equation =0, so the term of order 0 (first column) should integrate -1/R0
  
  out = sapply(1:nrow(lambda), function(j) polyroot(lambda[j,]) )
  
  out[abs(Arg(out))<=1e-10] = Re(out)[abs(Arg(out))<=1e-10]
  out[abs(Arg(out))>1e-10] = NA # here we have a resolution problem : in case of survey=5 and R0=5 (for instance), we have no solution in Real
  dt = dim(out)
  out = as.numeric(out)
  dim(out) = dt
  if (polydeg>2){out = apply(out,2,min,na.rm=T)}
  names(out)= paste("FOI", dat$adm0_adm1[adm],sep="_")
  return(out)
  # consider to check whether the first FOI is the one of the first survey...
}



#### for a given set of R0 values, calculate the burden (liste resu with s_at, n_inf_at)

fun.recurrence.seroprev.whole = function(adm, R0, t0_vac_adm, vac_eff_arg=1) { # 
  # adm specifies the adm where to calculate the burden: 27 adm if we want to calculate b, 479 adm if we calculate burden in whole africa
    
  # first create the 3d table to keep the results for all provinces
  s_atau1=rep(NA, length(adm)*length(dn2)*length(dn3))
  dim(s_atau1)= c(length(adm), length(dn2), length(dn3))
  n_inf_at = s_atau3 = s_atau2= s_atau1 
  Stot_at = Ninf_t = rep(NA, length(adm)*length(dn2))
  dim(Stot_at) = dim(Ninf_t) = c( length(adm),length(dn2) )
  
  # create the 4 compartiments accounting for natural infection i (+/-) and vaccination v (+/-)
  # I need these compartiments at the 3 tau times
  i_neg_v_neg_tau1 = i_neg_v_neg_tau2 = i_neg_v_neg_tau3 =  rep(NA, length(adm)*length(dn2)*length(dn3)) # i- v-
  i_pos_v_neg_tau1 = i_pos_v_neg_tau2 = i_pos_v_neg_tau3 =  rep(NA, length(adm)*length(dn2)*length(dn3)) # i+ v-
  i_neg_v_pos_tau1 = i_neg_v_pos_tau2 = i_neg_v_pos_tau3 =  rep(NA, length(adm)*length(dn2)*length(dn3)) # i- v+
  i_pos_v_pos_tau1 = i_pos_v_pos_tau2 = i_pos_v_pos_tau3 =  rep(NA, length(adm)*length(dn2)*length(dn3))  # i+ v+
  
  dim(i_neg_v_neg_tau1) = dim(i_neg_v_neg_tau2) = dim(i_neg_v_neg_tau3) = c( length(adm),length(dn2), length(dn3) )
  dim(i_pos_v_neg_tau1) = dim(i_pos_v_neg_tau2) = dim(i_pos_v_neg_tau3) = c( length(adm),length(dn2), length(dn3) )
  dim(i_neg_v_pos_tau1) = dim(i_neg_v_pos_tau2) = dim(i_neg_v_pos_tau3) = c( length(adm),length(dn2), length(dn3) )
  dim(i_pos_v_pos_tau1) = dim(i_pos_v_pos_tau2) = dim(i_pos_v_pos_tau3) = c( length(adm),length(dn2), length(dn3) )
    
  # I firstly need to calculate for each adm the FOI at t0_vac
  pop.mom=NULL
  for(j in 1:length(adm)){
    #pop.mom = rbind (pop.mom, pop.moments.whole[adm[j], t0_vac_adm[j]-1939,])
    pop.mom = rbind (pop.mom, pop.moments.whole[adm[j], ]) # changed the 09/02
  }
  pop.mom = cbind(dat$adm0_adm1[adm], pop.mom)
  
  foi = foi_whole_prevac(adm = adm, R0=R0, pop.moments = pop.mom, polydeg=6)
  
  for (i in 1:length(adm) ){
    
    time_limit = min(1981, t0_vac_adm[i])
    
    for (y in  1:length(dn2) ){
      if ( dn2[y] < time_limit){
        i_neg_v_neg_tau1[i,y,] = exp( - foi[i]*avec )
        i_pos_v_neg_tau1[i,y,] = 1- exp( - foi[i]*avec )
        i_neg_v_pos_tau1[i,y,] = 0 # no vaccination
        i_pos_v_pos_tau1[i,y,] = 0 # no vaccination
        
        # aging
        i_neg_v_neg_tau2[i,y,] = i_neg_v_neg_tau3[i,y,] = exp( - foi[i]*(avec+1) )
        i_pos_v_neg_tau2[i,y,] = i_pos_v_neg_tau3[i,y,] = 1- exp( - foi[i]*(avec+1) )
        i_neg_v_pos_tau2[i,y,] = i_neg_v_pos_tau3[i,y,] = 0
        i_pos_v_pos_tau2[i,y,] = i_pos_v_pos_tau3[i,y,] = 0
        
        #record infections
        n_inf_at[i,y,]= ( i_pos_v_neg_tau2[i,y,] -  i_pos_v_neg_tau1[i,y,])*p_prop3d[adm[i],y,]*P_tot_2d[adm[i],y]
        Ninf_t[i,y] = sum( n_inf_at[i,y,], na.rm=T)
        Stot_at[i,y] = sum(i_neg_v_neg_tau1[i,y,]*p_prop3d[adm[i],y,], na.rm=T)
      
      } else if (dn2[y] == time_limit){
        i_neg_v_neg_tau1[i,y,] = exp( - foi[i]*avec )
        i_pos_v_neg_tau1[i,y,] = 1- exp( - foi[i]*avec )
        i_neg_v_pos_tau1[i,y,] = 0 # no vaccination
        i_pos_v_pos_tau1[i,y,] = 0 # no vaccination
        
        # aging
        i_neg_v_neg_tau2[i,y,] = i_neg_v_neg_tau3[i,y,] = exp( - foi[i]*(avec+1) )
        i_pos_v_neg_tau2[i,y,] = i_pos_v_neg_tau3[i,y,] = 1- exp( - foi[i]*(avec+1) )
        i_neg_v_pos_tau2[i,y,] = i_neg_v_pos_tau3[i,y,] = 0
        i_pos_v_pos_tau2[i,y,] = i_pos_v_pos_tau3[i,y,] = 0  
        
        #record infections
        n_inf_at[i,y,]= ( i_pos_v_neg_tau2[i,y,] -  i_pos_v_neg_tau1[i,y,])*p_prop3d[adm[i],y,]*P_tot_2d[adm[i],y]
        Ninf_t[i,y] = sum( n_inf_at[i,y,], na.rm=T)
        Stot_at[i,y] = sum(i_neg_v_neg_tau1[i,y,]*p_prop3d[adm[i],y,], na.rm=T)
        
        #tau3: vaccination arise
        i_neg_v_neg_tau3[i,y,] = i_neg_v_neg_tau2[i,y,]*(1-vac_eff_arg*inc_v3d[adm[i],y,]) # eq (8)
        i_pos_v_neg_tau3[i,y,] = i_pos_v_neg_tau2[i,y,]*(1-vac_eff_arg*inc_v3d[adm[i],y,]) # eq (8)
        i_neg_v_pos_tau3[i,y,] = i_neg_v_pos_tau2[i,y,]+i_neg_v_neg_tau2[i,y,]*(vac_eff_arg*inc_v3d[adm[i],y,])
        i_pos_v_pos_tau3[i,y,] = i_pos_v_pos_tau2[i,y,]+i_pos_v_neg_tau2[i,y,]*(vac_eff_arg*inc_v3d[adm[i],y,])
        
      } else if (dn2[y] > time_limit){
        
        # children born susceptible and unvaccinated
        i_neg_v_neg_tau1[i,y,1] = 1
        i_pos_v_neg_tau1[i,y,1] = 0
        i_neg_v_pos_tau1[i,y,1] = 0 
        i_pos_v_pos_tau1[i,y,1] = 0
                
        for (a in 2:length(dn3)){
          # tau1: aging   eq(4)
          i_neg_v_neg_tau1[i,y,a] = i_neg_v_neg_tau3[i,y-1,a-1]
          i_pos_v_neg_tau1[i,y,a] = i_pos_v_neg_tau3[i,y-1,a-1]
          i_neg_v_pos_tau1[i,y,a] = i_neg_v_pos_tau3[i,y-1,a-1]
          i_pos_v_pos_tau1[i,y,a] = i_pos_v_pos_tau3[i,y-1,a-1]
        }
        
        if(dn2[y] < 1990) { 
          Stot_at[i,y] = sum(i_neg_v_neg_tau1[i,y,]*p_prop3d[adm[i],y,], na.rm=T)
        } else if (dn2[y] >= 1990){ # from 1990, we have ages btw 77 and 100y, what we don't have before
          Stot_at[i,y] = sum(i_neg_v_neg_tau1[i,y,!is.na(i_neg_v_neg_tau1[i,y,])]*p_prop3d[adm[i],y,!is.na(i_neg_v_neg_tau1[i,y,])])/
            sum(p_prop3d[adm[i],y,!is.na(i_neg_v_neg_tau1[i,y,])]) # total susceptibles are only those uninfected and unvaccinated
        }
        Ninf_t[i,y] = max (0,  P_tot_2d[adm[i],y]*(Stot_at[i,y] - 1/R0[i]) ) # eq (5) 
        n_inf_at[i,y,] = Ninf_t[i,y]*i_neg_v_neg_tau1[i,y,]*p_prop3d[adm[i],y,]/Stot_at[i,y] # eq (6) 
        
        #tau2: new infections arise
        i_neg_v_neg_tau2[i,y,] = i_neg_v_neg_tau1[i,y,] - n_inf_at[i,y,]/(p_prop3d[adm[i],y,]*P_tot_2d[adm[i],y] ) # eq (7)  
        i_pos_v_neg_tau2[i,y,] = i_pos_v_neg_tau1[i,y,] + n_inf_at[i,y,]/(p_prop3d[adm[i],y,]*P_tot_2d[adm[i],y] )
        i_neg_v_pos_tau2[i,y,] = i_neg_v_pos_tau1[i,y,]
        i_pos_v_pos_tau2[i,y,] = i_pos_v_pos_tau1[i,y,]
        
        #tau3: vaccination arise
        i_neg_v_neg_tau3[i,y,] = i_neg_v_neg_tau2[i,y,]*(1-vac_eff_arg*inc_v3d[adm[i],y,]) # eq (8)
        i_pos_v_neg_tau3[i,y,] = i_pos_v_neg_tau2[i,y,]*(1-vac_eff_arg*inc_v3d[adm[i],y,]) # eq (8)
        i_neg_v_pos_tau3[i,y,] = i_neg_v_pos_tau2[i,y,]+i_neg_v_neg_tau2[i,y,]*(vac_eff_arg*inc_v3d[adm[i],y,])
        i_pos_v_pos_tau3[i,y,] = i_pos_v_pos_tau2[i,y,]+i_pos_v_neg_tau2[i,y,]*(vac_eff_arg*inc_v3d[adm[i],y,])
        
       
      }
    }
  }
  
  liste_resu= list(i_neg_v_neg_tau3=i_neg_v_neg_tau3, i_pos_v_neg_tau3=i_pos_v_neg_tau3, i_neg_v_pos_tau3=i_neg_v_pos_tau3, 
                   i_pos_v_pos_tau3=i_pos_v_pos_tau3, Stot_at=Stot_at, Ninf_t_province=Ninf_t,
                   n_inf_at_province=n_inf_at)
  #return(Ninf_t)
  return(liste_resu)
}



##########################################################################
# calculating factor b
fun.calc.pdetect = function(params){#,mypreds,mypreds.nc) {
  # calculate the b terme from the left part of equation 6(plos Med paper)
  
  adms = match(unlist(adm1s),dat$adm0_adm1) 
  #adms = adms[-grep("CAF", dat$adm0_adm1[adms])]# this bit to change if I don't want CAF included in the estimates
  t0_vac_adm = t0_vac_africa[adms]
  
  # now I have to calculate the Number of infections over 30 years among each PROVINCE covered by a survey
  R0.vec = NULL
  for(k in c(1:max(n.serosurveys))) { # repeat the numbers of time of each adm1s
    R0.vec = c(R0.vec,rep(params[which(names(params)==paste("R0_",sero.studies[k],sep=""))],length(adm1s[[k]])) )
  }
  
  res = fun.recurrence.seroprev.whole (adms, R0.vec, t0_vac_adm, vac_eff_arg=params[1])
  
  # res returns the total Nb of infection on the whole period, I only need 1984 to 2013
  Ninf.30 = rowSums( res$Ninf_t_province[,which(dn2==1984):which(dn2==2013)] )
  Ninf.30 = ifelse(Ninf.30<1, 1, Ninf.30) #otherwise, we can have Ninf btw 0 & 1, thus log(Ninf)<0 and p.detect>0
  
  mypreds.nc = fun.calcPred(coefs = params[ii],newdata=x,type="link",varsin=varsin.nc) #ca correspond au terme X*beta - beta_c
  
  b = mypreds.nc[adms]-log(Ninf.30)
  
  return( mean(b) )
}

##########################################################################
# calculating factor b - same but returns the vector un-averaged
fun.calc.pdetect.multi = function(params){#,mypreds,mypreds.nc) {
  # calculate the b terme from the left part of equation 6(plos Med paper)
  
  adm = match(unlist(adm1s),dat$adm0_adm1) # mm= position of the provinces with serosurveys (27 provinces)
  t0_vac_adm = t0_vac_africa[adm]
  
  # now I have to calculate the Number of infections over 30 years among each PROVINCE covered by a survey
  R0.vec = NULL
  for(k in 1:n.serosurveys) { # repeat the numbers of time of each adm1s
    R0.vec = c(R0.vec,rep(params[which(names(params)==paste("R0_",sero.studies[k],sep=""))],length(adm1s[[k]])) )
  }
  
  res = fun.recurrence.seroprev.whole (adm, R0.vec, t0_vac_adm, vac_eff_arg=params[1])
  
  # res returns the total Nb of infection on the whole period, I only need 1984 to 2013
  Ninf.30 = rowSums( res$Ninf_t_province[,which(dn2==1984):which(dn2==2013)] )
  Ninf.30 = ifelse(Ninf.30<1, 1, Ninf.30)
  
  mypreds.nc = fun.calcPred(coefs = status$params[ii],newdata=x,type="link",varsin=varsin.nc) #ca correspond au terme X*beta - beta_c
  
  b = mypreds.nc[adm]-log(Ninf.30)
  names(b) = unlist(adm1s)
  
  return( b )
}


######################################################
## calculate R0 throughout Africa

fun.calcR0.Africa = function(status) {
  
  mypreds.nc  = fun.calcPred(coefs = status$params[ii],newdata=x,type="link",varsin=varsin.nc) #ca correspond au terme X*beta - beta_c
  p.detect.link = fun.calc.pdetect(params=status$params)
  
  Ninf.30_whole = exp( mypreds.nc - p.detect.link) ## ICI il manque une variable pays!!!!!
  
  
  R0_whole = rep(NA, length(Ninf.30_whole))
  for (i in 1: length(dn1) ){
    inf.bound = findInterval(Ninf.30_whole[i], R0_lookup[i,])
    sup.bound = inf.bound + 1
    
    if( sup.bound <= ncol(R0_lookup)){
      x.R0_1 = R0_lookup[i,inf.bound]# find inf bound
      x.R0_2 = R0_lookup[i,sup.bound] # sup bound
      
      y.R0_1 = as.numeric(colnames(R0_lookup)[inf.bound]) # find inf bound
      y.R0_2 = as.numeric(colnames(R0_lookup)[sup.bound]) # sup bound
      
      ## manualy calculate the linear interpolation
      a.lin.inter = (y.R0_2 - y.R0_1)/(x.R0_2-x.R0_1)
      b.lin.inter = y.R0_1 - a.lin.inter*x.R0_1
      
      R0_whole[i] = a.lin.inter*Ninf.30_whole[i] + b.lin.inter
      
    } else if (sup.bound  > ncol(R0_lookup)){
      R0_whole[i] = as.numeric(colnames(R0_lookup)[ncol(R0_lookup)])
    }
    
  }
  
  return(R0_whole)
}
#summary(R0_whole)
#hist(R0_whole)


######################################################
## Estimate yearly burden based on R0

fun.calcBurden = function(R0s,year,adm1 = NA, seed = NA, vac_eff_arg=1, total_run, current_run){
    # here, year as to be a single value
  # adm1 should be a adm0_adm1 code, not an index
  # total_run and current_run are just here to be sure that we use the same CFR with the same R0 estimates btw different scenarios
  # total_run : nb of posterior samples generated
  # current_run = index of the current sample we're estimating
  
  
  #if(length(year)>1) { stop("year has length>1")}
  # this function is defined for year having length = 1
  
  if(any(year<1984)) { stop("error: 1 year value is <1984, impossible to calculate burden if year <1984 (or you'll have to modify the fun.recurrence.seroprev.whole function)") } 
  
  if(is.na(adm1[1])) adm1=dat$adm0_adm1
  mm = match(adm1,dat$adm0_adm1) # mm is an INDEX
  
  if(length(R0s) != length(adm1)){ stop("error: R0s and adm1 have not the same length)") } 
  
  t0_vac_fun = t0_vac_africa[mm]
  
  burd = fun.recurrence.seroprev.whole(mm, R0s, t0_vac_adm=t0_vac_fun, vac_eff_arg=vac_eff_arg)
  
  burden.by.age = burd$n_inf_at_province[ ,match(year, dn2),]
  adm0 = dat$adm0[mm]
  
  cases = severe = deaths = life.years.left = NULL
  
  
  # generate CFR distributions
  if (!is.na(seed)) {set.seed(seed) 
  } else { set.seed(101)
  }
  prop.severe_all = rlnorm(total_run,meanlog=-2.222, sdlog=0.427)
  prop.severe_all[prop.severe_all>1]=1
  prop.death.all.cases_all = rlnorm(total_run, meanlog= -3.00, sdlog=0.442)
  prop.death.all.cases_all[prop.death.all.cases_all>1]=1
  
  prop.severe = prop.severe_all[current_run]
  prop.death.all.cases = prop.death.all.cases_all[current_run]
  
  for (yyyy in year){
    
    life1 = life0[life0$year %in% yyyy,]
    mm.adm0 = match(adm0,life1$country)
    life1 = as.matrix(life1[mm.adm0,-(1:2)])
    
    # here we have in ~10% of the cases prop.death.all.cases > prop.severe
    # this is not a technical problem, but may be a theoretical one
    cases_tmp = rowSums(burden.by.age[,which(year == yyyy),],na.rm=T)
    cases = cbind(cases, cases_tmp)
    
    severe_tmp = prop.severe*cases_tmp
    severe = cbind(severe, severe_tmp)
    
    deaths_tmp = prop.death.all.cases*cases_tmp
    deaths = cbind(deaths, deaths_tmp)
    
    YLL = prop.death.all.cases*rowSums(burden.by.age[,which(year == yyyy),]*life1,na.rm=TRUE)
    # times and disability weights are those from LaBeaud et al
    if(prop.death.all.cases>prop.severe)prop.death.all.cases = prop.severe # otherwise i can have negative value for (prop.severe-prop.death.all.cases)
    YLD = 17.8/365.25*0.172 *prop.severe*rowSums(burden.by.age[,which(year == yyyy),]*life1,na.rm=TRUE) + #acute phase among all severe cases
      28/365.25*0.024*(prop.severe-prop.death.all.cases)*rowSums(burden.by.age[,which(year == yyyy),]*life1,na.rm=TRUE)
    
    life.years.left_tmp = YLL + YLD
    life.years.left = cbind(life.years.left, life.years.left_tmp)
    
    # YLL_central = 0.0548*rowSums(burden.by.age[,which(year == yyyy),]*life1,na.rm=TRUE)
    # YLD_central = 17.8/365.25*0.172 *0.118*rowSums(burden.by.age[,which(year == yyyy),]*life1,na.rm=TRUE) + #acute phase among all severe cases
    #   28/365.25*0.024*(0.118-0.0548)*rowSums(burden.by.age[,which(year == yyyy),]*life1,na.rm=TRUE)
    # life.years.left_central_tmp = YLL_central + YLD_central # this is for GAVI central outputs - need a constant fraction of deaths
    # life.years.left_central = cbind(life.years.left_central, life.years.left_central_tmp)
    
   
  }
  
  colnames(cases) = colnames(severe) =colnames(life.years.left) = year
  
  return(list(cases=cases,severe=severe,deaths=deaths,life.years.left=life.years.left))
  #return(list(burd_tot=burd_tot, life.years.left=life.years.left))
}



################################################################################
fun.aggBurden = function(R0s,years,adm1 = NA, seed= NA, vac_eff_arg=1, total_run, current_run ) {
  if(is.na(adm1[1])) adm1=dat$adm0_adm1
  
  if(length(R0s) != length(adm1)){ stop("error: R0s and adm1 have not the same length)") }
  
  burden.adm1 = fun.calcBurden(R0s=R0s,year = years, adm1=adm1, seed = seed, vac_eff_arg=vac_eff_arg,
                               total_run = total_run, current_run=current_run)
  
  cases = aggregate(burden.adm1$cases,by=list(adm0 = dat$adm0),sum)
  severe = aggregate(burden.adm1$severe,by=list(adm0 = dat$adm0),sum)
  deaths = aggregate(burden.adm1$deaths,by=list(adm0 = dat$adm0),sum)
  lyleft = aggregate(burden.adm1$life.years.left,by=list(adm0 = dat$adm0),sum)
  #lyleft_central = aggregate(burden.adm1$life.years.left_central,by=list(adm0 = dat$adm0),sum)
  #names(cases)[2] = names(severe)[2] =names(deaths)[2] =names(lyleft)[2] = names(lyleft_central)[2] = years[1]
  
  return(list(cases=cases,severe.cases=severe, deaths=deaths, life.years.left=lyleft))
}











######## main functions for MCMC


## tweak MCMC
tweak.mcmc.R0.model = function(chain.length, bb.max, omit, 
                                sdj, sd.prior, transform,varsin.nc,
                                print_tweak = F){
  mcmc.params = matrix(NA,nrow=chain.length/omit, ncol=length(status$params)+length(status$lnL))
  lnLserosurvey=NULL
  for (i in 1:n.serosurveys) {lnLserosurvey[i] = paste("lnL serosurvey", i, sep="_")}
  colnames(mcmc.params) = c(names(status$params),"lnL GLM", lnLserosurvey, "lnL Prior")#names(status$lnL))
  rm(lnLserosurvey)
  
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
  ptm = proc.time()
  
  for(bb in 1:bb.max) { # tweak obtained after 13-15 batches 
    if(!grepl("NFS", homedir)) print(paste("tweak bb=", bb)) # means if we are not working on the cluster
    for(chain in 1:chain.length) {
      #print(i)
      for(index in 1:length(status$params)) { # k index les parametres
        status =fun.herd.mcmc.step(status, sdj, index)
        accept[bb,index] = accept[bb,index]+status$accept
      } # end for(k in 1:length(par))
      if(chain%%omit==0) {
        mcmc.params[chain/omit,] = c(status$params,status$lnL) # une fois que j'ai fait toutes les chaines pour 1 batch, j'enregistre les valeurs de param dans la tab mcmc.params
      } # end if(i%%omit==0)
    }  # end for(i in 1:chain.length)
    ## tweaking the AR:
    n.to.tweak[bb] = length(status$params) #-length(i.ex)
    n.non.static[bb] = 0
    for(k in 1:length(status$params)) {  # k index les parametres
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
    
    if(n.to.tweak[bb]==0 & n.non.static[bb]==0) break
    
  } ## end for(bb in 1:bb.max)

  ### print acept + params
  date=format(Sys.time(),"%Y%m%d")
  nb.runs = chain.length*bb.max
  name.dir = paste(currdir, "tweak_MCMC_R0_nb_runs=", nb.runs, "_", date, sep='')
  
  if(!dir.exists(name.dir)) dir.create(name.dir,  showWarnings = TRUE)
  setwd(name.dir)
  accept=accept/chain.length
  write.csv(accept, file="tweaking_accept.csv")
  
  write.csv(sdj, file=paste0(currdir, "sdj_after_tweaking_nb_run=", nb.runs, ".csv"), row.names=F)
  
  return(list(mcmc.params=mcmc.params, sd_mat=sd_mat, accept=accept,
              status = status, sdj=sdj))
}



##### run mcmc
run.mcmc.r0.model = function(chain.length, omit, bb.tot, nb_codes,index_code,
                             status,
                              sdj, sd.prior, transform,varsin.nc,
                              print_predictions, print_p.detect.nb_inf, print_burden){
  
  nb.runs = chain.length*bb.tot/omit
  batch_by_code=bb.tot/nb_codes
  bb.max= index_code*batch_by_code
  bb.min = bb.max-batch_by_code+1  
  
  date=format(Sys.time(),"%Y%m%d")
  name.dir = paste(currdir,"finalMCMC_R0_nb_runs=", nb.runs, "_", date, sep='')
  if(!dir.exists(name.dir)) dir.create(name.dir,  showWarnings = TRUE)
  setwd(name.dir)
  
  # write the last burn-in stauts
  to_write=c(mcmc.params[dim(mcmc.params)[1],], status$accept)
  names(to_write)[length(to_write)]="accept"
  write.csv(to_write, file=paste("last_mcmc_params_after_burnin-bbmax=",bb.max,"_nbchains=",chain.length, ".csv", sep="" ))
  
  
  # keep MCMC params
  mcmc.params = matrix(NA,nrow=chain.length/omit, ncol=length(status$params)+length(status$lnL)) # 1 table for each batch
  lnLserosurvey=NULL
  for (i in 1:n.serosurveys) {lnLserosurvey[i] = paste("lnL serosurvey", i, sep="_")}
  colnames(mcmc.params) = c(names(status$params),"lnL GLM", lnLserosurvey, "lnL Prior")#names(status$lnL))
  rm(lnLserosurvey)
  
  accept = matrix(0,nrow=(bb.max-bb.min+1),ncol=length(status$params))
  colnames(accept) = names(status$params)
  
  
  # keep FOI/R0
  R0_runs = matrix(NA, nrow=length(dat$adm0_adm1), ncol=chain.length/omit  )
  rownames(R0_runs) = dat$adm0_adm1
  
  # keep p.detect.link
  p.detect.link_runs= rep(NA, chain.length/omit)
  
  # keep p.detect.multi
  p.detect.multi_runs = rep(NA, (chain.length/omit)*length(unlist(adm1s)) )
  dim(p.detect.multi_runs) = c(length(unlist(adm1s)), chain.length/omit)
  
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
    if(!grepl("NFS", homedir)) print(paste("tweak bb=", bb)) # means if we are not working on the cluster
    for(i in 1:chain.length) {
      #print(paste("i chain=", i))
      for(index in (1:length(status$params))){
        #print(k)
        status =fun.herd.mcmc.step(status, sdj, index)
        accept[bb-bb.min+1,index] = accept[bb-bb.min+1,index]+status$accept
      }
      if(i%%omit==0) {
        mcmc.params[i/omit,] = c(status$params,status$lnL)
        R0_runs[,i/omit] = fun.calcR0.Africa(status)
        
        #mypreds.nc = fun.calcPred(coefs = params[ii],newdata=x,type="link",varsin=varsin.nc) 
        #p.detect.link_runs[i/omit] = fun.calc.pdetect(params= status$params)
        p.detect.multi_runs[,i/omit] = fun.calc.pdetect.multi(params= status$params)
        p.detect.link_runs[i/omit] = mean( p.detect.multi_runs[,i/omit])
        
        if(print_predictions==T){
          vac_eff = status$params[1]
          sero.out= NULL
          for (index_survey in 1:n.serosurveys){
            ## d abord les prediction serosurveys
            
            vcfac = ifelse(!is.na(vc.factor[index_survey]), vc.factor[index_survey], status$params[grep("vc.fac", names(status0$params))])
            #vc=vac_eff*(1-vcfact)*vc.agg3d[index_survey, study.years[[index_survey]]-1939,]
            
            R0surv = status$params[length(ii)+index_survey+1]
            res = fun.recurrence.seroprev.survey(R0=R0surv, t0_vac=t0_vac_sauv[index_survey], 
                                                 p_at_survey= p_at_survey_3d[index_survey,,],
                                                 P_tot_survey= P_tot_survey_2d[index_survey,], inc_v3d=inc_v3d.agg[index_survey,,], 
                                                 pop.moments=pop.moments.agg[index_survey,],
                                                 vac_eff_arg = vac_eff)
            if (vcfac==1){
              res_out = 1 - res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]
            } else if(vcfac==0){
              res_out = 1 - res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]/(res$i_pos_v_neg_tau3[study.years[[index_survey]]-1939,] +res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,])
            } else if(vcfac>0 & vcfac<1 ){# for CMRs
              res_out = vcfac*(1 - res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]) +
                (1-vcfac)*(1 - res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]/(res$i_pos_v_neg_tau3[study.years[[index_survey]]-1939,] +res$i_neg_v_neg_tau3[study.years[[index_survey]]-1939,]))
            }
            sero.out = rbind(sero.out,  res_out)
            
          }
          seroprev_predict_survey[ (1+((i/omit)-1)*n.serosurveys):(((i/omit))*n.serosurveys) ,] = sero.out
          
          
          ## puis les prediction modele complet
          
          sero.fit_glm_with_herd_adm1 = list()
          for(index_survey in 1:n.serosurveys){
            sero.fit_glm_with_herd = sero.fit_glm_with_herd_indiv_tmp= NULL
            
            for ( j in 1:length(adm1s[[index_survey]]) ){
              tadm=match(adm1s[[index_survey]][j], dat$adm0_adm1)
              t_survey = study.years[[index_survey]]
              t0_vac_pro = t0_vac_africa[tadm]
              R0_pred = R0_runs[tadm,i/omit]
              
              res = fun.recurrence.seroprev.whole(tadm , R0_pred, t0_vac_pro, vac_eff_arg=vac_eff)#$i_neg_v_neg_tau3[1, t_survey-1939,])
              vcfac = ifelse(!is.na(vc.factor[index_survey]), vc.factor[index_survey], status$params[grep("vc.fac", names(status0$params))])
              
              if (vcfac==1){
                res_out = 1 - res$i_neg_v_neg_tau3[1,study.years[[index_survey]]-1939,]
              } else if(vcfac==0){
                res_out = 1 - res$i_neg_v_neg_tau3[1,study.years[[index_survey]]-1939,]/(res$i_pos_v_neg_tau3[1,study.years[[index_survey]]-1939,] +res$i_neg_v_neg_tau3[1,study.years[[index_survey]]-1939,])
              } else if(vcfac>0 & vcfac<1 ){# for CMRs
                res_out = vcfac*(1 - res$i_neg_v_neg_tau3[1,study.years[[index_survey]]-1939,]) +
                  (1-vcfac)*(1 - res$i_neg_v_neg_tau3[1,study.years[[index_survey]]-1939,]/(res$i_pos_v_neg_tau3[1,study.years[[index_survey]]-1939,] +res$i_neg_v_neg_tau3[1,study.years[[index_survey]]-1939,]))
              }
              
              sero.fit_glm_with_herd_indiv_tmp = res_out
              #sero.fit_glm_with_herd_indiv_tmp = ifelse( sero.fit_glm_with_herd_indiv_tmp< 1e-10, 0, sero.fit_glm_with_herd_indiv_tmp)
              
              sero.fit_glm_with_herd = rbind(sero.fit_glm_with_herd,  sero.fit_glm_with_herd_indiv_tmp)
            }
            rownames(sero.fit_glm_with_herd) = adm1s[[index_survey]]
            sero.fit_glm_with_herd_adm1[[index_survey]] = sero.fit_glm_with_herd
          }
          
          sero.fit_glm_with_herd_agg = NULL 
          for (index_survey in 1:n.serosurveys){
            tadm=match(adm1s[[index_survey]], dat$adm0_adm1)
            vec_tmp = sero.fit_glm_with_herd_adm1[[index_survey]][,]*p_prop3d[tadm, study.years[[index_survey]]-1939,]*
              P_tot_2d[tadm, study.years[[index_survey]]-1939]
            
            if( length(adm1s[[index_survey]])>1 ) vec_tmp = colSums(vec_tmp)
            pop_tot_prov_surv = sum( P_tot_2d[tadm, study.years[[index_survey]]-1939])
            vec_tmp = vec_tmp/pop.agg3d[index_survey,study.years[[index_survey]]-1939,]
            sero.fit_glm_with_herd_agg = rbind(sero.fit_glm_with_herd_agg, vec_tmp)
          }
          sero.fit.glm[(1+((i/omit)-1)*n.serosurveys):(((i/omit))*n.serosurveys),] = sero.fit_glm_with_herd_agg
          
        }#if(print_predictions==T)
      }# if(i%%omit==0)
      
    }
    
    write.csv(mcmc.params,paste("mcmc_bb",bb,".csv",sep=""),row.names=FALSE)
    write.csv(R0_runs,paste("R0_bb",bb,".csv",sep=""),row.names=TRUE)
    write.csv(p.detect.link_runs,paste("p.detect.link_bb",bb,".csv",sep=""),row.names=FALSE)
    
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
  
