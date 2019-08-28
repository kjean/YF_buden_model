
create_status0 = function(x, beta0, n.serosurveys){

  parnames = c("vac_eff",paste("log",names(beta0),sep="."),paste("foi",sero.studies[1:n.serosurveys],sep="."), "vc.factor.CMRs")
  
  pars.ini = c(0.982,rep(NA,length(parnames)-1))
  names(pars.ini) =parnames
  ii = 2:(ncol(x)+1) 
 
  pars.ini[ii] = beta0
  pars.ini[(max(ii)+1):(max(ii)+6)] = rep(0.005, 6) # previous serosurveys
  
  #new serosurveys
  pars.ini[(max(ii)+7):(length(pars.ini)-1)] = 0.0005
  
  # vc.factor
  pars.ini[length(pars.ini)] = 0.6
  
  varsin.nc = c(1:ncol(x))[-grep("adm0",colnames(x))]# liste des var independantes moins les variables pays
  
  mu.prior = rep(NA, length(pars.ini)) # not really useful anymore
  mu.prior[grep("^log.adm0",names(pars.ini))] = 0 # not really useful anymore
  
  sd.prior=rep(2, length(pars.ini)) 
  
  transform = c(3,rep(0,length(beta0)),rep(1,n.serosurveys),2)
  
  sdj_vac_eff = 0.15
  
  sdj_glm = rep(1, length(beta0))
  
  # sdj_glm=c(0.32,0.24,1.88,5.37,4.88,2.69,
  #           3.84,2.63,2.69,4.44,2.22,4.88,
  #           2.89,3.76,1.88,3.84,0.39,0.0456,
  #           0.86,0.15)
  
  #sdj_foi= c(0.51,0.29,0.60,0.39,0.94,0.86) #old surveys
  sdj_foi=c(0.51,0.29,0.60,0.39,0.94,0.86,#old surveys
            rep(0.6, 25))
  sdj_CMRfac=0.61
  sdj=c(sdj_vac_eff,sdj_glm, sdj_foi, sdj_CMRfac)
  length(sdj)
  length(pars.ini)
  
  # Initialisating status
  status0 = list(params=pars.ini,lnL=rep(NA, 3+n.serosurveys),accept=0) # juste put an arbitrary high lnL for the very first round so that lnL is not NA
  
  lnLserosurvey=NULL
  for (i in 1:n.serosurveys) {lnLserosurvey[i] = paste("lnL", sero.studies[i], sep="_")}
  names(status0$lnL) = c( "lnL vac_eff","lnL GLM", lnLserosurvey, "lnL Prior")#names(status$lnL))
  status0$lnL[1:length(status0$lnL)] = -100000
  
  return(list(status0=status0, sdj=sdj, pars.ini=pars.ini,sd.prior=sd.prior, 
              mu.prior=mu.prior, transform=transform, varsin.nc=varsin.nc, ii =ii))
}