
create_status0 = function(x, beta0, n.serosurveys){

  parnames=  c("vac_eff", paste("log",names(beta0),sep="."), paste("R0", sero.studies, sep="_"), "vc.factor.CMRs")  
  pars.ini = rep(NA,length(parnames))
  names(pars.ini) = parnames
  
  pars.ini[1] = 0.98
  ii = grep("log", parnames)
  pars.ini[ii] = beta0 
  
  pars.ini[(max(grep("log", parnames))+1):(max(grep("log", parnames))+6)] = c(1.2,1.5, 1.4, 1.8, 1.3, 1.2)
  pars.ini[grep("zone", parnames)] = rep(1.01, 25)
  
  pars.ini[length(pars.ini)]=0.60 #vc.factor.CMRs
  
  status0 = list(params=pars.ini,lnL=rep(-1000, n.serosurveys+2),accept=0) 
  names(status0$lnL) = c("lnL_GLM" ,paste("lnL", sero.studies, sep="_"), "lnL_prior")
  
  varsin.nc = ii[-grep("adm0",colnames(x))] - 1 
  
  if(file.exists("sdj_tweaking.csv")){
    sdj = read.csv("sdj_tweaking.csv", h= T)
    sdj=as.numeric(sdj[nrow(sdj),-1])
  }
  else{
    # vc.factor
    sdj_vac_eff = 0.2
    sdj_glm= rep(1.5, length(beta0))
    sdj_glm[18]= 0.5
    sdj_r0= rep(0.03, n.serosurveys)
    sdj_r0[2]= 0.02
    
    sdj_r0=c(0.12, 0.02, 0.2, 0.11, 0.05, 0.08, #old serosurveys
             0.005, 0.03, 0.098, 0.073, 0.087, #UGA
             0.03,0.005, 0.014, 0.001, #RWA
             0.007,0.008, #ZMB
             0.041, 0.028, 0.055, 0.055,#SDN
             0.0001, 0.06, 0.0003, 0.02, 0.028,#KEN
             0.061, 0.003, 0.035, 0.0002, 0.01)#ETH
    sdj_CMRfac=1.5
    
    sdj = c(sdj_vac_eff, sdj_glm, sdj_r0, sdj_CMRfac)
    
  }
    
  transform = c(3,rep(0,length(ii)),rep(1,n.serosurveys), 2) # penser a ajouter ensuite vc.fact (transform =2)
  length(transform)
  sd.prior=rep(2, length(pars.ini)) 
  
  return(list(status0=status0, sdj=sdj, pars.ini=pars.ini,sd.prior=sd.prior, 
              transform=transform, varsin.nc=varsin.nc, ii =ii))
}