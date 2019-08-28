
###### hard code serological data

sero.studies = c("NGA","CAF","COD","COG","CMRn","CMRs",
                 "UGA_zone1","UGA_zone2", "UGA_zone3", "UGA_zone4", "UGA_zone5",
                 "RWA_zone1","RWA_zone2", "RWA_zone3", "RWA_zone4",
                 "ZMB_zone1", "ZMB_zone2",
                 "SDN_zone1", "SDN_zone2", "SDN_zone3", "SDN_zone4",
                 "KEN_zone1","KEN_zone2", "KEN_zone3", "KEN_zone4", "KEN_zone5",
                 "ETH_zone1","ETH_zone2", "ETH_zone3", "ETH_zone4", "ETH_zone5")
### NOTE that only data for "NGA","CAF","COD","COG","CMRn","CMRs" is publicly available
### the data for the further surveys is property of the WHO and are available from World Health Organization (contact: William Perea, pereaw@who.int). 
### data entered below is thus made up



n.serosurveys = length(sero.studies)

#study.years = list(1990,2009,1985,1985,1987,2000:2003)
# I fix study.year = 2002 pour CMRs - need to try if same results with 2001, etc
# I tried switching 2001 to 2002 and the result was really close
study.years = list(1990,2009,1985,1985,1987,2001,
                   2012,2012,2012,2012,2012,#UGA
                   2012,2012,2012,2012,# RWA
                   2013,2013,#ZMB
                   2012,2012,2012,2012, #SDN
                   2013,2013,2013,2013,2013,# KEN
                   2014,2014,2014,2014,2014) #ETH

length(study.years)
length(unlist(study.years))

#adm1s_gadm1 = list(2046,501:517,635,645,c(626,628),c(624,625,627,631,632)) # this is in gadm1

# in gadm2
adm1s_CAF = c(rep(NA, 17))
for (i in 1:17) {adm1s_CAF[i]<- paste("CAF", i, sep="_")}
adm1s_CAF =sort(adm1s_CAF)
adm1s = list("NGA_31", adm1s_CAF, "COD_3", "COG_4", c("CMR_4", "CMR_7"), c("CMR_2","CMR_3","CMR_5","CMR_10","CMR_9"),
             c("UGA_12","UGA_38") , c("UGA_6","UGA_40"),c("UGA_2","UGA_22"),c("UGA_17", "UGA_39"), "UGA_16",
             "RWA_3",c("RWA_8","RWA_10"),c("RWA_2","RWA_9"),"RWA_7",
             "ZMB_6", "ZMB_7",
             "SDN_6", "SDN_3", c("SDN_4", "SDN_6"), "SDN_1",
             c("KEN_3", "KEN_5"), c("KEN_5","KEN_7"),"KEN_8","KEN_7", "KEN_4",
             "ETH_9", c("ETH_4", "ETH_10"), "ETH_8", c("ETH_3", "ETH_8"), "ETH_8")
rm(adm1s_CAF)

length(adm1s)
length(unlist(adm1s))

sero.dat = list()
sero.dat[[1]] = matrix(c(3,50,105,26,0,2,24,11),ncol=2)
sero.dat[[2]] = matrix(c(91,109,107,453,149,29,3,12,10,59,34,7),ncol=2)
sero.dat[[3]] = matrix(c(14,42,32,26,26,3,18,15,13,14),ncol=2)
sero.dat[[4]] = matrix(c(0,360,0,0,72,0),ncol=2) 
sero.dat[[5]] = matrix(c(586+254,0,10+17,0),ncol=2) 
sero.dat[[6]] = matrix(c(0,73,51,40,28,64,0,3,4,9,14,39),ncol=2)


#### Warning:
### from here up to line 116, number of positive samples is made up
sero.dat[[7]] = matrix(c(19,65,23,14,
                         1,1,1,1),ncol=2)
sero.dat[[8]] = matrix(c(22,65,33,1,
                         1,1,1,1),ncol=2)
sero.dat[[9]] = matrix(c(35,55,24,9,
                         1,1,1,1),ncol=2)
sero.dat[[10]] = matrix(c(40,42,23,9,
                          1,1,1,1),ncol=2)
sero.dat[[11]] = matrix(c(11,83,9,2,
                          1,1,1,1),ncol=2)


sero.dat[[12]] = matrix(c(94,76,25,4,
                          1,1,1,1),ncol=2)
sero.dat[[13]] = matrix(c(113,79,57,11,
                          1,1,1,1),ncol=2)
sero.dat[[14]] = matrix(c(84,105,23,3,
                          1,1,1,1),ncol=2)
sero.dat[[15]] = matrix(c(270,227,85,30,
                          1,1,1,1),ncol=2)

sero.dat[[16]] = matrix(c(156,395,354,351,249,282,
                          1,1,1,1,1,1),ncol=2)
sero.dat[[17]] = matrix(c(210,489,396,253,218,326,
                          1,1,1,1,1,1),ncol=2)

sero.dat[[18]] = matrix(c(40,362,232,20,
                          1,1,1,1),ncol=2)
sero.dat[[19]] = matrix(c(33,322,243,57,
                          1,1,1,1),ncol=2)
sero.dat[[20]] = matrix(c(18,133,61,9,
                          1,1,1,1),ncol=2)
sero.dat[[21]] = matrix(c(81,133,61,9,
                          1,1,1,1),ncol=2)

sero.dat[[22]] = matrix(c(72,175,177,
                          1,1,1),ncol=2)
sero.dat[[23]] = matrix(c(66,159,238,
                          1,1,1),ncol=2)
sero.dat[[24]] = matrix(c(55,112,219,
                          1,1,1),ncol=2)
sero.dat[[25]] = matrix(c(46,149,236,
                          1,1,1),ncol=2)
sero.dat[[26]] = matrix(c(48,104,104,
                          1,1,1),ncol=2)

sero.dat[[27]] = matrix(c(46,67,31,8,
                          1,1,1,1),ncol=2)
sero.dat[[28]] = matrix(c(21,35,7,1,
                          1,1,1,1),ncol=2)
sero.dat[[29]] = matrix(c(191,337,103,25,
                          1,1,1,1),ncol=2)
sero.dat[[30]] = matrix(c(39,49,26,2,
                          1,1,1,1),ncol=2)
sero.dat[[31]] = matrix(c(197,342,107,11,
                          1,1,1,1),ncol=2)




n.age.groups=NULL

for(i in 1:n.serosurveys) {
  n.age.groups[i] = nrow(sero.dat[[i]])
}
age.min = list(c(0,10,20,40),c(0,5,10,15,40,65),c(0,15,30,40,50),c(0,6,16),c(0,14),c(0,16,26,36,46,56),
               c(0,15,40,65),c(0,15,40,65),c(0,15,40,65),c(0,15,40,65),c(0,15,40,65), #UGA
               c(0,15,40,65),c(0,15,40,65),c(0,15,40,65),c(0,15,40,65), # RWA
               c(0,5,15,25,35,45), c(0,5,15,25,35,45), # ZMB
               c(0,15,40,65),c(0,15,40,65),c(0,15,40,65),c(0,15,40,65), #SDN
               c(0,5,20), c(0,5,20), c(0,5,20), c(0,5,20), c(0,5,20), #KEN
               c(0,15,40,65),c(0,15,40,65),c(0,15,40,65),c(0,15,40,65),c(0,15,40,65))

vc.factor = c(0,0,1,1,1,NA,  rep(0,25)) # factor indicating if we should account for vaccination for estimating FOI from serosurveys

names(vc.factor) = paste("vc.factor.",sero.studies[1:n.serosurveys],sep="")


#il faut charger vc1 avant ca

if(FALSE){
  t0_vac = adm1s
  for (i in 1:n.serosurveys){
    print(i)
    for (j in 1:length(adm1s[[i]]) ){
      print(j)
      year_i = 1940
      sum_vac=0
      while(sum_vac == 0 & year_i<2051){
        sum_vac = sum(vc1[vc1$adm0_adm1==adm1s[[i]][j] & vc1$year==year_i, 4:104], na.rm=T)
        year_i = year_i + 1
      }
      t0_vac[[i]][j]=year_i-1
    }
  }
}



#hard coding the first year of vaccination
#t0_vac = list(1988,rep(1946,17),1972,1946, c(1946,1946),c(1946,1946,1946,1946,1946)) 
t0_vac = c(1988,1946,1972,1946,1946,1946,
           2050,2050,2050,2050,
           2050,2050,2050,2050,2050,
           2050,2050,2050,2050,2050,
           2015,2003,1994,2003,1994,
           2003,2050,2014,2050,2050,
           2050) # same t0 for all provinces within the same survey

names(t0_vac)=sero.studies
t0_vac_sauv = t0_vac

save(list=ls(), file = "serological_data.Rdata")



### calculate overall seroprevalence
library(Hmisc)

glob.seroprev = data.frame(matrix(rep(NA, 3*n.serosurveys), ncol=3))
for (i in 1:length(sero.dat)){
  total.eff = colSums(sero.dat[[i]])
  seroprev = binconf(total.eff[2], total.eff[1])
  glob.seroprev[i,] = seroprev
}
glob.seroprev = cbind(sero.studies,vc.factor,glob.seroprev)
summary(glob.seroprev[,1])
sero.studies[which(glob.seroprev[,3] == min(glob.seroprev[,3]))]
sero.studies[22]
glob.seroprev[22,]
study.years[22]

sero.studies[3]
glob.seroprev[3,]
study.years[3]
