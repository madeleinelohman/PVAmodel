library("popbio")

####finds relationship between density and seed production
#x is population/m^2 and y is the respective seedoutput
x <- data.frame(x=c(660,4000,178,4850,844.7,1372.67,1012.26,1772.87,692.51,1050.21,308.58,612.68),
                y=c(10.7,1.75,10,1.61,10.56,4.95,12.98,5,9.28,4.89,10.47,4.5))
reprodensityeq <- lm(x$y~log(x$x))
reprodensity <- reprodensityeq$coefficients
#highest population seed outputs, used if the population is above K
otherwise <- mean(c(1.75,1.61))
plot(x$y~log(x$x),xlab=expression("Plants per"~m^{2}),ylab="Seeds per plant",pch=19,ylim=c(0,20),
     main="Density Dependent Seed Production")
abline(reprodensityeq,xpd=F)


#unneeded to run if environment loaded
d <- read.csv("CONTROL.csv",header = T)

#####derives control values#####
controlf <- matrix(nrow=6,ncol=3)
colnames(controlf) <- c("parameter","mean","sd")
controlf <- as.data.frame(controlf)

controlf$parameter[1] <- "survival"
controlf$mean[1] <- mean(d$SURVIVAL,na.rm = T)/100
controlf$sd[1] <- sd(d$SURVIVAL,na.rm = T)/100

controlf$parameter[2] <- "adults"
controlf$mean[2] <- mean(d$ADULTS,na.rm = T)
controlf$sd[2] <- sd(d$ADULTS,na.rm = T)

controlf$parameter[3] <- "cover"
controlf$mean[3] <- mean(d$COVER,na.rm = T)
controlf$sd[3] <- sd(d$COVER,na.rm = T)

controlf$parameter[4] <- "seeds"
controlf$mean[4] <- mean(d$SEEDS,na.rm = T)
controlf$sd[4] <- sd(d$SEEDS,na.rm = T)

controlf$parameter[5] <- "seedbank"
controlf$mean[5] <- mean(d$SEEDBANK,na.rm = T)
controlf$sd[5] <- sd(d$SEEDBANK,na.rm = T)

controlf$parameter[6] <- "lambda"
controlf$mean[6] <- mean(d$LAMBDA,na.rm = T)
controlf$sd[6] <- sd(d$LAMBDA,na.rm = T)


###Needed to specify for throughout the model
K <- 4850
e <- exp(1)
nreps <- 500


####Function to determine mean and sd of parameters and set up certain parameters###

treatment <- function(nreps,dpop,dseed){
  
  data <- matrix(nrow=nreps,ncol=11)
  colnames(data) <- c("controlseeds","controlseedbank","controlN","controlsurvival","controllambda","DeltaPop",
                      "DeltaSeeds","survival","reproduction","Lambda","R")
  data <- as.data.frame(data)
  
  for(rep in 1:nreps){
    
    #for control values
    data$controlseeds[[rep]] <- max(0,rnorm(1,controlf$mean[4],controlf$sd[4])) 
    data$controlseedbank[[rep]] <- max(0,rnorm(1,controlf$mean[5],controlf$sd[5]))
    data$controlN[[rep]] <- max(0,rnorm(1,controlf$mean[2],controlf$sd[2]))
    data$controlsurvival[[rep]] <- max(0,rnorm(1,controlf$mean[1],controlf$sd[1]))
    data$controllambda[[rep]] <- max(0,rnorm(1,controlf$mean[6],controlf$sd[6]))
    
    ##This section gives uncertainty to parameters that describe change in population and seedbank
    #some parameters only have one value
    if(length(dpop)==1){
      data$DeltaPop[[rep]] <- max(0,rnorm(1,mean(dpop),0.15))
      } else{
        data$DeltaPop[[rep]] <- max(0,rnorm(1,mean(dpop),sd(dpop)))
    }
    
    if(length(dseed)==1){
      data$DeltaSeeds[[rep]] <- max(0,rnorm(1,mean(dseed),0.15))
    }
    else(data$DeltaSeeds[[rep]] <- max(0,rnorm(1,mean(dseed),sd(dseed))))
  }
  
  return(data)
}

#Proportion of population/seedbank present after treatment
grazepop <- c(0,0.25,0.36,0.72,0.21)
grazeseed <- c(0.21,0.52,0.12,0.19)

burnpop <- c(0.63,0,1.92,0.95,0.63,0.6,0.69,0.43,1.04,0.66)
burnseed <- c(1.13,5.46,2.63,0.78)

grazeburnpop <- 0.1
grazeburnseed <- c(0.55,0.5)

herbpop <- c(0.33,0.01)
herbseed <- c(1)

#Putting treatment function values into variables
control <- treatment(500,c(1,1),c(1,1))
csd <- mean(control$DeltaSeeds)
csdsd <- sd(control$DeltaSeeds)
cpd <- mean(control$DeltaPop)
cpdsd <- sd(control$DeltaPop)
cl <- mean(control$Lambda)
clsd <- sd(control$Lambda)
csurv <- mean(control$controlsurvival)
csurvsd <- sd(control$controlsurvival)
initseedbank <- mean(control$controlseedbank)
initseedbanksd <- sd(control$controlseedbank)
initabund <- mean(control$controlN)
initabundsd <- sd(control$controlN)
seeds <- mean(control$controlseeds)
seedssd <- sd(control$controlseeds)

grazing <- treatment(500,grazepop,grazeseed)
gsd <- mean(grazing$DeltaSeeds)
gsdsd <- sd(grazing$DeltaSeeds)
gpd <- mean(grazing$DeltaPop)
gpdsd <- sd(grazing$DeltaPop)

burning <- treatment(500,burnpop,burnseed)
bsd <- mean(burning$DeltaSeeds)
bsdsd <- sd(burning$DeltaSeeds)
bpd <- mean(burning$DeltaPop)
bpdsd <- sd(burning$DeltaPop)

grazeburn <- treatment(500,grazeburnpop,grazeburnseed)
gbsd <- mean(grazeburn$DeltaSeeds)
gbsdsd <- sd(grazeburn$DeltaSeeds)
gbpd <- mean(grazeburn$DeltaPop)
gbpdsd <- sd(grazeburn$DeltaPop)

herbicide <- treatment(500,herbpop,herbseed)
hsd <- mean(herbicide$DeltaSeeds)
hsdsd <- sd(herbicide$DeltaSeeds)
hpd <- mean(herbicide$DeltaPop)
hpdsd <- sd(herbicide$DeltaPop)

###########################################################
#PVA
###########################################################

#Describes population growth
Population <- function(lastyear,prev_abund,DeltaPop,DeltaPopsd){ 
  
  mylist3 <- list()
  
  #last year's seedbank calculated without year to year changes
  mylist3$y <- max(0,mylist3$z*(reprodensity[[2]]*log(prev_abund)+reprodensity[[1]]))
  #seedbank survival to adulthood to find population
  z <- max(0,lastyear*rnorm(1,csurv,csurvsd))
  #the adult population proportion present after treatment
  mylist3$x <- max(0,z*rnorm(1,DeltaPop,DeltaPopsd))
  
  if(mylist3$x > 4850){
    mylist3$x <- 4850
    mylist3$y <- lastyear
    return(mylist3)
  }
  
  if(mylist3$x==0){
    mylist3$x <- prev_abund
    mylist3$y <- lastyear
    return(mylist3)
  }
  
  return(mylist3)
}

Seedbank <- function(prev_abund,lastyear,DeltaSeeds,DeltaSeedssd){
  
  mylist2 <- list()
  
  if(prev_abund > 20000){
    mylist2$x <- otherwise*lastyear
    mylist2$z <- lastyear
    return(mylist2)
  }
  
  #find adult population
  mylist2$z <- max(0,rnorm(1,csurv,csurvsd)*prev_abund)
  
  #seedbank; seed output multiplied by population
  y <- max(0,mylist2$z*(reprodensity[[2]]*log(lastyear)+reprodensity[[1]]))
  
  #change in seedbank
  mylist2$x <- max(0,y*rnorm(1,DeltaSeeds,DeltaSeedssd))
  
  if(is.nan(mylist2$x)|is.infinite(mylist2$x)){mylist2$x <- 0}
  
  if(mylist2$x==0){
    mylist2$x <- otherwise*lastyear
  }
  
  return(mylist2)
}


mylist <- list()

PVA <- function(nreps,nyears,Init_N,Init_Nsd,DeltaSeeds,DeltaSeedssd,DeltaPop,DeltaPopsd){
  
  PopArray2 <- matrix(0,nrow=nyears+1,ncol=nreps)
  Seeds <- matrix(0,nrow=nyears+1,ncol=nreps)
  Seedpop <- matrix(0,nrow=nyears+1,ncol=nreps)
  Popseed <- matrix(0,nrow=nyears+1,ncol=nreps)
  
  for(rep in 1:nreps){
    
    PopArray2[1,rep] <- max(0,rnorm(1,Init_N,Init_Nsd))
    Seeds[1,rep] <- max(0,rnorm(1,initseedbank,initseedbanksd))
    Seedpop[1,rep] <- max(0,rnorm(1,Init_N,Init_Nsd))
    Popseed[1,rep] <- max(0,rnorm(1,initseedbank,initseedbanksd))
    
    
    for(y in 2:(nyears+1)){
      nextyear <- Population(Popseed[y-1,rep],PopArray2[y-1,rep],DeltaPop=DeltaPop,DeltaPopsd=DeltaPopsd)
      PopArray2[y,rep] <- nextyear$x
      Popseed[y,rep] <- nextyear$y
      nextyearseed <- Seedbank(Seeds[y-1,rep],Seedpop[y-1,rep],DeltaSeeds=DeltaSeeds,DeltaSeedssd=DeltaSeedssd)
      Seeds[y,rep] <- nextyearseed$x
      Seedpop[y,rep] <- nextyearseed$z
    }
  }
  mylist$poparray <- PopArray2
  mylist$seeds <- Seeds
  return(mylist)
}


####CONTROL####
EndControl <- PVA(nreps=500,nyears=5,Init_N=initabund,Init_Nsd=initabundsd,DeltaSeeds=csd,
                  DeltaSeedssd=0.15,DeltaPop=cpd,DeltaPopsd=0.15)



fakecontrol <- matrix(NA,nrow=6,ncol=2)
for(i in 1:nrow(EndControl$poparray)){
  mean <- mean(EndControl$poparray[i,])
  sd <- sd(EndControl$poparray[i,])
  n <- 500
  error <- qnorm(0.975)*sd[i]/sqrt(n)
  fakecontrol[i,1] <- mean-error
  fakecontrol[i,2] <- mean+error
}

controlnew <- NULL
for(i in 1:nrow(EndControl$poparray)){
  controlnew[i] <- mean(EndControl$poparray[i,])
}
time <- 1:6
plot(controlnew~time,ylim=c(700,4500),xaxp=c(1,6,5),xlab="Year",ylab=expression("Plants per"~m^{2}),pch=19,
     col="black",main="Control Population Projection",cex=1)
arrows(time,fakecontrol[,1],time,fakecontrol[,2],length=0)

endpop <- mean(EndControl$poparray[nrow(EndControl$poparray),])
endseed <- mean(EndControl$seeds[nrow(EndControl$seeds),])

acontrol <- endpop
scontrol <- sd(EndControl$poparray[nrow(EndControl$poparray),])
ncontrol <- nreps
errorcontrol <- qnorm(0.975)*scontrol/sqrt(ncontrol)
leftcontrol <- acontrol-errorcontrol
rightcontrol <- acontrol+errorcontrol
endci <- c(leftcontrol,rightcontrol)

acontrols <- endseed
scontrols <- sd(EndControl$seeds[nrow(EndControl$seeds),])
ncontrols <- nreps
errorcontrols <- qnorm(0.975)*scontrols/sqrt(ncontrols)
leftcontrols <- acontrols-errorcontrols
rightcontrols <- acontrols+errorcontrols
endcis <- c(leftcontrols,rightcontrols) 


############Grazing#############

GrazeDefault <- PVA(nreps=500,nyears=5,Init_N=initabund,Init_Nsd=initabundsd,DeltaSeeds=gsd,DeltaSeedssd=gsdsd,
                    DeltaPop=gpd,DeltaPopsd=gpdsd)

grazeendpop <- mean(GrazeDefault$poparray[nrow(GrazeDefault$poparray),])
grazeendseed <- mean(GrazeDefault$seeds[nrow(GrazeDefault$seeds),])


agraze <- grazeendpop
sgraze <- sd(GrazeDefault$poparray[nrow(GrazeDefault$poparray),])
ngraze <- nreps
errorgraze <- qnorm(0.975)*sgraze/sqrt(ngraze)
leftgraze <- agraze-errorgraze
rightgraze <- agraze+errorgraze
cigraze <- c(leftgraze,rightgraze)


agrazes <- grazeendseed
sgrazes <- sd(GrazeDefault$seeds[nrow(GrazeDefault$seeds),])
ngrazes <- nreps
errorgrazes <- qnorm(0.975)*sgrazes/sqrt(ngrazes)
leftgrazes <- agrazes-errorgrazes
rightgrazes <- agrazes+errorgrazes
cigrazes <- c(leftgrazes,rightgrazes)


############Burning#############

BurnDefault <- PVA(nreps=500,nyears=5,Init_N=initabund,Init_Nsd=initabundsd,DeltaSeeds=bsd,DeltaSeedssd=bsdsd,
                   DeltaPop=bpd,DeltaPopsd=bpdsd)


burnendpop <- mean(BurnDefault$poparray[nrow(BurnDefault$poparray),])
burnendseed <- mean(BurnDefault$seeds[nrow(BurnDefault$seeds),])


aburn <- burnendpop
sburn <- sd(BurnDefault$poparray[nrow(BurnDefault$poparray),])
nBurn <- nreps
errorburn <- qnorm(0.975)*sburn/sqrt(nBurn)
leftburn <- aburn-errorburn
rightburn <- aburn+errorburn
ciburn <- c(leftburn,rightburn)

aburns <- burnendseed
sburns <- sd(BurnDefault$seed[nrow(BurnDefault$seed),])
nBurns <- nreps
errorburns <- qnorm(0.975)*sburns/sqrt(nBurns)
leftburns <- aburns-errorburns
rightburns <- aburns+errorburns
ciburns <- c(leftburns,rightburns)


############Grazing/Burning#############

GrazeBurnDefault <- PVA(nreps=500,nyears=5,Init_N=initabund,Init_Nsd=initabundsd,DeltaSeeds=gbsd,
                        DeltaSeedssd=gbsdsd,DeltaPop=gbpd,DeltaPopsd=gbpdsd)

fakegb <- matrix(NA,nrow=6,ncol=2)
for(i in 1:nrow(GrazeBurnDefault$poparray)){
  mean <- mean(GrazeBurnDefault$poparray[i,])
  sd <- sd(GrazeBurnDefault$poparray[i,])
  n<- nreps
  error <- qnorm(0.975)*sd/sqrt(n)
  fakegb[i,1] <- mean-error
  fakegb[i,2] <- mean+error
}

gbv <- NULL
for(i in 1:nrow(GrazeBurnDefault$poparray)){
  gbv[i] <- mean(GrazeBurnDefault$poparray[i,])
}
plot(gbv~time,type="p",xaxp=c(1,6,5),xlab="Year",ylab="Plants per meter^2",pch=19,
     col="black",main="Graze/Burn Population Projection",ylim=c(-4,300))
arrows(time,fakegb[,1],time,fakegb[,2],length=0)

new <- GrazeBurnDefault$seeds[6,]>480
newefake <- which(GrazeBurnDefault$poparray[6,]<480)

grazeburnendpop <- mean(GrazeBurnDefault$poparray[nrow(GrazeBurnDefault$poparray),])
grazeburnendseed <- mean(GrazeBurnDefault$seeds[nrow(GrazeBurnDefault$poparray),])

agrazeburn <- grazeburnendpop
sgrazeburn <- sd(GrazeBurnDefault$poparray[nrow(GrazeBurnDefault$poparray),])
ngrazeburn <- nreps
errorgrazeburn <- qnorm(0.975)*sgrazeburn/sqrt(ngrazeburn)
leftgrazeburn <- agrazeburn-errorgrazeburn
rightgrazeburn <- agrazeburn+errorgrazeburn
cigrazeburn <- c(leftgrazeburn,rightgrazeburn)

agrazeburns <- grazeburnendseed
sgrazeburns <- sd(GrazeBurnDefault$seeds[nrow(GrazeBurnDefault$seeds),])
ngrazeburns <- nreps
errorgrazeburns <- qnorm(0.975)*sgrazeburns/sqrt(ngrazeburns)
leftgrazeburns <- agrazeburns-errorgrazeburns
rightgrazeburns <- agrazeburns+errorgrazeburns
cigrazeburns <- c(leftgrazeburns,rightgrazeburns)


############Herbicide#############

HerbicideDefault <- PVA(nreps=500,nyears=5,Init_N=initabund,Init_Nsd=initabundsd,DeltaSeeds=hsd,
                        DeltaSeedssd=hsdsd,DeltaPop=hpd,DeltaPopsd=hpdsd)


herbendpop <- mean(HerbicideDefault$poparray[nrow(HerbicideDefault$poparray),])
herbendseed <- mean(HerbicideDefault$seeds[nrow(HerbicideDefault$seeds),])


aherb <- herbendpop
sherb <- sd(HerbicideDefault$poparray[nrow(HerbicideDefault$poparray),])
nherb <- nreps
errorherb <- qnorm(0.975)*sherb/sqrt(nherb)
leftherb <- aherb-errorherb
rightherb <- aherb+errorherb
ciherb <- c(leftherb,rightherb)

aherbs <- herbendseed
sherbs <- sd(HerbicideDefault$seed[nrow(HerbicideDefault$seed),])
nherbs <- nreps
errorherbs <- qnorm(0.975)*sherbs/sqrt(nherbs)
leftherbs <- aherbs-errorherbs
rightherbs <- aherbs+errorherbs
ciherbs <- c(leftherbs,rightherbs)




mydfpop <- c(EndControl$poparray[nrow(EndControl$poparray),],
             BurnDefault$poparray[nrow(BurnDefault$poparray),],
             GrazeDefault$poparray[nrow(GrazeDefault$poparray),],
             GrazeBurnDefault$poparray[nrow(GrazeBurnDefault$poparray),],
             HerbicideDefault$poparray[nrow(HerbicideDefault$poparray),]
)
mylistpop <- c("Control","Graze","Burn","GrazeBurn","Herbicide")
mydf1pop <- data.frame(rep(NA,7500),rep(mylistpop,each=1500))
mydf1pop[,1] <- mydfpop
colnames(mydf1pop) <- c("values","treatment")
mynewpop <- aov(values ~ treatment, data = mydf1pop)
summary(mynewpop)


popci <- c(endci,cigraze,ciburn,cigrazeburn,ciherb)
popcomparepoints <- c(endpop,grazeendpop,burnendpop,grazeburnendpop,herbendpop)
pop <- data.frame(c(popci,popcomparepoints))
popcomparetreatment <- c("Control","Grazing","Burning",paste("Grazing","and","Burning",sep="\n"),"Herbicide")

barCenters <- barplot(popcomparepoints,names.arg = popcomparetreatment,beside = true, las = 2,cex.names = 1, 
                      main = "Population after 5 Years",ylab = expression("Plants per"~m^{2}),cex=0.9,border = "black",
                      ylim=c(0,max(popci)))

segments(barCenters, popci[c(2,4,6,8,10)], barCenters,
         popci[c(1,3,5,7,9)], lwd = 1.5)

arrows(barCenters, popci[c(2,4,6,8,10)], barCenters,
       popci[c(1,3,5,7,9)], lwd = 1.5, angle = 90,
       code = 3, length = 0.05)



mydfs <- c(EndControl$seeds[nrow(EndControl$seeds),],
           GrazeDefault$seeds[nrow(GrazeDefault$seeds),],
           BurnDefault$seeds[nrow(BurnDefault$seeds),],
           GrazeBurnDefault$seeds[nrow(GrazeBurnDefault$seeds),],
           HerbicideDefault$seeds[nrow(HerbicideDefault$seeds),]
)
mylist <- c("Control","Graze","Burn","GrazeBurn","Herbicide")
mydf1s <- data.frame(rep(NA,7500),rep(mylist,each=1500))
mydf1s[,1] <- mydfs
colnames(mydf1s) <- c("values","treatment")
mynews <- aov(values ~ treatment, data = mydf1)
summary(mynews)

3

seedci <- c(endcis,cigrazes,ciburns,cigrazeburns,ciherbs)
seedcomparepoints <- c(endseed,grazeendseed,burnendseed,grazeburnendseed,herbendseed)
seedcomparetreatment <- c("Control","Grazing","Burning",paste("Grazing","and","Burning",sep="\n"),"Herbicide")

barCenterss <- barplot(seedcomparepoints,names.arg=seedcomparetreatment,beside = true, las = 2,cex.names = 1, 
                       main = "Seedbank after 5 Years",ylab = expression("Seeds per"~m^{2}),cex=0.9,border = "black",
                       ylim=c(0,max(seedci)))

segments(barCenters, seedci[c(2,4,6,8,10)], barCenters,
         seedci[c(1,3,5,7,9)], lwd = 1.5)

arrows(barCenters, seedci[c(2,4,6,8,10)], barCenters,
       seedci[c(1,3,5,7,9)], lwd = 1.5, angle = 90,
       code = 3, length = 0.05)


