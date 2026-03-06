rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(ieggr)
library(reshape2)
library(boot)
library(lavaan)
library(tibble)
library(NST)
save.wd <- iwd(choose.dir())

### Calculation of beta diversity indices -------------------------------------------------------------------
com.file = "OTU-16S-09-16-Rarefied.csv" #Prokaryotes
com.file = "OTU-ITS-09-16-Rarefied.csv" #Fungi
com.file = "Geochip-09-16.csv" #Functional genes

treat.file = "Treatment-Env-09-16.csv"
comm = read.csv(com.file, row.names = 1, header = T, sep = ",")
treat = read.csv(treat.file, row.names = 1, header = T, sep = ",")

### After column 64, Listed as annotated information
comm = comm[,1:64]
comm = as.data.frame(t(comm))

dim(comm)
comm = comm[, colSums(comm) > 0]
dim(comm)

dim(comm)
comm = comm[rowSums(comm) > 0,]
dim(comm)

### Align OTU tables, and treatment information
samp = match.name(rn.list = list(treat = treat, comm = comm))
comm = samp$comm
treat = samp$treat
#comm <- comm / rowSums(comm)
#comm <- sqrt(comm)

### Bray–Curtis distance
dist.bray = vegdist(comm)
### Sorensen distance
dist.soren = vegdist(comm, binary = TRUE)

### Non-parametric test (Adonis) ------------------------------------------------------
### whole community ###
adonis = adonis2(dist.bray ~ Precipitation*Year + block, data = treat, permutations = 999)

### Each year ###
#result <- c()
prefix <- "GeoChip" #Prokaryotes,Fungi,GeoChip
for (Year in unique(treat$Year)){
  dist.bray.year = vegdist(comm[treat$Year==Year,])
  adonis = adonis2(dist.bray.year ~ Precipitation + block, data = treat[treat$Year==Year,], permutations = 999)
  result1 <- adonis["Precipitation",]
  row.names(result1) <- paste0(Year,"_",prefix) 
  result <- rbind(result,result1)
}

write.csv(result, file = "2-Adonis-year.csv")

### ANOVA test ------------------------------------------------------
### Environment ###
#result <- c()
for (i in c(11:27,ncol(treat))){
  test <- data.frame(Env=treat[,i],treat[,1:6])
  test <- na.omit(test)
  mean <- tapply(test$Env, test$Precipitation, mean)
  se <- tapply(test$Env, test$Precipitation, sd)/ sqrt(tapply(test$Env, test$Precipitation, length))
  
  if (length(unique(test$Year))>1){
    result1 <- summary(aov(Env ~ Precipitation*Year + block, data = test))
    output <- data.frame(Taxa=colnames(treat)[i],Year="Y09-Y16",
                         Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                         Item="Year",Fvalue=result1[[1]]$`F value`[2],Pvalue=result1[[1]]$`Pr(>F)`[2])
    result <- rbind(result,output)
    output <- data.frame(Taxa=colnames(treat)[i],Year="Y09-Y16",
                         Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                         Item="Treat",Fvalue=result1[[1]]$`F value`[1],Pvalue=result1[[1]]$`Pr(>F)`[1])
    result <- rbind(result,output)
    output <- data.frame(Taxa=colnames(treat)[i],Year="Y09-Y16",
                         Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                         Item="Year*Treat",Fvalue=result1[[1]]$`F value`[4],Pvalue=result1[[1]]$`Pr(>F)`[4])
    result <- rbind(result,output)
  }

  for (Year in unique(test$Year)){
    test1 <- test[test$Year==Year,]
    mean <- tapply(test1$Env, test1$Precipitation, mean)
    se <- tapply(test1$Env, test1$Precipitation, sd)/ sqrt(tapply(test1$Env, test1$Precipitation, length))
    
    result1 <- summary(aov(Env ~ Precipitation + block, data = test1))
    output <- data.frame(Taxa=colnames(treat)[i],Year=Year,
                         Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                         Item="Treat",Fvalue=result1[[1]]$`F value`[1],Pvalue=result1[[1]]$`Pr(>F)`[1])
    result <- rbind(result,output)
  }
  row.names(result) <- NULL
} 

write.csv(result, file = "1-ANOVA-Envmean.csv")

### Relative abundance of microbe ###
comm = comm[,c(1:64,ncol(comm))]
comm1 <- comm %>% group_by(Taxa) %>% summarise(across(1:64, sum))
comm1 <- column_to_rownames(comm1, var = names(comm1)[1])
comm1 <- comm1/colSums(comm1)
comm1 <- as.data.frame(t(comm1))

#result <- c()
for (i in 1:ncol(comm1)){
  test <- data.frame(RA=comm1[,i],treat)
  mean <- tapply(test$RA, test$Precipitation, mean)
  se <- tapply(test$RA, test$Precipitation, sd)/ sqrt(tapply(test$RA, test$Precipitation, length))
  
  result1 <- summary(aov(RA ~ Precipitation*Year + block, data = test))
  output <- data.frame(Taxa=colnames(comm1)[i],Year="Y09-Y16",
                       Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                       Item="Year",Fvalue=result1[[1]]$`F value`[2],Pvalue=result1[[1]]$`Pr(>F)`[2])
  result <- rbind(result,output)
  output <- data.frame(Taxa=colnames(comm1)[i],Year="Y09-Y16",
                        Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                        Item="Treat",Fvalue=result1[[1]]$`F value`[1],Pvalue=result1[[1]]$`Pr(>F)`[1])
  result <- rbind(result,output)
  
  for (Year in unique(test$Year)){
    test1 <- test[test$Year==Year,]
    mean <- tapply(test1$RA, test1$Precipitation, mean)
    se <- tapply(test1$RA, test1$Precipitation, sd)/ sqrt(tapply(test1$RA, test1$Precipitation, length))
    
    result1 <- summary(aov(RA ~ Precipitation + block, data = test1))
    output <- data.frame(Taxa=colnames(comm1)[i],Year=Year,
                         Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                         Item="Treat",Fvalue=result1[[1]]$`F value`[1],Pvalue=result1[[1]]$`Pr(>F)`[1])
    result <- rbind(result,output)
  }
  row.names(result) <- NULL
} 
  
write.csv(result, file = "2-ANOVA-RelAbund.csv")

### Relative abundance of genes ###
comm1 <- comm[comm$Subcategory1=="Carbon degradation",c(1:64,ncol(comm))]
comm1 <- comm[comm$Subcategory1=="Carbon fixation",c(1:64,ncol(comm))]
comm1 <- comm[comm$Gene_category=="Nitrogen",c(1:64,ncol(comm))]
comm1 <- comm[comm$Gene_category=="Phosphorus",c(1:64,ncol(comm))]
comm1 <- comm[comm$Gene_category=="Sulfur",c(1:64,ncol(comm))]

colnames(comm1)[ncol(comm1)] <- "Class"
comm1 <- na.omit(replace(comm1, comm1 == "", NA))
prefix = "Sulfur"

comm1 <- comm1 %>% group_by(Class) %>% summarise(across(1:64, sum))
comm1 <- column_to_rownames(comm1, var = names(comm1)[1])
comm1 <- as.data.frame(t(comm1))

#result <- c()
for (i in 1:ncol(comm1)){
  test <- data.frame(RA=comm1[,i],treat)
  mean <- tapply(test$RA, test$Precipitation, mean)
  se <- tapply(test$RA, test$Precipitation, sd)/ sqrt(tapply(test$RA, test$Precipitation, length))
  change <- test %>% group_by(block, Year) %>%
    summarise(pct = ifelse(RA[Precipitation == "Control"] == 0,NA,
        (RA[Precipitation == "Wet"] - RA[Precipitation == "Control"]) / 
          RA[Precipitation == "Control"] * 100))
  change.mean <- mean(change$pct,na.rm=TRUE)
  change.se <- sd(change$pct,na.rm=TRUE)/sqrt(nrow(change))
  
  result1 <- summary(aov(RA ~ Precipitation*Year + block, data = test))
  output <- data.frame(Category=prefix,Taxa=colnames(comm1)[i],Year="Y09-Y16",
                       Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                       Change.mean=change.mean,Change.se=change.se,
                       Item="Year",Fvalue=result1[[1]]$`F value`[2],Pvalue=result1[[1]]$`Pr(>F)`[2])
  result <- rbind(result,output)
  output <- data.frame(Category=prefix,Taxa=colnames(comm1)[i],Year="Y09-Y16",
                       Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                       Change.mean=change.mean,Change.se=change.se,
                       Item="Treat",Fvalue=result1[[1]]$`F value`[1],Pvalue=result1[[1]]$`Pr(>F)`[1])
  result <- rbind(result,output)
  
  for (Year in unique(test$Year)){
    test1 <- test[test$Year==Year,]
    mean <- tapply(test1$RA, test1$Precipitation, mean)
    se <- tapply(test1$RA, test1$Precipitation, sd)/ sqrt(tapply(test1$RA, test1$Precipitation, length))
    change <- test1 %>% group_by(block) %>%
      summarise(pct = ifelse(RA[Precipitation == "Control"] == 0,NA,
                             (RA[Precipitation == "Wet"] - RA[Precipitation == "Control"]) / 
                               RA[Precipitation == "Control"] * 100))
    change.mean <- mean(change$pct,na.rm=TRUE)
    change.se <- sd(change$pct,na.rm=TRUE)/sqrt(nrow(change))
    
    result1 <- summary(aov(RA ~ Precipitation + block, data = test1))
    output <- data.frame(Category=prefix,Taxa=colnames(comm1)[i],Year=Year,
                         Control.mean=mean[1],Control.se=se[1],Wet.mean=mean[2],Wet.se=se[2],
                         Change.mean=change.mean,Change.se=change.se,
                         Item="Treat",Fvalue=result1[[1]]$`F value`[1],Pvalue=result1[[1]]$`Pr(>F)`[1])
    result <- rbind(result,output)
  }
  row.names(result) <- NULL
} 

write.csv(result, file = "3-ANOVA-Gene.csv")

### Time-decay relationships (TDRs) ------------------------------------------------------
### tdc.lmm: function used to calculate the TDRs
tdc.lmm<-function(betai,treat,prefixi=NULL)
{
  wet=treat$Precipitation
  wet.lev=unique(wet)
  out=list()
  for(j in 1:length(wet.lev))
  {
    idj=which(wet==wet.lev[j])
    sampj=rownames(treat)[idj]
    betaij=betai[idj,idj]
    betaij3=dist.3col(betaij)
    dtj=abs(treat$year[match(betaij3[,1],rownames(treat))]-treat$year[match(betaij3[,2],rownames(treat))])
    plotj1=treat$plot[match(betaij3[,1],rownames(treat))]
    plotj2=treat$plot[match(betaij3[,2],rownames(treat))]
    idj.use=which(plotj1==plotj2)
    cij.use=1-betaij3[idj.use,3]
    dtj.use=dtj[idj.use]
    plotj.use=plotj1[idj.use]
    logCij=log(cij.use)
    logCij[logCij==-Inf]<-NA
    logdtj=log(dtj.use)
    logdtj[logdtj==-Inf]<-NA
    lmij=lmer(logCij~logdtj+((1+logdtj)|plotj.use))
    lmijsm=summary(lmij)
    AIC1=AIC(lmij)
    r2ij=rsquared.glmm(lmij)
    lmijCS=Anova(lmij,type = "II")
    out[[j]]=c(slop.fix=lmijsm$coefficients[2,1],slop.sd=lmijsm$coefficients[2,2],
               R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,
               P.typeII=lmijCS[[3]],Chisq=lmijCS[[1]])
  }
  outs=Reduce(cbind,out)
  colnames(outs)=paste0(prefixi,"_",method,"_",wet.lev) 
  outs
}

### Overall community ###
#result <- c()
#output <- c()
betai=as.matrix(dist.bray) 
betai=as.matrix(dist.soren) 

prefixi = "GeoChip" #Prokaryotes,Fungi,GeoChip
method = "Bray" #Bray,Sorenson
tdci=tdc.lmm(betai=betai,treat = treat,prefixi = prefixi)
slope.obs = as.vector(tdci[1,])
r2.obs=as.vector(tdci[3:4,])
aic.obs=as.vector(tdci[5:6,])
ds.obs=(-tdci[1,1])-(-tdci[1,2])
tdci=as.data.frame(t(tdci))

### Permutation tests
year.lev=unique(treat$year)
rand = 1000
year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
trace.seq=seq(from=1,to=rand,by = 100)
ind.rand=lapply(1:nrow(year.perm),
                function(k)
                {
                  if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                  out=list()
                  idi=year.perm[k,]
                  perm.treat=treat
                  perm.treat[,"year"]=year.lev[idi[match(treat$year,year.lev)]]
                  tdcr=tdc.lmm(betai = betai, treat = perm.treat,prefixi = NULL)
                  out$slop.fix=as.vector(tdcr[1,])
                  out$r2=as.vector(tdcr[3:4,])
                  out$aic=as.vector(tdcr[5:6,])
                  out$ds=(-tdcr[1,1])-(-tdcr[1,2])
                  out
                })

slop.ran =sapply(1:length(ind.rand),function(k){ind.rand[[k]]$slop.fix})
r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
ds.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$ds})
EPS <- sqrt(.Machine$double.eps)
p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS),na.rm=TRUE)+1)/(ncol(r2.ran)+1)
p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS),na.rm=TRUE)+1)/(ncol(aic.ran)+1)
p.ds=if (ds.obs < 0){mean(ds.ran < ds.obs, na.rm = TRUE)}else{
  mean(ds.ran > ds.obs, na.rm = TRUE)}
result1=data.frame(tdci,P.R2M=p.r2[c(1,3)],P.R2C=p.r2[c(2,4)],
                   P.AIC1=p.aic[c(1,3)],P.AIC2=p.aic[c(2,4)],P.diff=p.ds)
result <- rbind(result,result1)

betaij3=dist.3col(betai)
betaij3$dtj=abs(treat$year[match(betaij3[,1],rownames(treat))]-treat$year[match(betaij3[,2],rownames(treat))])
betaij3$plotj1=treat$plot[match(betaij3[,1],rownames(treat))]
betaij3$plotj2=treat$plot[match(betaij3[,2],rownames(treat))]
idj.use=which(betaij3$plotj1==betaij3$plotj2)
betaij3 = betaij3[idj.use,]
betaij3$Treatment=treat$Precipitation[match(betaij3[,1],rownames(treat))]
betaij3$logCij = log(1-betaij3[,3])
betaij3[betaij3$logCij==-Inf]$logCij<-NA
betaij3$logdtj = log(betaij3[,4])
colnames(betaij3)[(ncol(betaij3)-1):ncol(betaij3)] <- paste0(prefixi,"_",method,"_",colnames(betaij3)[(ncol(betaij3)-1):ncol(betaij3)])
if (is.null(output)) {output <- betaij3} else {
  output <- cbind(output, betaij3[, (ncol(betaij3)-1):ncol(betaij3)])}

write.csv(result,"2-TDR-Significance.csv")
write.csv(output,"2-TDR-data.csv")

### Microbial taxa ###
comm <- comm[comm$Subcategory1=="Carbon degradation",c(1:64,(ncol(comm)-1))]
comm <- comm[comm$Gene_category=="Nitrogen",c(1:64,(ncol(comm)-1))]
comm <- comm[comm$Gene_category=="Phosphorus",c(1:64,(ncol(comm)-1))]
comm <- comm[comm$Gene_category=="Sulfur",c(1:64,(ncol(comm)-1))]

colnames(comm)[ncol(comm)] <- "Taxa"
comm <- na.omit(replace(comm, comm == "", NA))

#result <- c()
method = "Sorenson" #Bray,Sorenson
for (taxa in unique(comm$Taxa)){
  treat = read.csv(treat.file, row.names = 1, header = T, sep = ",")
  comm1 <- comm[comm$Taxa==taxa,1:64]
  comm1 = as.data.frame(t(comm1))
  comm1 = comm1[, colSums(comm1) > 0]
  comm1 = comm1[rowSums(comm1) > 0,]
  
  samp = match.name(rn.list = list(treat = treat, comm1 = comm1))
  comm1 = samp$comm1
  treat = samp$treat
  
  ### Bray–Curtis distance
  #dist.bray = vegdist(comm1)
  #betai=as.matrix(dist.bray) 
  
  ### Sorensen distance
  dist.soren = vegdist(comm1, binary = TRUE)
  betai=as.matrix(dist.soren) 
  
  tdci=tdc.lmm(betai=betai,treat = treat,prefixi = taxa)
  slope.obs = as.vector(tdci[1,])
  r2.obs=as.vector(tdci[3:4,])
  aic.obs=as.vector(tdci[5:6,])
  ds.obs=(-tdci[1,1])-(-tdci[1,2])
  tdci=as.data.frame(t(tdci))
  
  ### Permutation tests
  year.lev=unique(treat$year)
  rand = 1000
  year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
  trace.seq=seq(from=1,to=rand,by = 100)
  ind.rand=lapply(1:nrow(year.perm),
                  function(k)
                  {
                    if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                    out=list()
                    idi=year.perm[k,]
                    perm.treat=treat
                    perm.treat[,"year"]=year.lev[idi[match(treat$year,year.lev)]]
                    tdcr=tdc.lmm(betai = betai, treat = perm.treat,prefixi = NULL)
                    out$slop.fix=as.vector(tdcr[1,])
                    out$r2=as.vector(tdcr[3:4,])
                    out$aic=as.vector(tdcr[5:6,])
                    out$ds=(-tdcr[1,1])-(-tdcr[1,2])
                    out
                  })
  
  slop.ran =sapply(1:length(ind.rand),function(k){ind.rand[[k]]$slop.fix})
  r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
  aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
  ds.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$ds})
  EPS <- sqrt(.Machine$double.eps)
  p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS),na.rm=TRUE)+1)/(ncol(r2.ran)+1)
  p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS),na.rm=TRUE)+1)/(ncol(aic.ran)+1)
  p.ds=if (ds.obs < 0){mean(ds.ran < ds.obs, na.rm = TRUE)}else{
    mean(ds.ran > ds.obs, na.rm = TRUE)}
  result1=data.frame(tdci,P.R2M=p.r2[c(1,3)],P.R2C=p.r2[c(2,4)],
                     P.AIC1=p.aic[c(1,3)],P.AIC2=p.aic[c(2,4)],P.diff=p.ds)
  result <- rbind(result,result1)
}

write.csv(result,"2-TDR-Significance.taxa.csv")

### Mantel test ------------------------------------------------------
env <- treat[,c(1:6,11:12,14:18,21:23,25:26)]

### Geochip - Category ###
#result <- c()
comm1 <- comm[comm$Subcategory1=="Carbon degradation",c(1:64,(ncol(comm)-1))]
comm1 <- comm[comm$Gene_category=="Nitrogen",c(1:64,(ncol(comm)-1))]
comm1 <- comm[comm$Gene_category=="Phosphorus",c(1:64,(ncol(comm)-1))]
comm1 <- comm[comm$Gene_category=="Sulfur",c(1:64,(ncol(comm)-1))]

colnames(comm1)[ncol(comm1)] <- "Class"
comm1 <- na.omit(replace(comm1, comm1 == "", NA))
  
for (taxa in unique(comm1$Class)){
  test <- comm1[comm1$Class==taxa,1:64]
  test = test[rowSums(test) > 0, ]
  test = as.data.frame(t(test))
  
  samp = match.name(rn.list = list(env = env, test = test))
  test = samp$test
  env = samp$env
  
  if(ncol(test)>2){
    dist.test <- vegdist(test, method = "bray")
    for (i in 7:ncol(env)){
      dist.env <- dist(env[,i], method = "euclidean")
      result1 <-  mantel(dist.test,dist.env,method = "pearson",permutations = 999,na.rm = TRUE)
      result1 <- data.frame(Taxa = taxa,Env = colnames(env)[i],
                            r=result1$statistic,p=result1$signif)
      result <- rbind(result,result1)
    }
  }
}

write.csv(result,"3-Geochip-MantelCategory.csv")

### Geochip - Gene ###
#result <- c()
comm$Class2 <- paste0(comm$Class1,"_",comm$Class2)

comm1 <- comm[comm$Subcategory1=="Carbon degradation",c(1:64,ncol(comm))]
comm1 <- comm[comm$Subcategory1=="Carbon fixation",c(1:64,ncol(comm))]
comm1 <- comm[comm$Gene_category=="Nitrogen",c(1:64,ncol(comm))]
comm1 <- comm[comm$Gene_category=="Phosphorus",c(1:64,ncol(comm))]
comm1 <- comm[comm$Gene_category=="Sulfur",c(1:64,ncol(comm))]

colnames(comm1)[ncol(comm1)] <- "Class"
comm1 <- na.omit(replace(comm1, comm1 == "", NA))
comm1 <- comm1[comm1$Class!="_",]
prefix = "Sulfur"

for (taxa in unique(comm1$Class)){
  test <- comm1[comm1$Class==taxa,1:64]
  test = test[rowSums(test) > 0, ]
  test = as.data.frame(t(test))
  
  samp = match.name(rn.list = list(env = env, test = test))
  test = samp$test
  env = samp$env
  
  if(ncol(test)>1){
    dist.test <- vegdist(test, method = "bray")
    for (i in 7:ncol(env)){
      dist.env <- dist(env[,i], method = "euclidean")
      result1 <-  mantel(dist.test,dist.env,method = "pearson",permutations = 999,na.rm = TRUE)
      result1 <- data.frame(Category=prefix,Taxa = taxa,Env = colnames(env)[i],
                            r=result1$statistic,p=result1$signif)
      result <- rbind(result,result1)
    }
  }
}

write.csv(result,"3-Geochip-MantelGene.csv")

### Overall community ###
#result <- c()
env <- treat[,c(4,11:18,21:27)]
prefixi = "GeoChip" #Prokaryotes,Fungi,GeoChip

for (i in 1:ncol(env)){
  if (i==1){
    dist.env <- dist(env, method = "euclidean")
    result1 <-  mantel(dist.bray,dist.env,method = "pearson",permutations = 999,na.rm = TRUE)
    result1 <- data.frame(Taxa = prefixi,Env = "All",
                          r=result1$statistic,p=result1$signif)
    result <- rbind(result,result1)
  }

  dist.env <- dist(env[,i], method = "euclidean")
  result1 <-  mantel(dist.bray,dist.env,method = "pearson",permutations = 999,na.rm = TRUE)
  result1 <- data.frame(Taxa = prefixi,Env = colnames(env)[i],
                        r=result1$statistic,p=result1$signif)
  result <- rbind(result,result1)
}

write.csv(result,"4-Mantel.csv")

### Temporal trend ------------------------------------------------------
### Environment ###
#result <- c()
for (i in c(11:27,ncol(treat))){
  test <- data.frame(Env=treat[,i],treat[,1:6])
  test <- na.omit(test)
  test$year <- test$year-2009
  
  for (Treat in unique(test$Precipitation)){
    test1 <- test[test$Precipitation==Treat,]
    lmij=lmer(Env~year+((1+year)|block),data=test1)
    lmijsm=summary(lmij)
    AIC1=AIC(lmij)
    r2ij=rsquared.glmm(lmij)
    lmijCS=Anova(lmij,type = "II")
    
    output <- data.frame(Env=colnames(treat)[i],Treat=Treat,
                         slop.fix=lmijsm$coefficients[2,1],slop.sd=lmijsm$coefficients[2,2],t=lmijsm$coefficients[2,3],
                         R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,P.typeII=lmijCS[[3]])
    result <- rbind(result,output)
  }
  row.names(result) <- NULL
} 

result$P.adjust <- p.adjust(result$P.typeII, method = "fdr")
write.csv(result, file = "1-Temporal-Env.csv")

### Relative abundance of microbe ###
comm = comm[,c(1:64,ncol(comm))]
comm1 <- comm %>% group_by(Taxa) %>% summarise(across(1:64, sum))
comm1 <- column_to_rownames(comm1, var = names(comm1)[1])
comm1 <- comm1/colSums(comm1)
comm1 <- as.data.frame(t(comm1))

prefix = "Fungi" #Prokaryotes,Fungi
#result <- c()
for (i in 1:ncol(comm1)){
  test <- data.frame(RA=comm1[,i],treat[,1:6])
  test <- na.omit(test)
  test$year <- test$year-2009
  
  for (Treat in unique(test$Precipitation)){
    test1 <- test[test$Precipitation==Treat,]
    lmij=lmer(RA~year+((1+year)|block),data=test1)
    lmijsm=summary(lmij)
    AIC1=AIC(lmij)
    r2ij=rsquared.glmm(lmij)
    lmijCS=Anova(lmij,type = "II")
    
    output <- data.frame(Kingdom=prefix,Taxa=colnames(comm1)[i],Treat=Treat,
                         slop.fix=lmijsm$coefficients[2,1],slop.sd=lmijsm$coefficients[2,2],t=lmijsm$coefficients[2,3],
                         R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,P.typeII=lmijCS[[3]])
    result <- rbind(result,output)
  }
  row.names(result) <- NULL
} 

result$P.adjust <- p.adjust(result$P.typeII, method = "fdr")
write.csv(result, file = "2-Temporal-RelAbund.csv")

### Relative abundance of genes ###
comm$Class2 <- paste0(comm$Class1,"_",comm$Class2)
comm1 <- comm[comm$Subcategory1=="Carbon degradation",c(1:64,ncol(comm))]
comm1 <- comm[comm$Subcategory1=="Carbon fixation",c(1:64,ncol(comm))]
comm1 <- comm[comm$Gene_category=="Nitrogen",c(1:64,ncol(comm))]
comm1 <- comm[comm$Gene_category=="Phosphorus",c(1:64,ncol(comm))]
comm1 <- comm[comm$Gene_category=="Sulfur",c(1:64,ncol(comm))]

colnames(comm1)[ncol(comm1)] <- "Class"
comm1 <- na.omit(replace(comm1, comm1 == "", NA))
comm1 <- comm1[comm1$Class!="_",]
prefix = "Sulfur"

comm1 <- comm1 %>% group_by(Class) %>% summarise(across(1:64, sum))
comm1 <- column_to_rownames(comm1, var = names(comm1)[1])
comm1 <- as.data.frame(t(comm1))

#result <- c()
for (i in 1:ncol(comm1)){
  test <- data.frame(RA=comm1[,i],treat[,1:6])
  test <- na.omit(test)
  test$year <- test$year-2009
  
  for (Treat in unique(test$Precipitation)){
    test1 <- test[test$Precipitation==Treat,]
    lmij=lmer(RA~year+((1+year)|block),data=test1)
    lmijsm=summary(lmij)
    AIC1=AIC(lmij)
    r2ij=rsquared.glmm(lmij)
    lmijCS=Anova(lmij,type = "II")
    
    output <- data.frame(Class=prefix,Gene=colnames(comm1)[i],Treat=Treat,
                         slop.fix=lmijsm$coefficients[2,1],slop.sd=lmijsm$coefficients[2,2],t=lmijsm$coefficients[2,3],
                         R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,P.typeII=lmijCS[[3]])
    result <- rbind(result,output)
  }
  row.names(result) <- NULL
} 

result$P.adjust <- p.adjust(result$P.typeII, method = "fdr")
write.csv(result, file = "3-Temporal-Gene.csv")

### CCA - VPA ------------------------------------------------------
### 1 ### CCA
env <- treat[,c(4,11:15,18,21:23)]

prefixi = "Fungi" #Prokaryotes,Fungi,GeoChip
comm <- decostand(comm, "hellinger") 
CCA <- cca(comm ~ ., data = env)
CCA.summary <- summary(CCA)

### Variation explained
CCA$CCA$eig[1]/CCA$tot.chi*100
CCA$CCA$eig[2]/CCA$tot.chi*100

result <- anova(CCA)
result1 <- anova(CCA, by = "term") 

#output <- c()
output1 <- data.frame(Type=prefixi,Index=c(row.names(result)[1],row.names(result1)[1:10]),
                      Fvalue=c(result$F[1],result1$F[1:10]),Pvalue=c(result$`Pr(>F)`[1],result1$`Pr(>F)`[1:10]))
output <- rbind(output,output1)
write.csv(output, file = "4-CCA.csv")

env.site <- as.data.frame(CCA$CCA$biplot[,1:2])
row.names(env.site) <- c("Time","Moisture","pH","Temperature",
                         "NO3−","NH4+","ANPP","PR","ER","GPP")
micro.site <- data.frame(CCA1 = CCA.summary$sites[,1], 
                         CCA2 = CCA.summary$sites[,2],
                         Taxa = "Functional genes") #Prokaryotes,Fungi,Functional genes
micro.site$Treatment <- treat$Precipitation[match(rownames(micro.site),rownames(treat))]
micro.site$Year <- treat$Year[match(rownames(micro.site),rownames(treat))]

### 2 ### partial CCA-based VPA
Time <- env[,1]
Soil <- env[,c(2:6)]
Plant <- env[,c(7:10)]
result <- varpart(comm, Soil, Plant, Time)
plot(result)

### NMDS - PLS ------------------------------------------------------
### NMDS
prefixi = "GeoChip" #Prokaryotes,Fungi,GeoChip
set.seed(1028)
result <- metaMDS(comm = comm,distance = "bray")
result <- metaMDS(comm = comm,distance = "bray", binary = TRUE)

result$stress
Sites <- as.data.frame(scores(result, display = "sites")) 
colnames(Sites) <- paste0(prefixi,"_bray_",colnames(Sites))

#output <- c()
output <- cbind(output,Sites)
output <- data.frame(output,treat[,c(1:5)])

write.csv(output, file = "0-NMDS.csv")

### PLS
############
#install.packages("OmicsPLS")
library(OmicsPLS)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ropls")
library(ropls)
#################
#################
R2ff<-function(xm,ym,o2p)
{
  c(R2x=mean(sapply(1:ncol(xm),function(i){1-sum((o2p$X_hat[,i]-xm[,i])^2)/sum((xm[,i]-mean(xm[,i]))^2)})),
    R2y=mean(sapply(1:ncol(ym),function(i){1-sum((o2p$Y_hat[,i]-ym[,i])^2)/sum((ym[,i]-mean(ym[,i]))^2)})))
}
R2sep<-function(xm,pls){(((getVipVn(pls))^2)/ncol(xm))*getSummaryDF(pls)[,'R2Y(cum)']}

plstest<-function(xm,ym,rand=100)
{
  nx=ncol(xm)
  combs=list()
  for(i in 1:nx)
  {
    message("i=",i," ",date())
    cbni=combn(nx,i)
    combs=c(combs,lapply(1:ncol(cbni),function(i){cbni[,i]}))
  }
  message("Total of ",length(combs)," combinations. ",date())
  #trac=seq(from=1,to=length(combs),by=50)
  id=t(rep(0,nx))
  colnames(id)=colnames(xm)
  
  dfna=t(rep(NA,8))
  colnames(dfna)=c('R2X(cum)','R2Y(cum)','Q2(cum)','RMSEE','pre','ort','pR2Y','pQ2')
  
  res=lapply(1:length(combs),
             function(i)
             {
               message("----- PLS i=",i," ",date())
               xmi=xm[,combs[[i]],drop=FALSE]
               pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
               if(class(pls)=="try-error")
               {
                 pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
               }
               out=dfna
               if(class(pls)!="try-error"){dfi=getSummaryDF(pls);out[1,match(colnames(dfi),colnames(dfna))]=as.vector(as.matrix(dfi[1,]))}
               idni=id
               idni[combs[[i]]]=1
               cbind(idni,out)
             })
  res
}

plsfw<-function(xm,ym,r2buf=0.98,Q2ck=FALSE,SEEck=FALSE,rand=100)
{
  nx=ncol(xm)
  fs<-list(integer(0))
  dfna=t(rep(NA,8))
  colnames(dfna)=c('R2X(cum)','R2Y(cum)','Q2(cum)','RMSEE','pre','ort','pR2Y','pQ2')
  id=t(rep(0,nx))
  colnames(id)=colnames(xm)
  R2Yn=list()
  fijn=list()
  outrc=list()
  if(Q2ck){Q2n=list()}
  if(SEEck){SEEn=list()}
  k=1
  kn=1
  for(fn in 1:nx)
  {
    for(i in 1:length(fs))
    {
      fsi=fs[[i]]
      fai=which(!((1:nx) %in% fsi))
      for(j in 1:length(fai))
      {
        fij=c(fsi,fai[j])
        
        message("----- PLS fn=",fn," in ",nx,", i=",i," in ",length(fs),", j=",j," in ",length(fai),". ",date())
        xmi=xm[,fij,drop=FALSE]
        pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
        if(class(pls)=="try-error")
        {
          pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
        }
        out=dfna
        if(class(pls)!="try-error"){dfi=getSummaryDF(pls);out[1,match(colnames(dfi),colnames(dfna))]=as.vector(as.matrix(dfi[1,]))}
        idni=id
        idni[fij]=1
        outrc[[k]]=cbind(idni,out)
        k=k+1
        R2Yn[[kn]]=out[,"R2Y(cum)"][[1]]
        if(Q2ck){Q2n[[kn]]=out[,"Q2(cum)"][[1]]}
        if(SEEck){SEEn[[kn]]=out[,"RMSEE"][[1]]}
        fijn[[kn]]=fij
        kn=kn+1
      }
    }
    maxR2n=max(unlist(R2Yn),na.rm = TRUE)
    kns=which(unlist(R2Yn)>=(maxR2n*r2buf))
    if(Q2ck)
    {
      if(length(kns)>1)
      {
        Q2nk=Q2n[kns]
        maxQ2nk=max(unlist(Q2nk),na.rm = TRUE)
        if(maxQ2nk>=0){maxQ2nkb=maxQ2nk*r2buf}else{maxQ2nkb=maxQ2nk*(1+(1-r2buf))}
        knsk1=which(unlist(Q2nk)>=maxQ2nkb)
        kns=kns[knsk1]
      }
    }
    if(SEEck)
    {
      if(length(kns)>1)
      {
        SEEnk=SEEn[kns]
        minSEE=min(unlist(SEEnk),na.rm = TRUE)
        knsk2=which(unlist(SEEnk)<=(minSEE*(1+(1-r2buf))))
        kns=kns[knsk2]
      }
    }
    EPS <- (.Machine$double.eps)
    if((max(kns)<=length(fs)) | sum(is.na(kns))>0 | maxR2n >= (1-EPS)){break}else{
      fs=fijn[kns[which(kns>length(fs))]]
      R2Yn=R2Yn[kns]
      fijn=fijn[kns]
      kn=length(R2Yn)+1
    }
  }
  outrcm=Reduce(rbind,outrc)
}

rp.pls<-function(xm,ym,rand=100)
{
  R2sepf<-function(xm,ym,predI)
  {
    pls=try(opls(x=xm,y=ym,predI=predI,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))
    if(class(pls)=="try-error"){pls=try(opls(x=xm,y=ym,predI=1,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))}
    if(class(pls)=="try-error"){out1=rep(NA,1+ncol(xm))}else{
      out1=c(R2Y=getSummaryDF(pls)[,"R2Y(cum)"],(((getVipVn(pls))^2)/ncol(xm))*getSummaryDF(pls)[,'R2Y(cum)'])
    }
    out1
  }
  pls=try(opls(x=xm,y=ym,predI=NA,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))
  predI=NA
  if(class(pls)=="try-error"){predI=1;pls=try(opls(x=xm,y=ym,predI=1,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))}
  if(class(pls)=="try-error"){out=rep(NA,2*(1+ncol(xm)))}else{
    R2obs=R2sepf(xm,ym,predI)
    perm=permute::shuffleSet(nrow(xm),nset = rand)
    tracs=seq(from=1,to=nrow(perm),by=20)
    R2rm=sapply(1:nrow(perm),
                function(i)
                {
                  if(i %in% tracs){message("i=",i,". ",date())}
                  ymri=ym[perm[i,],,drop=FALSE]
                  rownames(ymri)=rownames(ym)
                  R2sepf(xm,ymri,predI)
                })
    EPS <- (.Machine$double.eps)
    dR2=((R2rm-R2obs)>=(-EPS))
    out=c(R2obs,rowSums(dR2,na.rm = TRUE)/rand)
  }
  names(out)=c(paste0("R2.",c('Y',colnames(xm))),paste0("P.",c('R2Y',colnames(xm))))
  out
}

r2adj<-function(r2,n,p)
{
  idx=which((n-p-1)<0)
  out=1-((1-r2)*((n-1)/(n-p-1)))
  out[idx]=NA
  out
}

### data: To account for potential temporal autocorrelation, plot-level means were calculated by averaging microbial and environmental measurements across time points within each plot
data$NO3.N <- log10(data$NO3.N)
data[,c(4:ncol(data))] <- scale(data[,c(4:ncol(data))])

be=data[,'Heterotrophic',drop=FALSE]
ym=as.matrix(be)
Yname=paste(colnames(ym))

colnames(data)
#"Treat","annual_moisture","pH","NO3.N","plant.richness","Prokaryotes_NMDS1","Fungi_NMDS1","GeoChip_NMDS1","Heterotrophic"
data1=as.matrix(data[,c("Treat","annual_moisture","pH","NO3.N","plant.richness","Prokaryotes_NMDS1","Fungi_NMDS1","GeoChip_NMDS1"), drop=FALSE])
head(data1)
plstm=plsfw(data1,ym,r2buf=0.98,Q2ck=TRUE,SEEck=FALSE,rand=100) 
ieggr::save.file(plstm,prefix="PLS",filename = paste0(Yname,".PLSFWtest"))

# Optimuize result: R2Y larger than 98% of maximum R2Y, and P values for R2Y and Q2 should be <0.05; the least factor number; if multiple hits, use minimum RMSEE; if still multiple hits, choose the one with better biological sense.
sel.id=c(1) #annual_moisture
sel.id=c(2,3) #pH
sel.id=c(3) #NO3.N
sel.id=c(1,2,5,6) #Plant richness
sel.id=c(2,3,5) #Prokaryotes_NMDS1
sel.id=c(3,5,8) #Fungi_NMDS1
sel.id=c(2,4) #GeoChip_NMDS1
sel.id=c(7,8) #Heterotrophic
data2=data1[,sel.id,drop=FALSE]
head(data2)

pls=try(opls(x=data2,y=ym,predI=NA,orthoI=0,permI=1000,fig.pdfC="none"))
if(class(pls)=="try-error"){pls=try(opls(x=data2,y=ym,predI=1,orthoI=0,permI=1000,fig.pdfC="none"))}

vip=getVipVn(pls)
names(vip)=paste0("VIP.",names(vip))
rpi=R2sep(data2,pls)
names(rpi)=paste0("R2.",names(rpi))

Rsig<-Psig<-rep(NA,ncol(data2))
names(Rsig)<-paste0("r.Single.",colnames(data2))
names(Psig)<-paste0("P.Single.",colnames(data2))
for(j in 1:ncol(data2)){
  result <- cor.test(data2[,j],ym,method = "pearson")
  r <- result$estimate
  p <- result$p.value
  Rsig[j]=r
  Psig[j]=p
}

### Stochasticity (MST) ------------------------------------------------------
prefixi = "Prokaryotes" #Prokaryotes,Fungi,GeoChip

#result <- c()
for (year in unique(treat$Year)){
  group <- treat[treat$Year == year, 5, drop = FALSE]
  colnames(group) <- "treatment"
  comm1 <- comm[treat$Year==year,]
  comm1 = comm1[, colSums(comm1) > 0]

  set.seed(1028)
  tnst <-tNST(comm=comm1,group=group,dist.method="bray",
              abundance.weighted=TRUE)
  
  result1 <- data.frame(Taxa=prefixi,Year=year,tnst$index.pair.grp) 
  result <- rbind(result,result1)
}

write.csv(result, file = "5-MST.csv")
