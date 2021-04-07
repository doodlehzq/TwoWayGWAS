###################################################### 
#                                                    #
#  test genome-genome interaction with linear model  #
#                                                    #
######################################################

# pheno ~ geno_rice + geno_xoo + geno_rice:geno_xoo + xoo covariates + rice covariates

library(parallel)
library(reshape2)

#setwd("path_to_demo_data")



# read Xoo subset genotype data
genoXoo<-read.csv("xoo_geno.csv",header=T)
#Xoo SNP coordinates
XP<-as.numeric(sapply(colnames(genoXoo)[-1],function(x) substr(x,2,nchar(x))))

# read phenotype data
pheno<-read.csv("rice_xoo_pheno.csv",header=T)
# read rice subset genotype data
genoRice<-read.csv("rice_geno.csv",header=T)
#'0' represents missing values in genotype data
genoRice[genoRice==0]<-NA

#Xoo number
nXoo<-nrow(genoXoo)
#rice number
nRice<-ncol(genoRice)-2
#Xoo variant number
nXooVar<-ncol(genoXoo)-1
#Rice variant number
nRiceVar<-nrow(genoRice)


#xoo covariate vectors
xoo_c<-matrix(rep(0,nRice*nXoo*(nXoo-1)),ncol=nXoo-1)
for(i in 1:(nXoo-1)){
  xoo_c[(nRice*(i-1)+1):(nRice*i),i] <- 1
}

#rice covariate vectors
rice_c<-matrix(rep(0,nRice*nXoo*(nRice-1)),ncol=nRice-1)
for(i in 1:(nRice-1)){
  rice_c[nRice*c(0:(nXoo-1))+i,i]<-1
}

# prepare phenotype values
value=c()
for(i in 2:length(pheno)){
  value<-c(value,pheno[,i])
}

#transfer lession length to log scale
#value=log2(value)

# function to conduct interaction tests for the qth xoo variant (to all rice genotypes)
interA<-function(q){
  k1=genoRice
  k1$XP<-rep(XP[q],nrow(k1))
  k1$P=rep("NA",length(k1$BP))
  k=genoXoo[,q]
  xt=c()
  for(i in 1:length(k)){
    xt=c(xt,rep(k[i],nRice))
  }
  hh<-function(i){
    n=k1[i,3:(2+nRice)]
    rt=as.numeric(as.factor(rep(t(n),nXoo)))-1
    xt=as.numeric(as.factor(xt))-1
    inter= as.numeric(rt==xt)
    m = cbind(rt,xt,inter,xoo_c,rice_c)

    t=summary(lm(value~m))
    mm=t[[4]][,4]
    return(mm['minter'])
  }
  
  k1$P=sapply(1:length(k1$BP),hh)
  k1[,c(1,2,ncol(k1)-1,ncol(k1))]
}

#using 1 cpu core
r <- lapply(2:ncol(genoXoo),interA)

#using multiple cpu cores

#no_cores <- detectCores() - 1
#cl <- makeCluster(no_cores)
#clusterExport(cl=cl, varlist=c("genoXoo","pheno","genoRice","interA","value","xoo_c","rice_c","XP","nRice","nXoo"))
#r <- parLapply(cl,2:ncol(genoXoo),interA)

# merge results
res <- r[[1]]
for(i in 2: length(r)){
  res=rbind(res,r[[i]])
}

res$p.adj<-p.adjust(res$P,method="bonferroni")

#res
#CHR       BP    XP            P      p.adj
#1   1  2651540 56506 1.665282e-03 1.00000000
#2   1 10323754 56506 6.664005e-05 0.09996008
#3   1 32274685 56506 7.779066e-01 1.00000000
# ... ... ...
# ... ... ...



