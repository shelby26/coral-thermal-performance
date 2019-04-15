library(dplyr)
library(nlme)
library(multcomp)

data<-read.table(file.choose(),header=T,sep="\t")

## code from: http://rcompanion.org/rcompanion/d_07.html
model8<-lme(avg~species,random=~1|Month,data=Allmeans,method="REML")
anova.lme(model8,type="sequential",adjustSigma=F)
posthoc=glht(model8,linfct=mcp(species="Tukey"))
mcs=summary(posthoc,test=adjusted("single-step"))
mcs
cld(mcs,level=0.05,decreasing=T)
hist(residuals(model8),col="darkgray")
plot(fitted(model8),residuals(model8))


 
 