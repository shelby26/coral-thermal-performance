##Lipid Plots
##Use dataset "lipid"

##Symbiont lipids

lipidgroup<-group_by(lipid,Species)
sum_lipidgroup<-summarise(lipidgroup,StDev=sd(S_lipid),SLMean=mean(S_lipid))

symlipidaov<-aov(S_lipid~Species+Error(Rep),data=lipid)
summary(symlipidaov)
#TukeyHSD(symlipidaov)
pairwise.t.test(lipid$S_lipid,lipid$Species,p.adjust="bonferroni")

model_Slip<-lme(S_lipid~Species,random=~1|Month,data=lipid,method="REML")
anova.lme(model_Slip,type="sequential",adjustSigma=F)
posthoc=glht(model_Slip,linfct=mcp(Species="Tukey"))
mcs=summary(posthoc,test=adjusted("single-step"))
mcs
cld(mcs,level=0.05,decreasing=T)
hist(residuals(model_Slip),col="darkgray")
plot(fitted(model_Slip),residuals(model_Slip))

model_Hlip<-lme(H_lipid~Species,random=~1|Month,data=lipid,method="REML")
anova.lme(model_Hlip,type="sequential",adjustSigma=F)
posthoc=glht(model_Hlip,linfct=mcp(Species="Tukey"))
mcs=summary(posthoc,test=adjusted("single-step"))
mcs
cld(mcs,level=0.05,decreasing=T)
hist(residuals(model_Hlip),col="darkgray")
plot(fitted(model_Hlip),residuals(model_Hlip))


ggplot(lipid, aes(x=Species, y=S_lipid))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
  axis.ticks.length=unit(-0.20, "cm"), axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+ 
  coord_cartesian(ylim = c(0, 80))+
  ylab("Symbiont lipid content")+
  annotate("text",x=1,y=75,label="A")+
  annotate("text",x=2,y=70,label="A")+
  annotate("text",x=3,y=70,label="A,B")+
  annotate("text",x=4,y=75,label="B")+
  annotate("text",x=5,y=75,label="B")

