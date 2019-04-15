##Function for running any species and any response variable through TPC with 5 models (including Weibull)

##First load packages
{
library(minpack.lm)
library(pracma)
library(ggplot2)
library(devtools)
install_github('evpickett/evp')
library(evp)
library(MuMIn)
library(dplyr)
}

coral_data_pr<-read.table(file.choose(),header=T,sep="\t")  #Data file: PR PAM STATS.txt "Productivity" and FvFm "Photosynthetic Capacity"


TPC_modelfit<-function(datafile,SPECIES,var1,axismax)  ##SPECIES and var1 need to be in "", e.g. "ACRO"
{ 
  data1<-datafile
  data2<-subset(data1,data1$species==SPECIES)
  data2<-droplevels(data2)
  response<-(data2[,var1])
  temp<-data2$temp
  
gauss <- nlsLM(response ~ p1 * exp(-0.5*((abs(temp - p2)) / p3)^ 2),
                   start = list(p1 = 0, p2 = 1, p3 = 3), trace = T, control = nls.lm.control(maxiter = 100))

emg<-nlsLM(start=list(p1=0,p2=5,p3=10,p4=1),response~p1*p3*(sqrt(2*pi))*(1/2*p4)*exp((p3^2)/(2*(p4^2))+p2-temp/p4)*(p4/abs(p4)-erf((p2-temp)/sqrt(2*p3)+p3/sqrt(2*p4))),trace=T,control=nls.lm.control(maxiter=100))

weibull<-nlsLM(response~p1*((p3-1)/p3)^((1-p3)/p3)*(abs((temp-p4)/p2+((p3-1)/p3)^(1/p3))^(p3-1))*exp(-abs((temp-p4)/p2+((p3-1)/p3)^(1/p3))^p3+(p3-1)/p3),start=list(p1=1,p2=15,p3=1,p4=1),trace=T,control=nls.lm.control(maxiter=100))

lognor<-nlsLM(response~p1*exp(-0.5*(log((temp)/p3,base=exp(1))/p2)^2),start=list(p1=1,p2=1,p3=1),trace=T,control=nls.lm.control(maxiter=100))

linear<-nlsLM(start=list(p1=0,p2=0),response~p1*temp+p2,control=nls.lm.control(maxiter=100))

##SCORE MODELS

a<-(AIC(gauss)) 
b<-(AIC(emg))
c<-(AIC(weibull))
d<-(AIC(lognor))
e<-(AIC(linear))

AICvals<-c(a,b,c,d,e)
AICmin<-which.min(AICvals)
AICWeights<-Weights(AIC(gauss, emg, weibull, lognor, linear))
cat
#cat(AICmin)
cat("### AIC VALUES ###","\n","Gaussian=",a,"  (",AICWeights[1],")","\n","Exp mod Gaussian=",b,"  (",AICWeights[2],")","\n","Weibull=",c,"  (",AICWeights[3],")","\n","Log Normal=",d,"  (",AICWeights[4],")","\n","Linear=",e,"  (",AICWeights[5],")","\n")

dataplot=data.frame(temp,response)

if(AICmin=="1") {
  

    pred_gauss <- data.frame(pred = predict(gauss, newdata = data.frame(temp = c(-1:45))),
                            temp = c(-1:45))
    
    ggplot(pred_gauss, aes(y = pred, x = temp)) +
        geom_line() +
        theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        coord_cartesian(ylim = c(0, axismax),xlim=c(5,35))+
       # geom_point(data=dataplot,aes(temp,response))+
        ylab(var1)+
        ggtitle(bquote(atop(.(SPECIES))))
        #+geom_vline(aes(xintercept=),linetype=2,colour="black")

  }

else if(AICmin=="2") {

    pred_emg <- data.frame(pred = predict(emg, newdata = data.frame(temp = c(-1:45))),
                            temp = c(-1:45))
  
    ggplot(pred_emgpor, aes(y = pred, x = temp)) +
        geom_line() +
        theme_bw()+
        coord_cartesian(ylim = c(0, axismax))+
        geom_point(data=dataplot,aes(temp,response))+
        ylab(var1)+
      ggtitle(bquote(atop(.(SPECIES),atop(italic("Exp mod gaussian"),""))))
  }

else if(AICmin=="3") {
  
    pred_weibull <- data.frame(pred = predict(weibull, newdata = data.frame(temp = c(-1:45))),
                              temp = c(-1:45))

    ggplot(pred_weibull, aes(y = pred, x = temp)) +
        geom_line() +
        theme_bw()+
        coord_cartesian(ylim = c(0, axismax))+
        geom_point(data=dataplot,aes(temp,response))+
        ylab(var1)+
        ggtitle(bquote(atop(.(SPECIES),atop(italic("weibull"),""))))
  }

else if(AICmin=="4") {

    pred_lognor <- data.frame(pred = predict(lognor, newdata = data.frame(temp = c(-1:45))),
                             temp = c(-1:45))

    ggplot(pred_lognor, aes(y = pred, x = temp)) +
        geom_line() +
        theme_bw()+
        coord_cartesian(ylim = c(0, axismax))+
        geom_point(data=dataplot,aes(temp,response))+
        ylab(var1)+
        ggtitle(bquote(atop(.(SPECIES),atop(italic("log normal"),""))))
  }

else if(AICmin=="5") {

    pred_linear <- data.frame(pred = predict(linear, newdata = data.frame(temp = c(-1:45))),
                             temp = c(-1:45))

    ggplot(pred_linear, aes(y = pred, x = temp)) +
          geom_line() +
          theme_bw()+
          coord_cartesian(ylim = c(0, axismax))+
          geom_point(data=dataplot,aes(temp,response))+
          ylab(var1)+
          ggtitle(bquote(atop(.(SPECIES),atop(italic("linear"),""))))
  }

else(print("not working"))
}

  
  