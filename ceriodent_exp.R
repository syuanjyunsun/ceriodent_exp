####Statistical analyses for cerio dent experiment####

data=read.csv("ceriodent_exp.csv",header=T)
library(lme4)
library(car)
library(emmeans)

#parasite attack and physical barrier
#1. parasite attack
model=glmer(totalspore~currenthost*metshsource+(1|block),poisson,data=data)
Anova(model,type=3)

#post-hoc comparisons between different spore-sourced hosts while controlling for exposure hosts
a=emmeans (model,  ~  metshsource|currenthost , adjust="tukey")
pairs(a)

#post-hoc comparisons between different exposure hosts while controlling for spore-sourced hosts
b=emmeans (model,  ~  currenthost|metshsource , adjust="tukey")
pairs(b)

#2. successfully penetrating spores
model=glmer(penet~currenthost*metshsource+(1|block),poisson,data=data[data$totalspore!="0",])
Anova(model,type=3)

#post-hoc comparisons between different spore-sourced hosts while controlling for exposure hosts
a=emmeans (model,  ~  metshsource|currenthost , adjust="tukey")
pairs(a)

#post-hoc comparisons between different exposure hosts while controlling for spore-sourced hosts
b=emmeans (model,  ~  currenthost|metshsource , adjust="tukey")
pairs(b)

#3. proportion of penetrating spores out of attacking spores
model=glmer(penet~currenthost*metshsource+offset(log(totalspore+1))+(1|block),poisson,data=data[data$totalspore!=0,])
Anova(model,type=3)

#post-hoc comparisons between different spore-sourced hosts while controlling for exposure hosts
a=emmeans (model,  ~  metshsource|currenthost , adjust="tukey")
pairs(a)

#post-hoc comparisons between different exposure hosts while controlling for spore-sourced hosts
b=emmeans (model,  ~  currenthost|metshsource , adjust="tukey")
pairs(b)

#Host immune responses
#4. Total hemocyte number
model=glmer(hemocy~currenthost+metshsource+(1|block),poisson,data=data[data$totalspore!="0",])
Anova(model,type=3)

#5. Hemocytes per spore
model=glmer(log(hemocyperspore+1)~currenthost+metshsource+(1|block),gaussian,data=data[data$totalspore!="0",])
Anova(model,type=3)

####wilcox.test####
datac=data[data$currenthost=="cerio",]
datamd=data[data$metshsource=="dent",]

wilcox.test(hemocyperspore ~ metshsource, data=datac) 
wilcox.test(hemocyperspore ~ currenthost, data=datamd) 
wilcox.test(hemocy ~ metshsource, data=datac) 
wilcox.test(hemocy ~ currenthost, data=datamd) 

#Terminal infection
library(logistf)
####Implementing Firth's bias-Reduced penalized-likelihood logistic regression.
model=logistf(infection~currenthost+metshsource,firth = TRUE,pl=TRUE,data=data)
summary(model)

#focusing on D. dentifera-sourced parasites in dentifera, to look for the effects of immune responses on terminal infection
datadent<-data[data$currenthost=="dent",]
datadentdent<-datadent[datadent$metshsource=="dent",]

model=glmer(infection~hemocy+(1|block),binomial,data=datadentdent)
Anova(model,type=3)

model=glmer(infection~hemocyperspore+(1|block),data=datadentdent)
Anova(model,type=3)


