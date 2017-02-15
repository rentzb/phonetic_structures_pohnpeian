### Code by: Bradley Rentz
### R Script for Phonetic Structures of Pohnpeian, Bradley Rentz and Victoria Anderson, 2017
### Mac OS 10.10.5
### Base R version 3.2.2 "Fire Safety"
### lmerTest version 2.0-30 (lme4 version 1.1-12)
### dplyr version 0.4.3
### ggplot2 version 2.1.0


library(lmerTest)
library(ggplot2)
library(dplyr)
library(coefplot2)
library(tidyr)
library(phonR)
library(sjPlot)
library(sjmisc)
library(outliers)
library(multcomp)
library(LMERConvenienceFunctions)
library(scatterplot3d)
library(ggthemes)
library(rstanarm)
library(plotly)
library(yarrr)


options (mc.cores=parallel::detectCores ()) # Run on multiple cores

set.seed (3875)

#### outline of script ######
# 1. Geminates
# 2. Laminal-Apical Locus Equations
# 3. pw-p locus equations
# 4. laminal-apical VOT
# 5. pw-p VOT
# 6. Vowel allophony
# 7. Frication

## individual section outline
# a. import data
# b. subsetting
# c. initial plotting
# d. outlier removal
# e. optional data coding (ie., contrast coding)
# f. inferential statistics (LMERs)
# g. post hoc tests if applicable
# h. display plotting of data (as needed)

setwd("~/Documents/UH/QP/qp1_data/results") # change this to your wd
################################################ Begin Geminate Section #################################
#load geminate data (303 obs)
geminates <- read.csv("geminates.csv")

#subsets of geminate data by speaker
geminates.m01 <- geminates[geminates$speaker=="M01",]
geminates.m02 <- geminates[geminates$speaker=="M02",]
geminates.m03 <- geminates[geminates$speaker=="M03",]
geminates.f01 <- geminates[geminates$speaker=="F01",]
geminates.f02 <- geminates[geminates$speaker=="F02",]

# initial box plots
ggplot(geminates, aes(x=consonant,y=length,fill=c_length)) + geom_boxplot()+ facet_grid(speaker~word_place,scales="free")
#plot boxplots by speaker
ggplot(geminates.m01, aes(x=consonant, y=length, fill=c_length)) + geom_boxplot() + facet_grid(~word_place)
ggplot(geminates.m02, aes(x=consonant, y=length, fill=c_length)) + geom_boxplot() + facet_grid(~word_place)
ggplot(geminates.m03, aes(x=consonant, y=length, fill=c_length)) + geom_boxplot() + facet_grid(~word_place)
ggplot(geminates.f01, aes(x=consonant, y=length, fill=c_length)) + geom_boxplot() + facet_grid(~word_place)
ggplot(geminates.f02, aes(x=consonant, y=length, fill=c_length)) + geom_boxplot() + facet_grid(~word_place)

#remove 3 sd or greater for length grouped by speaker, consonant, and c_length (0 removed)
geminates.3sd <- geminates %>%
  group_by(speaker,consonant,word_place,c_length) %>%
  mutate(plusthreesdmean = mean(length)+3*sd(length),minusthreesd=mean(length)-3*sd(length)) 
geminates.clean <-geminates.3sd %>%
  filter(length < plusthreesdmean) %>%
  filter(length>minusthreesd)

# values over 0.35s are way too long so must be measurement error so remove
geminates.clean <- geminates.clean %>%
  filter(length< 0.35) # removes 3 out of 303 obs so 0.99%


gem3d <- with(geminates.clean, scatterplot3d(c_length, consonant, length, type="h" ))
fit <- with(geminates.clean,lm(length ~c_length+consonant) )
summary(fit)
gem3d$plane3d(fit)
# set up contrast coding


# change to contrast coding
contrasts(geminates.clean$c_length) <- c(-.5, .5)

relLik(geminates.lmer,geminates.lmer1)
mcp.fnc(geminates.lmer)

#lmers for geminates

# full model (doesn't run)
geminates.lmer = lmer(length ~ consonant * c_length * word_place * frame + (1+consonant*c_length * word_place*frame|speaker) + (1+consonant*c_length*word_place*frame|word),data=geminates.clean,REML=F)
# reduce model until converges by taking out correlation of random slopes, then random slopes, then fixed effects interactions as needed

# Full Bayesian model

geminates.blmer = stan_lmer(length ~ consonant * c_length * word_place * frame + (1+consonant*c_length * word_place*frame|speaker) + (1+consonant*c_length*word_place*frame|word),data=geminates.clean, 
                            prior_intercept = normal(0, 10),
                             prior = normal(0, 1),
                             prior_covariance = decov(regularization = 2),
                             
                             chains = 4,
                             iter = 2000)
#


geminates.lmer = lmer(length*1000 ~ consonant * c_length * word_place * frame + (1+consonant   |speaker) + (1+consonant|word),data=geminates.clean,REML=F)
summary(geminates.lmer)
sjt.lmer(geminates.lmer)
sjp.lmer(geminates.lmer)
sjp.lmer(geminates.lmer,type="pred",vars=c("c_length","word_place"),facet.grid=F)

# take out non-sig interaction of fixed effects
geminates.lmer1 = lmer(length ~ c_length + consonant + frame + word_place + consonant:word_place + consonant:frame  + word_place:frame + (1+consonant|speaker) + (1+consonant|word),data=geminates.clean,REML=F)
summary(geminates.lmer1) 

# not sig different but close (p=0.056) so go with lmer1
anova(geminates.lmer,geminates.lmer1)

# take out fixed effect frame
geminates.lmer2 = lmer(length ~ c_length + consonant +word_place + consonant:word_place + consonant:frame + word_place:frame + (1+consonant|speaker) + (1+consonant|word),data=geminates.clean,REML=F)
summary(geminates.lmer2)
anova(geminates.lmer1,geminates.lmer2) # sig different so keep lmer1

# take out consonant:frame
geminates.lmer3 = lmer(length ~ c_length + consonant + word_place + consonant:word_place + word_place:frame + (1+consonant|speaker) + (1+consonant|word),data=geminates.clean,REML=F)
summary(geminates.lmer3)
anova(geminates.lmer1,geminates.lmer3) # not sig so go with lmer3

# check random slopes for significant
geminates.lmer4 = lmer(length ~ c_length + consonant + word_place + consonant:word_place + word_place:frame + (1+consonant|speaker) + (1|word),data=geminates.clean,REML=F)
anova(geminates.lmer3,geminates.lmer4) # not sig so go with lmer4


geminates.lmer5 = lmer(length ~ c_length + consonant  + word_place +consonant:word_place + word_place:frame + (1|speaker) + (1|word), data=geminates.clean, REML=FALSE)
summary(geminates.lmer5)
anova(geminates.lmer4,geminates.lmer5) # not sig diff so go with lmer5

geminates.blmer5 = stan_lmer(length ~ c_length * consonant  + word_place +consonant:word_place + word_place:frame + frame+ (1+consonant+c_length|speaker) + (1+consonant+c_length|word), data=geminates.clean,
                             prior_intercept = normal(0, 1),
                              prior = normal(0, 2),
                            prior_covariance = decov(regularization = 2),
                            chains = 4,
                            iter = 2000)

plot(geminates.blmer5,pars=c("(Intercept)","c_length1","consonantm","consonantmw","consonantn","consonantr","word_placeinitial","word_placemedial","consonantm:word_placeinitial","consonantmw:word_placemedial"))

sink("sink-summary-blmer-geminate5.txt")
print(summary(geminates.blmer5),digits=4)
sink()
launch_shinystan(geminates.blmer5) # view diagnostic plots and others
#check residuals
qqnorm(residuals(geminates.lmer5)) # good, normal
plot(fitted(geminates.lmer5),residuals(geminates.lmer5)) # ok. mostly homoscedastic?

# post hoc comparisons 
summary(glht(geminates.lmer5, linfct=mcp(c_length = "Tukey")), test=adjusted(type="BH")) # c_length
summary(glht(geminates.lmer5,linfct=mcp(word_place = "Tukey")),test=adjusted(type="BH")) # word_place
# interesting word_place comparisons. what does that mean? also what about the covariate interaction warning? what should the contrast be




summary(glht(geminates.lmer5, linfct=mcp(c_length = "Tukey")), test=adjusted(type="BH"))
#data subset by place in word
geminates.clean.initial <- geminates.clean[geminates.clean$word_place=="initial",]
geminates.clean.medial <- geminates.clean[geminates.clean$word_place=="medial",]
geminates.clean.final <- geminates.clean[geminates.clean$word_place=="final",]

#lmer by place in word subgroup
#the full thing converges!
geminates.clean.initial.lmer <- lmer(length ~ consonant * c_length * frame + (1+consonant*c_length+frame|speaker), REML=FALSE, data=geminates.clean.initial)
summary(geminates.clean.initial.lmer)

geminates.clean.initial.blmer = stan_lmer(length ~ consonant * c_length * frame + (1+consonant*c_length+frame|speaker) + (1+consonant*c_length+frame|word), data=geminates.clean.initial,
                             prior_intercept = normal(0, 1),
                             prior = normal(0, 2),
                             prior_covariance = decov(regularization = 2),
                             chains = 4,
                             iter = 2000)

plot(geminates.clean.initial.blmer,pars=c("(Intercept)","c_length1","consonantmw"))

sink("sink-summary-blmer-geminate-initial.txt")
print(summary(geminates.clean.initial.blmer),digits=4)
sink()

#Full medial does not converge
geminates.clean.medial.lmer <- lmer(length ~ consonant * c_length * frame + (1+consonant*c_length*frame|speaker), REML=FALSE, data=geminates.clean.medial)

geminates.clean.medial.blmer = stan_lmer(length ~ consonant * c_length * frame + (1+consonant*c_length+frame|speaker) + (1+consonant*c_length+frame|word), data=geminates.clean.medial,
                                          prior_intercept = normal(0, 1),
                                          prior = normal(0, 2),
                                          prior_covariance = decov(regularization = 2),
                                          chains = 4,
                                          iter = 2000)

plot(geminates.clean.medial.blmer,pars=c("(Intercept)","c_length1","consonantm","consonantmw","consonantn","consonantr"))

sink("sink-summary-blmer-geminate-medial.txt")
print(summary(geminates.clean.medial.blmer),digits=4)
sink()
#smaller model converges, though no sig interactions, so remove them
geminates.clean.medial.lmer1 <- lmer(length ~ consonant * c_length * frame + (1+consonant+c_length|speaker), REML=FALSE, data=geminates.clean.medial)
summary(geminates.clean.medial.lmer1)

#use this model for medial:
geminates.clean.medial.lmer2 <- lmer(length ~ consonant + c_length + (1+consonant+c_length|speaker), REML=FALSE, data=geminates.clean.medial)
summary(geminates.clean.medial.lmer2)

##almost full model for final
#converges!
geminates.clean.final.lmer <- lmer(length ~  c_length * consonant * frame + (1+consonant*c_length+frame|speaker), REML=FALSE, data=geminates.clean.final)
summary(geminates.clean.final.lmer) # no length distinction for final c

geminates.clean.final.blmer = stan_lmer(length ~ consonant * c_length * frame + (1+consonant*c_length+frame|speaker) + (1+consonant*c_length+frame|word), data=geminates.clean.final,
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0, 2),
                                         prior_covariance = decov(regularization = 2),
                                         chains = 4,
                                         iter = 2000)

plot(geminates.clean.final.blmer,pars=c("(Intercept)","c_length1","consonantmw"))

sink("sink-summary-blmer-geminate-final.txt")
print(summary(geminates.clean.final.blmer),digits=4)
sink()


### plots

# change order of levels for word_place
geminates.clean$word_place=factor(geminates.clean$word_place,levels=c("initial","medial","final"))

plot.gem <- ggplot(geminates.clean,aes(x=consonant,y=1000*length,fill=c_length)) + 
  geom_boxplot() + 
  facet_grid(~word_place) +
  theme_bw() +
  ylab("duration (ms)") +
  scale_fill_manual(name="length",values=c("#ffffff","#bfbfbf")) +
  theme(legend.position=c(.95,0.1))

plot.gem

pirateplot(length*1000~c_length+consonant+word_place,data=geminates.clean,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="Duration (ms)",gl.col = "white",xlab="Consonant",cex.lab=0.4)

### means for geminate pairs

geminates.means <- geminates.clean %>%
  group_by(consonant, c_length,word_place) %>%
  summarise(cmean = mean(length), stderr = sd(length)/sqrt(length(length)))

geminates.meansoverall <- geminates.clean %>%
  group_by(consonant, c_length) %>%
  summarise(means = mean(length))

################################################ End Geminate Section #####################################

################################################ Begin Apical Laminal Locus Section #######################
t_locus <- read.csv("t_d_locus.csv")


# clean data (removed 0 obs)
t_locus_clean <- t_locus %>%
  group_by(speaker,consonant,vowel) %>%
  mutate(sdmeanf2ss = mean(f2ss)+3*sd(f2ss),negsdmeanf2ss=mean(f2ss)-3*sd(f2ss),
         sdmeanf2i=mean(f2i)+3*sd(f2i),negsdmeanf2i=mean(f2i)-3*sd(f2i)) %>%
  filter(f2i<sdmeanf2i)%>%
  filter(f2i>negsdmeanf2i)%>%
  filter(f2ss<sdmeanf2ss)%>%
  filter(f2ss>negsdmeanf2ss)


### lms together

locus.all <- merge(pw_locus,t_locus,all=TRUE) # have to load pw_locus too

locus.all.lm <- lm(f2i~f2ss + consonant,data=locus.all)
locus.all.lmer <- lmer(f2i~f2ss +(1|consonant),data=locus.all)
sjp.lmer(locus.all.lmer)
# Bayesian LM
locus.blm <- stan_lmer(f2i ~ f2ss | consonant, data=locus.all) # high R2 prior
launch_shinystan(locus.blm) # plots
locus.blm
plot(locus.blm)
stan_plot(locus.blm$fit,show_density=T)

locus.t.lm <- lm(f2i ~ f2ss + consonant, data=t_locus_clean)

summary(glht(locus.t.lm, linfct=mcp(consonant = "Tukey")), test=adjusted(type="BH"))
summary(locus.t.lm)

sjp.lm(locus.t.lm)

TukeyHSD(locus.t.lm, "consonant")

plot(locus.t.lm)
#subset by consonant t=laminal, d=apical (from PNI orthography)
t_locus_t <- t_locus[t_locus$consonant=="t̻",]
t_locus_d <- t_locus[t_locus$consonant=="t",]

### t
t_locus.blmer = stan_lmer(f2i ~ f2ss + (1+f2ss|speaker), data=t_locus_t, 
                              prior_intercept = normal(0, 50),
                              prior = normal(0, 2),
                              prior_covariance = decov(regularization = 2), 
                              chains = 4,
                              iter = 2000,adapt_delta=0.9999)

blmer_t_locus_plot<-plot(t_locus.blmer)
plotly_POST(blmer_t_locus_plot, filename = "blmer_t_locus_plot", sharing="public")



draws <- as.data.frame(t_locus.blmer)
colnames(draws)[1:2] <- c("a", "b")

write.csv(draws, "draws.csv")
 # laminal
base <- ggplot(t_locus_t, aes(x = f2ss, y = f2i)) + 
  geom_point(size = 1) + xlim(1000,3000) + ylim(1000,3000) + theme_bw() + xlab("F2 steady state (Hz)")+
  ylab("F2 Initial (Hz)") + ggtitle("Laminal Locus Equations")
base

t_plot <- base + 
  geom_abline(data = draws, aes(intercept = a, slope = b), color = "skyblue", size = 0.2, alpha = 0.1)

grid.arrange(d_plot, t_plot,p_plot,pw_plot,  ncol=2)

(t_plotly <- ggplotly(t_plot))
plotly_POST(t_plot, filename = "qp1/t_locus", sharing="public")
#### /d/
draws2 <- as.data.frame(d_locus.blmer)
colnames(draws2)[1:2] <- c("a", "b")
base2 <- ggplot(t_locus_d, aes(x = f2ss, y = f2i)) + 
  geom_point(size = 1) + xlim(1000,3000) + ylim(1000,3000) + theme_bw() + xlab("F2 steady state (Hz)")+
  ylab("F2 Initial (Hz)") + ggtitle("/t/ Locus Equations")
base2

d_plot <- base2 + 
  geom_abline(data = draws2, aes(intercept = a, slope = b), 
              color = "seagreen4", size = 0.2, alpha = 0.1)
(d_plotly <- ggplotly(d_plot))
plotly_POST(d_plotly, filename = "qp1/d_locus", world_readable=TRUE)

base3 <- ggplot(t_locus, aes(x = f2ss, y = f2i)) + 
  geom_point(size = 1)
base3 + 
  geom_abline(data = draws2, aes(intercept = a, slope = b), 
              color = "skyblue", size = 0.2, alpha = 0.25) +
  geom_abline(data = draws, aes(intercept = a, slope = b), 
              color = "green", size = 0.2, alpha = 0.25)
 
sink("sink-summary-blmer-t_locus.txt")
print(summary(t_locus.blmer),digits=4)
sink()

### d
d_locus.blmer = stan_lmer(f2i ~ f2ss + (1+f2ss|speaker), data=t_locus_d, 
                          prior_intercept = normal(0, 50),
                          prior = normal(0, 2),
                          prior_covariance = decov(regularization = 2), 
                          chains = 4,
                          iter = 2000, adapt_delta = 0.999 )

blmer_d_locus_plot<-plot(d_locus.blmer)
plotly_POST(blmer_d_locus_plot, filename = "blmer_d_locus_plot", sharing="public")


sink("sink-summary-blmer-d_locus.txt")
print(summary(d_locus.blmer),digits=4)
sink()


#subset by speaker
t_locus.m01 <- t_locus[t_locus$speaker=="M01",]
t_locus.m01.t <- t_locus.m01[t_locus.m01$consonant=="t̻",]
t_locus.m01.d <- t_locus.m01[t_locus.m01$consonant=="t",]
t_locus.m02 <- t_locus[t_locus$speaker=="M02",]
t_locus.m02.t <- t_locus.m02[t_locus.m02$consonant=="t̻",]
t_locus.m02.d <- t_locus.m02[t_locus.m02$consonant=="t",]
t_locus.m03 <- t_locus[t_locus$speaker=="M03",]
t_locus.m03.t <- t_locus.m03[t_locus.m03$consonant=="t̻",]
t_locus.m03.d <- t_locus.m03[t_locus.m03$consonant=="t",]
t_locus.f01 <- t_locus[t_locus$speaker=="F01",]
t_locus.f01.t <- t_locus.f01[t_locus.f01$consonant=="t̻",]
t_locus.f01.d <- t_locus.f01[t_locus.f01$consonant=="t",]
t_locus.f02 <- t_locus[t_locus$speaker=="F02",]
t_locus.f02.t <- t_locus.f02[t_locus.f02$consonant=="t̻",]
t_locus.f02.d <- t_locus.f02[t_locus.f02$consonant=="t",]

#linear regressions
#for laminal
t_locus_t.lm <- lm(f2i ~ f2ss, data=t_locus_t)
summary(t_locus_t.lm)
#for apical
t_locus_d.lm <- lm(f2i ~ f2ss, data=t_locus_d)
summary(t_locus_d.lm)

x <- c(0.5,1,1.5,2,1,2,3,4)
y<- c(1,2,3,4,1,2,3,4)
c <- c("t","t","t","t","d","d","d","d")
df <- data.frame(x,y,c)

x1 <- c(0.5,1,1.5,2)
y1<- c(1,2,3,4)

df2 <- data.frame(x1,y1)

lm1 <- lm(data=df2, y1 ~ x1)
lm1
summary(lm1)
       
ggplot(df,aes(x=x,y=y,color=c)) + geom_smooth(method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5)
#by speaker
#m01
t_locus.m01.t.lm <-lm(f2i ~ f2ss, data=t_locus.m01.t)
summary(t_locus.m01.t.lm)
t_locus.m01.d.lm <-lm(f2i ~ f2ss, data=t_locus.m01.d)
summary(t_locus.m01.d.lm)
#m02
t_locus.m02.t.lm <-lm(f2i ~ f2ss, data=t_locus.m02.t)
summary(t_locus.m02.t.lm)
t_locus.m02.d.lm <-lm(f2i ~ f2ss, data=t_locus.m02.d)
summary(t_locus.m02.d.lm)
#m03
t_locus.m03.t.lm <-lm(f2i ~ f2ss, data=t_locus.m03.t)
summary(t_locus.m03.t.lm)
t_locus.m03.d.lm <-lm(f2i ~ f2ss, data=t_locus.m03.d)
summary(t_locus.m03.d.lm)
#f01
t_locus.f01.t.lm <-lm(f2i ~ f2ss, data=t_locus.f01.t)
summary(t_locus.f01.t.lm)
t_locus.f01.d.lm <-lm(f2i ~ f2ss, data=t_locus.f01.d)
summary(t_locus.f01.d.lm)
#f02
t_locus.f02.t.lm <-lm(f2i ~ f2ss, data=t_locus.f02.t)
summary(t_locus.f02.t.lm)
t_locus.f02.d.lm <-lm(f2i ~ f2ss, data=t_locus.f02.d)
summary(t_locus.f02.d.lm)

# rename to laminal and apical since ipa symbols hard to discern difference
t_locus <- t_locus %>%
  mutate(consonant = ifelse(consonant=="t̻", "tt","t"))

#plotting the lms
t_d_locus_plot_speaker<-ggplot(t_locus, aes(y=f2i, x=f2ss, color=consonant, shape=consonant, linetype=consonant)) + 
  geom_point() +
  theme_bw() + 
  geom_smooth(method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5) + 
  facet_wrap(~speaker,scales="free") +
  ylab("F2 initial (Hz)") + 
  xlab("F2 steady state (Hz)") +
  scale_shape_manual(values=c(1,2))+
  scale_colour_manual(values=c("#000000","#900009"))+
                      #name="consonant", 
                      #breaks=c("t", "t̻"), 
                      #labels = c("apical", "laminal")) +
  #scale_linetype_discrete()+#name="consonant", 
                         # breaks=c("t", "t̻"), 
                          #labels = c("apical", "laminal")) +
  #guides(shape=F)+
  theme(legend.position=c(.8,0.3))

(t_d_locus_speaker_plotly <- ggplotly(t_d_locus_plot_speaker))
plotly_POST(t_d_locus_speaker_plotly, filename = "t_d_locus_speaker", sharing="public")

ggplot(t_locus, aes(y=f2ss, x=consonant, color=vowel)) + 
  geom_boxplot() +
  theme_bw() + 
 # geom_smooth(method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5) + 
 # facet_wrap(~speaker,scales="free") +
  ylab("F2 steady state (Hz)") + 
  xlab("Consonant") +
 # scale_shape_manual(values=c(1,20))+
  scale_colour_manual(values=c("#000000","#900009"))
  #name="consonant", 
  #breaks=c("t", "t̻"), 
  #labels = c("apical", "laminal")) +
  #scale_linetype_discrete()+#name="consonant", 
  # breaks=c("t", "t̻"), 
  #labels = c("apical", "laminal")) +
  #guides(shape=F)+
#  theme(legend.position=c(.8,0.3))

## merge locus together


locus_all_data <- merge(pw_locus,t_locus,all=T)
library(grid)
library(gridExtra)
pirate1 <- pirateplot(f2i~vowel+consonant,data=locus_all_data,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="F2 Initial (Hz)",gl.col = "white",ylim=c(600,3000))
pirate2<- pirateplot(f2ss~vowel+consonant,data=locus_all_data,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="F2 Steady State (Hz)",gl.col = "white",ylim=c(600,3000))

pirateplot(f2i~vowel+consonant,data=t_locus,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="F2 Initial (Hz)",gl.col = "white",ylim=c(1200,3200))
 pirateplot(f2ss~vowel+consonant,data=t_locus,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="F2 Steady State (Hz)",gl.col = "white",ylim=c(1200,3200))


locus_f2i_vowel<-ggplot(NULL, aes(y=f2i, x=consonant, fill=vowel)) + 
  geom_boxplot(data=t_locus) +
  geom_boxplot(data=pw_locus) +
  theme_bw() + 
  # geom_smooth(method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5) + 
  # facet_wrap(~speaker,scales="free") +
  ylab("F2 initial (Hz)") + 
  xlab("Consonant") +
  theme(legend.position=c(.95,0.09)) +
  # scale_shape_manual(values=c(1,20))+
   # scale_colour_manual(values=c("#000000","#900009","#666666","#145314"))
  scale_fill_manual(values=c("#A9A9A9","#ffffff")) + ylim(500,3000)



locus_f2ss_vowel<-ggplot(NULL, aes(y=f2ss, x=consonant, fill=vowel)) + 
  geom_boxplot(data=t_locus) +
  geom_boxplot(data=pw_locus) +
  theme_bw() + 
  # geom_smooth(method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5) + 
  # facet_wrap(~speaker,scales="free") +
  ylab("F2 steady state (Hz)") + 
  xlab("Consonant") +
  theme(legend.position=c(.95,0.09)) +
  # scale_shape_manual(values=c(1,20))+
  # scale_colour_manual(values=c("#000000","#900009","#666666","#145314"))
  scale_fill_manual(values=c("#A9A9A9","#ffffff")) + ylim(500,3000)

require(gridExtra)
grid.arrange(locus_f2i_vowel, locus_f2ss_vowel, ncol=2)

locus_f2i<-ggplot(NULL, aes(y=f2i, x=consonant)) + 
  geom_boxplot(data=t_locus) +
  geom_boxplot(data=pw_locus) +
  theme_bw() + 
  # geom_smooth(method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5) + 
  # facet_wrap(~speaker,scales="free") +
  ylab("F2 initial (Hz)") + 
  xlab("Consonant") +
  theme(legend.position=c(.95,0.09)) +
  # scale_shape_manual(values=c(1,20))+
  # scale_colour_manual(values=c("#000000","#900009","#666666","#145314"))
  scale_fill_manual(values=c("#A9A9A9","#ffffff")) + ylim(1000,3000)
locus_f2i

(locus_f2i_plotly <- ggplotly(locus_f2i))
plotly_POST(locus_f2i_plotly, filename = "locus_f2i", sharing="public")

ggplot(NULL, aes(y=f2ss, x=consonant, fill=vowel)) + 
  geom_boxplot(data=t_locus) +
  geom_boxplot(data=pw_locus) +
  theme_bw() + 
  #scale_fill_brewer(palette="RdBu") +
  # geom_smooth(method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5) + 
  # facet_wrap(~speaker,scales="free") +
  ylab("F2 steady state (Hz)") + 
  xlab("Consonant") +
  # scale_shape_manual(values=c(1,20))+
  scale_fill_manual(values=c("#A9A9A9","#ffffff")) +
#name="consonant", 
#breaks=c("t", "t̻"), 
#labels = c("apical", "laminal")) +
#scale_linetype_discrete()+#name="consonant", 
# breaks=c("t", "t̻"), 
#labels = c("apical", "laminal")) +
#guides(shape=F)+
  theme(legend.position=c(.95,0.1)) 

locus_f2ss <-ggplot(NULL, aes( x=consonant)) + 
  geom_boxplot(data=t_locus,aes(y=f2ss)) +

  geom_boxplot(data=pw_locus,aes(y=f2ss)) +

  theme_bw() + 
  #scale_fill_brewer(palette="RdBu") +
  # geom_smooth(method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5) + 
  # facet_wrap(~speaker,scales="free") +
  ylab("F2 steady state (Hz)") + 
  xlab("Consonant") +
  # scale_shape_manual(values=c(1,20))+
  scale_fill_manual(values=c("#A9A9A9","#ffffff")) +
  #name="consonant", 
  #breaks=c("t", "t̻"), 
  #labels = c("apical", "laminal")) +
  #scale_linetype_discrete()+#name="consonant", 
  # breaks=c("t", "t̻"), 
  #labels = c("apical", "laminal")) +
  #guides(shape=F)+
  theme(legend.position=c(.95,0.1)) +
  ylim(1000,3000)
locus_f2ss




#plotting the lms
locus_plots<-ggplot(NULL, aes(y=f2i, x=f2ss, color=consonant, shape=consonant, linetype=consonant)) + 
  geom_point(data=t_locus) + geom_point(data=pw_locus) +
  theme_bw() + 
  geom_smooth(data=t_locus,method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5)  +
    geom_smooth(data=pw_locus,method=lm,fullrange=T,se=F,alpha=.15,size=0.5)+
  ylab("F2 initial (Hz)") + 
  xlab("F2 steady state (Hz)") +
  scale_shape_manual(values=c(0,1,5,6))+
  scale_colour_manual(values=c("#000000","#900009","#666666","#145314"))+
  scale_linetype_manual(values=c("solid", "longdash", "dotted", "twodash")) +
  #name="consonant", 
  #breaks=c("t", "t̻"), 
  #labels = c("apical", "laminal")) +
  #scale_linetype_discrete()+#name="consonant", 
  # breaks=c("t", "t̻"), 
  #labels = c("apical", "laminal")) +
  #guides(shape=F)+
  theme(legend.position=c(.8,0.3))


t_locus_means <- t_locus_clean %>%
  group_by(consonant) %>%
  summarise(meanf2i = mean(f2i), meanf2ss = mean(f2ss))

t_locus_means_vowels <- t_locus_clean %>%
  group_by(consonant, vowel) %>%
  summarise(meanf2i = mean(f2i), meanf2ss = mean(f2ss))

################################################################## End Laminal Locus Section #################################

################################################################## Begin P-Pw Locus Section ##################################
pw_locus <- read.csv("p_pw_locus.csv")
# clean data (removed 0 obs)
pw_locus_clean <- pw_locus %>%
  group_by(speaker,consonant,vowel) %>%
  mutate(sdmeanf2ss = mean(f2ss)+3*sd(f2ss),negsdmeanf2ss=mean(f2ss)-3*sd(f2ss),
         sdmeanf2i=mean(f2i)+3*sd(f2i),negsdmeanf2i=mean(f2i)-3*sd(f2i)) %>%
  filter(f2i<sdmeanf2i)%>%
  filter(f2i>negsdmeanf2i)%>%
  filter(f2ss<sdmeanf2ss)%>%
  filter(f2ss>negsdmeanf2ss)

#subset by consonant p='plain p', pw=labialized-velarized p (from PNI orthography)
pw_locus_pw <- pw_locus[pw_locus$consonant=="pw",]
pw_locus_p <- pw_locus[pw_locus$consonant=="p",]

### pw
pw_locus.blmer = stan_lmer(f2i ~ f2ss + (1+f2ss|speaker), data=pw_locus_pw, 
                          prior_intercept = normal(0, 50),
                          prior = normal(0, 2),
                          prior_covariance = decov(regularization = 2), 
                          chains = 4,
                          iter = 2000, adapt_delta = 0.999 )

plot(pw_locus.blmer)

draws3 <- as.data.frame(pw_locus.blmer)
colnames(draws3)[1:2] <- c("a", "b")
base3 <- ggplot(pw_locus_pw, aes(x = f2ss, y = f2i)) + 
  geom_point(size = 1) + xlim(500,3000) + ylim(500,3000) + theme_bw() + xlab("F2 steady state (Hz)")+
  ylab("F2 Initial (Hz)") + ggtitle("/pw/ Locus Equations")
base3

base3 + 
  geom_abline(data = draws3, aes(intercept = a, slope = b), 
              color = "orchid4", size = 0.2, alpha = 0.1)


sink("sink-summary-blmer-pw_locus.txt")
print(summary(pw_locus.blmer),digits=4)
sink()

### p
p_locus.blmer = stan_lmer(f2i ~ f2ss + (1+f2ss|speaker), data=pw_locus_p, 
                           prior_intercept = normal(0, 50),
                           prior = normal(0, 2),
                           prior_covariance = decov(regularization = 2), 
                           chains = 4,
                           iter = 2000, adapt_delta = 0.999 )

plot(p_locus.blmer)

draws4 <- as.data.frame(p_locus.blmer)
colnames(draws4)[1:2] <- c("a", "b")
base4 <- ggplot(pw_locus_p, aes(x = f2ss, y = f2i)) + 
  geom_point(size = 1) + xlim(1000,3000) + ylim(1000,3000) + theme_bw() + xlab("F2 steady state (Hz)")+
  ylab("F2 Initial (Hz)") + ggtitle("/p/ Locus Equations")
base4

base4 + 
  geom_abline(data = draws4, aes(intercept = a, slope = b), 
              color = "darkred", size = 0.2, alpha = 0.1)

sink("sink-summary-blmer-p_locus.txt")
print(summary(p_locus.blmer),digits=4)
sink()






#subset by speaker
pw_locus.m01 <- pw_locus[pw_locus$speaker=="M01",]
pw_locus.m01.pw <- pw_locus.m01[pw_locus.m01$consonant=="pw",]
pw_locus.m01.p <- pw_locus.m01[pw_locus.m01$consonant=="p",]
pw_locus.m02 <- pw_locus[pw_locus$speaker=="M02",]
pw_locus.m02.pw <- pw_locus.m02[pw_locus.m02$consonant=="pw",]
pw_locus.m02.p <- pw_locus.m02[pw_locus.m02$consonant=="p",]
pw_locus.m03 <- pw_locus[pw_locus$speaker=="M03",]
pw_locus.m03.pw <- pw_locus.m03[pw_locus.m03$consonant=="pw",]
pw_locus.m03.p <- pw_locus.m03[pw_locus.m03$consonant=="p",]
pw_locus.f01 <- pw_locus[pw_locus$speaker=="F01",]
pw_locus.f01.pw <- pw_locus.f01[pw_locus.f01$consonant=="pw",]
pw_locus.f01.p <- pw_locus.f01[pw_locus.f01$consonant=="p",]
pw_locus.f02 <- pw_locus[pw_locus$speaker=="F02",]
pw_locus.f02.pw <- pw_locus.f02[pw_locus.f02$consonant=="pw",]
pw_locus.f02.p <- pw_locus.f02[pw_locus.f02$consonant=="p",]

#linear regressions
#for pw
pw_locus_pw.lm <- lm(f2i ~ f2ss, data=pw_locus_pw)
summary(pw_locus_pw.lm)
#for p
pw_locus_p.lm <- lm(f2i ~ f2ss, data=pw_locus_p)
summary(pw_locus_p.lm)

#by speaker
#m01
pw_locus.m01.pw.lm <-lm(f2i ~ f2ss, data=pw_locus.m01.pw)
summary(pw_locus.m01.pw.lm)
pw_locus.m01.p.lm <-lm(f2i ~ f2ss, data=pw_locus.m01.p)
summary(pw_locus.m01.p.lm)
#m02
pw_locus.m02.pw.lm <-lm(f2i ~ f2ss, data=pw_locus.m02.pw)
summary(pw_locus.m02.pw.lm)
pw_locus.m02.p.lm <-lm(f2i ~ f2ss, data=pw_locus.m02.p)
summary(pw_locus.m02.p.lm)
#m03
pw_locus.m03.pw.lm <-lm(f2i ~ f2ss, data=pw_locus.m03.pw)
summary(pw_locus.m03.pw.lm)
pw_locus.m03.p.lm <-lm(f2i ~ f2ss, data=pw_locus.m03.p)
summary(pw_locus.m03.p.lm)
#f01
pw_locus.f01.pw.lm <-lm(f2i ~ f2ss, data=pw_locus.f01.pw)
summary(pw_locus.f01.pw.lm)
pw_locus.f01.p.lm <-lm(f2i ~ f2ss, data=pw_locus.f01.p)
summary(pw_locus.f01.p.lm)
#f02
pw_locus.f02.pw.lm <-lm(f2i ~ f2ss, data=pw_locus.f02.pw)
summary(pw_locus.f02.pw.lm)
pw_locus.f02.p.lm <-lm(f2i ~ f2ss, data=pw_locus.f02.p)
summary(pw_locus.f02.p.lm)

#plotting the lms
ggplot(pw_locus_clean, aes(y=f2i, x=f2ss, color=consonant, shape=consonant, linetype=consonant)) + 
  geom_point() +
  theme_bw() + 
  geom_smooth(method=lm, fullrange=TRUE, se=FALSE,alpha = .15,size=0.5)  +
  ylab("F2 initial (Hz)") + 
  xlab("F2 steady state (Hz)") +
  scale_shape_manual(values=c(1,2))+
  scale_colour_manual(values=c("#000000","#900009"))+
  #name="consonant", 
  #breaks=c("t", "t̻"), 
  #labels = c("apical", "laminal")) +
  #scale_linetype_discrete()+#name="consonant", 
  # breaks=c("t", "t̻"), 
  #labels = c("apical", "laminal")) +
  #guides(shape=F)+
  theme(legend.position=c(.8,0.3))



################################################################ End P-Pw Section ##################################################

################################################################ Begin P-PW VOT ############################################
#import csv
p_vot <- read.csv("p_pw_vot.csv")

# remove outliers (1 removed)
p_vot_clean = p_vot %>%
  group_by(speaker_id,consonant2)%>%
  mutate(votsdmean = mean(vot) + 3*sd(vot),negvotsdmean = mean(vot)-3*sd(vot))%>%
  filter(vot<votsdmean)%>%
  filter(vot>negvotsdmean)
# still think there is one outlier for p and that it is measurment error or noise in recording error

# there is one outlier so readjust filter to 2.5sd+mean (total removed 2 including previous filter)
p_vot_clean = p_vot %>%
  group_by(speaker_id,consonant2)%>%
  mutate(votsdmean = mean(vot) + 2.5*sd(vot))%>%
  filter(vot<votsdmean)


#subset by speaker
p_vot.m01 <- p_vot[p_vot$speaker_id=="M01",]
p_vot.m02 <- p_vot[p_vot$speaker_id=="M02",]
p_vot.m03 <- p_vot[p_vot$speaker_id=="M03",]
p_vot.f01 <- p_vot[p_vot$speaker_id=="F01",]
p_vot.f02 <- p_vot[p_vot$speaker_id=="F02",]



#plot cleaned data
ggplot(p_vot_clean, aes(x=consonant2, y=vot)) + geom_boxplot() + theme_bw() + facet_grid(~speaker_id)


# nice plot
vot_plot<- ggplot(NULL, aes(x=consonant2,y=vot)) + 
  geom_boxplot(data=p_vot_clean) + geom_boxplot(data=t_vot) +
  theme_bw() + 
  ylab("VOT (s)") + xlab("consonant")

 

#lmers
#full model
p_vot_clean.lmer <- lmer(vot ~ consonant + (1+consonant|speaker_id) + (1+consonant|word), data=p_vot_clean, REML=FALSE)
summary(p_vot_clean.lmer)

p_vot.clean.blmer = stan_lmer(vot ~ consonant + (1+consonant|speaker_id) + (1+consonant|word), data=p_vot_clean, 
                              prior_intercept = normal(0, 5),
                              prior = normal(0, 2),
                              prior_covariance = decov(regularization = 2), 
                              chains = 4,
                              iter = 2000)

blmer_p_vot_plot<-plot(p_vot.clean.blmer)
plotly_POST(blmer_p_vot_plot, filename = "blmer_p_vot_plot", sharing="public")


sink("sink-summary-blmer-p_vot.txt")
print(summary(p_vot.clean.blmer),digits=4)
sink()

#check assumptions

qqnorm(residuals(p_vot_clean.lmer)) # normal
plot(fitted(p_vot_clean.lmer),residuals(p_vot_clean.lmer)) ## hmm. is that heteroscedastic? or just limited data points?

sjp.lmer(p_vot_clean.lmer, type="fe.pred") # that is interesting looking 

#model converges and no sig different in vot for the two segments

# mean VOT
p_means <- p_vot_clean %>%
  group_by(consonant2) %>%
  summarise(meanvot = mean(vot),stderr = sd(vot)/sqrt(length(vot)))
####################################################### End P-Pw VOT ############################################

####################################################### Begin Laminal-Apical VOT ################################
#import csv
t_vot <- read.csv("t_d_vot.csv")

#subset by speaker
t_vot.m01 <- t_vot[t_vot$speaker_id=="M01",]
t_vot.m02 <- t_vot[t_vot$speaker_id=="M02",]
t_vot.m03 <- t_vot[t_vot$speaker_id=="M03",]
t_vot.f01 <- t_vot[t_vot$speaker_id=="F01",]
t_vot.f02 <- t_vot[t_vot$speaker_id=="F02",]

# plot box plots to check for outliers
# change names of consonant for graphing
t_vot <- t_vot %>%
  mutate(consonant2 = ifelse(consonant2 == "t̻","tt","t"))

ggplot(t_vot, aes(x=consonant2, y=vot)) + geom_boxplot() + theme_bw() + facet_grid(~speaker_id)

ggplot(t_vot, aes(x=consonant2,y=vot)) + geom_boxplot() + theme_bw() + xlab("consonant") + ylab("VOT (s)")

#no outliers so good to go

#lmers
#full model
t_vot.lmer <- lmer(vot ~ consonant + (1+consonant|speaker_id) + (1+consonant|word), data=t_vot, REML=FALSE)
summary(t_vot.lmer)

t_vot.clean.blmer = stan_lmer(vot ~ consonant + (1+consonant|speaker_id) + (1+consonant|word), data=t_vot, 
                                      prior_intercept = normal(0, 5),
                                      prior = normal(0, 2),
                                      prior_covariance = decov(regularization = 2), 
                                      chains = 4,
                                      iter = 2000)

blmer_t_vot_plot<-plot(t_vot.clean.blmer)
#(blmer_t_vot_plot_plotly <- ggplotly(blmer_skew_plot))
plotly_POST(blmer_t_vot_plot, filename = "blmer_t_vot_plot", sharing="public")


sink("sink-summary-blmer-t_vot.txt")
print(summary(t_vot.clean.blmer),digits=4)
sink()

# some random diagnostic plots 
sjp.lmer(t_vot.lmer, type="fe.ri")
sjp.lmer(t_vot.lmer, type="fe.pred")
#it converges and there is a sig different in vot for the two segments

## mean vot

t_means <- t_vot %>%
  group_by(consonant2) %>% summarise(means = mean(vot), stderr= sd(vot)/sqrt(length(vot)))

## combined pirate plots

pirateplot(vot*1000~consonant2,data=t_vot,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="VOT (ms)",gl.col = "white",xlab="Consonant",pal="gray",bean.f.col=c("seagreen4","skyblue"))

vot_data_all<-merge(p_vot_clean,t_vot,all=T)
 pirateplot(vot*1000~consonant2,data=vot_data_all,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="VOT (ms)",gl.col = "white",xlab="Consonant",pal="gray")
################################################## End Laminal-Apical VOT #######################


##################################################### Start vowel allophony ##############################
# import csv
vowel <- read.csv("vowel_allophony.csv")

#plot boxplots and check for outliers
#f2
ggplot(vowel, aes(x=vowel, y=f2, fill=v_length)) + geom_boxplot() + theme_bw() + facet_grid(c_initial~speaker_id)

# clean data (removed 26 out of 771 obs)
vowel.clean = vowel %>%
  group_by(speaker_id,vowel,v_length,c_initial,c_final)%>%
  mutate(f2sdmean = mean(f2) + 3*sd(f2),negf2sdmean = mean(f2)-3*sd(f2),f1sdmean = mean(f1) + 3*sd(f1),negf1sdmean = mean(f1)-3*sd(f1))%>%
  filter(f2<f2sdmean)%>%
  filter(f2>negf2sdmean)%>%
  filter(f1<f1sdmean)%>%
  filter(f1>negf1sdmean) 



# there appear to be some measurement errors still with some very extreme values
# there are 2 F2 values for e for M01 that are really low (760Hz and below) so remove them (rows 172 and 174)
# also remove 1 F2 value for u for M02 above 2000 Hz, row 246
vowel.clean <- vowel.clean[-c(172,174,246),]
# but there are still within 3sigma
#plot cleaned data
ggplot(vowel.clean, aes(x=vowel, y=f2, fill=c_initial)) + geom_boxplot() + theme_bw() + facet_grid(v_length~speaker_id) 


# add contrast coding for v_length, c_initial, c_final, gender
# check original coding
contrasts(vowel.clean$v_length) <- c(-0.5,0.5)
contrasts(vowel.clean$gender) <- c(-0.5,0.5)
contrasts(vowel.clean$c_initial) <- c(-0.5,0.5)
contrasts(vowel.clean$c_final) <- c(-0.5,0.5)

#lmers
#f2 full model . It also fails to converge so will try smaller one
vowel.clean.lmer0 <- lmer(f2 ~ vowel * gender * c_initial * c_final * v_length + (1+vowel*gender*c_initial*c_final*v_length|speaker_id) + (1+vowel*gender*c_initial*c_final+v_length|word), data=vowel.clean, REML=FALSE)

options (mc.cores=parallel::detectCores ()) # Run on multiple cores

set.seed (3875)
# don't run
vowel.clean.blmer0 = stan_lmer(f2 ~ vowel * gender * c_initial * c_final * v_length + (1+vowel*gender*c_initial*c_final*v_length|speaker_id) + (1+vowel*gender*c_initial*c_final+v_length|word), data=vowel.clean, 
                            prior_intercept = normal(0, 40),
                            prior = normal(0, 20),
                            prior_covariance = decov(regularization = 2), 
                            chains = 4,
                            iter = 2000)
launch_shinystan(vowel.clean.blmer0)


loo_bglmer0 <- loo(vowel.clean.blmer0)
#

# maximal convergent model with only sig interactions 
vowel.clean.lmer <- lmer(f2 ~ vowel *v_length +gender +  c_initial + c_final+ (1+vowel|speaker_id) + (1|word), data=vowel.clean, REML=FALSE)
summary(vowel.clean.lmer)

vowel.clean.blmer = stan_lmer(f2 ~ vowel *v_length *gender +  c_initial * c_final+ (1+vowel|speaker_id) + (1+vowel|word), data=vowel.clean, 
                               prior_intercept = normal(0, 50),
                               prior = normal(0, 30),
                               prior_covariance = decov(regularization = 2), 
                               chains = 4,
                               iter = 2000)





loo.vowel.clean.blmer = loo(vowel.clean.blmer)

vowel.clean.blmer2 = stan_lmer(f2 ~ vowel *v_length +gender + vowel:gender +  c_initial * c_final+ (1+vowel|speaker_id) + (1+vowel|word), data=vowel.clean, 
                              prior_intercept = normal(0, 40),
                              prior = normal(0, 20),
                              prior_covariance = decov(regularization = 2), 
                              chains = 4,
                              iter = 2000,adapt_delta=0.9999)


sink("sink-summary-blmer-vowels.txt")
print(summary(vowel.clean.blmer2),digits=4)
sink()

loo.vowel.clean.blmer2 = loo(vowel.clean.blmer2)
compare(loo.vowel.clean.blmer2,loo.vowel.clean.blmer)

plot(vowel.clean.blmer2,pars=c("(Intercept)","v_length1","gender1","c_initial1","c_final1","c_initial1:c_final1"))

stan_hist(vowel.clean.blmer2,pars=c("c_initial1","c_final1","c_initial1:c_final1"))
plot(vowel.clean.blmer)
launch_shinystan(vowel.clean.blmer2)
stan_plot(vowel.clean.blmer2,point_est="median",pars=c("c_initial1","c_final1","c_initial1:c_final1"))

pp_check(vowel.clean.blmer2, check = "test", test="median", overlay = FALSE, nreps = 5)
stan_hist(vowel.clean.blmer2)

#visualize the bayes models
data(radon)
library(ggmcmc)
L.radon <- data.frame(
  Parameter=c(
    "alpha[1]",
    seq(1,446,1)),
  Label=rep(radon$counties$County, 2),
  Uranium=rep(radon$counties$uranium, 2),
  Location=rep(radon$counties$ns.location, 2),
  Coefficient=gl(2, length(radon$counties$id.county), 
                 labels=c("Intercept", "Slope")))

vowel.clean.blmer2.df <- as.data.frame(vowel.clean.blmer2)
vowel.bayes1<-ggs(vowel.clean.blmer2) # convert stan output to df
vowel.bayes2<-ggs(vowel.clean.blmer2,)
ggmcmc(vowel.bayes1,file="model_vowels.pdf") # plot several outputs as a pdf with default 5 plots per page
model <- mcmc.list(lapply(1:ncol(vowel.clean.blmer2), function(x) mcmc(as.array(vowel.clean.blmer2)[,x,])))
plot(vowel.clean.blmer2)
# check residuals 

qqnorm(residuals(vowel.clean.lmer))
plot(fitted(vowel.clean.lmer),residuals(vowel.clean.lmer))
## add the residuals from unbalanced.model3 to the dataframe
vowel.clean$resids <- residuals(vowel.clean.lmer)


## make a scatterplot to check the distribution of residuals
ggplot(vowel.clean, aes(x=speaker_id, y=resids))+geom_point()
# M01 and M02 are somewhat diff but still seem ok?

## this is cool...looks random which is good but see some  possible outliers
ggplot(vowel.clean, aes(x=speaker_id, y=resids, color=vowel))+geom_point()

## how about for c_final...pretty random
ggplot(vowel.clean,aes(x=speaker_id,y=resids,color=c_final))+geom_point()

### residuals look pretty random which is good overall



### how about for c_initial...pretty random
ggplot(vowel.clean, aes(x=speaker_id,y=resids,color=c_initial))+geom_point()
# now for v_length...pretty random
ggplot(vowel.clean,aes(x=speaker_id,y=resids,color=v_length))+geom_point()




# visualize lmer (interesting but can I used them?) [note they create several graphs, so click forward and backward in your plots area to see all the ones created]
sjp.lmer(vowel.clean.lmer, type="fe.pred", facet.grid=TRUE)
sjp.lmer(vowel.clean.lmer, type="fe.ri", facet.grid=TRUE)

# subset for short vowel
vowel.clean.short <- vowel.clean %>%
  filter(v_length=="short")
#subset for long vowels
vowel.clean.long <- vowel.clean%>%
  filter(v_length=="long")

### find group averages for front and back for both long and short

vowel.short.avg <- vowel.clean.short %>%
  group_by(vowel, c_initial)%>%
  summarise(mean.f2 = mean(f2),mean.f1 = mean(f1))

vowel.short.avg

vowel.long.avg <- vowel.clean.long %>%
  group_by(vowel, c_initial) %>%
  summarise(mean.f2 = mean(f2),mean.f1=mean(f1))
vowel.long.avg

vowel.mean <- vowel.clean %>%
  group_by(vowel, v_length)%>%
  summarise(mean.f1 = mean(f1), mean.f2 = mean(f2))

vowel.mean
### lmers for short vowels
# full model again (does not converge)
vowel.clean.short.lmer <- lmer(f2~c_initial*c_final*gender*vowel+(1+c_initial*c_final*gender|speaker_id)+(1*c_initial*c_final*gender|word),data=vowel.clean.short)

# remove non-sig interactions in predictors (doesn't converge at first so took out interaction in random slopes, the slopes, then predicotr interactions until converged)
vowel.clean.short.lmer1 <- lmer(f2 ~ vowel *   c_initial + c_final + gender+ (1+vowel|speaker_id) + (1|word), data=vowel.clean.short, REML=FALSE)
summary(vowel.clean.short.lmer1)

vowel.clean.short.lmer2 <- lmer(f2~vowel + c_initial + c_final + gender + (1+vowel|speaker_id) + (1|word),data=vowel.clean.short,REML=F)
anova(vowel.clean.short.lmer1,vowel.clean.short.lmer2) # not sig dif so go with 2
summary(vowel.clean.short.lmer2) # go with this model

vowel.clean.short.blmer2 = stan_lmer(f2~vowel + c_initial + c_final + gender + (1+vowel|speaker_id) + (1|word),data=vowel.clean.short, 
                               prior_intercept = normal(0, 40),
                               prior = normal(0, 20),
                               prior_covariance = decov(regularization = 2), 
                               chains = 4,
                               iter = 2000)

plot(vowel.clean.short.blmer2)

sink("sink-summary-blmer2-short.txt")
summary(vowel.clean.short.blmer2)
sink()

vowel.clean.short.blmer3 = stan_lmer(f2~vowel + c_initial * c_final + gender + vowel:gender + (1+vowel|speaker_id) + (1+vowel|word),data=vowel.clean.short, 
                                     prior_intercept = normal(0, 40),
                                     prior = normal(0, 20),
                                     prior_covariance = decov(regularization = 2), 
                                     chains = 4,
                                     iter = 2000)

loo.vowel.short3 = loo(vowel.clean.short.blmer3)
loo.vowel.short2 = loo(vowel.clean.short.blmer2)
compare(loo.vowel.short2,loo.vowel.short3)

plot(vowel.clean.short.blmer3)

sink("sink-summary-blmer3-short.txt")
summary(vowel.clean.short.blmer3)
sink()


#### lmers for long vowels

# full model for long (doesn't converge)
vowel.clean.long.lmer <- lmer(f2~vowel * c_initial * c_final * gender + (1+vowel*c_initial*c_final*gender|speaker_id)+(1+vowel*c_initial*c_final*gender|word),data=vowel.clean.long,REML=F)
summary(vowel.clean.lmer)

#reduce until converges and remove no sig predictors and interactions
vowel.clean.long.lmer1 <- lmer(f2~vowel * c_initial *gender * c_final + (1+vowel|speaker_id) + (1|word),data=vowel.clean.long,REML=F)
summary(vowel.clean.long.lmer1)
# c_initial is sig factor for both long and short!
vowel.clean.long.blmer3 = stan_lmer(f2~vowel + c_initial * c_final + gender + vowel:gender + (1+vowel|speaker_id) + (1+vowel|word),data=vowel.clean.long, 
                                     prior_intercept = normal(0, 40),
                                     prior = normal(0, 20),
                                     prior_covariance = decov(regularization = 2), 
                                     chains = 4,
                                     iter = 2000)

plot(vowel.clean.short.blmer3)

sink("sink-summary-blmer3-long.txt")
summary(vowel.clean.short.blmer3)
sink()

pirateplot(f2~c_initial+vowel+v_length,data=vowel.clean,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="F2 (Hz)",gl.col = "white",xlab="Vowel",cex.lab=0.4)
pirateplot(f2~c_final+vowel+v_length,data=vowel.clean,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="F2 (Hz)",gl.col = "white",xlab="Vowel",cex.lab=0.4)


ggplot(vowel.clean,aes(x=vowel,y=f2,fill=c_initial)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~v_length) +
  scale_fill_manual(name="Preceding C\n Class",values=c("#ffffff","#bfbfbf")) +
  ylab("F2 (Hz)")

ggplot(vowel.clean,aes(x=vowel,y=f2,fill=c_final)) + geom_boxplot() + theme_bw() +
  facet_grid(~v_length) + 
  scale_fill_manual(name="Following C\n Class",values=c("#ffffff","#bfbfbf")) + 
  ylab("F2 (Hz)")

pirateplot(f2~c_initial+vowel+v_length,data=vowel.clean,theme=2,inf="hdi",hdi.iter=5000,ylab="F2 (Hz)",cex.lab=0.5)
### means for geminate pairs

# plotting changes between vowels for front and back
mean_vowels = read.csv("mean_vowels.csv")

lobanov <- with(vowel.clean, normLobanov(cbind(f1, f2), group=speaker_id))

with(mean_vowels,plotVowels(cbind(f1front.norm,f1.back.norm),cbind(f2.front.norm,f2.back.norm), diph.arrows=T, group=vowel, diph.label.first.only=T,var.col.by = vowel, plot.tokens = T,diph.args.means=F, pch.tokens = vowel,pretty=T,cex.tokens = 1.2))

vowel_normed <- read.csv("vowel_norm.csv")
with(vowel_normed,plotVowels(f1norm,f2norm, group=c_intial,var.col.by = vowel, plot.tokens = T, pch.tokens = vowel,pretty=T,cex.tokens = 1))

### plotting vowels with ggplot

normvowels <- normLobanov(cbind(vowel.clean$f1,vowel.clean$f2),group=vowel.clean$speaker_id)
vowels <- vowel.clean$vowel
v_length <- vowel.clean$v_length
c_initial <- vowel.clean$c_initial
vowel_norm <- cbind(normvowels,vowels,v_length,c_initial)

vowel_norm_mean <- vowel_norm %>%
  group_by(vowels,v_length,c_initial) %>%
  summarise(meanf1 = mean(V1), meanf2 = mean(V2))



ggplot(vowel_norm_mean, aes(x=meanf2,y=meanf1,color=vowels,shape=v_length)) + 
  #geom_point() + 
  scale_y_reverse() + scale_x_reverse() + 
  #stat_ellipse(level=0.95) +
  theme_few() +
  geom_point()+
  scale_shape_manual(values=c(0,1))+
  xlab("F2") + ylab("F1") +
 # facet_grid(~v_length) +
  geom_text(aes(label=vowels),hjust=0, vjust=0)


##### skip for now. code not working
#transform data for vowel plots (combine two colums for vowel and v_length into 1)

vowels2 <- vowel.clean %>% unite(vowel2, vowel, v_length, sep=".")

with(vowel.clean, plotVowels(normLobanov(f1), normLobanov(f2), vowel=vowel,group=c_initial,plot.means=TRUE,pretty=T,var.col.by=vowel,plot.tokens=FALSE,pch.means = vowel))

with(vowel.clean, plotVowels(normLobanov(f1), normLobanov(f2), vowel=vowel, group=word, var.style.by=vowel, plot.means=TRUE, alpha.means=0.5, var.col.by=speaker_id, pretty=TRUE, plot.tokens=FALSE, pch.means=vowel, pch.tokens=word, legend.kwd='bottomright', ellipse.line = FALSE, cex.means=1, ellipse.fill = FALSE, fill.opacity = 0.1))
vowel.clean.short <- vowel.clean[vowel.clean$v_length=="short",]
vowel.clean.long <- vowel.clean[vowel.clean$v_length=="long",]
with(vowel.clean.short, plotVowels(normLobanov(f1), normLobanov(f2), vowel=vowel, group=vowel, plot.means=TRUE, var.col.by=vowel, pretty=TRUE, plot.tokens=FALSE, pch.means=vowel, legend.kwd='bottomright', ellipse.line = TRUE, cex.means=1, ellipse.fill = FALSE, fill.opacity = 0.1))
with(vowel.clean.long, plotVowels(normLobanov(f1), normLobanov(f2), vowel=vowel, group=vowel, plot.means=TRUE, var.col.by=vowel, pretty=TRUE, plot.tokens=FALSE, pch.means=vowel, legend.kwd='bottomright', ellipse.line = TRUE, cex.means=1, ellipse.fill = FALSE, fill.opacity = 0.1))


## heatmaps

force <- with(vowel.clean, repulsiveForce(f2, f1, vowel))
force.col <- plotrix::color.scale(force, cs1 = c(150, 330), cs2 = c(25, 75), cs3 = c(60, 
                                                                                     90), color.spec = "hcl")
coolcolors <- plotrix::color.scale(x = 0:100, cs1 = c(150, 330), cs2 = c(25, 75), 
                                   cs3 = c(60, 90), color.spec = "hcl")
par(mfrow = c(2, 2))

force.col <- plotrix::color.scale(force, cs1 = c(150, 330), cs2 = c(25, 75), cs3 = c(60, 
                                                                                     90), color.spec = "hcl")
with(vowel.clean, plotVowels(f1, f2, vowel, pch.tokens = vowel, cex.tokens = 1.2, col = force.col, 
                         var.col.by = NA, pretty = TRUE))
for (res in c(20, 40, 60)) {
  with(vowel.clean, plotVowels(f1, f2, vowel, plot.tokens = FALSE, heatmap = TRUE, 
                           heatmap.args = list(colormap = coolcolors, resolution = res, useRaster = TRUE), 
                           pretty = TRUE))
} 
                                                                                                    
                                                                                                               
                                                                                                                                                                                        resolution = 10, useRaster = TRUE), heatmap.legend = TRUE, pretty = TRUE))
########################################### End vowel allophony #################################################

########################################### Begin Frication #####################################################
#import csv
frication <- read.csv("t_d_frication.csv")

#plot and check for outliers
#cog (center of gravity)
ggplot(frication, aes(x=consonant, y=cog)) + geom_boxplot() + theme_bw() + facet_grid(~speaker)
#skewness
ggplot(frication, aes(x=consonant, y=skewness)) + geom_boxplot() +theme_bw() + facet_grid(~speaker)

# clean data (removes 3sigma outlier) [removed 0 obs.]
frication.clean <- frication %>%
  group_by(speaker,consonant)%>%
  mutate(skewsdmean = mean(skewness) + 3*sd(skewness),negskewsdmean = mean(skewness)-3*sd(skewness),cogsdmean = mean(cog) + 3*sd(cog),negcogsdmean = mean(cog)-3*sd(cog))%>%
  filter(cog<cogsdmean)%>%
  filter(cog>negcogsdmean) %>%
  filter(skewness<skewsdmean)%>%
  filter(skewness>negskewsdmean)




#lmers
#for cog (full model)
frication.clean.cog.lmer <- lmer(cog ~ consonant + (1+consonant|speaker) + (1+consonant|word), data=frication.clean, REML=FALSE)
summary(frication.clean.cog.lmer)

frication.clean.cog.blmer = stan_lmer(cog ~ consonant + (1+consonant|speaker) + (1+consonant|word), data=frication.clean, 
                                    prior_intercept = normal(0, 100),
                                    prior = normal(0, 50),
                                    prior_covariance = decov(regularization = 2), 
                                    chains = 4,
                                    iter = 2000)

plot(frication.clean.cog.blmer,pars=c("(Intercept)","consonantt","consonantt̻"))

sink("sink-summary-blmer-cog.txt")
print(summary(frication.clean.cog.blmer),digits=4)
sink()



# diagnostic plots
sjp.lmer(frication.clean.cog.lmer, type="rs.ri", sample.n=5)
ranef(frication.clean.cog.lmer)
# check residuals
qqnorm(residuals(frication.clean.cog.lmer)) # normal
plot(fitted(frication.clean.cog.lmer),residuals(frication.clean.cog.lmer)) # fairly homoscedastic, just few-ish data points

#for skewness (full model does not converge so used this one)
frication.clean.skewness.lmer <- lmer(skewness ~ consonant + (1|word) + (1+consonant|speaker), data=frication.clean, REML=FALSE)
summary(frication.clean.skewness.lmer)

frication.clean.skew.blmer = stan_lmer(skewness ~ consonant + (1+consonant|speaker) + (1+consonant|word), data=frication.clean, 
                                      prior_intercept = normal(0, 1),
                                      prior = normal(0, 2),
                                      prior_covariance = decov(regularization = 2), 
                                      chains = 4,
                                      iter = 2000)

blmer_skew_plot<- plot(frication.clean.skew.blmer,pars=c("(Intercept)","consonantt","consonantt̻"))

sink("sink-summary-blmer-skew.txt")
print(summary(frication.clean.skew.blmer),digits=3)

sink()

# check residuals
qqnorm(residuals(frication.clean.cog.lmer)) # normal
plot(fitted(frication.clean.cog.lmer),residuals(frication.clean.cog.lmer)) # fairly homoscedastic

# post hocs 

summary(glht(frication.clean.skewness.lmer, linfct=mcp(consonant = "Tukey")), test=adjusted(type="BH"))
# no laminal-apical distinction for skewness


summary(glht(frication.clean.cog.lmer, linfct=mcp(consonant = "Tukey")), test=adjusted(type="BH"))
# no laminal-apical distinction for cog

## nice plots
#cog
cog_plot<- ggplot(frication.clean,aes(x=consonant,y=cog)) + geom_boxplot() +
  ylab("Center of gravity (Hz)") + theme_bw()



# skewness
skewness_plot<-ggplot(frication.clean,aes(x=consonant,y=skewness)) + geom_boxplot() + theme_bw()


pirateplot(skewness~consonant,data=frication,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="Skewness",gl.col = "white",xlab="Consonant",pal="gray",bean.f.col=c("purple","seagreen4","skyblue"))
pirateplot(cog~consonant,data=frication,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="Center of Gravity (Hz)",gl.col = "white",xlab="Consonant",pal="gray",bean.f.col=c("purple","seagreen4","skyblue"))
