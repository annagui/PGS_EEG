### Statistical analyses for 
# "Polygenic contributions to face-sensitive cortical responses in infants with a family history of autism"
#author: "AnnaGui"
#date: "29/01/2021"

### Load packages
library(plyr)
library(ggplot2)
library(car)
library(dplyr)
library(ggpubr)

### Load files


WORKDIR="path/to/working/directory"
PHENOFILE="path/to/phenotype/csv/file"
IDMATCHFILE="path/to/file/for/matching/genetic/and/EEG/IDs"
PGSFILE_Aut="path/to/txt/file/with/all/thresholded/polygenic/scores"
PGSFILE_xDx="path/to/txt/file/with/all/thresholded/polygenic/scores"


#Phenotype and ID files

setwd(WORKDIR)
dataset<-read.csv(PHENOFILE)

ID_corr<-read.csv(IDMATCHFILE, header=T)
names(ID_corr)<-c("Family_ID","Genetic_ID","Phase_ID","tissue")


#Merge with Autism polygenic score (PGS) file (PT=0.01)
PRS_Aut<-read.table(PGSFILE_Aut, header=T)
PRS_Aut_id<-merge(ID_corr,PRS_Aut,by.x=c('Genetic_ID','Family_ID'), by.y=c('IID','FID'), all.y=T)
PRS_Aut_id$Phase_ID<-as.character(PRS_Aut_id$Phase_ID)

ds_prs_Aut<-merge(dataset, PRS_Aut_id[,c("Phase_ID","tissue","X0.01")], by.x='ID', by.y='Phase_ID',all = F)
ds_prs_Aut<-ds_prs_Aut[complete.cases(ds_prs_Aut$X0.01)&complete.cases(ds_prs_Aut$fn_n290_lat),]


#Merge with xDx PGS (PT=0.5)

PRS_xDx<-read.table(PGSFILE_xDx, header=T)
PRS_xDx_id<-merge(ID_corr,PRS_xDx,by.x=c('Genetic_ID','Family_ID'), by.y=c('IID','FID'), all.y=T)
PRS_xDx_id$Phase_ID<-as.character(PRS_xDx_id$Phase_ID)

ds_prs<-merge(ds_prs_Aut, PRS_xDx_id[,c("Phase_ID","X0.5")], by.x='ID', by.y='Phase_ID',all = F)


###Obtain summary info for N=104 infants with complete data

#*For entire sample:*
ddply(ds_prs, "sex", summarise,  N= length(unique(ID)))

mean(ds_prs_Aut$M_CA_8, na.rm = T)
sd(ds_prs$M_CA_8, na.rm = T)

mean(ds_prs$fn_n290_lat)
sd(ds_prs$fn_n290_lat)

mean(ds_prs$X0.01)
sd(ds_prs$X0.01)

mean(ds_prs$X0.5)
sd(ds_prs$X0.5)


#*By group:*

#rename outcome and recruitment group columns
ds_prs$Group<-ds_prs$HR_Outcome
ds_prs$HR_Outcome<-NULL

ds_prs$FH<-ds_prs$HR_LR
ds_prs$HR_LR<-NULL


ddply(ds_prs, c("sex","Group"), summarise,  N= length(unique(ID)))
ddply(ds_prs, "Group", summarise,  mean(M_CA_8, na.rm = T), sd(M_CA_8, na.rm = T))
ddply(ds_prs, "Group", summarise,  mean(fn_n290_lat), sd(fn_n290_lat))
ddply(ds_prs, "Group", summarise,  mean(X0.01), sd(X0.01))
ddply(ds_prs, "Group", summarise,  mean(X0.5), sd(X0.5))



###ANOVAs: 

##Age by Group

#check assumptions
ds_prs$Group<-as.factor(ds_prs$Group)
leveneTest(ds_prs$M_CA_8~ds_prs$Group)
shapiro.test(ds_prs[ds_prs$Group=='0',]$M_CA_8)
shapiro.test(ds_prs[ds_prs$Group=='1',]$M_CA_8)
shapiro.test(ds_prs[ds_prs$Group=='2',]$M_CA_8)
shapiro.test(ds_prs[ds_prs$Group=='3',]$M_CA_8)


#ANOVA
summary(aov(M_CA_8~Group, data=ds_prs))

#post-hoc
TukeyHSD(aov(M_CA_8~Group, data=ds_prs))


##F-N N290 latency by Group

#assumptions
leveneTest(ds_prs$fn_n290_lat~ds_prs$Group)
shapiro.test(ds_prs[ds_prs$Group=='0',]$fn_n290_lat)
shapiro.test(ds_prs[ds_prs$Group=='1',]$fn_n290_lat)
shapiro.test(ds_prs[ds_prs$Group=='2',]$fn_n290_lat)
shapiro.test(ds_prs[ds_prs$Group=='3',]$fn_n290_lat)

#ANOVA
summary(aov(fn_n290_lat~Group, data=ds_prs))

#post-hoc
TukeyHSD(aov(fn_n290_lat~Group, data=ds_prs))


#check with analyses by Familial History (noFH vs FH) and then by Outcome (FH-TD vs FH-Other vs FH-Aut)
summary(aov(ds_prs$fn_n290_lat~ds_prs$FH))
summary(aov(ds_prs[ds_prs$Group!='0',]$fn_n290_lat~as.factor(ds_prs[ds_prs$Group!='0',]$Group)))


#check with covariates:
summary(aov(fn_n290_lat~Group + sex + Phase + M_CA_8, data=ds_prs))


#Plot F-N N290 latency by Group

ds_prs_plot<-ds_prs[complete.cases(ds_prs$Group),]
ds_prs_plot$Group<-as.factor(ds_prs_plot$Group)
ds_prs_plot$Group<-droplevels(ds_prs_plot$Group)

#jpeg(file="N290_byGroup.jpeg",width = 5, height = 5, units = 'in', res = 300)
N290_byGroup<-ggplot(ds_prs_plot, aes(x=Group, y=fn_n290_lat, colour=Group, fill=Group, shape=Group)) + 
  geom_bar(position =position_dodge(.9), stat = 'summary', fun = 'mean', color="black",alpha = 0.5) +
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width=.2, color="black" ) +
  geom_jitter( aes(x=Group, y=fn_n290_lat, shape=Group), position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  ylab("F-N N290 latency (ms)\n Non-face > Face                                                 Face > Non-face") +
  xlab("") +
  scale_x_discrete(breaks=c("0", "1", "2","3"),labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_shape_manual(values=c("0"=1, "1"=13,"2"= 4,"3"=3), name=" ", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_fill_manual(values=c("0"='springgreen4', "1"='royalblue',"2"= 'darkorange1',"3"='firebrick'), name=" ", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_color_manual(values=c("0"='springgreen4', "1"='royalblue',"2"= 'darkorange1',"3"='firebrick'), name=" ", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust=-1, size=10, hjust=0.5), axis.line = element_line(colour = "black"), text = element_text(size=12),  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(fill=F, color=F, shape=F)  
print(N290_byGroup)
#dev.off()


##Autism PGS (PT=0.01) by Group

#assumptions
leveneTest(ds_prs$X0.01~as.factor(ds_prs$Group))
shapiro.test(ds_prs[ds_prs$Group=='0',]$X0.01)
shapiro.test(ds_prs[ds_prs$Group=='1',]$X0.01)
shapiro.test(ds_prs[ds_prs$Group=='2',]$X0.01)
shapiro.test(ds_prs[ds_prs$Group=='3',]$X0.01)

#ANOVA
summary(aov(ds_prs$X0.01~as.factor(ds_prs$Group)))

#post-hoc
TukeyHSD(aov(ds_prs$X0.01~as.factor(ds_prs$Group)))

#check with analyses by Familial History (noFH vs FH) and then by Outcome (FH-TD vs FH-Other vs FH-Aut)
summary(aov(ds_prs$X0.01~as.factor(ds_prs$FH)))
summary(aov(ds_prs[ds_prs$Group!='0',]$X0.01~as.factor(ds_prs[ds_prs$Group!='0',]$Group)))


#check results with covariates
summary(aov(X0.01~Group + sex + Phase + M_CA_8 + tissue, data=ds_prs))


#Plot Autism PGS by Group

#jpeg(file="ASD_byGroup.jpeg",width = 5, height = 5, units = 'in', res = 300)
ASD_byGroup<-ggplot(ds_prs_plot, aes(x=Group, y=X0.01, colour=Group, fill=Group, shape=Group)) + 
  geom_bar(position =position_dodge(.9), stat = 'summary', fun = 'mean', color="black",alpha = 0.5) +
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width=.2, color="black" ) +
  geom_jitter( aes(x=Group, y=X0.01, shape=Group), position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  ylab("Autism polygenic score (PT=0.01)") +
  xlab("") +
  scale_x_discrete(breaks=c("0", "1", "2","3"),labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_shape_manual(values=c("0"=1, "1"=13,"2"= 4,"3"=3), name=" ", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_fill_manual(values=c("0"='springgreen4', "1"='royalblue',"2"= 'darkorange1',"3"='firebrick'), name=" ", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_color_manual(values=c("0"='springgreen4', "1"='royalblue',"2"= 'darkorange1',"3"='firebrick'), name=" ", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust=-1, size=10, hjust=0.5), axis.line = element_line(colour = "black"), text = element_text(size=12),  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(fill=F, color=F, shape=F)  
print(ASD_byGroup)
#dev.off()



##xDx PGS (PT=0.5) by Group 

#assumptions
leveneTest(ds_prs$X0.5~as.factor(ds_prs$Group))
shapiro.test(ds_prs[ds_prs$Group=='0',]$X0.5)
shapiro.test(ds_prs[ds_prs$Group=='1',]$X0.5)
shapiro.test(ds_prs[ds_prs$Group=='2',]$X0.5)
shapiro.test(ds_prs[ds_prs$Group=='3',]$X0.5)

#ANOVA
summary(aov(ds_prs$X0.5~as.factor(ds_prs$Group)))

#post-hoc
TukeyHSD(aov(ds_prs$X0.5~as.factor(ds_prs$Group)))

#check with analyses by Familial History (noFH vs FH) and then by Outcome (FH-TD vs FH-Other vs FH-Aut)
summary(aov(ds_prs$X0.5~as.factor(ds_prs$FH)))
summary(aov(ds_prs[ds_prs$Group!='0',]$X0.5~as.factor(ds_prs[ds_prs$Group!='0',]$Group)))


#check results with covariates

summary(aov(X0.5~Group + sex + Phase + M_CA_8 + tissue, data=ds_prs))



#Plot xDx PGS by Group

#jpeg(file="xDx_byGroup.jpeg",width = 5, height = 5, units = 'in', res = 300)
xDx_byGroup<-ggplot(ds_prs_plot, aes(x=Group, y=X0.5, colour=Group, fill=Group, shape=Group)) + 
  geom_bar(position =position_dodge(.9), stat = 'summary', fun = 'mean', color="black",alpha = 0.5) +
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width=.2, color="black" ) +
  geom_jitter( aes(x=Group, y=X0.5, shape=Group), position=position_jitter(width=0.3, height=0.2), alpha=0.9) +
  ylab("xDx polygenic score (PT=0.5)") +
  xlab("") +
  scale_x_discrete(breaks=c("0", "1", "2","3"),labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_shape_manual(values=c("0"=1, "1"=13,"2"= 4,"3"=3), name=" ", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_fill_manual(values=c("0"='springgreen4', "1"='royalblue',"2"= 'darkorange1',"3"='firebrick'), name=" ", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_color_manual(values=c("0"='springgreen4', "1"='royalblue',"2"= 'darkorange1',"3"='firebrick'), name=" ", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust=-1, size=10, hjust=0.5), axis.line = element_line(colour = "black"), text = element_text(size=12),  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(fill=F, color=F, shape=F)  
print(xDx_byGroup)
dev.off()


### Regressions

##F-N N290 latency ~ Autism PGS (PT=0.01)

#linear regression with PGS more predictive of Autism
summary(lm(fn_n290_lat~X0.01, data=ds_prs))

#check correlation
cor.test(ds_prs$fn_n290_lat,ds_prs$X0.01)


#check results with covariates

summary(lm(fn_n290_lat~X0.01+ sex + Phase + M_CA_8 + tissue, data=ds_prs))


#Plot N290 ~ Autism PGS

#jpeg(file="N290_ASDpgs.jpeg",width = 7, height = 5, units = 'in', res = 300)
N290_ASDpgs<-ggplot(data=ds_prs_plot, aes(x=X0.01, y=fn_n290_lat)) + 
  ylab("F-N N290 latency (ms)\n Non-face > Face                                                 Face > Non-face") +
  xlab("Autism polygenic score (PT=0.01)") +
  xlim(-2.1,3.75)+
  geom_smooth(method = "lm", se = T, color="black") +
  geom_point(aes(color=as.factor(Group), shape=as.factor(Group))) + 
  scale_color_manual(values=c("0"='springgreen4', "1"='royalblue',"2"= 'darkorange1',"3"='firebrick'), name="", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_shape_manual(values=c("0"=1, "1"=13,"2"= 4,"3"=3), name="", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust=-1, size=11, hjust=0.5), plot.title = element_text(hjust = 0.5,size=12), axis.line = element_line(colour = "black"), text = element_text(size=12),  panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
print(N290_ASDpgs)
#dev.off()

##F-N N290 latency ~ xDx PGS (PT=0.5)

#linear regression with PGS more predictive of atypical development (FH-Other + FH-Aut)
summary(lm(fn_n290_lat~X0.5, data=ds_prs))

#check correlations
cor.test(ds_prs$fn_n290_lat,ds_prs$X0.5)


#check results with covariates
summary(lm(fn_n290_lat~X0.5+ sex + Phase + M_CA_8 + tissue, data=ds_prs))


#Plot N290 ~ xDx PGS

#jpeg(file="N290_xDxpgs.jpeg",width = 7, height = 5, units = 'in', res = 300)
N290_xDxpgs<-ggplot(data=ds_prs_plot, aes(x=X0.5, y=fn_n290_lat)) + 
  ylab("F-N N290 latency (ms)\n Non-face > Face                                                 Face > Non-face") +
  xlab("xDx polygenic score (PT=0.5)") +
  xlim(-2.1,3.75)+
  geom_smooth(method = "lm", se = T, color="black") +
  geom_point(aes(color=as.factor(Group), shape=as.factor(Group))) + 
  scale_color_manual(values=c("0"='springgreen4', "1"='royalblue',"2"= 'darkorange1',"3"='firebrick'), name="", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  scale_shape_manual(values=c("0"=1, "1"=13,"2"= 4,"3"=3), name="", labels=c("noFH", "FH-TD", "FH-Other","FH-Aut")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust=-1, size=11, hjust=0.5), plot.title = element_text(hjust = 0.5,size=12), axis.line = element_line(colour = "black"), text = element_text(size=12),  panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
print(N290_xDxpgs)
#dev.off()

## Figure 1 for paper
#jpeg('Figure1.jpg', units='in', width=12, height = 6, res = 300)
Figure1<-ggarrange(N290_ASDpgs, N290_xDxpgs, labels = c("A", "B"),
                   ncol = 2, nrow = 1,
                   common.legend = TRUE, legend = "bottom")
print(Figure1)
#dev.off()

## Figure 0 combining barplots
jpeg('Figure0.jpg', units='in', width=12, height = 6, res = 300)
Figure0<-ggarrange(ASD_byGroup, xDx_byGroup, labels = c("A", "B"),
                   ncol = 2, nrow = 1,
                  common.legend = TRUE, legend = "bottom")
print(Figure0)
dev.off()
