#this goes a bit further than ISQoL_anal_network
#1- it tests for independence between nodes
#2- it does a network analysis where potential confounding factors have been taken into account

#start similar
rm(list=ls())
library(dplyr)
library(networktools)
library(sjPlot)
library(stringr)
library(qgraph)
library(mice)
library(ggplot2)

load(file="ISQoL.RData") 
df_all1 <- filter(df_all,Dx=="SCZ"  & !is.na(A)) %>%
  select(-c(Dx,STORI_REBUILDING,STORI_GROWTH,WEMWBS_TOT,SERS_TOT, starts_with("ISMI"))) %>% #ISMI_Alien:ISMI_Resist,
  mutate(A=case_when(IS_TOT>=9 ~ 1,
                     IS_TOT<9 ~0))

df_all=df_all1 %>%
  filter_at(vars(starts_with("SQoL")), all_vars(!is.na(.)))
mat_Y_A_W=df_all
SEL=mat_Y_A_W$SQoL18_SEL;ROM=mat_Y_A_W$SQoL18_ROM;PHY=mat_Y_A_W$SQoL18_PHY;
PSY=mat_Y_A_W$SQoL18_PSY;RES=mat_Y_A_W$SQoL18_RES;AUT=mat_Y_A_W$SQoL18_AUT;
FRI=mat_Y_A_W$SQoL18_FRI;FAM=mat_Y_A_W$SQoL18_FAM;

Sym=mat_Y_A_W$IS_Sympt;Dis=mat_Y_A_W$IS_Disease;Trt=mat_Y_A_W$IS_Treatment;

Network_df_IS1=data.frame(
  Dis,Trt,Sym,
  SEL,RES,AUT,FAM,FRI,PHY,PSY,ROM)

#Insight measures are ordinal, but not QoL dimensions
Network= Network_df_IS1 %>% mutate(Trt=Trt*2)  %>% 
  mutate(across(c(Trt,Sym,Dis),#everything(),#
                ~factor(.,order = TRUE))) #TRUE  makes it an ordinal variable


#independence btw nodes
goldbricker(Network, threshold=0.1)


#ISQoL network confound
#our netw anal may be biased by confounders (clinical severity and cognition)

#############################""
#add MEMCHIF_MCI
#####
#cognitive functioning
#####
cog_fun=read.csv("data_enpsy_07072023.csv", sep=";", header = TRUE) %>%
  select(c(StudySubjectID,MEMCHIF_MCI))

colnames(cog_fun)
colSums(is.na(cog_fun))/1056*100
summary(cog_fun)

df_all=df_all %>%
  left_join(cog_fun,by=join_by(Study.Subject.ID==StudySubjectID)) 

#############################""
#mice procedure

pred_mat <- quickpred(df_all, mincor = 0.25)
df_multimp <- mice(df_all, m=5,meth='pmm', seed = 5, predictorMatrix = pred_mat)


df_mice<-mice::complete(df_multimp,"long")
#save(df_mice,file="df_mice_withMEMCHIF_MCI.RData")
#############################""


load("df_mice_withMEMCHIF_MCI.RData")
recalculate=function(outcome,i) {
  df=df_mice %>% filter(.imp==i) %>% select(-c(.id,.imp)) #%>% cbind(task$data %>% select(all_of(starts_with("delta"))))
  #df=task$data
  
  Y=df %>% select(outcome) %>% unlist();
  X=df %>% select(-c(starts_with("SQoL18_"),starts_with("IS_"),A,Study.Subject.ID,CENTRE))
  #X= df %>% select(c(CGI,Education__12_years))
  data=data.frame(Y=Y,X)
  
  if (sub("_.*","",outcome)=='IS') {
  glm_fit=glm(Y ~ ., data=data, family=quasipoisson())
  } else {glm_fit=glm(Y ~ ., data=data)}
  #output=predict(glm_fit,newdata=X)
    
  return(glm_fit$residuals)
}
colnames(df_all)

#redo analysis with this new network
fig2_3=function(dat, big_title) {
  ggplot(dat,aes(x=id,y=value, colour=var, group=var) )+
    geom_ribbon(aes(ymin = lbound, ymax = ubound), 
                alpha=0.2,fill = "lightgrey", colour="lightgrey") +
    geom_line(linewidth=1)+
    geom_point(size=2) +
    scale_colour_discrete(labels=c('Bootstrap mean', 'Sample'))+
    theme(#legend.position = "none",
      legend.text=element_blank(),#element_text(size=15),
      legend.position = "none",
      legend.title = element_blank(),
      axis.text = element_text(size=20),
      #axis.text.x = element_text(angle=45),
      axis.title = element_blank(),#element_text(size = 15),
      strip.text = element_text(size=30, face="bold"),
      panel.background = element_rect(fill="white",colour = "black"),
      panel.grid.major.y  = element_line(colour="grey"),
      strip.background = element_rect(fill = "white",colour = "black"),
      plot.title = element_text(colour = "black",face="bold", size=17)
      #,aspect.ratio = 0.85 
    ) +
    geom_hline(yintercept = 0, colour='black', linetype = 'da')+
    facet_wrap(type ~ ., labeller=labeller(type=c("edge" = big_title))) +
    coord_flip()
}

fig3=list()
dat_btw=list()
dat1=list()
dat2=list()
for (i in 1:5) {
Dis=recalculate(outcome="IS_Disease",i)
Sym=recalculate(outcome="IS_Sympt",i)
Trt=recalculate(outcome="IS_Treatment",i)
RES=recalculate(outcome="SQoL18_RES",i)
SEL=recalculate(outcome="SQoL18_SEL",i)
AUT=recalculate(outcome="SQoL18_AUT",i)
PHY=recalculate(outcome="SQoL18_PHY",i)
PSY=recalculate(outcome="SQoL18_PSY",i)
FRI=recalculate(outcome="SQoL18_FRI",i)
FAM=recalculate(outcome="SQoL18_FAM",i)
ROM=recalculate(outcome="SQoL18_ROM",i)

Network_df_IS1=data.frame(
  Dis,Trt,Sym,
  SEL,RES,AUT,FAM,FRI,PHY,PSY,ROM)
CorMat_IS1=cor_auto(Network_df_IS1, npn.SKEPTIC = TRUE, missing="pairwise")#so that 
groups=list(Insight=c(1,2,3),
            QoL=c(4,5,6,7,8,9,10,11))
# Q_IS1=qgraph(CorMat_IS1, graph = "glasso", sampleSize=nrow(Network_df_IS1),#glasso
#              tuning =0.5, layout = "spring", title = "", details = FALSE,#gamma
#              threshold=FALSE, groups=groups, palette="colorblind")




Network_<-bootnet::estimateNetwork(Network_df_IS1, default = "EBICglasso", tuning=0.5,corMethod=c("npn"),
                                   threshold=FALSE)#,missing="listwise"
#plot(Network_)

# ##########
# #CI of edge-weight and bridge centrality
# ##########
communities=c('1','1','1',
              '2','2','2','2','2','2','2','2'
)
#plot(boot2, "edge", plot = "difference", onlyNonZero = FALSE, order = "id")
boot1 <- bootnet::bootnet(Network_, nBoots =1000, caseN = 10, type = "nonparametric",#nonparametric or case
                          default = "EBICglasso", computeCentrality = TRUE,
                          model="GGM",nCores = 8,
                          statistics = c("edge","bridgeStrength",'bridgeExpectedInfluence')
                          #"Strength","Betweenness","Closeness","ExpectedInfluence")
                          , communities = communities
)

p=plot(boot1, statistics="edge",labels = TRUE, order = "sample", onlyNonZero = FALSE,
       sampleColor = "red", samplelwd = 1,
       meanColor = "blue", meanlwd = 0.5,#bootstrap mean
       bootColor = "black", bootlwd = 1 #CI
)
# p
dat_all=p$data
dat_QoL=p$data %>%
  filter(str_detect(node1, '[:upper:]') & !str_detect(node1, '[:lower:]')) %>%
  filter(str_detect(node2, '[:upper:]') & !str_detect(node2, '[:lower:]')) 

dat_Ins=p$data %>%
  filter(str_detect(node1, '[:lower:]') ) %>%
  filter(str_detect(node2, '[:lower:]'))
dat_btw[[i]]=dat_all %>% anti_join(dat_Ins) %>% anti_join(dat_QoL)

# fig2_3(dat_QoL,"QoL")    
# fig2_3(dat_Ins,"Insight") 
fig3[[i]]=fig2_3(dat_all %>% anti_join(dat_Ins) %>% anti_join(dat_QoL),
                 "Insight--QoL") 


p1=plot(boot1, statistics="bridgeExpectedInfluence",labels = TRUE, order = "sample", onlyNonZero = FALSE,
        sampleColor = "red", samplelwd = 1,
        meanColor = "blue", meanlwd = 0.5,#bootstrap mean
        bootColor = "black", bootlwd = 1 #CI
)

dat1[[i]]=p1$data

p2=plot(boot1, statistics="bridgeStrength",labels = TRUE, order = "sample", onlyNonZero = FALSE,
        sampleColor = "red", samplelwd = 1,
        meanColor = "blue", meanlwd = 0.5,#bootstrap mean
        bootColor = "black", bootlwd = 1 #CI
)

dat2[[i]]=p2$data


}#end of for loop through mice datasets
i=3;fig3[[i]]
save.image(file="anal_confound_all.RData")


#range of edge-weights for sample measures in dat_btw
#range of centrality indices for sample measures in dat1 and dat2
Dis_SEL=list()
Trt_RES=list()
Dis=list()
SEL=list()
Trt=list()
RES=list()
for (i in 1:5) {#loop through 5 imputed datasets
  
  Dis_SEL[[i]]=dat_btw[[i]] %>% 
    filter(var=="sample") %>%
    select(c(id,value, lbound, ubound)) %>%
    filter(id=="Dis--SEL") %>%
    select(-id)
  
  Trt_RES[[i]]=dat_btw[[i]] %>% 
    filter(var=="sample") %>%
    select(c(id,value, lbound, ubound, id)) %>%
    filter(id=="Trt--RES") %>%
    select(-id)
  
  Dis[[i]]=dat2[[i]] %>% 
    filter(var=="sample") %>%
    select(c(id,value, lbound, ubound, id)) %>%
    filter(id=="Dis") %>%
    select(-id)
  
  SEL[[i]]=dat2[[i]] %>% 
    filter(var=="sample") %>%
    select(c(id,value, lbound, ubound, id)) %>%
    filter(id=="SEL") %>%
    select(-id)
  
  Trt[[i]]=dat1[[i]] %>% 
    filter(var=="sample") %>%
    select(c(id,value, lbound, ubound, id)) %>%
    filter(id=="Trt") %>%
    select(-id)
  
  RES[[i]]=dat1[[i]] %>% 
    filter(var=="sample") %>%
    select(c(id,value, lbound, ubound, id)) %>%
    filter(id=="RES")%>%
    select(-id)
  
}
Dis_SEL=do.call(rbind,Dis_SEL)  ; tab_df(Dis_SEL, digits = 3, use.viewer = FALSE)
Trt_RES=do.call(rbind,Trt_RES) ; tab_df(Trt_RES, digits = 3, use.viewer = FALSE)
Dis=do.call(rbind,Dis) ; tab_df(Dis, digits = 3, use.viewer = FALSE)
SEL=do.call(rbind,SEL) ; tab_df(SEL, digits = 3, use.viewer = FALSE)
Trt=do.call(rbind,Trt) ; tab_df(Trt, digits = 3, use.viewer = FALSE)
RES=do.call(rbind,RES)  ; tab_df(RES, digits = 3, use.viewer = FALSE)
