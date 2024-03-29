# Clean data and Load usual packages
#first run Decipher_PEC.R script
close.screen(all=TRUE)
rm(list=ls())
header=1;
library(dplyr)

# load data, remove missing data
load(file="ISQoL.RData") 

df_all1 <- filter(df_all,Dx=="SCZ"  & !is.na(A)) %>%
  select(-c(Dx,STORI_REBUILDING,STORI_GROWTH,WEMWBS_TOT,SERS_TOT, starts_with("ISMI"))) %>% #ISMI_Alien:ISMI_Resist,
  mutate(A=case_when(IS_TOT>=9 ~ 1,
                     IS_TOT<9 ~0))

#  Overall correlation between insight and QoL
cor.test(df_all1$IS_TOT,df_all1$SQoL18_TOT)

#  Table 1
#remove missing values of outcome only
df_all=df_all1 %>%
  filter_at(vars(starts_with("SQoL")), all_vars(!is.na(.)))

#Table1_data <- mat_A_W_Table1[complete.cases(mat_A_W_Table1),]
Table1_data <- dplyr::select(df_all,-Study.Subject.ID) %>%
   mutate(across(c(Sex_Male:Addictions_Addictions), as.factor)
          ) %>%
  select(c("Age","Sex_Male","Education__12_years","Employment_UNEMPLOYED","Fam_Not_in_a_relationship", "RQTH_RQTH","Marginalisation_Marginalisation",
           "CGI","GAF",
           "First_Contact_5_to_10_years","First_Contact__10_years","N_Admissions__3",
           "Addictions_Addictions","Dx2_Oui","Dx_SOMA_Oui", "Forensic_Forensic",
           starts_with("SQoL"), starts_with("IS")
           ))

library(table1)
my.render.cont <- function(x) {
    with(stats.default(x), 
         sprintf("%0.2f (%0.1f)", MEAN, SD))
}

pvalue <- function(x, ...) {
    # Construct vectors of data y, and groups (strata) g
    y <- unlist(x)
    g <- factor(rep(1:length(x), times=sapply(x, length)))
    if (is.numeric(y)) {
        # For numeric variables, perform a standard 2-sample t-test
        p <- t.test(y ~ g)$p.value
    } else {
        # For categorical variables, perform a chi-squared test of independence
        p <- chisq.test(table(y, g))$p.value
    }
    # Format the p-value, using an HTML entity for the less-than sign.
    # The initial empty string places the output on the line below the variable label.
    c("", sub("<", "&lt;", format.pval(p, digits=2, eps=0.0001)))
}

table1(~ .
       ,
       data=Table1_data,  
       rowlabelhead = "Variables",
       overall=T,# extra.col=list(`P-value`=pvalue)
       #,render.continuous=my.render.cont
       ) -> Table1_final

#Figure 1 - general network plot

library(qgraph)
library(huge)
library(networktools)
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
  mutate(across(c(Trt,Sym,Dis),
                ~factor(.,order = TRUE))) #TRUE  makes it an ordinal variable

 CorMat_IS1=cor_auto(Network, npn.SKEPTIC = TRUE, detectOrdinal = TRUE,
                       ordinalLevelMax = 10, missing="pairwise")#so that 
 groups=list(Insight=c(1,2,3),
             QoL=c(4,5,6,7,8,9,10,11))#, ISMI=c(12,13,14,15,16))
 Q_IS1=qgraph(CorMat_IS1, graph = "glasso", sampleSize = nrow(Network),#glasso
             tuning =0.5, layout = "spring", title = "", details = FALSE,#gamma
             threshold=FALSE, groups=groups, palette="colorblind", edge.labels=TRUE)
 
 #centralityPlot(Q_IS1)

#can use networktools
 b <- bridge(Q_IS1, communities=c('1','1','1',
                                  '2','2','2','2','2','2','2','2'
                                  #,'3','3','3','3','3'#,'3'
                                  )
             )
 bb=b
 bb$`Bridge Betweenness`=NULL
 bb$`Bridge Closeness`=NULL
plot(bb, zscore = TRUE, color = FALSE, order="given")

#using bootnet for further analysis: fig 2-4 and Supp figures
library(bootnet)
library(stringr)

Network= Network_df_IS1 %>% mutate(Trt=Trt*2) %>% 
                 mutate(across(c(Sym,Dis,Trt),#Trt
                ~factor(.,order = TRUE))) #TRUE  makes it an ordinal variable

communities=c('1','1','1',
              '2','2','2','2','2','2','2','2'
              )

    Network_<-bootnet::estimateNetwork(Network, default = "EBICglasso", tuning=0.5,corMethod=c("cor_auto"), #QoL some variables are continuous, Insight others are ordered factors
                                       threshold=FALSE)#,missing="listwise"
    #detectOrdinal = TRUE, ordinalLevelMax = 10
    plot(Network_,edge.labels=TRUE)
    
    # ##########
    # boot 1 #CI of edge-weight and bridge centrality
    # ##########
    boot1 <- bootnet::bootnet(Network_, nBoots =1000, caseN = 10, type = "nonparametric",#nonparametric or case
                              default = "EBICglasso", computeCentrality = TRUE,
                              model="GGM",nCores = 8,
                   statistics = c("edge","bridgeStrength",'bridgeExpectedInfluence')
                   #"Strength","Betweenness","Closeness","ExpectedInfluence")
                   , communities = communities
                   #statistics = c("Strength","Betweenness","Closeness","ExpectedInfluence")
                   #
                   )
    
#Figure 2-3, within and between-community edges

    p=plot(boot1, statistics="edge",labels = TRUE, order = "sample", onlyNonZero = FALSE,
               sampleColor = "red", samplelwd = 1,
               meanColor = "blue", meanlwd = 0.5,#bootstrap mean
               bootColor = "black", bootlwd = 1 #CI
               )

    dat_all=p$data
    dat_QoL=p$data %>%
      filter(str_detect(node1, '[:upper:]') & !str_detect(node1, '[:lower:]')) %>%
      filter(str_detect(node2, '[:upper:]') & !str_detect(node2, '[:lower:]')) 

    dat_Ins=p$data %>%
      filter(str_detect(node1, '[:lower:]') ) %>%
      filter(str_detect(node2, '[:lower:]'))
    dat_btw=dat_all %>% anti_join(dat_Ins) %>% anti_join(dat_QoL)
    
fig2_3=function(dat, big_title) {
    ggplot(dat,aes(x=id,y=value, colour=var, group=var) )+
        geom_ribbon(aes(ymin = lbound, ymax = ubound), 
                    alpha=0.2,fill = "lightgrey", colour="lightgrey") +
        geom_line(linewidth=2)+
      geom_point(size=3) +
       scale_colour_discrete(labels=c('Bootstrap mean', 'Sample'))+
    #scale_y_discrete(position = "right") +

       theme(#legend.position = "none",
        legend.text=element_blank(),#element_text(size=15),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_text(size=25, face='bold'),
        #axis.text.x = element_text(angle=45),
        axis.ticks.length=unit(.25, "cm"),
        axis.title = element_blank(),#element_text(size = 15),
        strip.text = element_text(size=40, face="bold"),
        panel.background = element_rect(fill="white",colour = "black"),
        panel.grid.major.y  = element_line(colour="grey"),
        strip.background = element_rect(fill = "white",colour = "black"),
        plot.title = element_text(colour = "black",face="bold", size=40)
        #,aspect.ratio = 1.5
    ) +
    geom_hline(yintercept = 0, colour='black', linetype = 'da')+
    facet_wrap(type ~ ., labeller=labeller(type=c("edge" = big_title))) +
      coord_flip()
}
    fig2_3(dat_QoL,"QoL")    
    fig2_3(dat_Ins,"Insight")    
    fig2_3(dat_btw,"Insight--QoL") 
    

#Fig 4 Bridge centrality
#"bridgeStrength",'bridgeExpectedInfluence'
 p1=plot(boot1, statistics="bridgeExpectedInfluence",labels = TRUE, order = "sample", onlyNonZero = FALSE,
               sampleColor = "red", samplelwd = 1,
               meanColor = "blue", meanlwd = 0.5,#bootstrap mean
               bootColor = "black", bootlwd = 1 #CI
               )

    dat1=p1$data
    
p2=plot(boot1, statistics="bridgeStrength",labels = TRUE, order = "sample", onlyNonZero = FALSE,
               sampleColor = "red", samplelwd = 1,
               meanColor = "blue", meanlwd = 0.5,#bootstrap mean
               bootColor = "black", bootlwd = 1 #CI
               )

    dat2=p2$data
    

    fig4=function(dat, c_sentence) {
    ggplot(dat,aes(x=id,y=value, colour=var, group=var) )+
        geom_ribbon(aes(ymin = lbound, ymax = ubound), 
                    alpha=0.2,fill = "lightgrey", colour="lightgrey") +
        geom_line(linewidth=2)+
      geom_point(size=3) +
       scale_colour_discrete(labels=c(' Bootstrap mean', 'Sample'))+
       theme(#legend.position = "none",
        legend.text=element_blank(),#element_text(size=15),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_text(size=25, face='bold'),
        #axis.text.x = element_text(angle=45),
        axis.ticks.length=unit(.25, "cm"),
        axis.title = element_blank(),#element_text(size = 15),
        strip.text = element_text(size=40, face="bold"),
        panel.background = element_rect(fill="white",colour = "black"),
        panel.grid.major.y  = element_line(colour="grey"),
        strip.background = element_rect(fill = "white",colour = "black"),
        plot.title = element_text(colour = "black",face="bold", size=40)
        ,aspect.ratio = 1.5
    )  +
    geom_hline(yintercept = 0, colour='black', linetype = 'da')+
    facet_wrap(type ~ ., labeller=labeller(type=c_sentence)) +
      coord_flip()
    }
fig4(dat2, c("bridgeStrength" = "Strength"))   
fig4(dat1,c("bridgeExpectedInfluence" = "Expected influence"))    

#plot centrality (non bridge centrality measures)
#not in the paper
        
        plot="area"; 
        plot="difference"
        plot(boot1, statistics="ExpectedInfluence",plot=plot,labels = TRUE, order = "sample", onlyNonZero = FALSE,
               sampleColor = "red", samplelwd = 1,
               meanColor = "blue", meanlwd = 0.5,#bootstrap mean
               bootColor = "black", bootlwd = 1 #CI
               )
        plot(boot1, statistics="Strength",plot=plot,labels = TRUE, order = "sample", onlyNonZero = FALSE,
               sampleColor = "red", samplelwd = 1,
               meanColor = "blue", meanlwd = 0.5,#bootstrap mean
               bootColor = "black", bootlwd = 1 #CI
               )
        plot(boot1, statistics="Betweenness",plot=plot,labels = TRUE, order = "sample", onlyNonZero = FALSE,
               sampleColor = "red", samplelwd = 1,
               meanColor = "blue", meanlwd = 0.5,#bootstrap mean
               bootColor = "black", bootlwd = 1 #CI
               )
        plot(boot1, statistics="Closeness",plot=plot,labels = TRUE, order = "sample", onlyNonZero = FALSE,
               sampleColor = "red", samplelwd = 1,
               meanColor = "blue", meanlwd = 0.5,#bootstrap mean
               bootColor = "black", bootlwd = 1 #CI
               )

    ##########
    #stability of centrality indices - see Supp Figure 1 and 3
    ##########
    boot2 <- bootnet::bootnet(Network_, nBoots =1000, caseN = 10, type = "case",#nonparametric or case
                              default = "EBICglasso", computeCentrality = TRUE,
                              model="GGM",nCores = 8,
                   statistics = c("bridgeInDegree","bridgeOutDegree","bridgeStrength","bridgeCloseness","bridgeBetweenness",'bridgeExpectedInfluence',"Strength","Closeness","Betweenness","ExpectedInfluence"), 
                   communities = communities
                   )
plot(boot2, statistics = c("bridgeStrength","bridgeCloseness","bridgeBetweenness",'bridgeExpectedInfluence'), order = "id",
               onlyNonZero = FALSE)+
  theme(#legend.position = "none",
        legend.text=element_text(size=13),
        axis.text = element_text(size=13),
        axis.title = element_text(size = 15),
        strip.text = element_text(size=13, face="bold"),
        panel.background = element_rect(fill = "white",colour = "black"),
        strip.background = element_rect(fill = "white",colour = "black"),
        plot.title = element_text(colour = "black",face="bold", size=17),
        aspect.ratio = 0.85 
  )

plot(boot2, statistics = c("Strength","Closeness","Betweenness","ExpectedInfluence"), order = "id",
               onlyNonZero = FALSE)+
  theme(#legend.position = "none",
        legend.text=element_text(size=13),
        axis.text = element_text(size=13),
        axis.title = element_text(size = 15),
        strip.text = element_text(size=13, face="bold"),
        panel.background = element_rect(fill = "white",colour = "black"),
        strip.background = element_rect(fill = "white",colour = "black"),
        plot.title = element_text(colour = "black",face="bold", size=17),
        aspect.ratio = 0.85 
  )

print(bootnet::corStability(boot2))

    # ##########
    # #testing for significant differences between node centrality
    # ##########
    
    pref_order=c("Dis","Sym","Trt","AUT","FAM","FRI","PHY","PSY","RES","ROM","SEL")
    
    pref_order=c("SEL","ROM","RES","PSY","PHY","FRI","FAM","AUT","Dis","Sym","Trt")

    ##centrality
    y_S=c("Sym","Trt","AUT","FAM","FRI","PHY","PSY","RES","ROM")
    x_S=c("SEL","Dis")
    
    y_I=c("Dis","Sym","AUT","FAM","FRI","PHY","PSY","ROM","SEL")
    x_I=c("RES","Trt")
    DiffT_1=NULL
    DiffT_2=NULL

      for (j in 1:length(y_S)) {
        
      p1_1=bootnet::differenceTest(boot1,x= x_S[1], y=y_S[j], "bridgeStrength")
      p1_2=bootnet::differenceTest(boot1,x= x_S[2], y=y_S[j], "bridgeStrength")
      p1=cbind(p1_1[,c(1,2,4,5)],p1_2[,c(1,2,4,5)])
      
      p2_1=bootnet::differenceTest(boot1,x= x_I[1], y=y_I[j], "bridgeExpectedInfluence")
      p2_2=bootnet::differenceTest(boot1,x= x_I[2], y=y_I[j], "bridgeExpectedInfluence")
      p2=cbind(p2_1[,c(1,2,4,5)],p2_2[,c(1,2,4,5)])


      DiffT_1=rbind(DiffT_1,as.data.frame(p1))
      DiffT_2=rbind(DiffT_2,as.data.frame(p2))

        
      }
      
   
p1=bootnet::differenceTest(boot1,x= "SEL", y="Dis", "bridgeStrength")
p2=bootnet::differenceTest(boot1,x= "RES", y="Trt", "bridgeExpectedInfluence")
DiffT_1=rbind(DiffT_1,as.data.frame(p1[,c(1,2,4,5)]))
DiffT_2=rbind(DiffT_2,as.data.frame(p[,c(1,2,4,5)]))

   library(sjPlot)

##########
# TABLES
##########
      tab_df(DiffT_1)
      tab_df(DiffT_2)

##########
# Supplementary FIGURES 2, 4, 5 (Difference Tests)
##########
new_label_BS <- c("bridgeStrength" = "Bridge strength")
new_label_EI <- c("bridgeExpectedInfluence" = "Bridge expected influence")

#bridge strength - difference plot
p=plot(boot1, statistics="bridgeStrength",plot="difference",order="sample") +
theme(#legend.position = "none",
        legend.text=element_text(size=11),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13, colour = "black"),
        axis.title = element_text(size = 15),
        strip.text = element_text(size=13, face="bold"),
        #panel.background = element_rect(fill = "white",colour = "black"),
        strip.background = element_blank(),#element_rect(fill = "white",colour = "black"),
        #plot.title = element_text(colour = "black",face="bold", size=17)
  ) +
  facet_wrap(type ~ ., labeller=labeller(type=new_label_BS))

#bridge EI - difference plot
p=plot(boot1, statistics="bridgeExpectedInfluence",plot="difference",order="sample"
               ) +
theme(#legend.position = "none",
        legend.text=element_text(size=11),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13, colour = "black"),
        axis.title = element_text(size = 15),
        strip.text = element_text(size=13, face="bold"),
        #panel.background = element_rect(fill = "white",colour = "black"),
        strip.background = element_blank(),#element_rect(fill = "white",colour = "black"),
        #plot.title = element_text(colour = "black",face="bold", size=17)
  ) +
  facet_wrap(type ~ ., labeller=labeller(type=new_label_EI))

SOM_centrality_overall_Network = function(statistics,new_label) {
#bridge strength - difference plot
p=plot(boot1, statistics=statistics,plot="difference",order="sample")+
theme(#legend.position = "none",
        legend.text=element_text(size=11),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13, colour = "black"),
        axis.title = element_text(size = 15),
        strip.text = element_text(size=13, face="bold"),
        #panel.background = element_rect(fill = "white",colour = "black"),
        strip.background = element_blank(),#element_rect(fill = "white",colour = "black"),
        #plot.title = element_text(colour = "black",face="bold", size=17)
  ) +
  facet_wrap(type ~ ., labeller=labeller(type=new_label))
}
p1=SOM_centrality_overall_Network("betweenness",c("betweenness" = "Betweenness"))
p2=SOM_centrality_overall_Network("closeness",c("closeness" = "Closeness"))
p3=SOM_centrality_overall_Network("ExpectedInfluence",c("expectedInfluence" = "Expected Influence"))
p4=SOM_centrality_overall_Network("strength",c("strength" = "Strength"))

