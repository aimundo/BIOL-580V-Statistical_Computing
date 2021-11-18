library(here)
library(tidyverse)
library(broom)
library(ggforce)
library(scico)
library(viridis)
library(patchwork)
library(svglite) #to save svg figures)

filenames<-list.files(here("Final project/data"),pattern=".csv",full.names=TRUE)
data<-read.csv(here(paste(filenames)))
#remove rows where Cq=NaN
data<-data %>% filter(Cq!='NaN')
data<-data %>%filter(SAMPLE!=57) #to remove sample 57 that seems to be an outlier
data$GENE<-as.factor(data$GENE)
data$GENERAL_GROUP<-as.factor(data$GENERAL_GROUP)


#Calculating average Cq per sample per gene
data1<-data %>%
    filter(GENERAL_GROUP %in% c('CG','MG','MTD')) %>%
    group_by(ID,GENE,WEEK) %>%
    mutate(mean_Cq=mean(Cq),Cq_SD=sd(Cq)) %>%
    filter(duplicated(mean_Cq)==FALSE)

avg_effs<-tibble("GENE"=as.factor(c("VEGFA",
                                    "GAPDH",
                                    "HIF-1a",
                                    "RPTOR",
                                    "DEK",
                                    "ALDOA",
                                    "PGK1",
                                    "STAT3")),
                 "Eff"=c(1.97853947,
                         1.96252309,
                         2.0172954,
                         1.994233037,
                         1.98612065,
                         1.97435611,
                         1.91996324,
                         1.93909892)
)

data1<-data1 %>%
    left_join(avg_effs,data1,by="GENE")


#Sort the data1 dataframe alphabetically so it matches the order of the Av_Eff dataframe
data1<-data1[order(data1$GENE),]

#Calculate expression
data1 <-data1 %>%
    mutate(Expression=Eff^mean_Cq)


Calibrators<-data1 %>%
    group_by(GENE) %>%
    filter(WEEK==1,GENERAL_GROUP=='CG')%>%
    summarize(Baseline=mean(Expression))


#Values in new dataframes
#Data frames per gene
GAPDH<-data1 %>%
    filter(GENE %in% 'GAPDH')

VEGF<- data1 %>%
    filter(GENE %in% 'VEGFA')

HIF_1<- data1 %>%
    filter(GENE %in% 'HIF-1a')

RPTOR<- data1 %>%
    filter(GENE %in% 'RPTOR')

#Pasting the calibrator values in each dataframe

GAPDH<-GAPDH %>%
    mutate(Baseline=as.numeric(rep(Calibrators[3,2])))

VEGF<-VEGF%>%
    mutate(Baseline=as.numeric(rep(Calibrators[9,2])))


HIF_1<-HIF_1 %>%
    mutate(Baseline=as.numeric(rep(Calibrators[4,2])))

RPTOR<-RPTOR %>%
    mutate(Baseline=as.numeric(rep(Calibrators[6,2])))


#Pasting the GAPDH sample expression

VEGF$GAPDH_ex<-GAPDH$Expression
VEGF$GAPDH_Baseline<-GAPDH$Baseline

HIF_1$GAPDH_ex<-GAPDH$Expression
HIF_1$GAPDH_Baseline<-GAPDH$Baseline

RPTOR$GAPDH_ex<-GAPDH$Expression
RPTOR$GAPDH_Baseline<-GAPDH$Baseline

#Calculating Ratios

VEGF<-VEGF%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))

HIF_1<-HIF_1%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))

RPTOR<-RPTOR%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))

Norm<-bind_rows(VEGF,HIF_1,RPTOR)%>%
    group_by(GENE) %>%
    filter(WEEK==1,GENERAL_GROUP=='CG')%>%
    summarize(Norm_Ratio=mean(Ratio))

# Create the "Normalized ratio column by dividing each Ratio of expression by the CG week 1 values

VEGF<-VEGF%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[3])

HIF_1<-HIF_1%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[1])

RPTOR<-RPTOR%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[2])


pd<-position_dodge(0.8)
sz<-2
txt<-20
thm1<-scale_color_viridis_d(option="turbo",end=0.7)
thm1<-scale_color_scico_d(palette="vik",end=0.7)

my_theme = theme(
    axis.title.x = element_text(size = txt),
    axis.text.x = element_text(size = txt),
    axis.title.y = element_text(size = txt))

g1<-VEGF%>%
    group_by(GENERAL_GROUP,WEEK)%>%
    summarize(m=median(Norm_ratio),
              s=sd(Norm_ratio)) %>%
    ggplot(aes(x=WEEK,
               y=m,
               group=GENERAL_GROUP,
               color=GENERAL_GROUP,
    ))+
    geom_errorbar(aes(ymin=m-s,
                      ymax=m+s),
                  size=sz-1,
                  width=0.4,
                  position=pd,
                  show.legend=FALSE)+
    geom_line(size=sz,
              position=pd,
              show.legend=FALSE)+
    geom_point(data=VEGF,
               aes(x=WEEK,y=Ratio,color=GENERAL_GROUP),
               size=sz+2,
               position=pd,
               alpha=0.6,
               show.legend=FALSE)+
    #geom_label(data=VEGF,
    #aes(x=WEEK,y=Ratio,label=SAMPLE)
    # )+
    ggtitle('VEGF')+
    theme_classic()+
    my_theme+
    labs(x='Week',y='Normalized \n relative expression')+
    scale_y_continuous(breaks=seq(0,2,0.5))+
    thm1




g2<-HIF_1%>%
    group_by(GENERAL_GROUP,WEEK)%>%
    summarize(m=median(Norm_ratio),
              s=sd(Norm_ratio)) %>%
    ggplot(aes(x=WEEK,
               y=m,
               group=GENERAL_GROUP,
               color=GENERAL_GROUP,
    ))+
    geom_errorbar(aes(ymin=m-s,
                      ymax=m+s),
                  size=sz-1,
                  width=0.4,
                  position=pd,
                  show.legend=FALSE)+
    geom_line(size=sz,
              position=pd,
              show.legend=FALSE)+
    geom_point(data=HIF_1,
               aes(x=WEEK,y=Ratio,color=GENERAL_GROUP),
               size=sz+2,
               position=pd,
               alpha=0.6,
               show.legend=FALSE)+
    #geom_label(data=VEGF,
    #aes(x=WEEK,y=Ratio,label=SAMPLE)
    # )+
    ggtitle('HIF-1a')+
    theme_classic()+
    my_theme+
    labs(x='Week',y='Normalized \n relative expression')+
    scale_y_continuous(breaks=seq(0,2,0.5))+
    thm1



g3<-RPTOR%>%
    group_by(GENERAL_GROUP,WEEK)%>%
    summarize(m=median(Norm_ratio),
              s=sd(Norm_ratio)) %>%
    ggplot(aes(x=WEEK,
               y=m,
               group=GENERAL_GROUP,
               color=GENERAL_GROUP,
    ))+
    geom_errorbar(aes(ymin=m-s,
                      ymax=m+s),
                  size=sz-1,
                  width=0.4,
                  position=pd,
                  show.legend=FALSE)+
    geom_line(size=sz,
              position=pd)+
    geom_point(data=RPTOR,
               aes(x=WEEK,y=Ratio,color=GENERAL_GROUP),
               size=sz+2,
               position=pd,
               alpha=0.6,)+
    #geom_label(data=VEGF,
    #aes(x=WEEK,y=Ratio,label=SAMPLE)
    # )+
    ggtitle('RPTOR')+
    theme_classic()+
    my_theme+
    labs(x='Week',y='Normalized \n relative expression')+
    scale_y_continuous(breaks=seq(0,2,0.5))+
    thm1

g1+g2+g3
