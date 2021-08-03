###### Powell Center: Phenological patterns of mosquitoes #######

# Travis McDevitt-Galles
# 07/21/2021
# title: 04_GAM_Weather_explore

# Early data exploration: Fitting GAMS to seasonal mosquito abundance patterns

library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rstan)
library(rstanarm)

# loading in the seasonal abundance data
load( here("Data", "Mosquito_Data_Clean.Rda" ) )

load( here("Data", "DailyPrismMod.Rda" ) )
### names from the mosquito data frames
names(complete.df)

dim(complete.df) # 619,554 X 21

# Summarize data to start to make decisions on what to explore first

complete.df <- complete.df %>% group_by(SciName, Plot, Year) %>% 
  mutate( Sp_Yr_Plot_Count =  sum(Count_adj))



taxa.df <- readRDS(  here("Data", "pheno.rds" ))


taxa.wide <- select(taxa.df, -"sdPheno")

taxa.wide <- taxa.wide %>% tidyr::pivot_wider(
  names_from = Pheno, values_from= mPheno
)


taxa.wide %>%  ggplot(aes(x =First, y= Last, color=SciName))+
  geom_point()+ stat_smooth(method="glm")


pheno.pca <- FactoMineR::PCA(taxa.wide[,4:8])

taxa.wide$PCA1 <- pheno.pca$svd$U[,1]
taxa.wide$PCA2 <- pheno.pca$svd$U[,2]


taxa.wide$sDuration <- NA
taxa.wide$sFirst <- NA
taxa.wide$sLast <- NA
taxa.wide$sPeak <- NA
taxa.wide$sTotal <- NA


for( i in 1:length(unique(taxa.wide$Site))){
  
  focSite <- unique(taxa.wide$Site)[i]
  
  taxa.wide$sDuration[taxa.wide$Site== focSite] <- scale(taxa.wide$Duration[taxa.wide$Site== focSite])
  taxa.wide$sFirst[taxa.wide$Site== focSite] <- scale(taxa.wide$First[taxa.wide$Site== focSite])
  taxa.wide$sLast[taxa.wide$Site== focSite] <- scale(taxa.wide$Last[taxa.wide$Site== focSite])
  taxa.wide$sPeak[taxa.wide$Site== focSite] <- scale(taxa.wide$Peak[taxa.wide$Site== focSite])
  taxa.wide$sTotal[taxa.wide$Site== focSite] <- scale(taxa.wide$Total[taxa.wide$Site== focSite])
}





panel.cor <- function( x, y){
  par(usr=c(0,1,0,1))
  txt <- as.character(format(cor(x,y),digits = 2))
  text( 0.5,0.5,txt, cex=3)
}


pairs(taxa.wide[,4:8], upper.panel =panel.cor)

pairs(taxa.wide[,11:15], upper.panel =panel.cor)

pheno.pca <- FactoMineR::PCA(taxa.wide[,11:15])

taxa.wide$sPCA1 <- pheno.pca$svd$U[,1]
taxa.wide$sPCA2 <- pheno.pca$svd$U[,2]


ggplot(taxa.wide, aes(x= PCA1, y= PCA2, color= SciName))+
  geom_point()+stat_ellipse()


ggplot(taxa.wide, aes(x= sPCA1, y= sPCA2, color= SciName))+
  geom_point()+stat_ellipse()


ggplot(taxa.wide, aes(x= PCA1, y= PCA2, color= Site))+
  geom_point()+stat_ellipse()


join.df <- unique(select( ungroup(complete.df), c("Site","Plot" , "Domain")) )

contigus.df <- left_join( contigus.df, join.df, by='Plot')

contigus.df <- contigus.df %>%  filter( Site %in% taxa.df$Site)


head(contigus.df)

contigus.df %>% 
  ggplot( aes( x= Julian, y = Tmean7, color= Year) ) +
  geom_point(size=2,alpha=.5) + facet_wrap(~Site)

contigus.df %>% 
  ggplot( aes( x= Julian, y = PPT14, color= Year) ) +
  geom_point(size=2,alpha=.5) + facet_wrap(~Site)


contigus.df %>% 
  ggplot( aes( x= Julian, y = (CumGDD), color= Year) ) +
  geom_point(size=2,alpha=.5) + facet_wrap(~Site)


weather.14 <- filter( contigus.df, Year == "2014" & Site == "DELA")


gdd.df <- contigus.df %>% group_by(Site, Year) %>% 
  summarise(
      GDDFirst = DOY[min(which( CumGDD >= 100 ))],
      GDDTotal = max( CumGDD )
  )

at1 <- weather.14$DOY[min(which(weather.14$CumGDD >= 50 ))]

z <- 1

contigus.df <- filter( contigus.df, Year != "2013")

for( s in 1:length(unique(contigus.df$Site))){
  focSite <- unique(contigus.df$Site)[s] 
  foc.df <- filter( contigus.df, Site== focSite )
  for( y in 1:length(unique(foc.df$Year)) ){
    focYear <- unique(foc.df$Year)[y]
    year.df <- filter( foc.df, Year == focYear)
  
    GddFirst <- year.df$DOY[min(which(year.df$CumGDD >= 100 ))]
    GddTotal <- max(year.df$CumGDD, na.rm=T)
    
    dum.df <- data.frame( Site = as.character(focSite),
                          Year = as.character(focYear),
                          GddFirst = as.numeric(GddFirst),
                          GddTotal = as.numeric(GddTotal) )
    if( z == 1){
      gdd.df <- dum.df
      z <- z + 1
    }else{
      gdd.df <- rbind.data.frame(gdd.df, dum.df)
    }
    }
}

View(gdd.df)


at1 <- left_join(taxa.wide, gdd.df, by=c("Site","Year"))

ggplot(at1, aes(x=GddFirst, y=First,color=SciName ))+ stat_smooth(method="glm")+
  geom_point(size=2, alpha=.5) + facet_wrap(~SciName, scales="free")+
  theme(legend.position = "none")

ggplot(at1, aes(x=GddFirst, y=Peak, color=SciName ))+ stat_smooth(method="glm")+
  geom_point(size=2, alpha=.5) + facet_wrap(~SciName, scales="free")+
  theme(legend.position = "none")

ggplot(at1, aes(x=GddTotal, y=Last, color=SciName ))+ stat_smooth(method="glm")+
  geom_point(size=2, alpha=.5) + facet_wrap(~SciName, scales="free")+
  theme(legend.position = "none")

ggplot(at1, aes(x=GddFirst, y=Last, color=SciName ))+ stat_smooth(method="glm")+
  geom_point(size=2, alpha=.5) + facet_wrap(~SciName, scales="free")+
  theme(legend.position = "none")

ggplot(at1, aes(x=GddFirst, y=PCA1, color=SciName ))+ stat_smooth(method="glm")+
  geom_point(size=2, alpha=.5) + facet_wrap(~SciName, scales="free")+
  theme(legend.position = "none")

ggplot(at1, aes(x=GddTotal, y=PCA2, color=SciName ))+ stat_smooth(method="glm")+
  geom_point(size=2, alpha=.5) + facet_wrap(~SciName, scales="free")+
  theme(legend.position = "none")


ggplot(at1, aes(x=GddTotal, y=Last, color=SciName ))+ stat_smooth(method="glm")+
  geom_point(size=2, alpha=.5) + facet_wrap(~SciName, scales="free")+
  theme(legend.position = "none")

ggplot(at1, aes(x=GddTotal, y=GddFirst ))+ stat_smooth(method="glm")+
  geom_point(size=2, alpha=.5)+
  theme(legend.position = "none")

ggplot(at1, aes(x=GddTotal, y=Duration ,group=Site))+stat_smooth(method="lm")+
  geom_point(size=2, alpha=.5, color=at1$Year) + facet_wrap(~Site, scales = "free")+
  theme(legend.position = "none")

ggplot(at1, aes(x=GddFirst, y=First ,group=Site))+stat_smooth(method="lm")+
  geom_point(size=2, alpha=.5, color=at1$Year) + facet_wrap(~Site, scales = "free")+
  theme(legend.position = "none")

library(lme4)

first.m <- lmer( First ~ scale(GddFirst)*SciName + (1|Site)+(1|Year),
                 data=at1)

summary(first.m)

duration.m <- lmer( Duration ~ scale(GddTotal)*SciName + (1|Site)+(1|Year),
                 data=at1)

summary(duration.m)
