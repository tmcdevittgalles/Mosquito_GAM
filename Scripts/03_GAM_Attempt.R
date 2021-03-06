###### Powell Center: Phenological patterns of mosquitoes #######

# Travis McDevitt-Galles
# 07/21/2021
# title: 03_GAM_Attempt

# Early data exploration: Fitting GAMS to seasonal mosquito abundance patterns

library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rstan)
library(rstanarm)

# loading in the seasonal abundance data
load( here("Data", "Mosquito_Data_Clean.Rda" ) )

source( here("Scripts", "02_GAM_Functions.R" ))

### names from the mosquito data frames
names(complete.df)

dim(complete.df) # 619,554 X 21

# Summarize data to start to make decisions on what to explore first

complete.df <- complete.df %>% group_by(SciName, Plot, Year) %>% 
  mutate( Sp_Yr_Plot_Count =  sum(Count))

hist(log10(complete.df$Sp_Yr_Plot_Count))

# Filtering out Species year and plots that have less than 100 total catch
gam.df <- filter( complete.df, Sp_Yr_Plot_Count > 100)

hist(log10(gam.df$Sp_Yr_Plot_Count))

gam.df <- gam.df %>%
  group_by(SciName, Plot, Year) %>% 
  mutate( nSamp =  length( unique(TrapEvent)) )

hist((gam.df$nSamp))

#  Filtering out speices and years with less than 10 unique detections

gam.df <- filter( gam.df , nSamp >=10) ## use to be 10 but lets change it to 5 

dim( gam.df ) # 44807 X 23 # removed over 600,000 observations WHOA

  ### 

length(unique(gam.df$SciName)) ## 63 unique species

ggplot(gam.df , aes(x = SciName, y= log10( Sp_Yr_Plot_Count) ) )+
  geom_boxplot()+ facet_wrap(~Domain, scales="free_x")

## only selecting sites with at least 2 species 

# gam.df <- filter(gam.df , Domain == "D01" |
#                    Domain == "D02" |
#                    Domain == "D03" |
#                    Domain == "D04" |
#                    Domain == "D05" |
#                    Domain == "D06" |
#                    Domain == "D08" |
#                    Domain == "D09" |
#                    Domain == "D11"|
#                    Domain == "D03" )


gam.df <- filter(gam.df ,Domain == "D01" |
                         Domain == "D02" |
                         Domain == "D03" |
                         Domain == "D04" |
                         Domain == "D05" |
                         Domain == "D06" |
                         Domain == "D08" |
                         Domain == "D09" |
                         Domain == "D10"|
                         Domain == "D11" )



gam.df <- gam.df %>% 
  group_by(SciName, Site) %>% 
  mutate( nYear =  length( unique(Year)) )

hist(gam.df$nYear)

gam.df <- filter( gam.df , nYear > 2)

ggplot(gam.df , aes(x = Domain, y= log10( Sp_Yr_Plot_Count), fill= Year ) )+
  geom_boxplot()+ facet_wrap(~SciName, scales="free_x")

gam.df <- filter(gam.df, TotalWeight > 0)

gam.df %>%  filter( SciName== "Aedes vexans") %>% 
  ggplot(aes(x = DOY, y= Count_adj , color= as.factor(Year) ) )+
  geom_point(size=2) + facet_wrap(~Site, scales="free_y")



vex.df <- gam.df %>% filter(   Site == "KONZ" |
                               Site == "TALL" |
                               Site == "UNDE" | 
                               Site == "WOOD" )

at1 <- pheno_gam( vex.df, Species = "Aedes vexans")


 pheno.df <- at1[[1]]

pheno.df %>% ggplot(aes( x=Year,y = mPheno, color=Site ))+
  geom_point() + facet_wrap(~Pheno, scales="free")


at1[[2]]


gam_figure <- append(gam_figure, at1[[2]])

#taxa.df <- pheno.df

taxa.df <-  rbind.data.frame(taxa.df  ,pheno.df)




taxa.df %>% ggplot(aes( x=Year,y = mPheno, color=Site ))+
  geom_point() + facet_wrap(~Pheno, scales="free")

taxa.df %>% ggplot(aes( x= mPheno, fill=Site ))+
  geom_histogram() + facet_wrap(~Pheno, scales="free")


## number of unique combinations 
nrow(unique(taxa.df[,1:3])) ##177 

taxa.df <- unique(taxa.df)

saveRDS( taxa.df ,"pheno.rds")

saveRDS( gam_figure ,"Gam_Figs.rds")





#taxa.df <- readRDS(  here("Data", "pheno.rds" ))


taxa.wide <- select(taxa.df, -"sdPheno")

taxa.wide <- taxa.wide %>% tidyr::pivot_wider(
  names_from = Pheno, values_from= mPheno
)


taxa.wide %>%  ggplot(aes(x =First, y= Last, color=SciName))+
  geom_point()+ stat_smooth(method="glm") + facet_wrap(~SciName, scales="free")


pheno.pca <- FactoMineR::PCA(taxa.wide[,4:7])

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



ggplot(taxa.wide, aes(x= sPCA1, y= sPCA2, color= SciName))+
  geom_point()+stat_ellipse()


ggplot(taxa.wide, aes(x= PCA1, y= PCA2, color= Site))+
  geom_point()+stat_ellipse()


panel.cor <- function( x, y){
  par(usr=c(0,1,0,1))
  txt <- as.character(format(cor(x,y),digits = 2))
  text( 0.5,0.5,txt, cex=3)
  }


pairs(taxa.wide[,4:8], upper.panel =panel.cor)

pheno.pca <- FactoMineR::PCA(taxa.wide[,11:15])

taxa.wide$sPCA1 <- pheno.pca$svd$U[,1]
taxa.wide$sPCA2 <- pheno.pca$svd$U[,2]


ggplot(taxa.wide, aes(x= sPCA1, y= sPCA2, color= SciName))+
  geom_point()+stat_ellipse()
pairs(taxa.wide[,11:15], upper.panel =panel.cor)


library(lme4)


var.first <- lmer( First ~ (1|SciName) + (1|Site) + (1|Year),
                                                     data=taxa.wide )
summary(var.first)
    # Variance 
    # Site:  ~ 24%
    # Year:  ~ 28%
    # Taxa:  ~ 17%
    # Resid: ~ 30%


var.dur <- lmer( Duration ~ (1|SciName) + (1|Site) + (1|Year),
                 data=taxa.wide )
summary(var.dur)
# Variance 
  # Site:  ~ 40%
  # Year:  ~ 20%
  # Taxa:  ~ 17%
  # Resid: ~ 23%

var.peak <- lmer( Peak ~ (1|SciName) + (1|Site) + (1|Year),
  data=taxa.wide )

summary(var.peak) 
# Variance 
# Site:  ~ 0%
# Year:  ~ 0%
# Taxa:  ~ 30%
# Resid: ~ 70%

var.last <- lmer( Last ~ (1|SciName) + (1|Site) + (1|Year),
                  data=taxa.wide )

summary(var.last) 
# Variance 
# Site:  ~ 41%
# Year:  ~ 0%
# Taxa:  ~ 23%
# Resid: ~ 3cv%



var.tot <- lmer( Total ~ (1|SciName) + (1|Site) + (1|Year),
                  data=taxa.wide )

summary(var.tot) 
# Variance 
# Site:  ~ 41%
# Year:  ~ 0%
# Taxa:  ~ 23%
# Resid: ~ 36%

