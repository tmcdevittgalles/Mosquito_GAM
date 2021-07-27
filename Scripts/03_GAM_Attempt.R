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
  mutate( Sp_Yr_Plot_Count =  sum(Count_adj))

hist(log10(complete.df$Sp_Yr_Plot_Count))

# Filtering out Species year and plots that have less than 100 total catch
gam.df <- filter( complete.df, Sp_Yr_Plot_Count > 100)

hist(log10(gam.df$Sp_Yr_Plot_Count))

gam.df <- gam.df %>% filter( Count >0) %>% 
  group_by(SciName, Plot, Year) %>% 
  mutate( nSamp =  length( unique(TrapEvent)) )

hist((gam.df$nSamp))

#  Filtering out speices and years with less than 10 unique detections

gam.df <- filter( gam.df , nSamp >= 10)

dim( gam.df ) # 18,895 X 23 # removed over 600,000 observations WHOA

### 

length(unique(gam.df$SciName)) ## 34 unique species

ggplot(gam.df , aes(x = Domain, y= log10( Sp_Yr_Plot_Count) ) )+
  geom_boxplot()+ facet_wrap(~SciName, scales="free_x")

## only selecting sites with at least 2 species 

gam.df <- filter(gam.df , Domain == "D01" |
                   Domain == "D03" |
                   Domain == "D05" |
                   Domain == "D06" |
                   Domain == "D08" |
                   Domain == "D09" |
                   Domain == "D11")

gam.df <- gam.df %>% 
  group_by(SciName, Site) %>% 
  mutate( nYear =  length( unique(Year)) )

hist(gam.df$nYear)

gam.df <- filter( gam.df , nYear > 2)

ggplot(gam.df , aes(x = Domain, y= log10( Sp_Yr_Plot_Count), fill= Year ) )+
  geom_boxplot()+ facet_wrap(~SciName, scales="free_x")


gam.df %>%  filter( SciName== "Aedes canadensis") %>% 
  ggplot(aes(x = DOY, y= Count_adj/TrapHours , color= Year ) )+
  geom_point(size=2) + facet_wrap(~Plot)



at2 <- pheno_gam( gam.df, Species = "Psorophora columbiae")

at1 %>% ggplot(aes( x= Year, y = mPheno,color=Site ))+
  geom_point() + facet_wrap(~Pheno, scales="free")


taxa.df <-  rbind.data.frame(taxa.df ,at2)

taxa.df <- unique(taxa.df)

saveRDS( taxa.df ,"pheno.rds")
