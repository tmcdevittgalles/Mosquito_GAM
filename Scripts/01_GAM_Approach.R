###### Powell Center: Phenological patterns of mosquitoes #######

# Travis McDevitt-Galles
# 07/21/2021
# title: 01_GAM_Approach

# Early data exploration: Fitting GAMS to seasonal mosquito abundance patterns

library(here)
library(dplyr)
library(ggplot2)
library(patchwork)

# loading in the seasonal abundance data
load( here("Data", "Mosquito_Data_Clean.Rda" ) )

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

## trying to get the bayesian gam to work

library(rstanarm)

prac.df <- gam.df %>%  filter( SciName== "Aedes canadensis")


prac.df <- filter(prac.df, Year == "2018")


prac.df <- prac.df %>% group_by(Site, Plot, TrapEvent) %>% 
  summarise(
    Count =  round(sum(Count_adj),  0),
    TrapHours = sum(TrapHours),
    DOY = min(DOY)
  )

TrapHours <- prac.df$TrapHours

test <- stan_gamm4(Count ~ s(DOY) + offset(log(TrapHours)),
                   #random= ~(1|Plot) ,
                   data= prac.df , family = 'poisson',
                   chains=4, iter =2000)

## creating a new data frame
DOY <- seq(range(prac.df$DOY)[1],range(prac.df$DOY)[2], by =1 )

new.data <-  data.frame(DOY = DOY)

new.data$TrapHours <- 24

dfit  <- posterior_epred(test, newdata= new.data, offset = log(24))

matplot(t(dfit/log(24)), type="l")


pred.df <- as.data.frame(predictive_interval(dfit, prob = .9))

pred.df$DOY <- new.data$DOY

colnames(pred.df) <- c("Low","High", "DOY")

ggplot(pred.df)+ geom_ribbon(aes(x=DOY, ymin=Low/24, ymax=High/24),alpha=.5,
                             color= "grey")+ 
  geom_point(data=prac.df,aes(x=DOY, y=Count/TrapHours))

pheno_extract <- function(fit, species, site, year, DOY){
  
  nIter <- nrow(fit) ## number of iterations in the fitted data frame
  
  # creating a data frame to store the key phenometrics
  # 
  pheno.df <- data.frame( SciName = as.factor( rep( species, nIter)),
                          Site = as.factor( rep( site, nIter)),
                          Year = as.factor( rep( year, nIter)),
                          First = as.integer( rep( NA, nIter)),
                          Last = as.integer( rep(NA, nIter)),
                          Duration = as.integer( rep(NA, nIter)),
                          Peak = as.integer( rep(NA, nIter)),
                          Total = as.numeric( rep(NA, nIter))
  )
  
  
  ### looping over each predicted draw to determine key phenometrics
  
  for( i in 1:nIter){
    focus <- fit[i,]
    
    pheno.df$Total[i] <- sum(focus)
    
    foc_cum_sum <- cumsum(focus)
    
    pro_sum <- foc_cum_sum/pheno.df$Total[i]
    
    pheno.df$First[i] <- DOY[min( which( pro_sum > 0.05))]
    pheno.df$Last[i] <- DOY[max ( which( pro_sum < 0.95))]
    pheno.df$Duration[i] <- pheno.df$Last[i]-pheno.df$First[i]
    pheno.df$Peak[i] <-  DOY[which.max(focus)]
    
  }
  
  metric_1 <- c("First", "Last", "Duration", "Peak", "Total")
  
  pheno.df <- tidyr::gather( pheno.df, Pheno, Value, all_of(metric_1))
  
  pheno.df <- pheno.df %>% group_by( SciName, Site, Year, Pheno) %>% 
    summarize(
      mPheno = mean(Value),
      sdPheno = sd(Value)
    ) %>% ungroup()
  
  return(pheno.df)
}



test.df <- pheno_extract(fit = dfit, species = "Aedes canadensis",
                         site= "HARV", year="2018", DOY = DOY)

ggplot(test.df, aes(x=Peak)) + geom_density(adjust = 4)
