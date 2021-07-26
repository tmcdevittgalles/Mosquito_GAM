###### Powell Center: Phenological patterns of mosquitoes #######

# Travis McDevitt-Galles
# 07/21/2021
# title: 02_GAM_Functions

# Functions to run to iterate through species, year and sites to quantify
# phenometrics

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
  
  metrics <- c("First", "Last", "Duration", "Peak", "Total")
  
  pheno.df <- tidyr::gather( pheno.df, Pheno, Value, metrics)
  
  pheno.df <- pheno.df %>% group_by( SciName, Site, Year, Pheno) %>% 
              summarize(
                mPheno = mean(Value),
                sdPheno = sd(Value)
              ) %>% ungroup()
  
  return(pheno.df)
}


## iteratively work throught species gam fits


pheno_gam <- function( data, Species) {
  
  foc.df <- filter(data , SciName == Species)

    t <- 1
    
    for( i in 1:length(unqiue(foc.df$Site))){
    
      site.df <- filter(foc.df, Site == unique(foc.df$Site)[i])
      
    
      for( j in 1:length(unique(site.df$Year))){
      
        gam.df <- filter( site.df, Year == unique(site.df$Year)[i])
        
        
        gam.df <- gam.df %>% group_by(SciName, Site, Plot, TrapEvent) %>% 
          summarise(
            Count =  round(sum(Count_adj),  0),
            TrapHours = sum(TrapHours),
            DOY = min(DOY)
          )
        
        TrapHours <- gam.df$TrapHours
        
        bayGam <- stan_gamm4(Count ~ s(DOY) + offset(log(TrapHours)),
                           #random= ~(1|Plot) ,
                           data= gam.df , family = 'poisson',
                           chains=4, iter =2000)
        
        ## creating a new data frame
        DOY <- seq(range(gam.df$DOY)[1],range(gam.df$DOY)[2], by =1 )
        
        new.data <-  data.frame(DOY = DOY)
        
        new.data$TrapHours <- 24
        
        dfit  <- posterior_epred(bayGam, newdata= new.data, offset = log(24))
        
        
        pheno.df <- pheno_extract(fit = dfit,
                                  species = unique(gam.df$SciName)[1],
                                  site = unique(gam.df$Site)[1],
                                  year = unique(gam.df$Year)[1], 
                                  DOY = DOY)
        
        if( t == 1){
          
          out.df <- pheno.df
          
          t <- t + 1
          
        }else(
          
          out.df <- rbind.data.frame( out.df, pheno.df)
        )
      
    }
  }
return(out.df)
}