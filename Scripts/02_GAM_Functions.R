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
    #pheno.df$AUC[i] <-  bayestestR::AUC(x = DOY, y =focus, method = "spline")
    
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


## iteratively work throught species gam fits


pheno_gam <- function( data, Species) {
  
  foc.df <- filter(data , SciName == Species)
  
  plot_gam <- list() # empty list to store gam plots

    t <- 1
    
    for( i in 1:length(unique(foc.df$Site))){
    
      site.df <- filter(foc.df, Site == unique(foc.df$Site)[i])
      
    
      for( j in 1:length(unique(site.df$Year))){
      
        gam.df <- filter( site.df, Year == unique(site.df$Year)[j])
        

        gam.df <- gam.df %>% group_by(SciName, Site, Plot, DOY,Year) %>% 
          summarise(
            Count =  Count,
            TrapHours = sum(TrapHours),
            TotalWeight = sum(TotalWeight),
            SubsetWeight = sum(SubsetWeight)
          ) %>% ungroup()
        
        offset <- (gam.df$TrapHours/24) * gam.df$SubsetWeight/gam.df$TotalWeight
        if( nrow(gam.df) > 10){
          
        bayGam <- stan_gamm4(Count ~ s(DOY, bs="tp" 
                                       
                                       ) + offset(log(offset)),
                           #random= ~(1|Plot) ,
                           data= gam.df , family = 'poisson',
                           chains=4, iter =4000,
                           control = list(adapt_delta = 0.99))
        
        ## creating a new data frame
        DOY <- seq(range(gam.df$DOY)[1],range(gam.df$DOY)[2], by =1 )
        
        new.data <-  data.frame(DOY = DOY)
        
        new.data$offset <- mean(offset)
        
        dfit  <- posterior_epred(bayGam, newdata= new.data, offset = log(offset))
        
        plot.df <- as.data.frame(matrixStats::colMeans2(dfit))
        
        plot.df$DOY <- new.data$DOY
        
        hold <- matrixStats::colQuantiles(dfit, probs=c(0.025, 0.975))
        
        plot.df$Q_2.5 <- hold[,1]
        
        plot.df$Q_97.5 <- hold[,2]
        
        colnames(plot.df)[1] <- "Count"
        

        p1 <- ggplot()+
            geom_line(data=plot.df, aes(x=DOY, y= Count))+
            geom_ribbon( data = plot.df, aes(x=DOY, ymin=Q_2.5, ymax=Q_97.5),
                         alpha=.5, color="grey") +theme_classic()+
            geom_point(data=gam.df, aes(x=DOY,y=Count)) +
            xlab(paste("DOY — " ,unique(gam.df$Year[1]), "—", 
                       unique(gam.df$Site)[1] )) +
            ylab(paste("Count — " ,unique(gam.df$SciName[1])) ) +
            theme(
                   axis.line.x = element_line(color="black") ,
                   axis.ticks.y = element_line(color="black"),
                   axis.ticks.x = element_line(color="black"),
                   axis.title.x = element_text(size = rel(1.8)),
                   axis.text.x  = element_text(vjust=0.5, color = "black"),
                   axis.text.y  = element_text(vjust=0.5,color = "black"),
                   axis.title.y = element_text(size = rel(1.8), angle = 90) ,
                   strip.text.x = element_text(size=20) )
          
  
       plot_gam[[t]] <- p1  # add each plot into plot list
        
        
        pheno.df <- pheno_extract( fit = dfit,
                                  species = as.factor(unique(gam.df$SciName)[1]),
                                  site = as.factor(unique(gam.df$Site)[1]),
                                  year = as.factor(unique(gam.df$Year)[1]), 
                                  DOY = DOY)
        
        if( t == 1 ){
          
          out.df <- pheno.df
          
          
        }else(
          
          out.df <- rbind.data.frame( out.df , pheno.df )
          
        )
        t <- t + 1
        }
    }
    }
    
    
return(list(out.df, plot_gam))
}