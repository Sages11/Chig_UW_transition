#LUKAS DEFILIPPO 1/28/18
#Script for automated execution of hierarchical logistic model and inference
#load required packages
library(lubridate)
library(viridis)
library(rstan)
library(shinystan)
library(RColorBrewer)
#parallelize
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Useful functions for plotting and summarizing posteriors
#Function to calculate 95% quantiles
cred <- function(x){
  return(quantile(x, c(0.025, 0.975)))
  
  #Function to calculate 50% quantiles
}
cred.50 <- function(x){
  return(quantile(x, c(0.25, 0.75)))
}
#calculate 95% quantiles for columns
cred.2 <- function(x){
  return(apply( x, 2, cred))
}
#calculate 50% quantiles for columns
cred.3 <- function(x){
  return(apply( x, 2, cred.50))
}
#Function to calculate the medians of a column
colMedian <- function(x){
  return(apply(x, 2, median))
}

#Read in daily catch and escapement data
C_E_read.func <- function(file_name='2010-2017ChignikCEbyday_Schindler7-3-18.csv'){
  C_E_dat <- read.csv(file=file_name)
  C_E_dat$Date <-  as.Date(C_E_dat$Date, format='%m/%d/%y')
  C_E_dat$Year <- year(C_E_dat$Date)
  #Convert date to 'day of year' format
  C_E_dat$DOY <- yday(C_E_dat$Date)
  min(C_E_dat$DOY)
  #C_E_dat <- subset(C_E_dat, C_E_dat$DOY < C_E_dat$DOY[C_E_dat$Date=='2010-8-30'] &
  #                    C_E_dat$DOY >  C_E_dat$DOY[C_E_dat$Date=='2010-6-01'])
  
  #Subtract the minimum day of the year so the index begins at 1
  C_E_dat$DOY_index <- C_E_dat$DOY - (min(C_E_dat$DOY)-1)
  #min(C_E_dat$DOY_index)
  #max(C_E_dat$DOY_index)
  
  #Pad each year's escapement data to be the same length
  for (i in 1:length(unique(C_E_dat$Year))){
    #print(length(C_E_dat$Escapement[C_E_dat$Year==unique(C_E_dat$Year)[i]]))
    if(min(C_E_dat$DOY_index[C_E_dat$Year==unique(C_E_dat$Year)[i]]) > 1){
      x <- as.data.frame(matrix(ncol=ncol(C_E_dat), nrow=min(C_E_dat$DOY_index[C_E_dat$Year==unique(C_E_dat$Year)[i]])-1))
      colnames(x) <- colnames(C_E_dat)
      x[,1] <- min(C_E_dat$Date[C_E_dat$Year==unique(C_E_dat$Year)[i]]) - rev(c(1:nrow(x)))
      x[,2] <- 0
      x[,3] <- 0
      x[,4] <- unique(C_E_dat$Year)[i]
      x[,5] <- c(1:nrow(x)) + (min(C_E_dat$DOY)-1)
      x[,6] <- c(1:nrow(x))
      #print(x)
      C_E_dat <- rbind(x, C_E_dat)
      C_E_dat <- C_E_dat[order(C_E_dat$Year),]
    }
    if(max(C_E_dat$DOY_index[C_E_dat$Year==unique(C_E_dat$Year)[i]]) < 135){
      x <- as.data.frame(matrix(ncol=ncol(C_E_dat), nrow= (135- max(C_E_dat$DOY_index[C_E_dat$Year==unique(C_E_dat$Year)[i]]))))
      colnames(x) <- colnames(C_E_dat)
      x[,1] <- max(C_E_dat$Date[C_E_dat$Year==unique(C_E_dat$Year)[i]]) + (c(1:nrow(x)))
      x[,2] <- 0
      x[,3] <- 0
      x[,4] <- unique(C_E_dat$Year)[i]
      x[,5] <- c(1:nrow(x)) + max(C_E_dat$DOY[C_E_dat$Year==unique(C_E_dat$Year)[i]])
      x[,6] <- c(1:nrow(x)) + max(C_E_dat$DOY[C_E_dat$Year==unique(C_E_dat$Year)[i]]) - (min(C_E_dat$DOY)-1)
      #print(x)
      C_E_dat <- rbind(x, C_E_dat)
      C_E_dat <- C_E_dat[order(C_E_dat$Year, C_E_dat$DOY),]
    }
  }
  return(C_E_dat)
}
C_E_dat <- C_E_read.func()

#Read in escapement compositional data --> new genetic samples input via the excel file
#Notes in the excel file template, I have added 3 columns to the original:
#'Mixed' denotes whether a sampling stratum was spread across multiple days
#'Stratum' is the first date of a multi-day sampling event
#'#'Stratum_2' is the last date of a multi-day sampling event
#The two columns will be the same for single day sampling strata
#Experiment with how the model treats new data by adding rows to the spreadsheet in the same format
#as the others
Comp_dat_read.func <- function(file_name='2017 Chignik Estimates Summary.csv'){
  Comp_dat  <- read.csv(file= file_name)
  Comp_dat$Stratum <- as.Date(Comp_dat$Stratum, format= '%m/%d/%y')
  Comp_dat$Stratum_2 <- as.Date(Comp_dat$Stratum_2, format= '%m/%d/%y')
  
  
  C_E_dat <- C_E_read.func()
  #Convert date to 'day of year' format
  Comp_dat$DOY <- yday(Comp_dat$Stratum)
  Comp_dat$DOY_2 <- yday(Comp_dat$Stratum_2)
  
  #Subtract the minimum day of the year from the catch and escapement data (which begins earlier)
  #So that this index begins repective to the first day for the catch and escapement data (day 22 in this case)
  #min(Comp_dat$DOY)
  Comp_dat$DOY_index <- Comp_dat$DOY - (min(C_E_dat$DOY)-1)
  Comp_dat$DOY_index_2 <- Comp_dat$DOY_2 - (min(C_E_dat$DOY)-1)
  for (i in 1:nrow(Comp_dat)){
    if(Comp_dat$DOY_index[i]==Comp_dat$DOY_index_2[i]){
      Comp_dat$DOY_index_3[i] <- Comp_dat$DOY_index[i]
    }
    if(Comp_dat$DOY_index_2[i]==(Comp_dat$DOY_index[i]+1)){
      Comp_dat$DOY_index_3[i] <- Comp_dat$DOY_index[i]
    }
    if(Comp_dat$DOY_index_2[i] > (Comp_dat$DOY_index[i]+1)){
      Comp_dat$DOY_index_3[i] <- round((Comp_dat$DOY_index[i] + Comp_dat$DOY_index_2[i])/2)
    }
  }
  #min(Comp_dat$DOY_index)
  Comp_dat$yr_index <- as.numeric(Comp_dat$Year) - min(as.numeric(Comp_dat$Year)-1)
  
  #Create a total index of sample dates
  for (i in 1:nrow(Comp_dat)){
    Comp_dat$DOY_index_II[i] <- Comp_dat$DOY_index[i] + (Comp_dat$Year-(min(Comp_dat$Year)))[i]*max(C_E_dat$DOY_index)
    Comp_dat$DOY_index_II_2[i] <- Comp_dat$DOY_index_2[i] + (Comp_dat$Year-(min(Comp_dat$Year)))[i]*max(C_E_dat$DOY_index)
    Comp_dat$DOY_index_III[i] <- Comp_dat$DOY_index_3[i] + (Comp_dat$Year-(min(Comp_dat$Year)))[i]*max(C_E_dat$DOY_index)
  }
  
  #Check that indices from the daily catch and escapement data align with those from the compositional data
  #Note: 2012, 2016 are leap years so the doy indices won't be the same for these years
  #print(C_E_dat$Date[C_E_dat$DOY_index==22])
  #print(Comp_dat$Stratum[Comp_dat$DOY_index==22])
  return(Comp_dat)
}
Comp_dat <- Comp_dat_read.func(file_name='2017 Chignik Estimates Summary.csv')

#number of years
n_year <- length(unique(Comp_dat$Year))

#season length 
n_day <- 135

#get the number of genetic samples for each year
length.vec <- c()
for (i in 1:n_year){
  length.vec[i] <- length(Comp_dat$DOY_index[Comp_dat$Year==unique(Comp_dat$Year)[i]])
}

#Function to Execute models in stan
exec_func <- function(script, exec=F, shiny=F, n_chain=4, n_iter=60000, adapt=0.999, tree=10, thin=5, file_name,
                      dat = list(s = length.vec, n_samples = sum(length.vec), comp_dat_x = round(Comp_dat$Proportion_Chignik*Comp_dat$n),comp_dat_index = Comp_dat$DOY_index_3,
                                 n_day = max(C_E_dat$DOY_index), comp_dat_index_2 = Comp_dat$DOY_index_III, comp_dat_N = Comp_dat$n),
                      init.list = list(x0_ind = rep(52,8), steep_ind=rep(0.17,8), mu_steep=0.17, mu_x0=52, k=20, sigma_x0=5, sigma_steep=0.04, steep_Z = rep(0.01,8))){
  
  init = replicate(n_chain, list())
  for(i in 1:n_chain){
    init[[i]] <- init.list
  }
  
  if(exec==T){
    mod <- stan(file = script, data = dat,
                iter = n_iter, chains = n_chain, control=list(adapt_delta=adapt, max_treedepth=tree), thin=thin, seed=666,
                #init=list(x1 <- init.list, x2 <- init.list, x3 <- init.list, x4 <- init.list))
                #init = replicate(n_chain, init.list))
                init = init)
    
    
    saveRDS(mod, paste(file_name, '.rds'))
    #print(mod, pars=c('steep_ind', 'x0_ind'))
    #print(mod, pars=c('mu_steep', 'mu_x0', 'sigma_steep', 'sigma_x0'))
    #print(mod, pars=c('pred_comp'))
    #print(mod, pars=c('ppc'))
  }
  
  #Shinystan diagnostics
  if (shiny==T){
    fit1.shiny <- as.shinystan(mod)  
    launch_shinystan(fit1.shiny)
  }
}

#Run Beta-Binomial model (select exec=T to run model, fit to all years of data - not for in-season management)
exec_func(script='Chignik_hierarchical_transition_beta_bin_A.stan', exec=T, shiny=F, n_chain=5, n_iter=5000, adapt=0.999, tree=10, thin=5,
          dat = list(s = length.vec, n_samples = sum(length.vec), comp_dat_x = round(Comp_dat$Proportion_Chignik*Comp_dat$n),comp_dat_index = Comp_dat$DOY_index_3,
                     n_day = max(C_E_dat$DOY_index), comp_dat_index_2 = Comp_dat$DOY_index_III, comp_dat_N = Comp_dat$n), file_name='Chignik_est_beta_bin',
          init.list = list(x0_ind = rep(52,n_year), steep_ind=rep(0.17,n_year), mu_steep=0.17, mu_x0=52, k=20, sigma_x0=5, sigma_steep=0.04, steep_Z = rep(0.01,n_year)))

#Execute and plot --> select preseason=T if forecasting for a given season before any samples have been collected
#otherwise preseason=F and the function will run the model iteratively adding each sample for the most recent year
#and plot. This function is meant to be the guts of what a manager would run
plot_func <- function(exec=T, iter=1000, preseason=F){
  #Read in genetic sampling data
  Comp_dat <- Comp_dat_read.func(file_name='2017 Chignik Estimates Summary.csv')
  #Read in escapement data
  C_E_dat <- C_E_read.func()
  #Total number of years evaluated
  n_year_fixed <- length(unique(Comp_dat$Year))
  #Number of days within a season
  n_day = max(C_E_dat$DOY_index)
  #Empty list for storing the pre-season curves
  hyper_list <- replicate(n_year_fixed, list())
  #Empty list for storing the pre-season midpoints (x0) terms
  x0_list <- replicate(n_year_fixed, list())
  #Empty list for storing the iterative curves
  pred_list <- replicate(n_year_fixed, list())
  #Empty lists for storing the assigned escapement to the early and late runs under the iterative curve
  in_esc_mat_late <- replicate(n_year_fixed,list())
  in_esc_mat_early <- replicate(n_year_fixed,list())
  #Empty lists for storing the assignment errors to the early and late runs under the iterative curve
  in_esc_error_late <- replicate(n_year_fixed,list())
  in_esc_error_early <- replicate(n_year_fixed,list())
  #Colors for plotting --> unique color for each GSI sample
  col.vec <- magma(table(Comp_dat$Year)[n_year_fixed], alpha = 1, begin = 0, end = 1, direction = -1)
  #Read in escapement goals (based on 2018 values) 
  #esc_goals <- read.csv(file='Esc_goals.csv')
  #Reformat dates to be in DOY format
  #esc_goals$Date <-  as.Date(esc_goals$Date, format='%m/%d/%y')
  #esc_goals$DOY <- yday(esc_goals$Date)
  #Remove last row
  #esc_goals <- esc_goals[1:(nrow(esc_goals)-1),]
  #Create DOY index compatible with the genetic samples
  #esc_goals$DOY_index <- esc_goals$DOY - (min(C_E_dat$DOY)-1)
  #Add entry for day 135 for plotting purposes
  #esc_goals$DOY_index[length(esc_goals$DOY_index)] <- n_day
  #Standardize the recorded escapement data to a fixed season length (n_day), filling in zeros for any given year's missing data
  for (i in 1:length(unique(C_E_dat$Year))){
    C_E_dat$Escapement_cumulative_late[C_E_dat$Year==unique(C_E_dat$Year)[i]][1] <- C_E_dat$Escapement[C_E_dat$Year==unique(C_E_dat$Year)[i]][1]
    for(j in 2:length(C_E_dat$Escapement[C_E_dat$Year==unique(C_E_dat$Year)[i]])){
      C_E_dat$Escapement_cumulative_late[C_E_dat$Year==unique(C_E_dat$Year)[i]][j] <- C_E_dat$Escapement[C_E_dat$Year==unique(C_E_dat$Year)[i]][j] + C_E_dat$Escapement_cumulative_late[C_E_dat$Year==unique(C_E_dat$Year)[i]][j-1]
    }
  }
  #Create plot showing the evolution of the transition curves
  pdf(file='loo_iter_curve_in_season.pdf')
  #Create lists (curve, cumulative assignments, and assignment errors) for each year with dimensions equal to the total number of samples for that year
  #plus 1 (to include the early period informed by the hyper curve, only for the assignment based lists)
  for(i in 1:n_year_fixed){
    pred_list[[i]] <- replicate(as.vector(table(Comp_dat$Year))[i], list())
    in_esc_mat_early[[i]] <- replicate(as.vector(table(Comp_dat$Year))[i]+1, list())
    in_esc_mat_late[[i]] <- replicate(as.vector(table(Comp_dat$Year))[i]+1, list())
    in_esc_error_early[[i]] <- replicate(as.vector(table(Comp_dat$Year))[i]+1, list())
    in_esc_error_late[[i]] <- replicate(as.vector(table(Comp_dat$Year))[i]+1, list())
  }
  
  #Open loop to iteratively exclude one year at a time from the compositional data --> the hyper-curves
  #from models fitted to these data wil be used as the preseason curve for each year before genetic samples
  #are collected
    y <- n_year_fixed
#exclude given year from data
    if(preseason==F){
     Comp_dat_2 <- Comp_dat[-which(Comp_dat$Year==unique(Comp_dat$Year)[y]),]
    }else{
    Comp_dat_2 <- Comp_dat
    }
    #redefine number of years
    n_year <- length(unique(Comp_dat_2$Year))
    #re-calculate the sample length vector
    length.vec <- as.vector(table(Comp_dat_2$Year))
    #Redefine cumulative indices based on years that were dropped
    for (i in 1:nrow(Comp_dat_2)){
      Comp_dat_2$DOY_index_II[i] <- Comp_dat_2$DOY_index[i] + (as.numeric(as.factor(Comp_dat_2$Year))-1)[i]*max(C_E_dat$DOY_index)
      Comp_dat_2$DOY_index_II_2[i] <- Comp_dat_2$DOY_index_2[i] + (as.numeric(as.factor(Comp_dat_2$Year))-1)[i]*max(C_E_dat$DOY_index)
      Comp_dat_2$DOY_index_III[i] <- Comp_dat_2$DOY_index_3[i] + (as.numeric(as.factor(Comp_dat_2$Year))-1)[i]*max(C_E_dat$DOY_index)
    }
    #Run the model if specificed
    if(exec==T){
      exec_func(script='Chignik_hierarchical_transition_beta_bin_A.stan', exec=T, shiny=F, n_chain=5, n_iter=iter, adapt=0.999, tree=10, thin=5,
                dat = list(s = length.vec, n_samples = sum(length.vec), comp_dat_x = round(Comp_dat_2$Proportion_Chignik*Comp_dat_2$n),comp_dat_index = Comp_dat_2$DOY_index_3,
                           n_day = max(C_E_dat$DOY_index), comp_dat_index_2 = Comp_dat_2$DOY_index_III, comp_dat_N = Comp_dat_2$n), file_name=paste('Chignik_beta_bin_year=',unique(Comp_dat$Year)[y]),
                init.list = list(x0_ind = rep(52,n_year), steep_ind=rep(0.17,n_year), mu_steep=0.17, mu_x0=52, k=20, sigma_x0=5, sigma_steep=0.04, steep_Z = rep(0.01, n_year)))
    }
    #otherwise load in preexisting RDA object
    df <- as.data.frame(readRDS(paste('Chignik_beta_bin_year=',unique(Comp_dat$Year)[y],".rds")))
    
    #Save the hyper-curve from the model fit to a list  
    hyper_list[[y]]  <- df[grep('hyper_curve', colnames(df))]
    #Save the estimated midpoint term from the model fit to a list  
    x0_list[[y]] <- df[grep('mu_x0', colnames(df))]
   if(preseason==F){
    #Loop through the number of samples for each year
    for(j in 1:length(which(Comp_dat$Year==unique(Comp_dat$Year)[y]))){
      #Incrementally add the rows containing each sample into the data frame, re-estimating the model with
      #the new data
      Comp_dat_3 <- rbind(Comp_dat_2, Comp_dat[which(Comp_dat$Year==unique(Comp_dat$Year)[y])[1:j],])
      #re-define number of years based on new data
      n_year <- length(unique(Comp_dat_3$Year))
      #redefine sample length vector based on new data
      length.vec <- as.vector(table(Comp_dat_3$Year))
      #Redefine cumulative indices based on years that were dropped
      for (i in 1:nrow(Comp_dat_3)){
        Comp_dat_3$DOY_index_II[i] <- Comp_dat_3$DOY_index[i] + (as.numeric(as.factor(Comp_dat_3$Year))-1)[i]*max(C_E_dat$DOY_index)
        Comp_dat_3$DOY_index_II_2[i] <- Comp_dat_3$DOY_index_2[i] + (as.numeric(as.factor(Comp_dat_3$Year))-1)[i]*max(C_E_dat$DOY_index)
        Comp_dat_3$DOY_index_III[i] <- Comp_dat_3$DOY_index_3[i] + (as.numeric(as.factor(Comp_dat_3$Year))-1)[i]*max(C_E_dat$DOY_index)
      }
      #Sort the new data frame by year
      Comp_dat_3 <- Comp_dat_3[order(Comp_dat_3$Year),]
      
      #Run the model if specified
      if(exec==T){
        exec_func(script='Chignik_hierarchical_transition_beta_bin_A.stan', exec=T, shiny=F, n_chain=5, n_iter=iter, adapt=0.999, tree=10, thin=5,
                  dat = list(s = length.vec, n_samples = sum(length.vec), comp_dat_x = round(Comp_dat_3$Proportion_Chignik*Comp_dat_3$n),comp_dat_index = Comp_dat_3$DOY_index_3,
                             n_day = max(C_E_dat$DOY_index), comp_dat_index_2 = Comp_dat_3$DOY_index_III, comp_dat_N = Comp_dat_3$n), file_name=paste('Chignik_beta_bin_year=',unique(Comp_dat$Year)[y],'sample',j),
                  init.list = list(x0_ind = rep(52,n_year), steep_ind=rep(0.17,n_year), mu_steep=0.17, mu_x0=52, k=20, sigma_x0=5, sigma_steep=0.04, steep_Z = rep(0.01, n_year)))
      }
      #otherwise load in the saved model object
      df <- as.data.frame(readRDS(paste('Chignik_beta_bin_year=',unique(Comp_dat$Year)[y],'sample',j,".rds")))
      pred <- df[grep('pred_comp', colnames(df))]
      
      #Slice up the predictions according to season length and year
      cut_start <- seq(from=1, to=ncol(pred), by=n_day) 
      cut_end <- seq(from=n_day, to=ncol(pred), by=n_day)
      
      #Save the selected curve to a list
      pred_list[[y]][[j]]  <- pred[,c(cut_start[y]:cut_end[y])]
    }
   }
    #Plot the hyper-curves
    par(oma=c(5,5,5,5))
    plot(c(1:n_day), colMedian(hyper_list[[y]]), type='l', lty=3, xaxt='n', yaxt='n', xlim=c(20,90), lwd=2, axes=F, ylab=NA, xlab=NA)
    #vertical line for the midpoint term  
    abline(v=colMedian(x0_list[[y]]), lty=3)
    #Polygon for shading the background
    #polygon(c(15:95, rev(15:95)), c(rep(-0.1, length(15:95)), rep(1.1, length(15:95))), col=rgb(1,1,1, max=255, alpha=0))
    if(preseason==F){
    #Plot the iterative fits --> loop through the total number of samples for a given year
    for(j in 1:length(which(Comp_dat$Year==unique(Comp_dat$Year)[y]))){
      #Allow plotting the final curve a different style from all previous curves
      if(j < length(which(Comp_dat$Year==unique(Comp_dat$Year)[y]))){
        points(colMedian(pred_list[[y]][[j]]),type='l', lty=1, col=rgb(col2rgb(col.vec)[1,j], col2rgb(col.vec)[2,j], col2rgb(col.vec)[3,j], max=255, alpha=255), lwd=1.75)
      }
      if(j==length(which(Comp_dat$Year==unique(Comp_dat$Year)[y]))){
        points(colMedian(pred_list[[y]][[j]]),type='l', lty=1, col=col.vec[j], lwd=2.5)
      }
     }
    
    #Plot the actual data points and the confidence bounds
    points(Comp_dat$DOY_index_3[Comp_dat$Year==unique(Comp_dat$Year)[y]], Comp_dat$Proportion_Chignik[Comp_dat$Year==unique(Comp_dat$Year)[y]], col='black', pch=21, bg=col.vec[1:length(Comp_dat$DOY_index_3[Comp_dat$Year==unique(Comp_dat$Year)[y]])], cex=1.75, lwd=0.25)
    arrows(Comp_dat$DOY_index_3[Comp_dat$Year==unique(Comp_dat$Year)[y]], Comp_dat$Lower_Chignik[Comp_dat$Year==unique(Comp_dat$Year)[y]], Comp_dat$DOY_index_3[Comp_dat$Year==unique(Comp_dat$Year)[y]], Comp_dat$Upper_Chignik[Comp_dat$Year==unique(Comp_dat$Year)[y]], angle=90, code=3, length=0.01, lwd=1.25, lty=1, col=col.vec[1:length(Comp_dat$DOY_index_3[Comp_dat$Year==unique(Comp_dat$Year)[y]])])
    
    legend(x=17.5, y=1.065, legend= c('Preseason', paste('Sample', 1:length(col.vec))),col=c('black', col.vec), pch=c(NA, rep(16, length(col.vec))), bty='n', cex=0.95, lty=c(3, rep(1, length(col.vec))), y.intersp=0.875, x.intersp = 0.75, pt.cex=1.25)
    legend(x=17.5, y=1.065, legend= c('Preseason', paste('Sample', 1:length(col.vec))),col='black', pch=c(NA, rep(1, length(col.vec))), bty='n', lwd=0.25, lty=rep(NA, length(col.vec)), cex=0.95, y.intersp=0.875, x.intersp = 0.75, pt.cex=1.25)
    text(x=85, y=0.01, unique(Comp_dat$Year)[y])
    }else{
      text(x=85, y=0.01, unique(Comp_dat$Year)[y]+1)
      
    }
    #add axes, legends and text
      axis(side=1, at=seq(15, 85, 10), labels= (seq(15, 85, 10) +  min(C_E_dat$DOY)) , cex.axis=1)
      axis(side=2, cex.axis=1)
      mtext(side=2, 'Proportion late run (Chignik)', line=2.5)
      mtext(side=1, 'Day of year', line=2.5)
      abline(h=0.5, lty=2)
 dev.off()
 return(pred_list)
}
x <- plot_func(exec=T, iter=1000, preseason=F)

#Isolate the median estimates of the current year's transition and write to a csv
mat_store <- matrix(nrow=n_day, ncol=length.vec[length(length.vec)])
for(i in 1:length.vec[length(length.vec)]){
  mat_store[,i] <- unlist(colMedian(x[[n_year]][[i]]))
}
write.csv(mat_store, file='current_year_curves.csv')
#plot(mat_store[,1])
#points(mat_store[,2])
#points(mat_store[,3])
#points(mat_store[,4])
#points(mat_store[,5])


