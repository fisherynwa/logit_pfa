################################################################################
#
#   Filename    :	case_study.R    												  
#                                                                                     
#   Project     :       BiomJ article "Multple two-sample testing under
#                                      arbitrary covariance dependency 
#                                      with an application in 
#                                      imaging mass spectrometry"        
#
#   Authors     :       V. Vutov and T. Dickhaus                                                                
#   Date        :       25.04.2022
#   Purpose     :       Approximating the FDP(t) and
#                       identifying the most assoc. m/z for cancer association
#
#		Source data : https://gitlab.informatik.uni-bremen.de/digipath/Supervised_NMF_Methods_for_MALDI.git. -
#                 Matlab file: (L1-8 tic-redc-spdn-adrs.mat)		
#
#   For an in-depth description of this dataset, see - https://pubmed.ncbi.nlm.nih.gov/30395171/
#   
#   Data Objects -  
#
#    data.tic: data.tic contains all (m/z) covariates. This frame has been standardized by tic and binned at 0.4 DA.
#    subtypes: (the target variable) A vector that contains information about the cancerous status per each spectrum.
#              The status '1' (y_i := 1) corresponds to Sqcc, otherwise ADC (i.e. y_i := 0).
#
#    mz_vector: This object maps the columns with the m/z values.
#   
#   R Version   :       R-4.1.2                                                                
#   
#
#   Input data files  :    task_ad_sq.RData                                                           
#   Output data files :    Table8.tex, Table9.tex, Fig3.pdf, Fig2.pdf, Fig3.pdf, Fig4.pdf
#
#   Required R packages :  ggplot2, lqr, gridExtra, dplyr, xtable, quantreg
#
#
################################################################################
  
  rm(list = ls())
  
  
###############################################################################
### Check for missing packages and install them
###############################################################################

  list.of.packages <- c("ggplot2", "lqr", 'gridExtra')
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  
  if(length(new.packages)) install.packages(new.packages)
  
###########################################
# Loading R packages, MALDI data & custom functions
###########################################
  
  library(xtable)
  library(dplyr)
  library(quantreg)
  library(ggplot2)
  require(gridExtra)
  
  load('data/task_ad_sq.RData')
  source('codes/MainFunctions.R')
  source('codes/pfa_custom.R')
  
  
######################################################
## Executing multiple marginal models for each j 
######################################################
  
  z_values = estimating_psi(X = data.tic, y = subtypes)
  
  covariance = marginal_covarinace_matrix(psi = z_values$Psi,
                                          n = z_values$num_obs, 
                                          p = z_values$num_feat)
  
######################################################
## Efron Empirical Correction, for more computation details, see
## 'source('codes/Ahat-Efron.R')' and the mentioned references in the paper
# We have saved the genuine Z-values - save(genuine_z_values, file='./results/z_values.rda')
######################################################
  
  Std <- sqrt(diag(covariance))
  
  estimates = z_values$betas
  
  genuine_z_values <- estimates/Std 
  
  Ahat(x0 = 1.0, z = genuine_z_values)
  
  zval_correction = genuine_z_values/2.63
  
######################################################
## Running PFA, a chunk of code taken from pfa::pfa.test()
######################################################
  
  pfa_k_10 <-  pfa(Z = zval_correction, Sigma = covariance,
                   reg = 'L1', K = 10, Kmax = 1699,  t=exp(-seq(1, 7.5, 0.2)))
  
  pfa_k_9 <-  pfa(Z = zval_correction, Sigma = covariance,
                  reg = 'L1', K = 9, Kmax = 1699,  t=exp(-seq(1, 7.5, 0.2)))
  
  pfa_k_8 <-  pfa(Z = zval_correction, Sigma = covariance,
                  reg = 'L1', K = 8, Kmax = 1699,  t=exp(-seq(1, 7.5, 0.2)))
  
  pfa_k_7 <-  pfa(Z = zval_correction, Sigma = covariance,
                  reg = 'L1', K = 7, Kmax = 1699,  t=exp(-seq(1, 7.5, 0.2)))
  
  pfa_k_6 <-  pfa(Z = zval_correction, Sigma = covariance,
                  reg = 'L1', K = 6, Kmax = 1699,  t=exp(-seq(1, 7.5, 0.2)))
  
  
######################################################
#### Number of R(t) and \widehat{FDP(t)} over a grid of plausible t - Table 8
######################################################  

  table_8 <-  pfa_k_6$FDP[c(20, 24, 26, 29, 31, 33), ]
  
  table_8$false.rejects <- NULL
  table_8$t <- format(table_8$t, scientific = TRUE, digits = 3) 
  
  table_8$FDP <- round(table_8$FDP, 4)
  
  print(xtable(table_8, type = "latex", digits = c(0, 3, 0, 4)), file = "./results/Table8.tex",
        include.rownames = FALSE, floating.environment = TRUE)
  
######################################################
#### Reporting Stat. Significant m/z Values - Table 9
######################################################
  
  ass_mz = cbind('mz_val' = mz_vector, 
                       'z_val' = round(zval_correction, 2)) %>% as.data.frame()

  
  pfa_z_SqCC <- sort(ass_mz$z_val, index.return = TRUE)
  
  pfa_z_ADC  <- sort(ass_mz$z_val, decreasing = TRUE ,index.return = TRUE)
  
  sqcc_table_9 <- ass_mz[c(pfa_z_SqCC$ix[1:20]),]
  
  adc_table_9 <- ass_mz[c(pfa_z_ADC$ix[1:20]),]
    
  table_9 <- cbind(sqcc_table_9, adc_table_9)
  
  colnames(table_9) <- c(rep(c("m/z values", "Z-scores"), 2))
   
  print(xtable(table_9, type = "latex"), file = "./results/Table9.tex",
        include.rownames = FALSE)
  
############################################
## Figure 1 - two examples of mass spectra ##
############################################
  
  pdf(file = 'results/Fig1.pdf')
  
  par(mfrow = c(1, 2))
  plot(x  = mz_vector, y = X[2050,], type = 'l', col = 'blue', 
       ylab = 'Intensity', xlab = 'm/z values')
  
  plot(x  = mz_vector, y = X[520,], type = 'l', col = 'blue', 
      ylab = 'Intensity', xlab = 'm/z values')
  
  dev.off()
  
##################################################
## Figure 2 - The Empirical z-values and p-values ##
##################################################
  
  P <- 2 * (1 - pnorm(abs(genuine_z_values)))
  df_z = data.frame(PF = genuine_z_values)
  
  df <- data.frame(PF = P)
  
  plot1 = ggplot(df_z, aes(x = PF)) + 
    geom_histogram(aes(y =..density..),
                   breaks = seq(-32, 32, by = 2), 
                   colour = "black", 
                   fill = "white") + ylab('Density') + xlab('Z-Values') +
    stat_function(fun = dnorm, args = list(mean = mean(df_z$PF), sd = sd(df_z$PF)), color = "red") 
  
  plot2 = ggplot(df, aes(x = PF)) + 
    geom_histogram(aes(y = ..density..),
                   breaks = seq(0, 1, by = 0.08), 
                   colour = "black", 
                   fill = "white") + xlab('P-values') + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab('Frequency') 
  
  
  pdf(file = './results/Fig2.pdf')
  
  
  grid.arrange(plot1, plot2, ncol=2)
  
  
  dev.off()

  
  
##################################################
## Figure 3 - The approximated number of V(t) and FDP(t) given different k
##################################################

  logit_pfa_6 = pfa_k_6$FDP
  logit_pfa_7 =  pfa_k_7$FDP
  logit_pfa_8 = pfa_k_8$FDP
  logit_pfa_9 = pfa_k_9$FDP
  logit_pfa_10 = pfa_k_10$FDP

  pdf(file = './results/Fig3.pdf')
  
  
  par(mfrow = c(1, 2))
  
  matplot(x = logit_pfa_6[,2][21:33], y = logit_pfa_6[,3][21:33], 
          ylim = range(c(10, 150)), xlim = range(c(250, 400)), type = c("b"),
          pch=1, col = 'black', xlab = '# of total rejections',
          ylab = 'Approximated # of False Rejections') #plot
  
  par(new=TRUE)
  
  lines(x = logit_pfa_7[,2][21:33], y = logit_pfa_7[,3][21:33], type = c("b"), pch= 2, col = 'green')
  lines(x = logit_pfa_8[,2][21:33], y = logit_pfa_8[,3][21:33], type = c("b"), pch= 3, col = 'blue')
  lines(x = logit_pfa_9[,2][21:33], y = logit_pfa_9[,3][21:33], type = c("b"), pch= 4, col = 'orange')
  lines(x = logit_pfa_10[,2][21:33], y = logit_pfa_10[,3][21:33], type = c("b"), pch= 5, col = 'red')
  
  legend("topleft", legend=c("k = 6", "k = 7", 'k = 8', 'k = 9', ' k = 10'), 
         col=c("black", "green", 'blue', 'orange', 'red'), 
         pch = c(1, 2, 3, 4, 5),
         bty = "n", 
         pt.cex = 0.8, 
         cex = 1, 
         text.col = "black", 
         horiz = F , 
         inset = c(0.1, 0.1))
  
  
  matplot(x = logit_pfa_6[,2][21:33], y = logit_pfa_6[,4][21:33], 
          ylim = range(c(0, 0.5)), xlim = range(c(250, 400)), type = c("b"),
          pch=1, col = 'black', xlab = '# of total rejections', 
          ylab = 'Approximation of FDP(t)') #plot
  
  par(new=TRUE)
  lines(x = logit_pfa_7[,2][21:33],  y = logit_pfa_7[,4][21:33], type = c("b"), pch= 2, col = 'green')
  lines(x = logit_pfa_8[,2][21:33],  y = logit_pfa_8[,4][21:33], type = c("b"), pch= 3, col = 'blue')
  lines(x = logit_pfa_9[,2][21:33],  y = logit_pfa_9[,4][21:33], type = c("b"), pch= 4, col = 'orange')
  lines(x = logit_pfa_10[,2][21:33], y = logit_pfa_10[,4][21:33], type = c("b"), pch= 5, col = 'red')
  
  
  
  legend("topleft", legend=c("k = 6", "k = 7", 'k = 8', 'k = 9', ' k = 10'), 
         col=c("black", "green", 'blue', 'orange', 'red'), 
         pch = c(1, 2, 3, 4, 5),
         bty = "n", 
         pt.cex = 0.8, 
         cex = 1, 
         text.col = "black", 
         horiz = F , 
         inset = c(0.1, 0.1))
  
  
  dev.off()
  
##################################################
## Figure 4 - The main results of the analysis for k = 6
##################################################
  
##################################################
# Results for FDP(t) and V(t) based on k = 6
##################################################
  Rt = pfa_k_6$FDP$rejects
  t = pfa_k_6$FDP$t
  Vt = pfa_k_6$FDP$false.rejects
  FDPt = pfa_k_6$FDP$FDP
  FDP_6 = pfa_k_6$FDP

  rt = ggplot(FDP_6, aes(-log(t, base = 10), Rt)) + 
    geom_point()  + ylim(250, 1200) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylab("R(t)") + xlab('-log(t)') 
  
  
  vt = ggplot(FDP_6, aes(-log(t, base = 10), Vt)) + 
    geom_point()  + ylim(15, 950) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylab("Estimated V(t)") + xlab('-log(t)')
  
  fdpgg = ggplot(FDP_6, aes(-log(t, base = 10), FDPt)) + 
    geom_point()  + ylim(0.05, 0.8) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylab("Approximation of FDP(t)") + xlab('-log(t)')
  
  
  pdf(file = './results/Fig4.pdf')
  
  grid.arrange(rt, vt, fdpgg, ncol=3)
  
  dev.off()
  
  

  
  
  
