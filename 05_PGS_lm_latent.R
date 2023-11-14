# ===========================
# 05_PGS_lm_latent 
# Run linear regression with IV = PGS, DV = latent baseline and change scores 
# Summmarise output in table format
# ===========================

# set up ====================
## libs
library(ggplot2)
library(dplyr)
library(ggpubr)
library(regclass)

## data 
load("MS_V2_26-09-2023/processed_data/pheno_wide.Rdata")
load("MS_V2_26-09-2023/processed_data/PGSs_28082023.Rdata")

Full_data <- left_join(dat_final, PGSs,
                       by = c( "GeneID"="IID"))
Full_data <- distinct(Full_data)
Full_data$DEP_PRS_standardized <- item.reverse(Full_data$DEP_PRS_standardized, min = -3, max = 4)


## Predictors
PGSes <- c("WBS_PRS_standardized",
           "General_NEU_PRS_standardized",
           "H_PRS_standardized")

## Outcomes
outcomes <- c("LatV_T1",
              "v_lcs1", 
              "v_lcs2",
              "LatV_T4",
              "S")

## outcome DF
results <- data.frame(matrix(nrow = 0, ncol = 9))

# Run Linear Regression for latent scores ================================
for(o in 1:length(outcomes)){
  outcome <- outcomes[[o]]
  for(i in 1:length(PGSes)){
    Score <- PGSes[[i]]
    form <- as.formula(paste(outcome,"~ Age + Sex + 
                             Income"))
    M0 <- lm(form,
             data = Full_data)
    M0_adjR <- summary(M0)$adj.r.squared
    
    M1 <- update(M0, paste(".~. + C1 + C2 + " ,Score))
    if(outcome == "LatV_T4"){print(summary(M1))}
    sumM1 <- round(as.data.frame(coef(summary(M1))),3)
    M1_adjR <- round(summary(M1)$adj.r.squared,3)
    Rdiff <- round(M1_adjR - M0_adjR, 3)
    aov <- anova(M0,M1)
    aovF <- round(aov$F[2],3)
    aovp <- round(aov$`Pr(>F)`[2],3)
    
    results[paste(Score, outcome), c(1:9)] <- c(sumM1[Score,], nobs(M1), 
                                M1_adjR, Rdiff,
                                aovF, aovp)
    #print(summary(M1))
    #print(aovp)
    
  }
}

names(results) <- c("beta",
                    "SE",
                    "t-statistic",
                    "P",
                    "n",
                    "adj. R^2",
                    "added adj. R^2",
                    "F",
                    "Pr(>F)")
results <- rownames_to_column(results, "PGS")
results$PGS <- gsub("_PRS_standardized", " " , results$PGS)
latent_res <- results
writexl::write_xlsx(latent_res, "MS_V2_26-09-2023/Tables/05_PGS_lm_latent.xlsx")

#done.