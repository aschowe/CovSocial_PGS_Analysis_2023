# ================================
# 03_ANOVA 
# Run anova to examine mean differences in PGS between symptom trajectories
# ================================

# set up  ========================
## libs 
library(dplyr)
library(rstatix)
library(writexl)

##func
source("MS_V2_26-09-2023/00_helper_functions.R")

## data
load("MS_V2_26-09-2023/processed_data/pheno_wide.Rdata")
load("MS_V2_26-09-2023/processed_data/PGSs_28082023.Rdata")

##merge
Full_data <- left_join(dat_final, PGSs,
                       by = c( "GeneID"="IID"))
Full_data <- distinct(Full_data)
Full_data$DEP_PRS_standardized <- item.reverse(Full_data$DEP_PRS_standardized, min = -3, max = 4)
Full_data$NEU_PRS_standardized <- item.reverse(Full_data$NEU_PRS_standardized, min = -4, max = 3)

## Predictors 
PGSes <- c(colnames(Full_data)[grep("_standardized", colnames(Full_data))])

## Output df 
aov_tab <- as.data.frame(matrix(nrow = 0,
                                ncol = 8))

# Run ANOVA ===================================

for(i in 1:length(PGSes)){
  PGS <- PGSes[[i]]
  form <- as.formula(paste(PGS, "~ class"))
  
  aov_tab[PGS,c(1:4)] <-c(mean(Full_data[[PGS]][Full_data$class == 1]), 
                          mean(Full_data[[PGS]][Full_data$class == 2]),
                          mean(Full_data[[PGS]][Full_data$class == 3]),
                          mean(Full_data[[PGS]][Full_data$class == 4]))
  
  res.aov <- aov(form, 
                 data = Full_data)
  res.aov <- anova_summary(res.aov)
  aov_tab[PGS,c(5:7)] <- res.aov[1,c(2,4,5)]
  aov_tab[PGS,c(8)] <- p.adjust(res.aov[1,c(5)], 
                                method = "bonferroni",
                                n = 19)
}
names(aov_tab) <- c("Most Vulnerable (n = 128)",
                    "More Vulnerable (n = 296)",
                    "More Resilient (n = 575)",
                    "Most Resilient (n = 317)",
                    "df",
                    "F",
                    "p",
                   "p-adj")


# Save output ======================
aov_tab <- aov_tab[order(aov_tab$p, decreasing=FALSE),]
aov_tab <- round(aov_tab,3)
aov_tab <- rownames_to_column(aov_tab, "PGS")
aov_tab$PGS <- gsub("_PRS_standardized", "", aov_tab$PGS)
write_xlsx(aov_tab, 
           "MS_V2_26-09-2023/Tables/03_PGS_aov_tab.xlsx")

#done.
