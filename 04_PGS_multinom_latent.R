# ===========================
# 04_PGSs_multinom_classes 
# Run multinomial regression of PGS predicting latent symptom trajectory
# Plot output 
# ===========================

# set up ===================
## libs 
library(nnet) #multinomial function
library(mcca) #rsquare for multinomial
library(ggplot2)
library(dplyr)
library(ggpubr)
library(regclass)
library(tibble)
library(patchwork)

##funct
source("MS_V2_26-09-2023/00_helper_functions.R")

## data 
load("MS_V2_26-09-2023/processed_data/pheno_wide.Rdata")
load("MS_V2_26-09-2023/processed_data/PGSs_28082023.Rdata")

Full_data <- left_join(dat_final, PGSs,
                       by = c( "GeneID"="IID"))
Full_data <- distinct(Full_data)

## Predictors
PGSes <- c("WBS_PRS_standardized",
           "General_NEU_PRS_standardized",
           "H_PRS_standardized")

## Outcome
Full_data$class4 <- relevel(Full_data$class, ref = "4") #most resilient as ref

## outcome DF 
results <- data.frame(matrix(nrow = 0, ncol = 7))

# run multinom ==================================
# ref = class 4 
for(i in 1:length(PGSes)){
  Score <- PGSes[[i]]
  
  #run M0
  M0 <- multinom(formula = paste("class4 ~ Age + Sex + Income"),
                 #), 
                 data = Full_data)     
  Rsq0 <- rsq(y = Full_data$class4,
              d = Full_data[,c("Age",
                               "Sex",
                               "Income")], 
              method = "multinom")
  #Adjusted R-squared value is calculated using the formula: 
  # 1 - (1 - R-squared) * ((n - 1)/(n - p - 1)), n=samples, p=predictors
  Rsq0_adj <- 1 - (1-Rsq0$measure)*((1316-1)/(1316-4-1))
  
  #run M1 for each score
  M1 <- multinom(formula = paste("class4 ~ Age + Sex + Income  + C1 + C2 + ", 
                                 Score
  ), 
  data = Full_data)   
  #compute p-value for output 
  Rsq1 <- rsq(y = Full_data$class4,
              d = Full_data[,c("Age",
                               "Sex",
                               "Income",
                               Score)], 
              method = "multinom")
  Rsq1_adj <- round(1 - (1-Rsq1$measure)*((1316-1)/(1316-5-1)),3)
  
  odds <- round(exp(coef(M1)),2) #get odds 
  ci <- exp(confint(M1, level=0.95))
  ci1 <- round(ci[,,1],2)
  ci2 <- round(ci[,,2],2)
  ci3 <- round(ci[,,3],2)
  z <- summary(M1)$coefficients/summary(M1)$standard.errors #derive z 
  p <- (1 - pnorm(abs(z), 0, 1)) * 2      
  odds[1,Score] <- ifelse(p[1,Score] < 0.05 & p[1,Score] >= 0.01, 
                          paste0(odds[1,Score],"*"), odds[1,Score])
  odds[2,Score] <- ifelse(p[2,Score] < 0.05 & p[2,Score] >= 0.01, 
                          paste0(odds[2,Score],"*"), odds[2,Score])
  odds[3,Score] <- ifelse(p[3,Score] < 0.05 & p[3,Score] >= 0.01, 
                          paste0(odds[3,Score],"*"), odds[3,Score])
  odds[1,Score] <- ifelse(p[1,Score] < 0.01 & p[1,Score] >= 0.001, 
                          paste0(odds[1,Score],"**"), odds[1,Score])
  odds[2,Score] <- ifelse(p[2,Score] < 0.01 & p[2,Score] >= 0.001, 
                          paste0(odds[2,Score],"**"), odds[2,Score])
  odds[3,Score] <- ifelse(p[3,Score] < 0.01 & p[3,Score] >= 0.001, 
                          paste0(odds[3,Score],"**"), odds[3,Score])
  odds[1,Score] <- ifelse(p[1,Score] < 0.001, 
                          paste0(odds[1,Score],"***"), odds[1,Score])
  odds[2,Score] <- ifelse(p[2,Score] < 0.001, 
                          paste0(odds[2,Score],"***"), odds[2,Score])
  odds[3,Score] <- ifelse(p[3,Score] < 0.001, 
                          paste0(odds[3,Score],"***"), odds[3,Score])
    # print("#############################")
    # print(Score)
    # print(p)
    # print("#############################")
  
  Rdiff <- round(Rsq1_adj - Rsq0_adj,3)
  aov <- anova(M0,M1)
  LR <- round(aov$`LR stat.`[2],3)
  LR_p <- round(aov$`Pr(Chi)`[2],3)
  
  results[Score,c(1:ncol(results))] <- c(paste0(odds[1,Score], "(",
                                         ci1[Score,1],",", ci1[Score,2], ")"),
                                         paste0(odds[2,Score], "(",
                                               ci2[Score,1],",", ci2[Score,2], ")"),
                                         p[2,Score], 
                                         paste0(odds[3,Score], "(",
                                         ci3[Score,1],",", ci3[Score,2], ")"),
                                         Rsq1_adj, 
                                         Rdiff,
                                         LR,
                                         LR_p)
   
}

names(results) <- c(paste0(c("OR"),c(1:3)), 
                    "adj_Rsqr",
                    "added_Rsqr",
                    "LRT",
                    "chiˆ2 p")
#results <- results[order(results$p1, decreasing=FALSE),]
#results <- round(results,3)
results <- rownames_to_column(results, "PGS")
results$PGS <- gsub("_PRS_standardized", " " , results$PGS)

# save output in table ====================
writexl::write_xlsx(results, "MS_V2_26-09-2023/Tables/04_PGS_multinom.xlsx")

# Boxplot WBS ====================================
levels(Full_data$class) <- c("Most Vulnerable",
                             "More Vulnerable",
                             "More Resilient",
                             "Most Resilient")

# setEPS()
# postscript("MS_V2_26-09-2023/Figures/04_PGS_multiom_boxplot_WeB.eps")

p1 <- ggplot(Full_data, aes(x = class, y = Full_data[, "WBS_PRS_standardized"], fill = class)) +
  geom_boxplot() +
  theme_classic() +
  labs(
    y = "WBS PGS",
    x = "") + 
  theme(legend.position = "none",
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.x=element_blank()) +
  scale_fill_manual(values = c("#aa2339", "#cc7b88", "#7b88cc", "#2339aa")) +
  scale_y_continuous(limits = c(NA,4.5)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_signif(comparisons = list(c("Most Resilient", "Most Vulnerable")),
              y_position = 1,
              annotations = "***",
              textsize = 8,
              vjust = 0.5,
              extend_line = - 0.02) +
  geom_signif(comparisons = list(c("Most Resilient", "More Vulnerable")),
              y_position = 1.8,
              annotations = "**",
              textsize = 8,
              vjust = 0.5,
              extend_line = - 0.02) +
    geom_signif(comparisons = list(c("Most Resilient", "More Resilient")),
                y_position = 2.6,
                annotations = "*",
                textsize = 8,
                vjust = 0.5,
                extend_line = - 0.02)
#dev.off()

# Boxplot NEU ====================================

p2 <- ggplot(Full_data, aes(x = class, y = Full_data[, "General_NEU_PRS_standardized"], fill = class)) +
  geom_boxplot() +
  #stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.3, size = 5) +
  theme_classic() +
  labs(
    y = "General NEU PGS",
    x = "") + 
  theme(legend.position = "none",
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.x=element_text(size = 12)) +
  scale_fill_manual(values = c("#aa2339", "#cc7b88", "#7b88cc", "#2339aa")) +
  scale_y_continuous(limits = c(NA,4.5)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_signif(comparisons = list(c("Most Resilient", "Most Vulnerable")),
              y_position = 1,
              annotations = "***",
              textsize = 8,
              vjust = 0.5,
              extend_line = - 0.02) +
  geom_signif(comparisons = list(c("Most Resilient", "More Vulnerable")),
              y_position = 1.8,
              annotations = "**",
              textsize = 8,
              vjust = 0.5,
              extend_line = - 0.02) +
  geom_signif(comparisons = list(c("Most Resilient", "More Resilient")),
              y_position = 2.6,
              annotations = "**",
              textsize = 8,
              vjust = 0.5,
              extend_line = - 0.02)

setEPS()
postscript("MS_V2_26-09-2023/Figures/04_PGS_multiom_boxplot.eps")
p1/p2 + plot_annotation(tag_levels = 'a')
dev.off()

# #done.
