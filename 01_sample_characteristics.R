# ============================
# Sample descriptives inc. demographic table and trajectory plot
# =============================

# Set up ========================
##libs ===============
library(dplyr)
library(gtsummary)
library(tidyverse)
library(tidyr)
library(forestplot)
library(ggplot2)
library(scales)

## data 
load("MS_V2_26-09-2023/processed_data/pheno_wide.Rdata")

# Table 1: Sample Characteristics =======================
demo <- as.data.frame(dat_final) %>% select(
  colnames(dat_final[,c("Age",
                  "Sex",
                  "EducationYears",
                  "Income",
                  "Diagnosis",
                  "class",
                  "LatV_T1",
                  "LatV_T4",
                  "v_lcs1",
                  "v_lcs2",
                  "S")]))


Tbl1 <- 
  tbl_summary(data = demo, 
              by = class,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"), 
              digits = list(all_continuous() ~ c(2, 2, 2)),           
              missing = "ifany",
              missing_text = "Missing"
  )
# #add p-value 
Paper.Table.1 <- add_p(Tbl1,
                       pvalue_fun = purrr::partial(style_pvalue, digits = 3)) %>% bold_p() %>% add_overall()


# Supplementary Table 1 =======================================
load("MS_V2_26-09-2023/processed_data/all_phenos.Rdata") 
demo <- as.data.frame(full_dat) %>% select(
  colnames(full_dat[,c("Age",
                        "Sex",
                        "EducationYears",
                        "Income",
                        "Diagnosis",
                        "class",
                       "genotypes",
                       "LatV_T1",
                       "LatV_T4",
                       "v_lcs1",
                       "v_lcs2",
                       "S")]))


Tbl1 <- 
  tbl_summary(data = demo, 
              by = genotypes,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"), 
              digits = list(all_continuous() ~ c(2, 2, 2)),           
              missing = "ifany",
              missing_text = "Missing"
  )
# #add p-value 
Paper.Table.1 <- add_p(Tbl1,
                       pvalue_fun = purrr::partial(style_pvalue, digits = 3)) %>% bold_p() 

# Plot latent classes ===================================
plot_dat <- dat_final[,c("id", "class",
                         colnames(dat_final)[grep("LatV", colnames(dat_final))])]

plot_dat_long <- pivot_longer(plot_dat,
                              cols = LatV_T1:LatV_T7,
                              names_to = "timepoint",
                              values_to = "latent_score"
                              )
plot_dat_long$timepoint <- gsub("LatV_", "", plot_dat_long$timepoint)
plot_dat_long$timepoint <- as.factor(plot_dat_long$timepoint)
levels(plot_dat_long$timepoint) <- c("T1 Jan 2020", "T2 Mar/Apr 2020", 
                                     "T3 Jun 2020", "T4 Nov 2020", "T5 Dec 2020", 
                                     "T6 Jan 2021", "T7 Mar/Apr 2021")
levels(plot_dat_long$class) <- c("Most Vulnerable (n = 128)",
                                 "More Vulnerable (n = 296)",
                                 "More Resilient (n = 575)",
                                 "Most Resilient (n = 317)")
  
png(file = "MS_V2_26-09-2023/Figures/01_RV_trajectories.png",
    width = 10, height = 8, units = 'in',res=300)

plot <- plot_dat_long %>% 
  group_by(timepoint, class) %>% 
  summarise(
    n = n(),
    score = mean(latent_score, na.rm = T),
    SE = sd(latent_score, na.rm = T)/sqrt(n))

ggplot(na.omit(plot), 
               aes(x = timepoint,
                   y = score,
                   colour = class)) +
  geom_line(aes(group = class)) + 
  geom_line(data = na.omit(plot)) + 
  geom_point(aes(x = timepoint), size = 2) +
  scale_color_manual(values = c("#aa2339", "#cc7b88", "#7b88cc", "#2339aa")) + 
  geom_errorbar(aes(ymin=score-SE, ymax=score+SE), width=.1) +
  ylab(expression(paste(symbol('\254'), " Resilience ", 
                        "Vulnerability ", symbol('\256')))) +
  scale_x_discrete(
                   labels = label_wrap(3)) + 
  theme_classic() + 
   theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.x=element_text(size = 14),
    axis.title.x = element_blank()) 
dev.off()

#done.