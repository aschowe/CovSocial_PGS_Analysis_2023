# =============================
# 2 PGS correlations & validation using pearson's r and linear regression
# ============================

# set up =====================
##libs
library(writexl)
library(misty)

##func
source("MS_V2_26-09-2023/00_helper_functions.R")

## data 
load("Phenotypic_data/validate_PRS.Rdata")
load("MS_V2_26-09-2023/processed_data/PGSs_28082023.Rdata")
load("MS_V2_26-09-2023/processed_data/pheno_wide.Rdata")

##merge
dat <- left_join(validate_PRS, PGSs,
                 by = c( "GeneID"="IID"))
dat <- subset(dat, GeneID %in% dat_final$GeneID)
dat$DEP_PRS_standardized <- item.reverse(dat$DEP_PRS_standardized, min = -3, max = 4)
dat$NEU_PRS_standardized <- item.reverse(dat$NEU_PRS_standardized, min = -4, max = 3)

# Pearson's Correlations  =================================
PRSes <- colnames(dat)[grep("_standardized", colnames(dat))]

cor.mat <- corstars(dat[,PRSes])
rownames(cor.mat) <- gsub("_PRS_standardized", "", rownames(cor.mat))
colnames(cor.mat) <- gsub("_PRS_standardized", "", colnames(cor.mat))
cor.mat <- cbind(" "=rownames(cor.mat), cor.mat)
writexl::write_xlsx(cor.mat, "MS_V2_26-09-2023/Tables/02_validity_cormat.xlsx")

# Linear regression ======================================
## Depressive symptoms ------------------------------------
depress_PhQ9_differences <- as.data.frame(matrix(nrow = 0,
                                                 ncol = 6))

for(i in 1:length(PRSes)){
  form <- as.formula(paste("mean_depression ~ C1 + C2 +", 
                           PRSes[[i]]))
  M1 <- lm(form,
           data = dat)
  sumM1 <- as.data.frame(coef(summary(M1)))
  depress_PhQ9_differences[PRSes[[i]],c(1:8)] <- c(sumM1[PRSes[[i]],c(1:4)], 
                                                   nobs(M1),
                                                   summary(M1)$adj.r.squared,
                                                   confint(M1, PRSes[[i]],
                                                           level = 0.95))

  names(depress_PhQ9_differences) <- c("beta", "SE",
                                       "t-statistic", "P",
                                       "n", "Rsqr", "CI_low","CI_high" )
}

depress_PhQ9_differences <- depress_PhQ9_differences[order(depress_PhQ9_differences$Rsqr, decreasing=T),]
depress_PhQ9_differences <- round(depress_PhQ9_differences, 3)
rownames(depress_PhQ9_differences) <- gsub("_PRS_standardized", "", rownames(depress_PhQ9_differences))
depress_PhQ9_differences <- rownames_to_column(depress_PhQ9_differences, "PGS")

## save table
writexl::write_xlsx(depress_PhQ9_differences, "MS_V2_26-09-2023/Tables/02_validity_depress.xlsx")

## Plot

setEPS()
postscript("MS_V2_26-09-2023/Figures/02_PGS_validity_forest_plots_depress.eps")
depress_PhQ9_differences |>
  forestplot(labeltext = c(PGS, beta, Rsqr),
             mean = beta,
             lower = CI_low,
             upper = CI_high,
             zero = 0,
             cex = 2,
             lineheight = unit(2, "cm"),
             line.margin = 10,             
             xlav = "label here",
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             xticks=c(-0.15,0.15),
  ) |>
  fp_add_header("PGS", "beta", "adj. R-square") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_lines(h_2 = gpar(lty = 1, columns = 1:4), 
               h_21 = gpar(lwd = 1,
                           columns = 1:4)) |>
  fp_decorate_graph(box = gpar(lty = 2, col = "black"),
                    graph.pos = 3)
dev.off()

## Anxiety  ------------------------------------
anx_STAI_differences <- as.data.frame(matrix(nrow = 0,
                                             ncol = 6))

for(i in 1:length(PRSes)){
  form <- as.formula(paste("STAI ~ C1 + C2 +", 
                           PRSes[[i]]))
  M1 <- lm(form,
           data = dat)
  sumM1 <- as.data.frame(coef(summary(M1)))
  anx_STAI_differences[PRSes[[i]],c(1:8)] <- c(sumM1[PRSes[[i]],c(1:4)], 
                                                   nobs(M1),
                                                   summary(M1)$adj.r.squared,
                                               confint(M1, PRSes[[i]],
                                                       level = 0.95))
  
  names(anx_STAI_differences) <- c("beta", "SE",
                                       "t-statistic", "P",
                                       "n", "Rsqr", "CI_low","CI_high" )
}

anx_STAI_differences <- anx_STAI_differences[order(anx_STAI_differences$Rsqr, decreasing=T),]
anx_STAI_differences <- round(anx_STAI_differences, 3)
rownames(anx_STAI_differences) <- gsub("_PRS_standardized", "", rownames(anx_STAI_differences))
anx_STAI_differences <- rownames_to_column(anx_STAI_differences, "PGS")

writexl::write_xlsx(anx_STAI_differences, "MS_V2_26-09-2023/Tables/02_validity_traitanx.xlsx")

## Plot
setEPS()
postscript("MS_V2_26-09-2023/Figures/02_PGS_validity_forest_plots_STAI.eps")
anx_STAI_differences |>
  forestplot(labeltext = c(PGS, beta, Rsqr),
             mean = beta,
             lower = CI_low,
             upper = CI_high,
             zero = 0,
             cex = 2,
             lineheight = unit(2, "cm"),
             line.margin = 10,             
             xlav = "label here",
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             xticks=c(-0.5, -1, -2,-3, 0.5, 1,2,3),
  ) |>
  fp_add_header("PGS", "beta", "adj. R-square") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_lines(h_2 = gpar(lty = 1, columns = 1:4), 
               h_21 = gpar(lwd = 1,
                           columns = 1:4)) |>
  fp_decorate_graph(box = gpar(lty = 2, col = "black"),
                    graph.pos = 3)
dev.off()


## Loneliness  ------------------------------------
loneli_UCLA_differences <- as.data.frame(matrix(nrow = 0,
                                                ncol = 6))
for(i in 1:length(PRSes)){
  form <- as.formula(paste("UCLA ~ C1 + C2 +", 
                           PRSes[[i]]))
  M1 <- lm(form,
           data = dat)
  sumM1 <- as.data.frame(coef(summary(M1)))
  loneli_UCLA_differences[PRSes[[i]],c(1:8)] <- c(sumM1[PRSes[[i]],c(1:4)], 
                                               nobs(M1),
                                               summary(M1)$adj.r.squared,
                                               confint(M1, PRSes[[i]],
                                                       level = 0.95))
  
  names(loneli_UCLA_differences) <- c("beta", "SE",
                                       "t-statistic", "P",
                                       "n", "Rsqr", "CI_low","CI_high" )
}

loneli_UCLA_differences <- loneli_UCLA_differences[order(loneli_UCLA_differences$Rsqr, decreasing=T),]
loneli_UCLA_differences <- round(loneli_UCLA_differences, 3)
rownames(loneli_UCLA_differences) <- gsub("_PRS_standardized", "", rownames(loneli_UCLA_differences))
loneli_UCLA_differences <- rownames_to_column(loneli_UCLA_differences, "PGS")
writexl::write_xlsx(loneli_UCLA_differences, "MS_V2_26-09-2023/Tables/02_validity_traitlone.xlsx")

## Plot
setEPS()
postscript("MS_V2_26-09-2023/Figures/02_PGS_validity_forest_plots_UCLA.eps")
loneli_UCLA_differences |>
  forestplot(labeltext = c(PGS, beta, Rsqr),
             mean = beta,
             lower = CI_low,
             upper = CI_high,
             zero = 0,
             cex = 2,
             lineheight = unit(2, "cm"),
             line.margin = 10,             
             xlav = "label here",
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             xticks=c(-0.2,0.2),
  ) |>
  fp_add_header("PGS", "beta", "adj. R-square") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_lines(h_2 = gpar(lty = 1, columns = 1:4), 
               h_21 = gpar(lwd = 1,
                           columns = 1:4)) |>
  fp_decorate_graph(box = gpar(lty = 2, col = "black"),
                    graph.pos = 3)
dev.off()


## Neuroticism  ------------------------------------
neurotic_NEO_differences <- as.data.frame(matrix(nrow = 0,
                                                 ncol = 6))

for(i in 1:length(PRSes)){
  form <- as.formula(paste("NEO ~ C1 + C2 +", 
                           PRSes[[i]]))
  M1 <- lm(form,
           data = dat)
  sumM1 <- as.data.frame(coef(summary(M1)))
  neurotic_NEO_differences[PRSes[[i]],c(1:8)] <- c(sumM1[PRSes[[i]],c(1:4)], 
                                                  nobs(M1),
                                                  summary(M1)$adj.r.squared,
                                                  confint(M1, PRSes[[i]],
                                                          level = 0.95))
  
  names(neurotic_NEO_differences) <- c("beta", "SE",
                                       "t-statistic", "P",
                                       "n", "Rsqr", "CI_low","CI_high" )
}

neurotic_NEO_differences <- neurotic_NEO_differences[order(neurotic_NEO_differences$Rsqr, decreasing=T),]
neurotic_NEO_differences <- round(neurotic_NEO_differences, 3)
rownames(neurotic_NEO_differences) <- gsub("_PRS_standardized", "", rownames(neurotic_NEO_differences))
neurotic_NEO_differences <- rownames_to_column(neurotic_NEO_differences, "PGS")

writexl::write_xlsx(neurotic_NEO_differences, "MS_V2_26-09-2023/Tables/02_validity_traitNEO.xlsx")

## Plot
setEPS()
postscript("MS_V2_26-09-2023/Figures/02_PGS_validity_forest_plots_NEO.eps")
neurotic_NEO_differences |>
  forestplot(labeltext = c(PGS, beta, Rsqr),
             mean = beta,
             lower = CI_low,
             upper = CI_high,
             zero = 0,
             cex = 2,
             lineheight = unit(2, "cm"),
             line.margin = 10,             
             xlav = "label here",
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             xticks=c(-2.5, -0.5, -1, -2, 0.5, 1,2, 2.5),
             #title = "Most Resilient vs More Resilient"
             #mar = unit(c(0, 0 , 0 ,2), "mm")
  ) |>
  fp_add_header("PGS", "beta", "adj. R-square") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_lines(h_2 = gpar(lty = 1, columns = 1:4), 
               h_21 = gpar(lwd = 1,
                           columns = 1:4)) |>
  fp_decorate_graph(box = gpar(lty = 2, col = "black"),
                    graph.pos = 3)
dev.off()

## Life Satisfaction  ------------------------------------
LifeSatis_SWLS_differences <- as.data.frame(matrix(nrow = 0,
                                                   ncol = 5))

for(i in 1:length(PRSes)){
  form <- as.formula(paste("SWLS ~ C1 + C2 +", 
                           PRSes[[i]]))
  M1 <- lm(form,
           data = dat)
  sumM1 <- as.data.frame(coef(summary(M1)))
  LifeSatis_SWLS_differences[PRSes[[i]],c(1:8)] <- c(sumM1[PRSes[[i]],c(1:4)], 
                                                   nobs(M1),
                                                   summary(M1)$adj.r.squared,
                                                   confint(M1, PRSes[[i]],
                                                           level = 0.95))
  
  names(LifeSatis_SWLS_differences) <- c("beta", "SE",
                                       "t-statistic", "P",
                                       "n", "Rsqr", "CI_low","CI_high" )
}

LifeSatis_SWLS_differences <- LifeSatis_SWLS_differences[order(LifeSatis_SWLS_differences$Rsqr, decreasing=T),]
LifeSatis_SWLS_differences <- round(LifeSatis_SWLS_differences, 3)
rownames(LifeSatis_SWLS_differences) <- gsub("_PRS_standardized", "", rownames(LifeSatis_SWLS_differences))
LifeSatis_SWLS_differences <- rownames_to_column(LifeSatis_SWLS_differences, "PGS")

writexl::write_xlsx(LifeSatis_SWLS_differences, 
                    "MS_V2_26-09-2023/Tables/02_validity_traitLiSat.xlsx")

## Plot
setEPS()
postscript("MS_V2_26-09-2023/Figures/02_PGS_validity_forest_plots_SWLS.eps")
LifeSatis_SWLS_differences |>
  forestplot(labeltext = c(PGS, beta, Rsqr),
             mean = beta,
             lower = CI_low,
             upper = CI_high,
             zero = 0,
             cex = 2,
             lineheight = unit(2, "cm"),
             line.margin = 10,             
             xlav = "label here",
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             xticks=c(-1.5, -1, -0.5, -0.2, 0.2, 0.5, 1,1.5),
             #title = "Most Resilient vs More Resilient"
             #mar = unit(c(0, 0 , 0 ,2), "mm")
  ) |>
  fp_add_header("PGS", "beta", "adj. R-square") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_lines(h_2 = gpar(lty = 1, columns = 1:4), 
               h_21 = gpar(lwd = 1,
                           columns = 1:4)) |>
  fp_decorate_graph(box = gpar(lty = 2, col = "black"),
                    graph.pos = 3)
dev.off()

## Educational attainment  ------------------------------------
edu_years_differences <- as.data.frame(matrix(nrow = 0,
                                                   ncol = 5))

for(i in 1:length(PRSes)){
  form <- as.formula(paste("EducationYears~ C1 + C2 +", 
                           PRSes[[i]]))
  M1 <- lm(form,
           data = dat)
  sumM1 <- as.data.frame(coef(summary(M1)))
  edu_years_differences[PRSes[[i]],c(1:8)] <- c(sumM1[PRSes[[i]],c(1:4)], 
                                                     nobs(M1),
                                                     summary(M1)$adj.r.squared,
                                                     confint(M1, PRSes[[i]],
                                                             level = 0.95))
  
  names(edu_years_differences) <- c("beta", "SE",
                                         "t-statistic", "P",
                                         "n", "Rsqr", "CI_low","CI_high" )
}

edu_years_differences <- edu_years_differences[order(edu_years_differences$Rsqr, decreasing=T),]
edu_years_differences <- round(edu_years_differences, 3)
rownames(edu_years_differences) <- gsub("_PRS_standardized", "", rownames(edu_years_differences))
edu_years_differences <- rownames_to_column(edu_years_differences, "PGS")

writexl::write_xlsx(edu_years_differences, 
                    "MS_V2_26-09-2023/Tables/02_validity_EDUyears.xlsx")

## Plot
setEPS()
postscript("MS_V2_26-09-2023/Figures/02_PGS_validity_forest_plots_EDUyears.eps")
edu_years_differences |>
  forestplot(labeltext = c(PGS, beta, Rsqr),
             mean = beta,
             lower = CI_low,
             upper = CI_high,
             zero = 0,
             cex = 2,
             lineheight = unit(2, "cm"),
             line.margin = 10,             
             xlav = "label here",
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             xticks=c(-1.5, -0.2, -0.5, -1, 0.2, 0.5, 1,1.5),
  ) |>
  fp_add_header("PGS", "beta", "adj. R-square") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_lines(h_2 = gpar(lty = 1, columns = 1:4), 
               h_21 = gpar(lwd = 1,
                           columns = 1:4)) |>
  fp_decorate_graph(box = gpar(lty = 2, col = "black"),
                    graph.pos = 3)
dev.off()

#done.
