#######################################
# libraries
#######################################
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(afex)
library(ggpubr)
library(car)
library(emmeans)
library(glmulti)
library(factoextra)
library(vegan)
library(sjPlot)

##########################
# wrapper function for glmulti
mixed_glmulti <- function (formula, data, random1 = "", random2 = "",...){
  lmer(paste(deparse(formula), random1, random2), data = data, REML = F, ...)
}

# functions to remove outliers using quantile method

remove_outliers <- function (data, col = "x") {
  
  Q1 <- quantile(pull(data, col), .25)
  Q3 <- quantile(pull(data, col), .75)
  IQR <- IQR(pull(data, col))
  new_data <- subset(data, 
                     pull(data, col) > (Q1 - 1.5*IQR) & pull(data, col) < (Q3 + 1.5*IQR))
return(new_data)

}


# custom theme for plots
theme <- theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill = NA, linewidth = 1),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.border = element_rect(linewidth  = 2))

#######################################
# Data import and cleaning
#####################################
# importing data 
nutrients <- read.csv("~/Desktop/GCS/GCS_git/Data/nutrients_final_updated.csv")


# removing mbc data with NAs 
mbc_final <- nutrients %>% 
  na.omit()

# removing negative mbc values from the mbc_final dataframe
mbc_final <- mbc_final %>% 
  filter(mbc > 0)


# averaging mbc data by field and sampling period
mbc_final <- mbc_final %>% 
  group_by(field, sampling.period, tillage, irrigation, fertilizer, crop_rotation, 
           cover_crop, residue) %>% 
  summarise(mbc = mean(mbc, na.rm = T),
            clay = mean(clay, na.rm = T),
            gwc = mean(gwc, na.rm = T),
            P = mean(P, na.rm = T),
            K = mean(K, na.rm = T),
            Mg = mean(Mg, na.rm = T),
            Ca = mean(Ca, na.rm = T),
            CEC = mean(CEC, na.rm = T),
            pH = mean(pH, na.rm = T),
            S = mean(S, na.rm = T),
            B = mean(B, na.rm = T),
            Zn = mean(Zn, na.rm = T),
            Mn = mean(Mn, na.rm = T),
            Fe = mean(Fe, na.rm = T),
            Cu = mean(Cu, na.rm = T),
            OM = mean(OM, na.rm = T),
            total.N = mean(total.N, na.rm = T))


mbc_final[, c(1:8)] <- lapply(mbc_final[, c(1:8)], factor)

mbc_final$tillage <- factor(mbc_final$tillage, 
                            levels = c("no","minimal","yes"))
mbc_final$residue <- factor(mbc_final$residue, 
                            levels = c("none-low", "medium-high"))
sapply(mbc_final, class)

###########################################
# fame data
fame <- read.csv("~/Desktop/GCS/GCS_git/Data/fame_final_updated.csv")


# converting into a longer format
fame_final <- fame%>% 
  pivot_longer(c(total.fame, BactSum, FungSum,FBRatio, GPosSum, GNegSum, ActSum,
                 SapSum, AMF, protozoa),
               names_to = "microbial.groups",
               values_to = "fame",
               values_drop_na = TRUE) %>% 
  relocate(microbial.groups, .before = "tillage")

fame_final[, c(2, 4:11)] <- lapply(fame_final[,c(2, 4:11)], factor)
fame_final$tillage <- factor(fame_final$tillage, 
                            levels = c("no","minimal","yes"))
fame_final$residue <- factor(fame_final$residue, 
                            levels = c("none-low", "medium-high"))

sapply(fame_final, class)


############################
# GWC
############################

boxplot(mbc_final$gwc) # no outliers

# GWC model
gwc_model <- glmulti(gwc ~ tillage + irrigation + fertilizer + crop_rotation +
                      residue + cover_crop, random1 = "+(1|field)", 
                    random2 = "+ (1|sampling.period)", data = mbc_final,
                    fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(gwc_model)$bestmodel
plot(gwc_model)
# the models below the red line has the aic difference of less than 2.

weightable(gwc_model)[1:2,] # checking the relative weight of top models
plot(gwc_model, type = "s") # checking the importance of the predictors
best_gwc_model <- gwc_model@objects[[1]] # extracting best model

Anova(best_gwc_model)
plot(best_gwc_model)
hist(residuals(best_gwc_model))
summary(best_gwc_model)


# gwc model with clay

gwc_model_2 <- glmulti(gwc ~ clay, random1 = "+(1|field)", 
                      random2 = "+ (1|sampling.period)", data = mbc_final,
                      fitfunction = mixed_glmulti, level = 1, crit = "aicc")

best_gwc_model_2 <- gwc_model_2@objects[[1]]
Anova(best_gwc_model_2)
plot(best_gwc_model_2)
summary(best_gwc_model_2)

# post hoc analysis to check the interaction effect between crop rotation and clay

gwc_int <- lmer(gwc ~ crop_rotation *clay + (1|field) + (1|sampling.period),
     data = mbc_final)
Anova(gwc_int)
summary(gwc_int)

predict_plot_gwc <-plot_model(gwc_int, type = "pred", term = c("clay", "crop_rotation"))

predict_data_gwc <- as.data.frame(predict_plot_gwc$data) %>%
  rename(crop_rotation = group)


# GWC-clay graph

(gwc_clay <- ggplot(mbc_final, aes(clay, gwc, color = crop_rotation)) +
    geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
    geom_line(data = predict_data_gwc, aes(x = x, y = predicted),
              linewidth = 1.5) +
    geom_ribbon(data = predict_data_gwc, aes(x = x, ymin = conf.low, ymax = conf.high,
                                             fill = crop_rotation), inherit.aes = F,
                alpha = 0.3  ) +
    scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                       values = c(15, 16, 17)) +
    scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                       values = c("#606060", "#FF9933")) +
    scale_fill_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                       values = c("#606060", "#FF9933")) +
    labs(y = expression("Gravimetric water content (g" *~g^-1* ")"),
         x = "Clay (%)" ) +
    theme +
    theme(legend.position = "top"))



#####################################

# averaging nutrients data by field 
nutrients_final <- nutrients %>% 
  group_by(field, sampling.period, tillage, irrigation, fertilizer, crop_rotation, 
           cover_crop, residue) %>% 
  summarise(clay = mean(clay, na.rm = T),
            gwc = mean(gwc, na.rm = T),
            P = mean(P, na.rm = T),
            K = mean(K, na.rm = T),
            Mg = mean(Mg, na.rm = T),
            Ca = mean(Ca, na.rm = T),
            CEC = mean(CEC, na.rm = T),
            pH = mean(pH, na.rm = T),
            S = mean(S, na.rm = T),
            B = mean(B, na.rm = T),
            Zn = mean(Zn, na.rm = T),
            Mn = mean(Mn, na.rm = T),
            Fe = mean(Fe, na.rm = T),
            Cu = mean(Cu, na.rm = T),
            OM = mean(OM, na.rm = T),
            total.N = mean(total.N, na.rm = T))

nutrients_final[, c(1:8)] <- lapply(nutrients_final[, c(1:8)], factor)

nutrients_final$tillage <- factor(nutrients_final$tillage, 
                                   levels = c("no","minimal","yes"))
nutrients_final$residue <- factor(nutrients_final$residue, 
                                   levels = c("none-low", "medium-high"))

sapply(nutrients_final, class)

###################
# N & P model
###################
boxplot(nutrients_final$total.N)

# removing outliers

nitrogen <- remove_outliers(nutrients_final, "total.N")
boxplot(nitrogen$total.N)
 
# N model 
N_model <- glmulti(log(total.N) ~ tillage + irrigation + fertilizer + crop_rotation +
                       residue + cover_crop, random1 = "+(1|field)", 
                     random2 = "+ (1|sampling.period)", data = nitrogen,
                     fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(N_model)$bestmodel
plot(N_model)
# the models below the red line has the aic difference of less than 2.

weightable(N_model)[1:2,] # checking the relative weight of top models
plot(N_model, type = "s") # checking the importance of the predictors
best_N_model <- N_model@objects[[1]] # extracting the best model

Anova(best_N_model)
plot(best_N_model)
hist(residuals(best_N_model))
summary(best_N_model)


# Total N model with continuous predictors

N_model_2 <- glmulti(log(total.N) ~ clay + OM + gwc, random1 = "+(1|field)", 
                       random2 = "+ (1|sampling.period)", data = nitrogen,
                       fitfunction = mixed_glmulti, level = 1, crit = "aicc")

best_N_model_2 <- N_model_2@objects[[1]]
Anova(best_N_model_2)
plot(best_N_model_2)
summary(best_N_model_2)


# post hoc analysis to check the interaction effect between residue and OM

N_int <- lmer(total.N ~ residue * OM + (1|field) + (1|sampling.period),
                data = nitrogen)
Anova(N_int)
summary(N_int)

(predict_plot_N <-plot_model(N_int, type = "pred", term = c("OM", "residue")))

predict_data_N <- as.data.frame(predict_plot_N$data) %>%
  rename(residue = group)


# N vs OM graph

(N_om <- ggplot(nitrogen, aes(OM, total.N, color = residue)) +
    geom_point(size = 3, alpha = 1, aes(shape = fertilizer)) +
    geom_line(data = predict_data_N, aes(x = x, y = predicted),
              linewidth = 1.5) +
    geom_ribbon(data = predict_data_N, aes(x = x, ymin = conf.low, ymax = conf.high,
                                             fill = residue), inherit.aes = F,
                alpha = 0.3  ) +
    scale_shape_manual(labels = c("Not Fertilized", "Fertilized"),
                       values = c(16, 17)) +
    scale_color_manual(labels = c("None-Low Residue","Medium-High Residue"),
                       values = c("#131E3A", "#980019")) +
    scale_fill_manual(labels = c("None-Low Residue","Medium-High Residue"),
                       values = c("#131E3A", "#980019")) + 
    labs(y = "Total Nitrogen (ppm)",
         x = "Organic Matter (%)" ) +
    theme +
    theme(legend.position = "top"))



#####################
# P
#####################
boxplot(nutrients_final$P)
# removing outliers

phosphorus <- remove_outliers(nutrients_final, "P")

boxplot(phosphorus$P)

# converting into kg/ha
phosphorus <- mutate(phosphorus,
                     P = P *1.121)


# P model
P_model <- glmulti(log(P) ~ tillage + irrigation + fertilizer + crop_rotation +
                     residue + cover_crop, random1 = "+(1|field)", 
                   random2 = "+ (1|sampling.period)", data = phosphorus,
                   fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(P_model)$bestmodel

plot(P_model)
# the models below the red line has the aic difference of less than 2.

weightable(P_model)[1:3,] # checking the relative weight of top models
plot(P_model, type = "s") # checking the importance of the predictors


best_P_model <- P_model@objects[[1]] # extracting the best model

Anova(best_P_model)
plot(best_P_model)
hist(residuals(best_P_model))
summary(best_P_model)


# total P model with continuous predictors

P_model_2 <- glmulti(log(P) ~ clay + OM + gwc, random1 = "+(1|field)", 
                     random2 = "+ (1|sampling.period)", data = phosphorus,
                     fitfunction = mixed_glmulti, level = 1, crit = "aicc")

best_P_model_2 <- P_model_2@objects[[1]]
Anova(best_P_model_2)
plot(best_P_model_2)
summary(best_P_model_2)

# post hoc analysis to check the interaction effect between irrigation and OM

P_int <- lmer(P ~ irrigation * OM + (1|field) + (1|sampling.period),
              data = phosphorus)
Anova(P_int)
summary(P_int)

(predict_plot_P <-plot_model(P_int, type = "pred", term = c("OM", "irrigation")))

predict_data_P <- as.data.frame(predict_plot_P$data) %>%
  rename(irrigation = group)


# total P vs OM graph

(P_om <- ggplot(phosphorus, aes(OM, P, color = irrigation)) +
    geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
    geom_line(data = predict_data_P, aes(x = x, y = predicted),
              linewidth = 1.5) +
    geom_ribbon(data = predict_data_P, aes(x = x, ymin = conf.low, ymax = conf.high,
                                           fill = irrigation), inherit.aes = F,
                alpha = 0.3  ) +
    scale_color_manual(labels = c("Dryland", "Irrigated"),
                       values = c("red", "blue")) +
    scale_fill_manual(labels = c("Dryland", "Irrigated"),
                      values = c("red", "blue")) +
    scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                       values = c(15, 16, 17)) +
    labs(y = "Total Phosphorus (kg/ha)",
         x = "Organic Matter (%)" ) +
    theme +
    theme(legend.position = "top"))


(NP_plot <- ggarrange(N_om + rremove("xlab") + rremove("x.text"),
          P_om,
          common.legend = FALSE, labels = "auto", align = "hv", ncol = 1,
          legend = "right"))


######################################################
# Organic matter (OM) 
#####################################################

boxplot(nutrients_final$OM) # no outliers

# OM model
om_model <- glmulti(OM ~ tillage + irrigation + fertilizer + crop_rotation +
                      residue + cover_crop, random1 = "+(1|field)", 
                    random2 = "+ (1|sampling.period)", data = nutrients_final,
                    fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(om_model)$bestmodel
plot(om_model) 
# the models below the red line has the aic difference of less than 2.

weightable(om_model)[1:2,] # checking the relative weight of top models
plot(om_model, type = "s") # checking the importance of the predictors
best_om_model <- om_model@objects[[1]] # extracting the best model

Anova(best_om_model)
plot(best_om_model)
hist(residuals(best_om_model))
summary(best_om_model)



# pairwise comparison
emmm_om_tillage <- emmeans(best_om_model, ~ tillage)
contrast(emmm_om_tillage, interaction = "pairwise")


# OM model with continuous predictors


om_model_2 <- glmulti(OM ~ clay + gwc, random1 = "+(1|field)", 
                    random2 = "+ (1|sampling.period)", data = nutrients_final,
                    fitfunction = mixed_glmulti, level = 1, crit = "aicc")

best_om_model_2 <- om_model_2@objects[[1]]
Anova(best_om_model_2)
plot(best_om_model_2)
hist(residuals(best_om_model_2))
summary(best_om_model_2)

# post hoc analysis to check the interaction effect between crop rotation and clay

OM_int <- lmer(OM ~ crop_rotation * clay + (1|field) + (1|sampling.period),
                data = nutrients_final)
Anova(OM_int)
summary(OM_int)

predict_plot_OM_clay <-plot_model(OM_int, type = "pred", term = c("clay", "crop_rotation"))

predict_data_OM_clay <- as.data.frame(predict_plot_OM_clay$data) %>%
  rename(crop_rotation = group)


# OM-clay graph

(OM_clay <- ggplot(nutrients_final, aes(clay, OM, color = crop_rotation)) +
    geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
    geom_line(data = predict_data_OM_clay, aes(x = x, y = predicted),
              linewidth = 1.5) +
    geom_ribbon(data = predict_data_OM_clay, aes(x = x, ymin = conf.low, ymax = conf.high,
                                             fill = crop_rotation), inherit.aes = F,
                alpha = 0.3  ) +
    scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                       values = c(15, 16, 17)) +
    scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                       values = c("#606060", "#FF9933")) +
    scale_fill_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                      values = c("#606060", "#FF9933")) +
    labs(y = "Organic matter (%)",
         x = "Clay (%)" ) +
    scale_y_continuous(breaks = seq(0, 2.5, 0.5)) +
    theme +
    theme(legend.position = "top"))



# post hoc analysis to check the interaction effect between crop rotation and GWC

OM_int_gwc <- lmer(OM ~ crop_rotation *gwc + (1|field) + (1|sampling.period),
                data = nutrients_final)
Anova(OM_int_gwc)
summary(OM_int_gwc)

(predict_plot_OM_gwc <-plot_model(OM_int_gwc, type = "pred", term = c("gwc", "crop_rotation")))

predict_data_OM_gwc <- as.data.frame(predict_plot_OM_gwc$data) %>%
  rename(crop_rotation = group)


# GWC-clay graph

(OM_gwc <- ggplot(nutrients_final, aes(gwc, OM, color = crop_rotation)) +
    geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
    geom_line(data = predict_data_OM_gwc, aes(x = x, y = predicted),
              linewidth = 1.5) +
    geom_ribbon(data = predict_data_OM_gwc, aes(x = x, ymin = conf.low, ymax = conf.high,
                                             fill = crop_rotation), inherit.aes = F,
                alpha = 0.3  ) +
    scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                       values = c(15, 16, 17)) +
    scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                       values = c("#606060", "#FF9933")) +
    scale_fill_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                      values = c("#606060", "#FF9933")) +
    labs(x = expression("Gravimetric water content (g" *~g^-1* ")"),
         y = "Organic matter (%)" ) +
    scale_y_continuous(breaks = seq(0,2.5, 0.5)) +
    theme +
    theme(legend.position = "top"))

# combining both graphs
ggarrange(OM_clay, OM_gwc + rremove("ylab") + rremove("y.text"),
          common.legend = TRUE, labels = "auto", align = "hv")


# Grand-average OM by different management practices

mean_om_tillage <- as.data.frame(emmeans(best_om_model, ~ tillage))

(om_tillage <- ggplot(nutrients_final, aes(tillage, OM)) +
  geom_violin() +
  geom_point(data = mean_om_tillage,
             aes(x = tillage, y = emmean), size = 4) +
  geom_errorbar(data = mean_om_tillage,
                aes(x = tillage, ymin = lower.CL, 
                    ymax = upper.CL), inherit.aes = FALSE,
                width = 0.1, linewidth = 1) +
  scale_x_discrete(labels = c("No-till", "Minimal-till", "Till")) +
  labs(x = "Tillage", y = "Organic matter (%)" ) +
  theme)

mean_om_irrigation <- as.data.frame(emmeans(best_om_model, ~ irrigation))

(om_irrigation <- ggplot(nutrients_final, aes(irrigation, OM)) +
  geom_violin() +
  geom_point(data = mean_om_irrigation,
             aes(x = irrigation, y = emmean), size = 4) +
  geom_errorbar(data = mean_om_irrigation,
                aes(x = irrigation, ymin = lower.CL, 
                    ymax = upper.CL), inherit.aes = FALSE,
                width = 0.1, linewidth = 1) +
  scale_x_discrete(labels = c("Dryland", "Irrigated")) +
  labs(x = "Irrigation", y = "Organic matter (%)" ) +
  theme)

mean_om_cr <- as.data.frame(emmeans(best_om_model, ~ crop_rotation))

(om_cr <- ggplot(nutrients_final, aes(crop_rotation, OM)) +
    geom_violin() +
    geom_point(data = mean_om_cr,
               aes(x = crop_rotation, y = emmean), size = 4) +
    geom_errorbar(data = mean_om_cr,
                  aes(x = crop_rotation, ymin = lower.CL, 
                      ymax = upper.CL), inherit.aes = FALSE,
                  width = 0.1, linewidth = 1) +
    scale_x_discrete(labels = c("Continuous Cotton", "Crop Rotation")) +
    labs(x = "Crop Rotation", y = "Organic matter (%)" ) +
    theme)

mean_om_cc <- as.data.frame(emmeans(best_om_model, ~ cover_crop))

(om_cc <- ggplot(nutrients_final, aes(cover_crop, OM)) +
    geom_violin() +
    geom_point(data = mean_om_cc,
               aes(x = cover_crop, y = emmean), size = 4) +
    geom_errorbar(data = mean_om_cc,
                  aes(x = cover_crop, ymin = lower.CL, 
                      ymax = upper.CL), inherit.aes = FALSE,
                  width = 0.1, linewidth = 1) +
    scale_x_discrete(labels = c("Fallow", "Winter Cover Crop")) +
    labs(x = "Cover Crops", y = "Organic matter (%)" ) +
    theme)

# combining all graphs
ggarrange(om_tillage + rremove("xlab"),
          om_cr + rremove("ylab") + rremove("y.text") + rremove("xlab"),
          om_irrigation + rremove("xlab"),
          om_cc + rremove("ylab") + rremove("y.text") + rremove("xlab"),
          common.legend = TRUE, labels = "auto", align = "hv")


#########################
# PCA of nutrients
###########################

nutrient_var <- c("P", "K", "Mg", "Ca", "S", "B", "Mn", "Zn", "Fe", "Cu",
                  "CEC")
# running PCA
pca_nutrients <- prcomp(mbc_final[, nutrient_var], scale = TRUE)
biplot(pca_nutrients)

#extracting pca coordinates
pca_coordinates <- data.frame(pca_nutrients$x)
pca_coordinates$tillage <- mbc_final$tillage
pca_coordinates$residue <- mbc_final$residue

# extracting pca loadings
pca_loadings <- pca_nutrients$rotation

# scaling loadings with eigen values
pca_loadings <- pca_loadings %*% diag(pca_nutrients$sdev^2)

pca_loadings <- as.data.frame(pca_loadings)

# calculating proportion of variance explained by PC1 and PC2
prop.expl <- summary(pca_nutrients)$importance
pc1 <- paste0("PC1 (", round(prop.expl[2] * 100,1), "%)")

pc2 <- paste0("PC2 (", round(prop.expl[5] * 100,1), "%)")

# pca_plot

(nutrient_PCA <-
  ggplot(pca_coordinates, aes(PC1, PC2)) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  geom_point(size = 3, alpha = 0.7, aes(shape = tillage, color = residue)) +
  geom_segment(data = 2* pca_loadings, aes( x = 0, xend = V1, y = 0, yend = V2),
               arrow = arrow(length = unit(0.15, "inches")), color = "black",
               linewidth = 1.2) +
  geom_text(data = 2* pca_loadings, aes(V1, V2), label = rownames(pca_loadings), 
            vjust = "outward",  hjust = "outward", fontface = "bold", color = "black") +
  scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                     values = c(15, 16, 17))+
  scale_color_manual(labels = c("None-Low","Medium-High"),
                     values = c("#131E3A", "#980019")) +
  labs(x = pc1, y = pc2) +
  theme +
  theme(legend.position = "top"))


###########################################
# adding PC1 and PC2 to original data frame
mbc_final$PC1<- pca_nutrients$x[, 1]
mbc_final$PC2<- pca_nutrients$x[, 2]

# contribution of variables to PC1
var <- get_pca_var(pca_nutrients)
fviz_contrib(pca_nutrients, choice = "var",axes = 1)


##########################
# mbc
##########################

boxplot(mbc_final$mbc)


# removing outliers
mbc_new <- remove_outliers(mbc_final, "mbc")
boxplot(mbc_new$mbc)

# mbc model
mbc_new$log.mbc <- log((mbc_new$mbc))

mbc_model <- glmulti(log.mbc ~ tillage + irrigation + fertilizer + crop_rotation +
                      residue + cover_crop, random1 = "+(1|field)", 
                    random2 = "+ (1|sampling.period)", data = mbc_new,
                    fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(mbc_model)$bestmodel
plot(mbc_model) 
# the models below the red line has the aic difference of less than 2.

weightable(mbc_model)[1:3,] # checking the relative weight of top models
plot(mbc_model, type = "s") # checking the importance of the predictors
best_mbc_model <- mbc_model@objects[[1]] # extracting the best model

Anova(best_mbc_model)
plot(best_mbc_model)
hist(residuals(best_mbc_model))
summary(best_mbc_model)

# pairwise comparison
contrast(emmeans(best_mbc_model, ~ tillage), interaction = "pairwise")
emmeans(best_mbc_model, ~ tillage)

######################################################

# mbc model with continuous predictors 
mbc_model_2 <- glmulti(log.mbc ~ clay + pH + OM + total.N + PC1 + PC2 + gwc,
                     random1 = "+(1|field)", 
                     random2 = "+ (1|sampling.period)", data = mbc_new,
                     fitfunction = mixed_glmulti, level = 1, crit = "aicc")


summary(mbc_model_2)$bestmodel
plot(mbc_model_2) 

best_mbc_model_2 <- mbc_model_2@objects[[1]]
Anova(best_mbc_model_2)
plot(best_mbc_model_2)
hist(residuals(best_mbc_model_2))
summary(best_mbc_model_2)

##########################################################
# post hoc analysis to check the interaction effect between residue and gwc

mbc_int_gwc <- lmer(mbc ~ residue * gwc + (1|field) + (1|sampling.period),
               data = mbc_new)
Anova(mbc_int_gwc)
summary(mbc_int_gwc)

(predict_plot_mbc_gwc <-plot_model(mbc_int_gwc, type = "pred",
                                  term = c("gwc", "residue")))

predict_data_mbc_gwc <- as.data.frame(predict_plot_mbc_gwc$data) %>%
  rename(residue = group)

(mbc_gwc <- ggplot(mbc_new, aes(gwc, mbc, color = residue)) +
  geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
  geom_line(data = predict_data_mbc_gwc, aes(x = x, y = predicted),
              linewidth = 1.5) +
  geom_ribbon(data = predict_data_mbc_gwc, aes(x = x, ymin = conf.low, ymax = conf.high,
                                                fill = residue), inherit.aes = F,
                alpha = 0.3  ) +
  scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                     values = c(15, 16, 17))+
  scale_color_manual(labels = c("None-Low","Medium-High"),
                     values = c("#131E3A", "#980019")) +
  scale_fill_manual(labels = c("None-Low","Medium-High"),
                       values = c("#131E3A", "#980019")) +
  labs(y = expression("Microbial Biomass Carbon (mg" *~kg^-1*")"),
       x = expression("Gravimetric water content (g" *~g^-1*")"))+
  theme +
  theme(legend.position = "top"))


# post hoc analysis to check the interaction effect between residue and clay

mbc_int_clay <- lmer(mbc ~ residue * clay + (1|field) + (1|sampling.period),
                    data = mbc_new)
Anova(mbc_int_clay)
summary(mbc_int_clay)

(predict_plot_mbc_clay <-plot_model(mbc_int_clay, type = "pred",
                                   term = c("clay", "residue")))

predict_data_mbc_clay <- as.data.frame(predict_plot_mbc_clay$data) %>%
  rename(residue = group)


(mbc_clay <- 
  ggplot(mbc_new, aes(clay, mbc, color = residue)) +
  geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
  geom_line(data = predict_data_mbc_clay, aes(x = x, y = predicted),
            linewidth = 1.5) +
  geom_ribbon(data = predict_data_mbc_clay, aes(x = x, ymin = conf.low, ymax = conf.high,
                                               fill = residue), inherit.aes = F,
              alpha = 0.3  ) +
  scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                     values = c(15, 16, 17))+
  scale_color_manual(labels = c("None-Low","Medium-High"),
                     values = c("#131E3A", "#980019")) +
  scale_fill_manual(labels = c("None-Low","Medium-High"),
                     values = c("#131E3A", "#980019")) +
  labs(y = expression("Microbial Biomass Carbon (mg" *~kg^-1*")"),
       x = "Clay (%)" ) +
  theme +
  theme(legend.position = "top"))

# combining both plots
 ggarrange(mbc_gwc,
           mbc_clay + rremove("ylab") + rremove("y.text"),
           common.legend = TRUE,
           labels = "auto", align = "hv")



################################################################################
# PCA of nutrients on fame dataset

# averaging FAME data by field and sampling period

fame_average <- fame_final %>% 
  group_by(field, sampling.period,microbial.groups, tillage, irrigation, fertilizer, crop_rotation, 
           cover_crop, residue) %>% 
  summarise(fame = mean(fame, na.rm = T),
            clay = mean(clay, na.rm = T),
            gwc = mean(gwc, na.rm = T),
            P = mean(P, na.rm = T),
            K = mean(K, na.rm = T),
            Mg = mean(Mg, na.rm = T),
            Ca = mean(Ca, na.rm = T),
            CEC = mean(CEC, na.rm = T),
            pH = mean(pH, na.rm = T),
            S = mean(S, na.rm = T),
            B = mean(B, na.rm = T),
            Zn = mean(Zn, na.rm = T),
            Mn = mean(Mn, na.rm = T),
            Fe = mean(Fe, na.rm = T),
            Cu = mean(Cu, na.rm = T),
            OM = mean(OM, na.rm = T),
            total.N = mean(total.N, na.rm = T))



pca_fnutrients <- prcomp(fame_average[, nutrient_var], scale = TRUE)
biplot(pca_fnutrients)
summary(pca_fnutrients)

# extracting PC1 & PC2
fame_average$PC1 <- pca_fnutrients$x[, 1]
fame_average$PC2 <- pca_fnutrients$x[, 2]


# total fame markers model
# subset total fame

total_fame <- fame_average %>% 
  filter(microbial.groups == "total.fame")

boxplot(total_fame$fame)
total_fame <- remove_outliers(total_fame, "fame")

boxplot(total_fame$fame)

# total fame model

total_fame_model <- glmulti(log(fame) ~ tillage + irrigation + fertilizer + crop_rotation +
                      residue + cover_crop, random1 = "+(1|field)", 
                    random2 = "+ (1|sampling.period)", data = total_fame,
                    fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(total_fame_model)$bestmodel
plot(total_fame_model) 
best_fame_model <- total_fame_model@objects[[1]]
Anova(best_fame_model)
plot(best_fame_model)
hist(residuals(best_fame_model))
summary(best_fame_model)

# pairwise comparison
contrast(emmeans(best_fame_model, ~ tillage), interaction = "pairwise")

# total fame with continuous predictors



total_fame %>%
  group_by(crop_rotation) %>% 
  summarize(mean = mean(fame))

total_fame_model_2 <- glmulti(log(fame) ~ clay + OM + pH + total.N + gwc + PC1 + PC2,
                              random1 = "+(1|field)", 
                              random2 = "+ (1|sampling.period)", data = total_fame,
                            fitfunction = mixed_glmulti, level = 1, crit = "aicc")



summary(total_fame_model_2)$bestmodel
best_fame_model_2 <- total_fame_model_2@objects[[1]]
Anova(best_fame_model_2)
plot(best_fame_model_2)
hist(residuals(best_fame_model_2))
summary(best_fame_model_2)

# post hoc analysis to check the interaction effect between crop rotation and OM

tfame_int_OM <- lmer(fame ~ crop_rotation * OM + (1|field) + (1|sampling.period),
                     data = total_fame)
Anova(tfame_int_OM)
summary(tfame_int_OM)

(predict_plot_tfame_OM <-plot_model(tfame_int_OM, type = "pred",
                                    term = c("OM", "crop_rotation")))

predict_data_tfame_OM <- as.data.frame(predict_plot_tfame_OM$data) %>%
  rename(crop_rotation = group)



# OM and fame graph

(tfame_om <- ggplot(total_fame, aes(OM, fame, color = crop_rotation)) +
  geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
  geom_line(data = predict_data_tfame_OM, aes(x = x, y = predicted),
              linewidth = 1.5) +
  geom_ribbon(data = predict_data_tfame_OM, aes(x = x, ymin = conf.low, ymax = conf.high,
                                                  fill = crop_rotation), inherit.aes = F,
                alpha = 0.3  ) +
  scale_fill_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                       values = c("#606060", "#FF9933")) +
  scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                     values = c("#606060", "#FF9933"))+
  scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                       values = c(15, 16, 17)) +
  labs(y = expression("Total FAME (nmol" *~g^-1*")"),
       x = "Organic Matter (%)" ) +
  theme +
  theme(legend.position = "top"))

# post hoc analysis to check the interaction effect between crop rotation and gwc

tfame_int_gwc <- lmer(fame ~ crop_rotation * gwc + (1|field) + (1|sampling.period),
                     data = total_fame)
Anova(tfame_int_gwc)
summary(tfame_int_gwc)

(predict_plot_tfame_gwc <-plot_model(tfame_int_gwc, type = "pred",
                                    term = c("gwc", "crop_rotation")))

predict_data_tfame_gwc <- as.data.frame(predict_plot_tfame_gwc$data) %>%
  rename(crop_rotation = group)




# total fame and gwc
 (tfame_gwc <- ggplot(total_fame, aes(gwc, fame, color = crop_rotation)) +
    geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
   
    geom_line(data = predict_data_tfame_gwc, aes(x = x, y = predicted),
               linewidth = 1.5) +
    geom_ribbon(data = predict_data_tfame_gwc, aes(x = x, ymin = conf.low, ymax = conf.high,
                                                   fill = crop_rotation), inherit.aes = F,
                 alpha = 0.3  ) +
    scale_fill_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                       values = c("#606060", "#FF9933")) +     
    scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                        values = c("#606060", "#FF9933"))+
    scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                        values = c(15, 16, 17)) +
    labs(y = expression("Total FAME (nmol" *~g^-1*")"),
         x = expression("Gravimetric water content (g" *~g^-1*")")) +
    theme +
    theme(legend.position = "top"))

ggarrange(tfame_om, tfame_gwc + rremove("ylab") + rremove("y.text"), labels = "auto",
           common.legend = T ,align = "hv")
 
#####################################################
# bacteria fame
# subset total fame

bacteria_fame <- fame_average %>% 
  filter(microbial.groups == "BactSum")

boxplot(bacteria_fame$fame)
# removing outliers
bacteria_fame <- remove_outliers(bacteria_fame, "fame")

boxplot(bacteria_fame$fame)

#bacterial fame model

bacteria_fame_model <- glmulti(log(fame) ~ tillage + irrigation + fertilizer + crop_rotation +
                              residue + cover_crop, random1 = "+(1|field)", 
                            random2 = "+ (1|sampling.period)", data = bacteria_fame,
                            fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(bacteria_fame_model)$bestmodel
plot(bacteria_fame_model) 
best_bacteria_model <- bacteria_fame_model@objects[[1]]
Anova(best_bacteria_model)
plot(best_bacteria_model)
hist(residuals(best_bacteria_model))
summary(best_bacteria_model)

# pairwise comparison
contrast(emmeans(best_bacteria_model, ~ tillage), interaction = "pairwise")


# bacteria fame with continuous predictors


bacteria_fame_model_2 <- glmulti(log(fame) ~ clay + OM + pH + gwc + total.N + PC1 + PC2, 
                              random1 = "+(1|field)", 
                              random2 = "+ (1|sampling.period)", data = bacteria_fame,
                              fitfunction = mixed_glmulti, level = 1, crit = "aicc")



summary(bacteria_fame_model_2)$bestmodel
best_bacteria_model_2 <- bacteria_fame_model_2@objects[[1]]
Anova(best_bacteria_model_2)
plot(best_bacteria_model_2)
hist(residuals(best_bacteria_model_2))
summary(best_bacteria_model_2)



# post hoc analysis to check the interaction effect between crop rotation and clay

bacteria_int_OM <- lmer(fame ~ crop_rotation * OM + (1|field) + (1|sampling.period),
                     data = bacteria_fame)
Anova(bacteria_int_OM)
summary(bacteria_int_OM)

(predict_plot_bacteria_OM <-plot_model(bacteria_int_OM, type = "pred",
                                    term = c("OM", "crop_rotation")))

predict_data_bacteria_OM <- as.data.frame(predict_plot_bacteria_OM$data) %>%
  rename(crop_rotation = group)

# OM and bacterial fame graph
(bfame_om <- ggplot(bacteria_fame, aes(OM, fame, color = crop_rotation)) +
  geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
  geom_line(data = predict_data_bacteria_OM, aes(x = x, y = predicted),
              linewidth = 1.5) +
  geom_ribbon(data = predict_data_bacteria_OM, aes(x = x, ymin = conf.low, ymax = conf.high,
                                                  fill = crop_rotation), inherit.aes = F,
                alpha = 0.3  ) +
  scale_fill_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                      values = c("#606060", "#FF9933")) +
  scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                     values = c("#606060", "#FF9933"))+
  scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                       values = c(15, 16, 17)) +
  labs(y = expression("Bacterial FAME (nmol" *~g^-1*")"),
       x = "Organic Matter (%)" ) +
  theme +
  theme(legend.position = "top"))



#################################################
# fungi fame
# subset fungi fame

fungi_fame <- fame_average %>% 
  filter(microbial.groups == "FungSum")

boxplot(fungi_fame$fame)
# removing outliers
fungi_fame <- remove_outliers(fungi_fame, "fame")

boxplot(fungi_fame$fame)

# fungi model 
fungi_fame_model <- glmulti(log(fame) ~ tillage + irrigation + fertilizer + crop_rotation +
                                 residue + cover_crop, random1 = "+(1|field)", 
                               random2 = "+ (1|sampling.period)", data = fungi_fame,
                               fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(fungi_fame_model)$bestmodel
plot(fungi_fame_model) 
best_fungi_model <- fungi_fame_model@objects[[1]]
Anova(best_fungi_model)
plot(best_fungi_model)
hist(residuals(best_fungi_model))
summary(best_fungi_model)

# pairwise comparison between tillage types
contrast(emmeans(best_fungi_model, ~ tillage), interaction = "pairwise")


# fungal fame with continuous predictors


fungi_fame_model_2 <- glmulti(log(fame) ~ clay + OM + pH + gwc + total.N + PC1 +PC2,
                              random1 = "+(1|field)", 
                              random2 = "+ (1|sampling.period)", data = fungi_fame,
                              fitfunction = mixed_glmulti, level = 1, crit = "aicc")



summary(fungi_fame_model_2)$bestmodel
best_fungi_model_2 <- fungi_fame_model_2@objects[[1]]
Anova(best_fungi_model_2)
plot(best_fungi_model_2)
hist(residuals(best_fungi_model_2))
summary(best_fungi_model_2)


# post hoc analysis to check the interaction effect between crop_rotation and OM

fungi_int_OM <- lmer(fame ~ crop_rotation* OM + (1|field) + (1|sampling.period),
                     data = fungi_fame)
Anova(fungi_int_OM)
summary(fungi_int_OM)

(predict_plot_fungi_OM <-plot_model(fungi_int_OM, type = "pred",
                                    term = c("OM", "crop_rotation")))

predict_data_fungi_OM <- as.data.frame(predict_plot_fungi_OM$data) %>%
  rename(crop_rotation = group)


# OM and fungi fame graph

(ffungi_om <- ggplot(fungi_fame, aes(OM, fame, color = crop_rotation)) +
  geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
  geom_line(data = predict_data_fungi_OM, aes(x = x, y = predicted),
              linewidth = 1.5) +
  geom_ribbon(data = predict_data_fungi_OM, aes(x = x, ymin = conf.low, ymax = conf.high,
                                                  fill = crop_rotation), inherit.aes = F,
                alpha = 0.3  ) +
  scale_fill_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                      values = c("#606060", "#FF9933")) +
  scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                       values = c("#606060", "#FF9933"))+
  scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                       values = c(15, 16, 17)) +
  labs(y = expression("Fungal FAME (nmol" *~g^-1*")"),
       x = "Organic Matter (%)" ) +
  theme +
  theme(legend.position = "top"))


# post hoc analysis to check the interaction effect between crop rotation and gwc

fungi_int_gwc <- lmer(fame ~ crop_rotation * gwc + (1|field) + (1|sampling.period),
                     data = fungi_fame)
Anova(fungi_int_gwc)
summary(fungi_int_gwc)

(predict_plot_fungi_gwc <-plot_model(fungi_int_gwc, type = "pred",
                                    term = c("gwc", "crop_rotation")))

predict_data_fungi_gwc <- as.data.frame(predict_plot_fungi_gwc$data) %>%
  rename(crop_rotation = group)


# gwc and fungi fame graph
(ffungi_gwc <- ggplot(fungi_fame, aes(gwc, fame, color = crop_rotation)) +
  geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
  geom_line(data = predict_data_fungi_gwc, aes(x = x, y = predicted),
              linewidth = 1.5) +
  geom_ribbon(data = predict_data_fungi_gwc, aes(x = x, ymin = conf.low, ymax = conf.high,
                                                  fill = crop_rotation), inherit.aes = F,
                alpha = 0.3  ) +
  scale_fill_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                      values = c("#606060", "#FF9933")) +
  scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                       values = c("#606060", "#FF9933")) +
  scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                       values = c(15, 16, 17)) +
  labs(y = expression("Fungal FAME (nmol" *~g^-1*")"),
       x = expression("Gravimetric water content (g" *~g^-1*")")) +
  theme +
  theme(legend.position = "top"))


ggarrange(ffungi_om, ffungi_gwc + rremove("ylab") + rremove("y.text"), labels = "auto",
          common.legend = T ,align = "hv")



###############################################
# microbial biomass by tillage & irrigation graphs
fame_new <- fame_final %>%
  group_by(tillage, irrigation, microbial.groups) %>% 
  summarize(mean.fame = mean(fame),
            se = sd(fame)/sqrt(n()),
            uci = mean.fame + 1.96 * se,
            lci = mean.fame - 1.96 * se)

fame_new$microbial.groups <- factor(fame_new$microbial.groups,
                                    levels = c("BactSum", "GPosSum", "GNegSum",
                                               "ActSum", "FungSum", "SapSum",
                                               "AMF", "protozoa",
                                               "FBRatio","total.fame"
                                               ))
# renaming the levels                                   
levels(fame_new$microbial.groups) <-  c("Bacteria", "G+ Bacteria", 
                                        "G- Bacteria","Actinomycetes", "Fungi",
                                        "Saprophytic Fungi", "AMF", "Protozoa",
                                        "FB Ratio", "Total FAME")

(microbes_barplot<-
  ggplot(fame_new, aes(tillage, mean.fame, fill = irrigation)) +
  geom_col(position = position_dodge(preserve = "single"),
            width = 0.9) +
  scale_fill_brewer(palette = "Set1",
                    labels = c("Dryland", "Irrigated")) +
  geom_errorbar(aes(ymin = lci, ymax = uci),
                width = 0.5,
                position = position_dodge(preserve = "single", width = 0.9)) +
  facet_wrap(~ microbial.groups, scales = "free_y") +
  labs(y = expression("FAME (nmol" *~g^-1*")"), x = NULL) +
  scale_x_discrete(labels = c("No-till", "Minimal-till", "Till")) +
  theme +
  theme(legend.position = "top"))


# AMF
AMF_fame <- fame_average%>% 
  filter(microbial.groups == "AMF")

boxplot(AMF_fame$fame)

# removing outliers
AMF_fame <- remove_outliers(AMF_fame, "fame")

boxplot(AMF_fame$fame)

AMF_model <- glmulti(log(fame) ~ tillage + irrigation + fertilizer + crop_rotation +
                       residue + cover_crop, random1 = "+(1|field)", 
                     random2 = "+ (1|sampling.period)", data = AMF_fame,
                     fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(AMF_model)$bestmodel
plot(AMF_model)
best_AMF_model <-AMF_model@objects[[1]]
Anova(best_AMF_model)
plot(best_AMF_model)
hist(residuals(best_AMF_model))
summary(best_AMF_model)


AMF_model_2 <- glmulti(log(fame) ~ clay + OM + pH + gwc + total.N + PC1 +PC2,
                              random1 = "+(1|field)", 
                              random2 = "+ (1|sampling.period)", data = AMF_fame,
                              fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(AMF_model_2)$bestmodel
plot(AMF_model_2)
best_AMF_model_2 <-AMF_model_2@objects[[1]]
Anova(best_AMF_model_2)
plot(best_AMF_model_2)
hist(residuals(best_AMF_model_2))
summary(best_AMF_model_2)

# post hoc analysis to check the interaction effect between crop_rotation and gwc

AMF_int_gwc <- lmer(fame ~ crop_rotation * gwc + (1|field) + (1|sampling.period),
                     data = AMF_fame)
Anova(AMF_int_gwc)
summary(AMF_int_gwc)

(predict_plot_AMF_gwc <-plot_model(AMF_int_gwc, type = "pred",
                                    term = c("gwc", "crop_rotation")))

predict_data_AMF_gwc <- as.data.frame(predict_plot_AMF_gwc$data) %>%
  rename(crop_rotation = group)


# AMF and gwc graph

(AMF_gwc <- ggplot(AMF_fame, aes(gwc, fame, color = crop_rotation)) +
    geom_point(size = 3, alpha = 1, aes(shape = tillage)) +
    geom_line(data = predict_data_AMF_gwc, aes(x = x, y = predicted),
              linewidth = 1.5) +
    geom_ribbon(data = predict_data_AMF_gwc, aes(x = x, ymin = conf.low, ymax = conf.high,
                                                  fill = crop_rotation), inherit.aes = F,
                alpha = 0.3  ) +
    scale_fill_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                      values = c("#606060", "#FF9933")) +
    scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                       values = c("#606060", "#FF9933")) +
    scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                       values = c(15, 16, 17)) +
    labs(y = expression("AMF (nmol" *~g^-1*")"),
         x = expression("Gravimetric water content (g" *~g^-1*")")) +
    theme +
    theme(legend.position = "top"))


# FB Ratio

FBR_fame <- fame_average%>% 
  filter(microbial.groups == "FBRatio")

boxplot(FBR_fame$fame)

# removing outliers
FBR_fame <- remove_outliers(FBR_fame, "fame")

boxplot(FBR_fame$fame)

fbr_model <- glmulti(log(fame) ~ tillage + irrigation + fertilizer + crop_rotation +
                                 residue + cover_crop, random1 = "+(1|field)", 
                               random2 = "+ (1|sampling.period)", data = FBR_fame,
                               fitfunction = mixed_glmulti, level = 1, crit = "aicc")

summary(fbr_model)$bestmodel
plot(fbr_model) 
best_fbr_model <-fbr_model@objects[[1]]
Anova(best_fbr_model)
plot(best_fbr_model)
hist(residuals(best_fbr_model))
summary(best_fbr_model)

# pairwise comparison
contrast(emmeans(best_fbr_model, ~ tillage), interaction = "pairwise")

FBR_fame %>% 
  group_by(irrigation) %>% 
  summarise(x = mean(fame))

# fbr fame with continuous predictors


fbr_model_2 <- glmulti(fame ~ clay + OM + pH + gwc + total.N + PC1 + PC2, 
                                 random1 = "+(1|field)", 
                                 random2 = "+ (1|sampling.period)", data = FBR_fame,
                                 fitfunction = mixed_glmulti, level = 1, crit = "aicc")



summary(fbr_model_2)$bestmodel
best_fbr_model_2 <-fbr_model_2@objects[[1]]
Anova(best_fbr_model_2)
plot(best_fbr_model_2)
hist(residuals(best_fbr_model_2))
summary(fbr_model)


##############################################
# NMDS on different microbial groups
fame_new <- fame %>% 
  na.omit()

nmds_data <- fame_new [, c(13:15, 17:22)]

# converting data into matrix
nmds_matrix <- as.matrix(nmds_data)

set.seed(100)
nmds_result <- metaMDS(nmds_matrix, distance = "bray")
plot(nmds_result)


plt_irrigation <- gg_ordiplot(nmds_result, groups = fame_new$irrigation,
                            plot = FALSE)
irrigation_centroid <- plt_irrigation$df_mean.ord
irrigation_ellipse <- plt_irrigation$df_ellipse



nmds_scores <- data.frame(scores(nmds_result)$sites)

nmds_scores$tillage <- fame_new$tillage
nmds_scores$crop_rotation <- fame_new$crop_rotation
nmds_scores$irrigation <- fame_new$irrigation
nmds_scores$cover_crop <- fame_new$cover_crop
nmds_scores$residue <- fame_new$residue

# extracting loadings vectors

nmds_loadings <- as.data.frame(scores(nmds_result)$species)
nmds_loadings$groups <- c("Protozoa", "Bacteria", "Fungi", "G+","G-", "Actinomycetes",
                          "Saprophytes", "AMF", "Total FAME")


nmds_irri <-
  ggplot() +
    geom_point(data = irrigation_centroid, aes(x, y, color = Group), size = 6,) +
    # geom_hline(yintercept = 0, lty = "dashed") +
    # geom_vline(xintercept = 0, lty = "dashed") +
    geom_segment(data = nmds_loadings, aes( x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                 arrow = arrow(length = unit(0.05, "inches")), color = "black",
                 linewidth = 1.0) +
    stat_ellipse(data = irrigation_ellipse, aes(x, y, color = Group), linewidth = 1) +
    geom_text(data = nmds_loadings, aes(NMDS1, NMDS2, label = groups), 
              vjust = "outward",  hjust = "outward", color = "black") +
    scale_color_manual(labels = c("Dryland","Irrigated"),
                       values = c("red", "blue"))  +
    scale_y_continuous(limits = c(-0.15, 0.15)) +
    scale_x_continuous(limits = c(-0.15, 0.15)) +
    labs(x = "NMDS1", y = "NMDS2") +
    theme + 
    theme(legend.position = "top");nmds_irri



plt_tillage <- gg_ordiplot(nmds_result, groups = fame_new$tillage,
                              plot = FALSE)
tillage_centroid <- plt_tillage$df_mean.ord
tillage_ellipse <- plt_tillage$df_ellipse

nmds_tillage <-
  ggplot() +
  geom_point(data = tillage_centroid, aes(x, y, color = Group), size = 6) +
  # geom_hline(yintercept = 0, lty = "dashed") +
  # geom_vline(xintercept = 0, lty = "dashed") +
  geom_segment(data = nmds_loadings, aes( x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.05, "inches")), color = "black",
               linewidth = 1.0) +
  stat_ellipse(data = tillage_ellipse, aes(x, y, color = Group), linewidth = 1) +
  geom_text(data = nmds_loadings, aes(NMDS1, NMDS2, label = groups), 
            vjust = "outward",  hjust = "outward", color = "black") +
  scale_color_manual(labels =  c("Minimal-till", "No-till","Till"),
                     values = c("#000000","#660099", "#880000")) +
  scale_y_continuous(limits = c(-0.15, 0.15)) +
  scale_x_continuous(limits = c(-0.15, 0.15)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme + 
  theme(legend.position = "top");nmds_tillage


# # NMDS ordination plot
# 
# nmds_irri <- ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
#   geom_point(size = 3, aes(color = irrigation), alpha = 0.6) +
#   geom_hline(yintercept = 0, lty = "dashed") +
#   geom_vline(xintercept = 0, lty = "dashed") +
#   geom_segment(data = nmds_loadings, aes( x = 0, xend = NMDS1, y = 0, yend = NMDS2),
#                arrow = arrow(length = unit(0.15, "inches")), color = "black",
#                linewidth = 1.2) +
#   stat_ellipse(aes(color = irrigation), linewidth = 1) +
#   geom_text(data = nmds_loadings, aes(NMDS1, NMDS2, label = groups), 
#             vjust = "outward",  hjust = "outward", fontface = "bold", color = "black") +
#   scale_color_manual(labels = c("Dryland","Irrigated"),
#                     values = c("red", "blue"))  +
#   scale_y_continuous(limits = c(-0.15, 0.15)) +
#   scale_x_continuous(limits = c(-0.15, 0.15)) +
#   theme + 
#   theme(legend.position = "top")
#   
# 
# 
# nmds_till <- ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
#   geom_point(size = 3, aes(color = tillage), alpha = 0.6) +
#   geom_hline(yintercept = 0, lty = "dashed") +
#   geom_vline(xintercept = 0, lty = "dashed") +
#   geom_segment(data = nmds_loadings, aes( x = 0, xend = NMDS1, y = 0, yend = NMDS2),
#                arrow = arrow(length = unit(0.15, "inches")), color = "black",
#                linewidth = 1.2) +
#   stat_ellipse(aes(color = tillage), linewidth = 1) +
#   geom_text(data = nmds_loadings, aes(NMDS1, NMDS2, label = groups), 
#             vjust = "outward",  hjust = "outward", fontface = "bold", color = "black") +
#   scale_color_manual(labels =  c("Minimal-till", "No-till","Till"),
#                      values = c("#000000","#660099", "#880000")) +
#   scale_y_continuous(limits = c(-0.15, 0.15)) +
#   scale_x_continuous(limits = c(-0.15, 0.15)) +
#   theme + 
#   theme(legend.position = "top")


(nmds_graphs <- ggarrange(nmds_irri,
          nmds_tillage + rremove("ylab") + rremove("y.text"),
          common.legend = F,
          ncol = 2,
          labels = "auto",
          legend = "top"))

# permanova test to check the mutivariate effects of management practices on 
# microbial abundances
adonis_result <- adonis2(nmds_data ~ fame_new$tillage + fame_new$irrigation + 
                           fame_new$crop_rotation + fame_new$cover_crop +
                           fame_new$residue, method = "bray")

adonis_result



# Graphs for the paper
setwd("~/Desktop/GCS/graphs/")

Fig_4 <- ggsave("Fig4.pdf", gwc_clay, dpi = 300, width = 8, height =6 ,
                units = "in", device = "pdf")

(OM_combined <- ggarrange(OM_clay, OM_gwc + rremove("ylab") + rremove("y.text"),
                common.legend = TRUE, labels = "auto", align = "hv"))

Fig_5 <- ggsave("Fig5.pdf", OM_combined, dpi = 300, width = 12, height = 6 ,
                units = "in", device = "pdf")

(OM_management_combined <- ggarrange(om_tillage + rremove("xlab"),
          om_cr + rremove("ylab") + rremove("y.text") + rremove("xlab"),
          om_irrigation + rremove("xlab"),
          om_cc + rremove("ylab") + rremove("y.text") + rremove("xlab"),
          common.legend = TRUE, labels = "auto", align = "hv"))

Fig_6 <- ggsave("Fig6.pdf", OM_management_combined, dpi = 300, width = 8,
                height = 6 , units = "in", device = "pdf")


Fig_7 <- ggsave("Fig7.pdf", microbes_barplot, dpi = 300, width = 15, height = 8,
                units = "in", device = "pdf")

(fame_graphs <- ggarrange(tfame_gwc,
                          ffungi_gwc,
                          AMF_gwc,
                          tfame_om,
                          ffungi_om,
                          bfame_om,
                          common.legend = T,
                          align = "hv",
                          labels = "auto"))

Fig_8 <- ggsave("Fig8.pdf", fame_graphs, dpi = 300, width = 15, height =8 ,
                units = "in", device = "pdf")


Fig_9 <- ggsave("Fig9.pdf", nmds_graphs, dpi = 300, width = 12, height =6 ,
                units = "in", device = "pdf")

# Supplemental Graphs

S2 <- ggsave("S2.pdf", NP_plot, dpi = 300, width = 10, height =8 ,
             units = "in", device = "pdf")
S3 <- ggsave("S3.pdf", nutrient_PCA, dpi = 300, width = 8, height = 6 ,
             units = "in", device = "pdf")

mbc_combined <- 
  ggarrange(mbc_gwc,
          mbc_clay + rremove("ylab") + rremove("y.text"),
          common.legend = TRUE,
          labels = "auto", align = "hv")

S4 <- ggsave("S4.pdf", mbc_combined, dpi = 300, width = 12, height = 6 ,
             units = "in", device = "pdf")




# ##############################
# # extra plots
# 
# # OM-clay graph by tillage
# 
# (om_clay <- ggplot(nutrients_final, aes(clay, OM, color = tillage)) +
#    geom_point(size = 3, alpha = 1, aes(shape = crop_rotation)) +
#    geom_smooth(method = "lm", linewidth = 1.5) +
#    scale_color_manual(labels =  c("No-till", "Minimal-till","Till"),
#                       values = c("#000000","#660099", "#880000"))+
#    scale_shape_manual(labels = c("Continuous Cotton", "Crop Rotation"),
#                       values = c(15, 16)) +
#    labs(y = "Organic matter (%)", x = "Clay (%)" ) +
#    theme +
#    theme(legend.position = "top"))
# 
# # OM_gwc
# 
# (om_gwc <- ggplot(nutrients_final, aes(gwc, OM, color = tillage)) +
#     geom_point(size = 3, alpha = 1, aes(shape = crop_rotation)) +
#     geom_smooth(method = "lm", linewidth = 1.5) +
#     scale_color_manual(labels =  c("No-till", "Minimal-till","Till"),
#                        values = c("#000000","#660099", "#880000"))+
#     scale_shape_manual(labels = c("Continuous Cotton", "Crop Rotation"),
#                        values = c(15, 16)) +
#     labs(y = "Organic matter (%)", x = "Gravimetric water content (g/g)" ) +
#     theme +
#     theme(legend.position = "top"))
# 
# 
# ggarrange(om_clay, om_gwc + rremove("ylab") + rremove("y.text"),
#           common.legend = TRUE, labels = "auto", align = "hv")
# 
# 
# ggarrange(om_gwc + rremove("xlab") + rremove("x.text"),
#           tfame_gwc + rremove("xlab") + rremove("x.text"),
#           ffungi_gwc,
#           AMF_gwc, 
#           common.legend = T,
#           align = "hv",
#           labels = "auto")
