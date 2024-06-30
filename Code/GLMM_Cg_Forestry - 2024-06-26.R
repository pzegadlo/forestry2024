library(ggpubr)
library(plyr)
library(tidyverse)
library(lme4)
library(AICcmodavg)
library(RVAideMemoire)

data <- read_csv("GLMM_Input_Cg_2024-06-25.csv", col_types = cols(
      VegH = col_factor(levels=c("1","2","3")),
      UndD = col_factor(levels=c("1","2","3")),
      I_VegH2_UndD = col_factor(),
      I_VegH3_UndD = col_factor(),
      I_PopHi_UndD = col_factor(),
      I_Yr_StS_R = col_factor(levels=c("0", "2017", "2018", "2019")),
      I_Yr_StS_S = col_factor(levels=c("0", "2017", "2018", "2019")),
      I_Yr_StS_T = col_factor(levels=c("0", "2018", "2019")),
      SquareID = col_factor()
    )
)

# Rescaling due to numerical issues during GLMM modelling
data$PopAf <- data$PopAf/100
data$PopCg <- data$PopCg/100

# Baseline GLMM
model_1_baseline <- glmer(formula = captures ~ 1 + VegH + UndD + SMI_4 + CWDv + CWDd45 +
                                 Plnt_1 + Plnt_2 + Plnt_3 + Plnt_4 + Plnt_5 + Plnt_8 + Plnt_11 + Plnt_12 +
                                 PopAf + PopCg +
                                 Year_2018 + Year_2019 +
                                 I_VegH2_UndD + I_VegH3_UndD +
                                 I_VegH2_CWDv + I_VegH3_CWDv + I_UndD2_CWDv + I_UndD3_CWDv + 
                                 I_Yr_StS_T +I_Yr_StS_R + I_Yr_StS_S + 
                                 I_PopHi_VegH_2 + I_PopHi_VegH_3 + I_PopHi_UndD + I_PopHi_CWDv +
                                 (1 | SquareID), data = data, family = poisson,
                                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=7e5)))  

summary(model_1_baseline)   
AICc(model_1_baseline)

overdisp.glmer(model_1_baseline)

# No humidity
model_2 <- glmer(formula = captures ~ 1 + VegH + UndD + CWDv + CWDd45 +
                   Plnt_1 + Plnt_2 + Plnt_3 + Plnt_4 + Plnt_5 + Plnt_8 + Plnt_11 + Plnt_12 +
                   PopAf + PopCg +
                   Year_2018 + Year_2019 +
                   I_VegH2_UndD + I_VegH3_UndD +
                   I_VegH2_CWDv + I_VegH3_CWDv + I_UndD2_CWDv + I_UndD3_CWDv + 
                   I_Yr_StS_T + I_Yr_StS_R + I_Yr_StS_S + 
                   I_PopHi_VegH_2 + I_PopHi_VegH_3 + I_PopHi_UndD + I_PopHi_CWDv +
                   (1 | SquareID), data = data, family = poisson,
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=7e5)))  

summary(model_2)   
AICc(model_2)

# No humidity + only chosen plant species groups
model_3 <- glmer(formula = captures ~ 1 + VegH + UndD + CWDv + CWDd45 +
                   Plnt_1 + Plnt_3 +
                   PopAf + PopCg +
                   Year_2018 + Year_2019 +
                   I_VegH2_UndD + I_VegH3_UndD +
                   I_VegH2_CWDv + I_VegH3_CWDv + I_UndD2_CWDv + I_UndD3_CWDv + 
                   I_Yr_StS_T + I_Yr_StS_R + I_Yr_StS_S + 
                   I_PopHi_VegH_2 + I_PopHi_VegH_3 + I_PopHi_UndD + I_PopHi_CWDv +
                   (1 | SquareID), data = data, family = poisson,
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=7e5)))  

summary(model_3)   
AICc(model_3)

# No humidity + only chosen plant species groups + only significant controls and interactions
model_4 <- glmer(formula = captures ~ 1 + VegH + UndD + CWDv + CWDd45 +
                   Plnt_1 + Plnt_3 +
                   PopCg + data$I_Yr_2019_StS_T +
                   I_UndD3_CWDv + I_PopHi_VegH_2 + I_PopHi_VegH_3 +
                   (1 | SquareID), data = data, family = poisson,
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=7e5)))   

summary(model_4)   
AICc(model_4)

# Residuals (actual - predicted), if they are negative, we overestimate, if they are positive, we underestimate
res <- residuals(model_4, type="pearson", scaled=TRUE)
data_tab <- data.frame(cbind(res, data$cwd_vol))
colnames(data_tab) <- c("res", "cwd_vol")
data_tab$cwd_vol <- round_any(data_tab$cwd_vol, 0.25)
data_chart <- data_tab %>% 
  group_by(cwd_vol) %>% 
  summarize(mean = mean(res), cnt = n())
ch_preprog = ggplot(data_chart, aes(x=cwd_vol, y=mean)) + 
  geom_point(aes(size = cnt)) + geom_hline(yintercept=0) + theme_bw() + 
  scale_x_continuous(limits=c(0,4), breaks=seq(0, 4, 0.5)) +
  scale_y_continuous(limits=c(-0.45,0.45), breaks=seq(-0.4, 0.4, 0.1)) +
  labs(x = bquote(atop('Coarse woody debris volume'~(m^3),'No CWD threshold')), y = "Mean Pearson residual", size = "Observation count") +
  theme(legend.position="bottom") + scale_size_area() +
  geom_rect(aes(xmin = 0.65, xmax = 1.10, ymin = -0.1, ymax = 0.2), 
            fill = "gray", alpha = 0.01, color = "black", linetype = 'dotted')


# Final choice

data$CWDv_T <- (data$cwd_vol > 0.75)

# Model 4 with added threshold
model_final <- glmer(formula = captures ~ 1 + VegH + UndD + CWDv + CWDd45 + CWDv_T +
                               Plnt_1 + Plnt_3 +
                               PopCg + data$I_Yr_2019_StS_T +
                               I_UndD3_CWDv + I_PopHi_VegH_2 + I_PopHi_VegH_3 +
                               (1 | SquareID), data = data, family = poisson,
                             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=7e5)))   

summary(model_final)   
AICc(model_final)

overdisp.glmer(model_final)

# Residuals plot for the final model
res <- residuals(model_final, type="pearson", scaled=TRUE)
data_tab <- data.frame(cbind(res, data$cwd_vol))
colnames(data_tab) <- c("res", "cwd_vol")
data_tab$cwd_vol <- round_any(data_tab$cwd_vol, 0.25)
data_chart <- data_tab %>% 
  group_by(cwd_vol) %>% 
  summarize(mean = mean(res), cnt = n())
ch_postprog = ggplot(data_chart, aes(x=cwd_vol, y=mean)) + 
  geom_point(aes(size = cnt)) + geom_hline(yintercept=0) + theme_bw() + 
  scale_x_continuous(limits=c(0,4), breaks=seq(0, 4, 0.5)) +
  scale_y_continuous(limits=c(-0.45,0.45), breaks=seq(-0.4, 0.4, 0.1)) +
  labs(x = bquote(atop('Coarse woody debris volume'~(m^3),'No CWD threshold')), y = "Mean Pearson residual", size = "Observation count") +
  theme(legend.position="bottom") + scale_size_area()

p_prog <- ggarrange(ch_preprog, ch_postprog, nrow = 1, common.legend = TRUE, legend="bottom")
#ggsave("Cg_Threshold.png", plot = p_prog, device = "png", width =9, height =6)

# R-squared calculation - code from Nakagawa and Schielzeth (2013)

# Creating a dummy variable that allows estimating additive dispersion in lmer 
Unit <- factor(1:length(data$captures))

# Fit null model without fixed effects (but including all random effects)
m0 <- glmer(captures ~ 1 + (1 | SquareID) + (1 | Unit), family = "poisson", data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=8e5)))

# Fit alternative model including fixed and all random effects
mF <- update(model_final, . ~ . + (1 | Unit))

# Calculation of the variance in fitted values
VarF <- var(as.vector(fixef(mF) %*% t(mF@pp$X)))

# R2GLMM(c) - conditional R2GLMM for full model
(VarF + VarCorr(mF)$SquareID[1])/(VarF + VarCorr(mF)$SquareID[1] + VarCorr(mF)$Unit[1] + log(1 + 1/exp(as.numeric(fixef(m0)))))
