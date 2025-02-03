########## infection prevalence
summarytools::freq(child$cut_BK2)
summarytools::freq(child$cut_JC2)
summarytools::freq(child$cut_KI2)
summarytools::freq(child$cut_WU2)
summarytools::freq(child$cut_MCV2)
summarytools::freq(child$EBV_class2)
summarytools::freq(child$CMV_class2)
summarytools::freq(child$HP_class2)
summarytools::freq(child$cut_Avd362)
summarytools::freq(child$cut_Toxo2)

## infection prevalence age 4 and 6 
age6 <- baseline[, c("cut_BK3","cut_JC3","cut_KI3","cut_WU3",
                     "cut_MCV3","EBV_class3","CMV_class3","HP_class3",
                     "cut_Avd363","cut_Toxo3", "childid")]

age11 <- baseline[, c("cut_BK4","cut_JC4","cut_KI4","cut_WU4",
                      "cut_MCV4","EBV_class4","CMV_class4","HP_class4",
                      "cut_Avd364","cut_Toxo4", "childid")]

age6_inf<-rename(child)
age11_inf<-rename(child)

age6_inf_1<-merge(age6_inf,age6, by ="childid")
age11_inf_1<-merge(age11_inf,age11, by ="childid")

age6_inf_1 <- age6_inf_1 %>%
  drop_na(cut_BK3, cut_JC3, cut_KI3, cut_WU3, cut_MCV3,
          EBV_class3, CMV_class3, HP_class3, cut_Avd363, cut_Toxo3)

age11_inf_1 <- age11_inf_1 %>%
  drop_na(cut_BK4, cut_JC4, cut_KI4, cut_WU4, cut_MCV4,
          EBV_class4, CMV_class4, HP_class4, cut_Avd364, cut_Toxo4)

summarytools::freq(age6_inf_1$cut_BK3)
summarytools::freq(age6_inf_1$cut_JC3)
summarytools::freq(age6_inf_1$cut_KI3)
summarytools::freq(age6_inf_1$cut_WU3)
summarytools::freq(age6_inf_1$cut_MCV3)
summarytools::freq(age6_inf_1$EBV_class3)
summarytools::freq(age6_inf_1$CMV_class3)
summarytools::freq(age6_inf_1$HP_class3)
summarytools::freq(age6_inf_1$cut_Avd363)
summarytools::freq(age6_inf_1$cut_Toxo3)

summarytools::freq(age11_inf_1$cut_BK4)
summarytools::freq(age11_inf_1$cut_JC4)
summarytools::freq(age11_inf_1$cut_KI4)
summarytools::freq(age11_inf_1$cut_WU4)
summarytools::freq(age11_inf_1$cut_MCV4)
summarytools::freq(age11_inf_1$EBV_class4)
summarytools::freq(age11_inf_1$CMV_class4)
summarytools::freq(age11_inf_1$HP_class4)
summarytools::freq(age11_inf_1$cut_Avd364)
summarytools::freq(age11_inf_1$cut_Toxo4)

### average time at risk
risk<-rename(all_antigens_survival_vars_child_20240805)

maternal$cut_BK
child$time_BK
risk<-rename(child)

time_risk_BK<-mean(risk$time_BK*risk$cut_BK)
time_risk_BK

time_risk_JC<-mean(risk$time_JC*risk$cut_JC)
time_risk_JC

time_risk_KI<-mean(risk$time_KI*risk$cut_KI)
time_risk_KI

time_risk_WU<-mean(risk$time_WU*risk$cut_WU)
time_risk_WU

time_risk_MCV<-mean(risk$time_MCV*risk$cut_MCV)
time_risk_MCV

time_risk_EBV_class<-mean(risk$time_EBV_class*risk$cut_EBV_class)
time_risk_EBV_class

time_risk_CMV_class<-mean(risk$time_CMV_class*risk$cut_CMV_class)
time_risk_CMV_class

time_risk_HP_class<-mean(risk$time_HP_class*risk$cut_HP_class)
time_risk_HP_class

time_risk_AVD<-mean(risk$time_Avd36*risk$cut_Avd36)
time_risk_AVD

time_risk_toxo<-mean(risk$time_Toxo*risk$cut_Toxo)
time_risk_toxo
###### crosstables
chisq.test(child$y4_totalsleep,child$y4_Nightmares)
chisq.test(child$y4_totalsleep,child$y4_Snoring)
chisq.test(child$y4_Nightmares,child$y4_Snoring)

#######################################################################  child
child<-rename(child_27122024)
colnames(child)
child$y4_totalsleep <- factor(child$y4_totalsleep)
child$y4_totalsleep<-relevel(child$y4_totalsleep, ref = "0")
##bk
bk<-coxph(Surv(time_BK, cut_BK)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y +sex , data = child)
summary(bk)
##jc
jc<-coxph(Surv(time_JC, cut_JC)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(jc)
###ki
ki<-coxph(Surv(time_KI, cut_KI)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(ki)
#### WU
wu<-coxph(Surv(time_WU, cut_WU)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(wu)
### mcv
mcv<-coxph(Surv(time_MCV, cut_MCV)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(mcv)
################# herpes
ebv<-coxph(Surv(time_EBV_class, cut_EBV_class)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(ebv)
cmv<-coxph(Surv(time_CMV_class, cut_CMV_class)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(cmv)
######## hp 
hp<-coxph(Surv(time_HP_class, cut_HP_class)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(hp)
######### avd 
avd<-coxph(Surv(time_Avd36, cut_Avd36)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(avd)
######### toxo
toxo<-coxph(Surv(time_Toxo, cut_Toxo)~ as.factor(y4_totalsleep) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(toxo)

extract_cox_results <- function(model, model_name) {
  tidy_model <- tidy(model, exponentiate = TRUE)
  conf_intervals <- confint(model)
  exp_conf_intervals <- exp(conf_intervals)
  results <- tidy_model %>%
    mutate(
      `lower.95` = exp_conf_intervals[, 1],
      `upper.95` = exp_conf_intervals[, 2],
      model = model_name
    ) %>%
    select(model, term, estimate, lower.95, upper.95, p.value) %>%
    rename(
      `exp(coef)` = estimate,
      `Pr(>|z|)` = p.value
    )
  return(results)
}

results_bk <- extract_cox_results(bk, "bk")
results_jc <- extract_cox_results(jc, "jc")
results_ki <- extract_cox_results(ki, "ki")
results_wu <- extract_cox_results(wu, "wu")
results_mcv <- extract_cox_results(mcv, "mcv")
results_ebv <- extract_cox_results(ebv, "ebv")
results_cmv <- extract_cox_results(cmv, "cmv")
results_hp <- extract_cox_results(hp, "hp")
results_avd <- extract_cox_results(avd, "avd")
results_toxo <- extract_cox_results(toxo, "toxo")

combined_results <- bind_rows(results_bk, results_jc, results_ki, results_wu,results_mcv,results_ebv,
                              results_cmv,results_hp,results_avd,results_toxo)

write.xlsx(combined_results, file = "C:/Users/rgalan/Desktop/y4_totalsleep_cox.xlsx", rowNames = FALSE)
####################################################################################
#################################################################################
######################################################################################
child$y4_Snoring <- factor(child$y4_Snoring)
child$y4_Snoring<-relevel(child$y4_Snoring, ref = "0")
##bk
bk<-coxph(Surv(time_BK, cut_BK)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y +sex , data = child)
summary(bk)
##jc
jc<-coxph(Surv(time_JC, cut_JC)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(jc)
###ki
ki<-coxph(Surv(time_KI, cut_KI)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(ki)
#### WU
wu<-coxph(Surv(time_WU, cut_WU)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(wu)
### mcv
mcv<-coxph(Surv(time_MCV, cut_MCV)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(mcv)
################# herpes
ebv<-coxph(Surv(time_EBV_class, cut_EBV_class)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(ebv)
cmv<-coxph(Surv(time_CMV_class, cut_CMV_class)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(cmv)
######## hp 
hp<-coxph(Surv(time_HP_class, cut_HP_class)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(hp)
######### avd 
avd<-coxph(Surv(time_Avd36, cut_Avd36)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(avd)
######### toxo
toxo<-coxph(Surv(time_Toxo, cut_Toxo)~ as.factor(y4_Snoring) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(toxo)


extract_cox_results <- function(model, model_name) {
  tidy_model <- tidy(model, exponentiate = TRUE)
  conf_intervals <- confint(model)
  exp_conf_intervals <- exp(conf_intervals)
  results <- tidy_model %>%
    mutate(
      `lower.95` = exp_conf_intervals[, 1],
      `upper.95` = exp_conf_intervals[, 2],
      model = model_name
    ) %>%
    select(model, term, estimate, lower.95, upper.95, p.value) %>%
    rename(
      `exp(coef)` = estimate,
      `Pr(>|z|)` = p.value
    )
  return(results)
}

results_bk <- extract_cox_results(bk, "bk")
results_jc <- extract_cox_results(jc, "jc")
results_ki <- extract_cox_results(ki, "ki")
results_wu <- extract_cox_results(wu, "wu")
results_mcv <- extract_cox_results(mcv, "mcv")
results_ebv <- extract_cox_results(ebv, "ebv")
results_cmv <- extract_cox_results(cmv, "cmv")
results_hp <- extract_cox_results(hp, "hp")
results_avd <- extract_cox_results(avd, "avd")
results_toxo <- extract_cox_results(toxo, "toxo")

combined_results <- bind_rows(results_bk, results_jc, results_ki, results_wu,results_mcv,results_ebv,
                              results_cmv,results_hp,results_avd,results_toxo)

write.xlsx(combined_results, file = "C:/Users/rgalan/Desktop/y4_Snoring_cox.xlsx", rowNames = FALSE)
###############################################
###############################################
#################################################
child$y4_Nightmares <- factor(child$y4_Nightmares)
child$y4_Nightmares<-relevel(child$y4_Nightmares, ref = "0")
##bk
bk<-coxph(Surv(time_BK, cut_BK)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y +sex , data = child)
summary(bk)
##jc
jc<-coxph(Surv(time_JC, cut_JC)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(jc)
###ki
ki<-coxph(Surv(time_KI, cut_KI)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(ki)
#### WU
wu<-coxph(Surv(time_WU, cut_WU)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(wu)
### mcv
mcv<-coxph(Surv(time_MCV, cut_MCV)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(mcv)
################# herpes
ebv<-coxph(Surv(time_EBV_class, cut_EBV_class)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(ebv)
cmv<-coxph(Surv(time_CMV_class, cut_CMV_class)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(cmv)
######## hp 
hp<-coxph(Surv(time_HP_class, cut_HP_class)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(hp)
######### avd 
avd<-coxph(Surv(time_Avd36, cut_Avd36)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(avd)
######### toxo
toxo<-coxph(Surv(time_Toxo, cut_Toxo)~ as.factor(y4_Nightmares) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(toxo)


extract_cox_results <- function(model, model_name) {
  tidy_model <- tidy(model, exponentiate = TRUE)
  conf_intervals <- confint(model)
  exp_conf_intervals <- exp(conf_intervals)
  results <- tidy_model %>%
    mutate(
      `lower.95` = exp_conf_intervals[, 1],
      `upper.95` = exp_conf_intervals[, 2],
      model = model_name
    ) %>%
    select(model, term, estimate, lower.95, upper.95, p.value) %>%
    rename(
      `exp(coef)` = estimate,
      `Pr(>|z|)` = p.value
    )
  return(results)
}

results_bk <- extract_cox_results(bk, "bk")
results_jc <- extract_cox_results(jc, "jc")
results_ki <- extract_cox_results(ki, "ki")
results_wu <- extract_cox_results(wu, "wu")
results_mcv <- extract_cox_results(mcv, "mcv")
results_ebv <- extract_cox_results(ebv, "ebv")
results_cmv <- extract_cox_results(cmv, "cmv")
results_hp <- extract_cox_results(hp, "hp")
results_avd <- extract_cox_results(avd, "avd")
results_toxo <- extract_cox_results(toxo, "toxo")

combined_results <- bind_rows(results_bk, results_jc, results_ki, results_wu,results_mcv,results_ebv,
                              results_cmv,results_hp,results_avd,results_toxo)

write.xlsx(combined_results, file = "C:/Users/rgalan/Desktop/y4_Nightmares_cox.xlsx", rowNames = FALSE)
###############################################
###############################################
#################################################
child$sleepscore <- factor(child$sleepscore)
child$sleepscore<-relevel(child$sleepscore, ref = "0")
##bk
bk<-coxph(Surv(time_BK, cut_BK)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y +sex , data = child)
summary(bk)
##jc
jc<-coxph(Surv(time_JC, cut_JC)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(jc)
###ki
ki<-coxph(Surv(time_KI, cut_KI)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(ki)
#### WU
wu<-coxph(Surv(time_WU, cut_WU)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(wu)
### mcv
mcv<-coxph(Surv(time_MCV, cut_MCV)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(mcv)
################# herpes
ebv<-coxph(Surv(time_EBV_class, cut_EBV_class)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(ebv)
cmv<-coxph(Surv(time_CMV_class, cut_CMV_class)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(cmv)
######## hp 
hp<-coxph(Surv(time_HP_class, cut_HP_class)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(hp)
######### avd 
avd<-coxph(Surv(time_Avd36, cut_Avd36)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(avd)
######### toxo
toxo<-coxph(Surv(time_Toxo, cut_Toxo)~ as.factor(sleepscore) + y4_zbmi +asthmadd4y +room_sharing4y +num_children4y+sex , data = child)
summary(toxo)


extract_cox_results <- function(model, model_name) {
  tidy_model <- tidy(model, exponentiate = TRUE)
  conf_intervals <- confint(model)
  exp_conf_intervals <- exp(conf_intervals)
  results <- tidy_model %>%
    mutate(
      `lower.95` = exp_conf_intervals[, 1],
      `upper.95` = exp_conf_intervals[, 2],
      model = model_name
    ) %>%
    select(model, term, estimate, lower.95, upper.95, p.value) %>%
    rename(
      `exp(coef)` = estimate,
      `Pr(>|z|)` = p.value
    )
  return(results)
}

results_bk <- extract_cox_results(bk, "bk")
results_jc <- extract_cox_results(jc, "jc")
results_ki <- extract_cox_results(ki, "ki")
results_wu <- extract_cox_results(wu, "wu")
results_mcv <- extract_cox_results(mcv, "mcv")
results_ebv <- extract_cox_results(ebv, "ebv")
results_cmv <- extract_cox_results(cmv, "cmv")
results_hp <- extract_cox_results(hp, "hp")
results_avd <- extract_cox_results(avd, "avd")
results_toxo <- extract_cox_results(toxo, "toxo")

combined_results <- bind_rows(results_bk, results_jc, results_ki, results_wu,results_mcv,results_ebv,
                              results_cmv,results_hp,results_avd,results_toxo)

write.xlsx(combined_results, file = "C:/Users/rgalan/Desktop/sleepscore_cox.xlsx", rowNames = FALSE)
##############################################################
#########################################################
############################################################
################################################################ forest plots 
####
install.packages("gridExtra")
library(gridExtra)
totalsleep<-rename(y4_totalsleep_cox)
totalsleep$color <- ifelse(totalsleep$p_value < 0.05, "red", "black")
totalsleep$Ag <- factor(totalsleep$Ag, levels = totalsleep$Ag)
totalsleep$Ag<-factor(totalsleep$Ag, levels = rev(totalsleep$Ag))
forest_plot_totalsleep <- ggplot(totalsleep, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange(aes(color = color), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  coord_flip() +
  labs(x = "Ag", y = "HR") +
  ggtitle("a")+
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0) 
  )
print(forest_plot_totalsleep)
totalsleep_table <- totalsleep[, c("HR", "low", "high", "p_value")]
totalsleep_table <- totalsleep_table %>%
  mutate(
    HR = round(HR, 2),
    low = round(low, 2),
    high = round(high, 2),
    p_value = round(p_value, 2)
  )
tabla_totalsleep <- tableGrob(totalsleep_table, rows = NULL, theme = ttheme_minimal(base_size = 10))
empty_plot <- ggplot() + 
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
combined_plot_totalsleep <- plot_grid(forest_plot_totalsleep, tabla_totalsleep, ncol = 2, align = 'h', axis = 'tb', rel_widths = c(2/3, 1/3))
print(combined_plot_totalsleep)
######
snoring<-rename(y4_Snoring_cox)
snoring$color <- ifelse(snoring$p_value < 0.05, "red", "black")
snoring$Ag <- factor(snoring$Ag, levels = snoring$Ag)
snoring$Ag<-factor(snoring$Ag, levels = rev(snoring$Ag))
forest_plot_snoring <- ggplot(snoring, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange(aes(color = color), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  coord_flip() +
  labs(x = "Ag", y = "HR") +
  ggtitle("b")+
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0) 
  )
print(forest_plot_snoring)
snoring_table <- snoring[, c("HR", "low", "high", "p_value")]
snoring_table <- snoring_table %>%
  mutate(
    HR = round(HR, 2),
    low = round(low, 2),
    high = round(high, 2),
    p_value = round(p_value, 2)
  )
tabla_snoring <- tableGrob(snoring_table, rows = NULL, theme = ttheme_minimal(base_size = 10))
empty_plot <- ggplot() + 
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
combined_plot_snoring <- plot_grid(forest_plot_snoring, tabla_snoring, ncol = 2, align = 'h', axis = 'tb', rel_widths = c(2/3, 1/3))
print(combined_plot_snoring)
##############
nightmares<-rename(y4_Nightmares_cox)
nightmares$color <- ifelse(nightmares$p_value < 0.05, "red", "black")
nightmares$Ag <- factor(nightmares$Ag, levels = nightmares$Ag)
nightmares$Ag<-factor(nightmares$Ag, levels = rev(nightmares$Ag))
forest_plot_nightmares <- ggplot(nightmares, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange(aes(color = color), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  coord_flip() +
  labs(x = "Ag", y = "HR") +
  ggtitle("c")+
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0) 
  )
print(forest_plot_nightmares)
nightmares_table <- nightmares[, c("HR", "low", "high", "p_value")]
nightmares_table <- nightmares_table %>%
  mutate(
    HR = round(HR, 2),
    low = round(low, 2),
    high = round(high, 2),
    p_value = round(p_value, 2)
  )
tabla_nightmares <- tableGrob(nightmares_table, rows = NULL, theme = ttheme_minimal(base_size = 10))
empty_plot <- ggplot() + 
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
combined_plot_nightmares <- plot_grid(forest_plot_nightmares, tabla_nightmares, ncol = 2, align = 'h', axis = 'tb', rel_widths = c(2/3, 1/3))
print(combined_plot_nightmares)


combined_outcomes_child<- plot_grid(
  combined_plot_totalsleep, 
  combined_plot_snoring, 
  combined_plot_nightmares, 
  ncol = 1, 
  nrow = 3, 
  align = 'v',
  rel_widths = c(1, 1)
)
print(combined_outcomes_child)
########
sleepscore<-rename(sleep_score_daysleep_cat)
sleepscore$color <- ifelse(sleepscore$p_value < 0.05, "red", "black")
sleepscore$Ag <- factor(sleepscore$Ag, levels = sleepscore$Ag)
sleepscore$Ag<-factor(sleepscore$Ag, levels = rev(sleepscore$Ag))
forest_plot_sleepscore <- ggplot(sleepscore, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange(aes(color = color), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  coord_flip() +
  labs(x = "Ag", y = "HR") +
  ggtitle("")+
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0) 
  )
print(forest_plot_sleepscore)
sleep_score_daysleep_cat$`p_value for trend`
sleepscore_table <- sleepscore[, c("HR", "low", "high", "p_value","p_value for trend")]
sleepscore_table <- sleepscore_table %>%
  mutate(
    HR = round(HR, 2),
    low = round(low, 2),
    high = round(high, 2),
    p_value = round(p_value, 2),
    `p_value for trend` = round(`p_value for trend`,2)
  )
library(gridExtra)
tabla_sleepscore <- tableGrob(sleepscore_table, rows = NULL, theme = ttheme_minimal(base_size = 10))
empty_plot <- ggplot() + 
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
combined_plot_sleepscore <- plot_grid(forest_plot_sleepscore, tabla_sleepscore, ncol = 2, align = 'h', axis = 'tb', rel_widths = c(2/3, 1/3))
print(combined_plot_sleepscore)
######################################### maternal 
sleepscore_maternal<-rename(sleep_score_cat_adjusted_smoking)
sleepscore_maternal$color <- ifelse(sleepscore_maternal$p_value < 0.05, "red", "black")
sleepscore_maternal$Ag <- factor(sleepscore_maternal$Ag, levels = sleepscore_maternal$Ag)
sleepscore_maternal$Ag<-factor(sleepscore_maternal$Ag, levels = rev(sleepscore_maternal$Ag))
forest_plot_sleepscore_maternal <- ggplot(sleepscore_maternal, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange(aes(color = color), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  coord_flip() +
  labs(x = "Ag", y = "HR") +
  ggtitle("")+
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0) 
  )
print(forest_plot_sleepscore_maternal)
sleepscore_maternal$`p_value for trend`
sleepscore_maternal_table <- sleepscore_maternal[, c("HR", "low", "high", "p_value","p_value for trend")]
sleepscore_maternal_table <- sleepscore_maternal_table %>%
  mutate(
    HR = round(HR, 2),
    low = round(low, 2),
    high = round(high, 2),
    p_value = round(p_value, 2),
    `p_value for trend` = round(`p_value for trend`,2)
  )
tabla_sleepscore_maternal <- tableGrob(sleepscore_maternal_table, rows = NULL, theme = ttheme_minimal(base_size = 10))
empty_plot <- ggplot() + 
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
combined_plot_sleepscore_maternal <- plot_grid(forest_plot_sleepscore_maternal, tabla_sleepscore_maternal, ncol = 2, align = 'h', axis = 'tb', rel_widths = c(2/3, 1/3))
print(combined_plot_sleepscore_maternal)

maternal$sleep_score<-factor(maternal$sleep_score)
maternal$sleep_score<-relevel(maternal$sleep_score, ref = "0")

##bk
bk<-coxph(Surv(time_BK, cut_BK)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(bk)
##jc
jc<-coxph(Surv(time_JC, cut_JC)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(jc)
###ki
ki<-coxph(Surv(time_KI, cut_KI)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(ki)
#### WU
wu<-coxph(Surv(time_WU, cut_WU)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(wu)
### mcv
mcv<-coxph(Surv(time_MCV, cut_MCV)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(mcv)

################# herpes
ebv<-coxph(Surv(time_EBV_class, cut_EBV_class)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(ebv)


cmv<-coxph(Surv(time_CMV_class, cut_CMV_class)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(cmv)

######## hp 
hp<-coxph(Surv(time_HP_class, cut_HP_class)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre +m_smk30, data = maternal)
summary(hp)

######### avd 
avd<-coxph(Surv(time_Avd36, cut_Avd36)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(avd)
######### toxo
toxo<-coxph(Surv(time_Toxo, cut_Toxo)~ sleep_score + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(toxo)

######################################################
maternal$w30_Epworth <- factor(maternal$w30_Epworth)
maternal$w30_Epworth<-relevel(maternal$w30_Epworth, ref = "0")


##bk
bk<-coxph(Surv(time_BK, cut_BK)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(bk)
##jc
jc<-coxph(Surv(time_JC, cut_JC)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(jc)
###ki
ki<-coxph(Surv(time_KI, cut_KI)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre +m_smk30, data = maternal)
summary(ki)
#### WU
wu<-coxph(Surv(time_WU, cut_WU)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(wu)
### mcv
mcv<-coxph(Surv(time_MCV, cut_MCV)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(mcv)

################# herpes
ebv<-coxph(Surv(time_EBV_class, cut_EBV_class)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(ebv)

cmv<-coxph(Surv(time_CMV_class, cut_CMV_class)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(cmv)

######## hp 
hp<-coxph(Surv(time_HP_class, cut_HP_class)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(hp)

######### avd 
avd<-coxph(Surv(time_Avd36, cut_Avd36)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre  +m_smk30, data = maternal)
summary(avd)
######### toxo
toxo<-coxph(Surv(time_Toxo, cut_Toxo)~ as.factor(w30_Epworth) + mage + m_educ+marital_pgn+bmi_pre  +m_smk30,data = maternal)
summary(toxo)


############################################################################################################## sleep disturbances
##############################################################################################################
maternal$w30_sleepdisturb <- factor(maternal$w30_sleepdisturb)
maternal$w30_sleepdisturb<-relevel(maternal$w30_sleepdisturb, ref = "0")

##bk
bk<-coxph(Surv(time_BK, cut_BK)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(bk)
##jc
jc<-coxph(Surv(time_JC, cut_JC)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(jc)
###ki
ki<-coxph(Surv(time_KI, cut_KI)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(ki)
#### WU
wu<-coxph(Surv(time_WU, cut_WU)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(wu)
### mcv
mcv<-coxph(Surv(time_MCV, cut_MCV)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(mcv)

################# herpes
ebv<-coxph(Surv(time_EBV_class, cut_EBV_class)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(ebv)


cmv<-coxph(Surv(time_CMV_class, cut_CMV_class)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre +m_smk30, data = maternal)
summary(cmv)

######## hp 
hp<-coxph(Surv(time_HP_class, cut_HP_class)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre +m_smk30, data = maternal)
summary(hp)

######### avd 
avd<-coxph(Surv(time_Avd36, cut_Avd36)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(avd)
######### toxo
toxo<-coxph(Surv(time_Toxo, cut_Toxo)~ as.factor(w30_sleepdisturb) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(toxo)

###################################################
####################################################### sleep duration 
maternal$w30_HoursSleeping<-factor(maternal$w30_HoursSleeping)
maternal$w30_HoursSleeping<-relevel(maternal$w30_HoursSleeping, ref = "0")

##bk
bk<-coxph(Surv(time_BK, cut_BK)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(bk)
##jc
jc<-coxph(Surv(time_JC, cut_JC)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(jc)
###ki
ki<-coxph(Surv(time_KI, cut_KI)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(ki)
#### WU
wu<-coxph(Surv(time_WU, cut_WU)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(wu)
### mcv
mcv<-coxph(Surv(time_MCV, cut_MCV)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(mcv)

################# herpes
ebv<-coxph(Surv(time_EBV_class, cut_EBV_class)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(ebv)


cmv<-coxph(Surv(time_CMV_class, cut_CMV_class)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(cmv)

######## hp 
hp<-coxph(Surv(time_HP_class, cut_HP_class)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(hp)

######### avd 
avd<-coxph(Surv(time_Avd36, cut_Avd36)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(avd)
######### toxo
toxo<-coxph(Surv(time_Toxo, cut_Toxo)~ as.factor(w30_HoursSleeping) + mage + m_educ+marital_pgn+bmi_pre +m_smk30, data = maternal)
summary(toxo)

###################################################
####################################################### snoring freq
maternal$Snoring<-factor(maternal$Snoring)
maternal$Snoring<-relevel(maternal$Snoring, ref = "0")

##bk
bk<-coxph(Surv(time_BK, cut_BK)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(bk)
##jc
jc<-coxph(Surv(time_JC, cut_JC)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(jc)
###ki
ki<-coxph(Surv(time_KI, cut_KI)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(ki)
#### WU
wu<-coxph(Surv(time_WU, cut_WU)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(wu)
### mcv
mcv<-coxph(Surv(time_MCV, cut_MCV)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(mcv)

################# herpes
ebv<-coxph(Surv(time_EBV_class, cut_EBV_class)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(ebv)


cmv<-coxph(Surv(time_CMV_class, cut_CMV_class)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(cmv)

######## hp 
hp<-coxph(Surv(time_HP_class, cut_HP_class)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(hp)

######### avd 
avd<-coxph(Surv(time_Avd36, cut_Avd36)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(avd)
######### toxo
toxo<-coxph(Surv(time_Toxo, cut_Toxo)~ as.factor(Snoring) + mage + m_educ+marital_pgn+bmi_pre+m_smk30 , data = maternal)
summary(toxo)


####### rename
w30<-rename(w30_adjusted_smoking)
disturb<-rename(sleep_disturb_adjusted_smoking)
duration<-rename(sleep_duration_adjusted_smoking)
snoring<-rename(snoring_adjusted_smoking)

######################################## w30

w30$Ag <- factor(w30$Ag, levels = w30$Ag)
w30$Ag<-factor(w30$Ag, levels = rev(w30$Ag))

forest_plot_w30 <- ggplot(w30, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "solid") +
  coord_flip() +
  labs(x = "Ag", y = "HR") +
  ggtitle("a")+
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0) 
  )

print(forest_plot_w30)


################################################ sleep disturb


disturb$Ag <- factor(disturb$Ag, levels = disturb$Ag)
disturb$Ag<-factor(disturb$Ag, levels = rev(disturb$Ag))
disturb$color <- ifelse(disturb$`p-value` < 0.05, "red", "black")

forest_plot_disturb <- ggplot(disturb, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange(aes(color = color), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  coord_flip() +
  labs(x = "Ag", y = "HR") +
  ggtitle("b")+
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0) 
  )

print(forest_plot_disturb)

################################################## sleep duration 

duration$Ag <- factor(duration$Ag, levels = duration$Ag)
duration$Ag<-factor(duration$Ag, levels = rev(duration$Ag))
duration$color <- ifelse(duration$`p-value` < 0.05, "red", "black")

forest_plot_duration <- ggplot(duration, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange(aes(color = color), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  coord_flip() +
  labs(x = "Ag", y = "HR") +
  ggtitle("c")+
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0) 
  )

print(forest_plot_duration)

###################################################### snoring
snoring<-rename(snoring_results)
snoring$Ag <- factor(snoring$Ag, levels = snoring$Ag)
snoring$Ag<-factor(snoring$Ag, levels = rev(snoring$Ag))
snoring$color <- ifelse(snoring$`p-value` < 0.05, "red", "black")

forest_plot_snoring <- ggplot(snoring, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange(aes(color = color), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  coord_flip() +
  labs(x = "Ag", y = "HR") +
  ggtitle("d")+
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0) 
  )

print(forest_plot_snoring)

####################################################### 
sleep_cat$Ag <- factor(sleep_cat$Ag, levels = sleep_cat$Ag)
sleep_cat$Ag<-factor(sleep_cat$Ag, levels = rev(sleep_cat$Ag))
sleep_cat$color <- ifelse(sleep_cat$p_value < 0.05, "red", "black")

forest_plot_sleep_cat <- ggplot(sleep_cat, aes(x = Ag, y = HR, ymin = low, ymax = high)) +
  geom_pointrange(aes(color = color), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  coord_flip() +
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5, lineheight = 0.9),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    plot.title = element_text(hjust = 0)  
  )

print(forest_plot_sleep_cat)

#########################################3
maternal<-rename(maternal_29122024)
summarytools::freq(maternal$infection_score)

excessive_lm <- lm(infection_score ~ w30_Epworth + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(excessive_lm)
confint(excessive_lm)

disturb_lm <- lm(infection_score ~ w30_sleepdisturb + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(disturb_lm)
confint(disturb_lm)

hours_lm <- lm(infection_score ~ w30_HoursSleeping + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(hours_lm)
confint(hours_lm)

snoring_lm <- lm(infection_score ~ Snoring + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(snoring_lm)
confint(snoring_lm)

sleep_lm <- lm(infection_score ~ sleep_score + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(sleep_lm)
confint(sleep_lm)
###########################
################################# diagnostics cox
bk_test<- cox.zph(bk)
bk_test
ggcoxzph(bk_test)


jc_test<- cox.zph(jc)
jc_test
ggcoxzph(jc_test)

ki_test<- cox.zph(ki)
ki_test
ggcoxzph(ki_test)

wu_test<- cox.zph(wu)
wu_test
ggcoxzph(wu_test)

mcv_test<- cox.zph(mcv)
mcv_test
ggcoxzph(mcv_test)

ebv_test<- cox.zph(ebv)
ebv_test
ggcoxzph(ebv_test)

cmv_test<- cox.zph(cmv)
cmv_test
ggcoxzph(cmv_test)

hp_test<- cox.zph(hp)
hp_test
ggcoxzph(hp_test)

avd_test<- cox.zph(avd)
avd_test
ggcoxzph(avd_test)

toxo_test<- cox.zph(toxo)
toxo_test
ggcoxzph(toxo_test)

child$y4_NightSleepDuration
sleep_lm <- lm(infection_score ~ y4_NightSleepDuration + y4_zbmi + asthmadd4y + 
                 room_sharing4y + num_children4y + sex, data = child)
summary(totalsleep_lm)
confint(totalsleep_lm)

snoring_lm <- lm(infection_score ~ y4_Snoring + y4_zbmi + asthmadd4y + 
                   room_sharing4y + num_children4y + sex, data = child)
summary(snoring_lm)
confint(snoring_lm)

nightmares_lm <- lm(infection_score ~ y4_Nightmares + y4_zbmi + asthmadd4y + 
                      room_sharing4y + num_children4y + sex, data = child)
summary(nightmares_lm)
confint(nightmares_lm)

sleepscore_lm <- lm(infection_score ~ sleep_score + y4_zbmi + asthmadd4y + 
                      room_sharing4y + num_children4y + sex, data = child)
summary(sleepscore_lm)
confint(sleepscore_lm)

maternal<-rename(maternal_29122024)
summarytools::freq(maternal$infection_score)

excessive_lm <- lm(infection_score ~ w30_Epworth + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(excessive_lm)
confint(excessive_lm)

disturb_lm <- lm(infection_score ~ w30_sleepdisturb + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(disturb_lm)
confint(disturb_lm)

hours_lm <- lm(infection_score ~ w30_HoursSleeping + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(hours_lm)
confint(hours_lm)

snoring_lm <- lm(infection_score ~ Snoring + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(snoring_lm)
confint(snoring_lm)

sleep_lm <- lm(infection_score ~ sleep_score + mage + m_educ + marital_pgn + bmi_pre, data = maternal)
summary(sleep_lm)
confint(sleep_lm)


######################## pvalue for trend sleep child and maternal 
maternal$cut_HP_class
maternal$sleep_score<-as.numeric(maternal$sleep_score)
model <- glm(cut_HP_class ~ sleep_score, data = maternal, family = binomial)
summary(model)
p_value_trend <- summary(model)$coefficients["sleep_score", "Pr(>|z|)"]
cat("P-value for trend:", p_value_trend, "\n")

child<-rename(child_03072024)
child<-merge(child, all_antigens_survival_vars_child_20240805, by = "childid")
child$sleep_score_daysleep
child$cut_Toxo

model <- glm( cut_Toxo~ sleep_score_daysleep, data = child, family = binomial)
summary(model)
p_value_trend <- summary(model)$coefficients["sleep_score_daysleep", "Pr(>|z|)"]
cat("P-value for trend:", p_value_trend, "\n")
