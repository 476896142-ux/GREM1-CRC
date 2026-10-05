# 先看生存表的列名
names(surv)
head(surv)

# ===== 1. 合并生存数据 + 亚型标签 =====
surv_df <- surv %>%
  dplyr::select(sample, OS_time = OS.time, OS_status = OS) %>%
  mutate(OS_time = as.numeric(OS_time),
         OS_status = as.numeric(OS_status)) %>%
  filter(!is.na(OS_time), !is.na(OS_status)) %>%
  inner_join(subtype, by = "sample")

cat("有生存数据的肿瘤样本数:", nrow(surv_df), "\n")
table(surv_df$Subtype, surv_df$OS_status)

# ===== 2. KM 生存曲线（4 亚型） =====
library(survival)
library(survminer)

fit_km <- survfit(Surv(OS_time, OS_status) ~ Subtype, data = surv_df)

p_km <- ggsurvplot(fit_km, data = surv_df,
                   pval = TRUE, conf.int = TRUE,
                   risk.table = TRUE,
                   palette = "jco",
                   title = "OS of 4 consensus subtypes",
                   xlab = "Time (days)",
                   legend.title = "Subtype")
print(p_km)
ggsave("MR/results/KM_4subtypes.png", p_km$plot, width = 8, height = 6, dpi = 200)

# ===== 3. 单因素 Cox（亚型整体预后差异） =====
cox_uni <- coxph(Surv(OS_time, OS_status) ~ Subtype, data = surv_df)
summary(cox_uni)

# ===== 4. 以亚型3（免疫热型）为参照，看各亚型 HR =====
surv_df$Subtype_f <- factor(surv_df$Subtype, levels = c("3","1","2","4"))
cox_ref <- coxph(Surv(OS_time, OS_status) ~ Subtype_f, data = surv_df)
summary(cox_ref)

# ===== 5. 多因素 Cox（亚型 + 年龄 + 性别 + 分期） =====
names(clin)   # 确认临床表列名
clin_df <- clin %>%
  dplyr::select(sample = sampleID,
                age = age_at_initial_pathologic_diagnosis,
                gender = gender,
                stage = pathologic_stage) %>%
  mutate(age = as.numeric(age),
         stage = factor(stage))

surv_clin <- surv_df %>%
  left_join(clin_df, by = "sample") %>%
  filter(!is.na(age), !is.na(gender), !is.na(stage))

cox_multi <- coxph(Surv(OS_time, OS_status) ~ Subtype_f + age + gender + stage,
                   data = surv_clin)
summary(cox_multi)

# ===== 6. 列线图（Nomogram） =====
library(rms)
dd <- datadist(surv_clin)
options(datadist = "dd")

fit_nom <- cph(Surv(OS_time, OS_status) ~ Subtype_f + age + gender + stage,
               data = surv_clin, x = TRUE, y = TRUE, surv = TRUE)

surv_fun <- Survival(fit_nom)
nom <- nomogram(fit_nom,
                fun = list(function(x) surv_fun(1095, x),   # 3年
                           function(x) surv_fun(1825, x)),  # 5年
                funlabel = c("3-year OS", "5-year OS"),
                lp = FALSE)

png("MR/results/nomogram.png", width = 1600, height = 1000, res = 200)
plot(nom)
dev.off()

cat("\n已保存: KM_4subtypes.png, nomogram.png\n")
##############
# ===== 合并分期为 I/II/III/IV 四大类 =====
table(surv_clin$stage, useNA = "ifany")

# ===== 合并分期为 I/II/III/IV 四大类 =====
surv_clin$stage_simple <- factor(
  dplyr::case_when(
    surv_clin$stage %in% c("Stage I", "Stage IA") ~ "I",
    surv_clin$stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "II",
    surv_clin$stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "III",
    surv_clin$stage %in% c("Stage IV", "Stage IVA", "Stage IVB") ~ "IV",
    TRUE ~ NA_character_
  ),
  levels = c("I", "II", "III", "IV"))

table(surv_clin$stage_simple, useNA = "ifany")

# 去掉分期不明的样本
surv_clin2 <- surv_clin %>% filter(!is.na(stage_simple))
cat("最终多因素分析样本数:", nrow(surv_clin2), "\n")

# ===== 重新多因素 Cox =====
cox_multi2 <- coxph(Surv(OS_time, OS_status) ~ Subtype_f + age + gender + stage_simple,
                    data = surv_clin2)
summary(cox_multi2)

# ===== 列线图 =====
library(rms)
dd <- datadist(surv_clin2)
options(datadist = "dd")

fit_nom <- cph(Surv(OS_time, OS_status) ~ Subtype_f + age + gender + stage_simple,
               data = surv_clin2, x = TRUE, y = TRUE, surv = TRUE)

surv_fun <- Survival(fit_nom)
nom <- nomogram(fit_nom,
                fun = list(function(x) surv_fun(1095, x),
                           function(x) surv_fun(1825, x)),
                funlabel = c("3-year OS", "5-year OS"),
                lp = FALSE)
summary(nom)
#png("MR/results/nomogram.png", width = 1600, height = 1000, res = 200)
png("MR/results/nomogram.png", width = 2400, height = 1400, res = 300)
par(mar = c(4, 5, 3, 2))
plot(nom, cex.axis = 0.9, cex.var = 1.0)
dev.off()



cat("已保存: nomogram.png\n")


#############----------------PCA曲线------------------------########################
# 3 年 OS 校准
cal3 <- calibrate(fit_nom, cmethod = "KM", method = "boot",
                  u = 3 * 365.25, B = 1000)
cal3

# 5 年 OS 校准
cal5 <- calibrate(fit_nom, cmethod = "KM", method = "boot",
                  u = 5 * 365.25, B = 1000)
cal5

png("MR/results/calibration_3yr.png", width=1200, height=1200, res=300)
plot(cal3, xlab = "Predicted 3-year OS", ylab = "Observed survival",
     main = "Calibration at 3 years")
dev.off()

png("MR/results/calibration_5yr.png", width=1200, height=1200, res=300)
plot(cal5, xlab = "Predicted 5-year OS", ylab = "Observed survival",
     main = "Calibration at 5 years")
dev.off()


