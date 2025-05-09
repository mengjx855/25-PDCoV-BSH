#### Jinxin Meng, 20250414, 20250415 ####
setwd("F:/project/20250313_PDCoV_BAC_Bile_Zhangyq/data_available/statistics/Experimental-Data//")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, rstatix, openxlsx)
library(lme4)
library(lmerTest)
library(emmeans)

#### TCID50 ####
data <- read.xlsx('data.xlsx', sheet = 1)

data_x <- data %>% gather('loc', 'value', -animal, -group)

model <- lmer(value ~ group * loc + (1 | animal), data = data_x)
summary(model)

shapiro.test(residuals(model))

emmeans(model, pairwise ~ group | loc)

#### 病毒载量 ####
data <- read.xlsx('data.xlsx', sheet = 2)

data_x <- data %>% gather('loc', 'value', -animal, -group)

model <- lmer(value ~ group * loc + (1 | animal), data = data_x)
summary(model)

shapiro.test(residuals(model))

emmeans(model, pairwise ~ group | loc)

#### 体重 ####
data <- read.xlsx('data.xlsx', sheet = 3)

data_x <- data %>% gather('day', 'value', -animal, -group)

model <- lmer(value ~ group * day + (day | animal), data = data_x)
summary(model)

shapiro.test(residuals(model))

emmeans(model, pairwise ~ group | day)

#### 粪便排毒 ####
data <- read.xlsx('data.xlsx', sheet = 4)

data_x <- data %>% gather('day', 'value', -animal, -group) %>%
  mutate(rank = rank(value),
         DPI = as.numeric(gsub('D', '', day)))

# 数据不符合正态分布
data_x$value %>% shapiro.test()
data_x$value %>% log %>% shapiro.test()

# 基于Huber的稳健回归
library(robustlmm)
robust_model <- rlmer(value ~ group * day + (1|animal), 
                      data = data_x, method = "DASvar")

residuals(robust_model) %>% shapiro.test()

# Gamma分布广义线性模型
library(glmmTMB)
gamma_model <- glmmTMB(value ~ group * day + (1 | animal), data = data_x,
                       family = Gamma(link = "log"))
residuals(gamma_model) %>% shapiro.test()

# 普通混合效应模型
model <- lmer(value ~ group * day + (1 | animal), data = data_x)
summary(model)

shapiro.test(residuals(model))
qqnorm(residuals(model))
qqline(residuals(model))
plot(model, which = 1)

emmeans(model, pairwise ~ group | day)

# 重复测量多因素方差检验
model <- aov_ez(id = "animal",   # 受试者编号
                dv = "value",        # 因变量（血压）
                within = "day",  # 时间点（重复测量因素）
                between = "group", # 药物类型（组间因素）
                data = data_x)
summary(model)

shapiro.test(residuals(model))

emmeans(model, pairwise ~ group | day, adjust = 'bonferroni')
