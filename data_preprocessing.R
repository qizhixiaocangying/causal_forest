# 清理环境
rm(list = ls())

# 包管理（自动安装并加载所需包）
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, survey, mice, grf, ggplot2, car, dplyr)

# 读取数据
load("nhanes_merged.rda")
#load("NHANES_20172018.rda")

# --- 数据预处理与变量构建 ---
nhanes_processed <- nhanes_merged %>%
  # 1. 仅保留20岁及以上，排除孕妇
  filter(RIDAGEYR >= 20, RIDEXPRG != 1 | is.na(RIDEXPRG)) %>%
  # 2. 构建权重
  mutate(WGT4YR = WTMEC2YR / 2) %>%
  filter(!is.na(WGT4YR) & WGT4YR > 0) %>%
  # 3. 构建代谢综合征各组分风险变量
  mutate(
    waist_risk   = as.integer((RIAGENDR == 'Male' & BMXWAIST >= 102) | (RIAGENDR == 'Female' & BMXWAIST >= 88)),
    trigly_risk  = as.integer(LBXTR >= 150),
    hdl_risk     = as.integer((RIAGENDR == 1 & LBDHDD < 40) | (RIAGENDR == 2 & LBDHDD < 50)),
    bp_risk      = as.integer(BPXSY1 >= 130 | BPXDI1 >= 85),
    glucose_risk = as.integer(LBXGLU >= 100)
  ) %>%
  # 4. 计算代谢综合征组分数及状态
  mutate(
    metS_components = rowSums(across(c(waist_risk, trigly_risk, hdl_risk, bp_risk, glucose_risk)), na.rm = FALSE),
    metS_status = case_when(
      metS_components >= 3 ~ 1,
      metS_components < 3 ~ 0,
      TRUE ~ NA_real_
    )) %>% filter(!is.na(metS_status)) %>%
  # 5. 睡眠变量处理
  mutate(
    # 将睡眠时间字符串转为小时（如"23:30"）
    bedtime_hours = suppressWarnings(
      ifelse(
        is.na(SLQ300) | SLQ300 == "Don't know",
        NA_real_,
        as.numeric(substr(SLQ300, 1, 2)) + as.numeric(substr(SLQ300, 4, 5))/60
      )
    ),
    waketime_hours = suppressWarnings(
      ifelse(
        is.na(SLQ310) | SLQ310 == "Don't know",
        NA_real_,
        as.numeric(substr(SLQ310, 1, 2)) + as.numeric(substr(SLQ310, 4, 5))/60
      )
    ),
    # 计算睡眠时长，跨午夜自动处理
    sleep_duration = case_when(
      is.na(bedtime_hours) | is.na(waketime_hours) ~ NA_real_,
      bedtime_hours > waketime_hours ~ waketime_hours + 24 - bedtime_hours,
      TRUE ~ waketime_hours - bedtime_hours
    ),
    # 睡眠时长分组
    sleep_cat = case_when(
      is.na(sleep_duration) ~ NA_character_,
      sleep_duration < 6 ~ "Short sleep (<6h)",
      sleep_duration <= 8 ~ "Optimal sleep (6-8h)",
      sleep_duration > 8 ~ "Long sleep (>8h)"
    )
  )%>%filter(!is.na(sleep_cat)) %>%
  
  
  # --- 构建调节变量 PA  ---
  
  # 步骤1：将Yes/No转换为逻辑值，NA保持为NA
  mutate(
    PAQ605_logic = case_when(
      PAQ605 == "Yes" ~ TRUE,
      PAQ605 == "No" ~ FALSE,
      TRUE ~ NA
    ),
    PAQ620_logic = case_when(
      PAQ620 == "Yes" ~ TRUE,
      PAQ620 == "No" ~ FALSE,
      TRUE ~ NA
    ),
    PAQ650_logic = case_when(
      PAQ650 == "Yes" ~ TRUE,
      PAQ650 == "No" ~ FALSE,
      TRUE ~ NA
    ),
    PAQ665_logic = case_when(
      PAQ665 == "Yes" ~ TRUE,
      PAQ665 == "No" ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  
  # 步骤2：处理数值型变量，将NA保持为NA
  mutate(
    PAQ610_num = as.numeric(ifelse(PAQ610 == "NA", NA, PAQ610)),
    PAQ625_num = as.numeric(ifelse(PAQ625 == "NA", NA, PAQ625)),
    PAQ655_num = as.numeric(ifelse(PAQ655 == "NA", NA, PAQ655)),
    PAQ670_num = as.numeric(ifelse(PAQ670 == "NA", NA, PAQ670)),
    PAD680_num = as.numeric(ifelse(PAD680 == "NA", NA, PAD680))
  ) %>%
  
  # 步骤3：计算每项活动的总时长（分钟/周）
  # 注意：这里我们假设每次活动持续60分钟，如果有实际时长数据请替换
  mutate(
    # 剧烈工作活动分钟/周 (假设每次60分钟)
    vigorous_work_min = ifelse(PAQ605_logic == TRUE, PAQ610_num * 60, 0),
    # 中等工作活动分钟/周 (假设每次60分钟)
    moderate_work_min = ifelse(PAQ620_logic == TRUE, PAQ625_num * 60, 0),
    # 剧烈休闲活动分钟/周
    vigorous_rec_min = ifelse(PAQ650_logic == TRUE, PAQ655_num * 60, 0),
    # 中等休闲活动分钟/周
    moderate_rec_min = ifelse(PAQ665_logic == TRUE, PAQ670_num * 60, 0)
  ) %>%
  
  # 步骤4：计算MET值
  mutate(
    # 剧烈活动MET-min/week
    vigorous_MET = (vigorous_work_min + vigorous_rec_min) * 8.0,
    # 中等活动MET-min/week
    moderate_MET = (moderate_work_min + moderate_rec_min) * 4.0,
    # 总体力活动MET-min/week
    total_PA_MET_min = vigorous_MET + moderate_MET,
    # 总体力活动MET-h/week
    total_PA_MET_hour = total_PA_MET_min / 60,
    # 久坐行为（小时/周）
    sedentary_hour_week = PAD680_num / 60,
    
    PA_level = factor(case_when(
      total_PA_MET_min < 600 ~ "Low",
      total_PA_MET_min >= 600 & total_PA_MET_min < 3000 ~ "Medium",
      total_PA_MET_min >= 3000 ~ "High",
      TRUE ~ NA_character_ # 可能存在所有PA变量都缺失的情况
    ), levels = c("Low", "Medium", "High"))
  ) %>%
  
  
  # --- 构建调节变量 DQ (极简化代理评分) ---
  mutate(
    DR1TKCAL = ifelse(is.na(DR1TKCAL) | DR1TKCAL <= 0, NA, DR1TKCAL), # 处理无效能量值
    K_Na_ratio = DR1TPOTA / DR1TSODI,
    Fiber_density = DR1TFIBE / (DR1TKCAL / 1000),
    # 处理比率中的Inf或NaN
    K_Na_ratio = ifelse(!is.finite(K_Na_ratio), NA, K_Na_ratio),
    Fiber_density = ifelse(!is.finite(Fiber_density), NA, Fiber_density),
    # 简化组合 (可能有NA)
    DQ_proxy_score = scale(K_Na_ratio)[,1] + scale(Fiber_density)[,1] # 使用 scale 并取第一列确保是向量
  ) %>%
  
  
  
  # --- 构建其他协变量 X ---
  mutate(
    BMI = BMXBMI,
    Age = RIDAGEYR,
    Income = INDFMPIR,
    Gender = factor(ifelse(RIAGENDR == "Female", "Female", "Male")),
    RaceEthnicity = factor(RIDRETH3),
    Education = factor(DMDEDUC2),
    SmokingStatus = factor(case_when(
      SMQ020 == "No" ~ "Never Smoker",
      SMQ020 == "Yes" & SMQ040 == "Not at all" ~ "Former Smoker",
      SMQ020 == "Yes" & SMQ040 %in% c("Every day", "Some days") ~ "Current Smoker",
      TRUE ~ NA_character_
    ), levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
    AlcoholAvgDrinks = ALQ130,
    AlcoholAvgDrinks = ifelse(AlcoholAvgDrinks %in% c(777, 999), NA, AlcoholAvgDrinks),
    AlcoholCat = factor(case_when(
      AlcoholAvgDrinks == 0 ~ "None", # 更清晰的None定义
      AlcoholAvgDrinks > 0 & AlcoholAvgDrinks <= 1 ~ "Light",
      AlcoholAvgDrinks > 1 & AlcoholAvgDrinks <= 3 ~ "Moderate", # 调整阈值示例
      AlcoholAvgDrinks > 3 ~ "Heavy",
      TRUE ~ NA_character_
    ), levels = c("Light", "Moderate", "Heavy")))


phq9_map <- c(
  "Not at all" = 0,
  "Several days" = 1,
  "More than half the days" = 2,
  "Nearly every day" = 3,
  "Don't know" = NA,
  "Refused" = NA
)

nhanes_processed <- nhanes_processed  %>%
  mutate(
    across(
      # 只处理PHQ-9原始题目（如DPQ030, DPQ040, ...），不影响已有_num变量
      matches("^DPQ\\d{3}$"),
      ~ unname(phq9_map[.]),
      .names = "{.col}_num"
    )
  ) %>%
  mutate(
    # 只统计新建的DPQ_num变量
    phq9_miss_count = rowSums(is.na(select(., matches("^DPQ\\d{3}_num$")))),
    PHQ9_score = ifelse(
      phq9_miss_count <= 2,
      rowSums(select(., matches("^DPQ\\d{3}_num$")), na.rm = TRUE),
      NA
    )
  ) %>%
  
  mutate(
    kappa = ifelse(RIAGENDR == 2, 0.7, 0.9),
    alpha = ifelse(RIAGENDR == 2, -0.241, -0.302),
    eGFR = ifelse(
      !is.na(LBXSCR) & LBXSCR > 0 & !is.na(RIDAGEYR),
      142 *
        (pmin(LBXSCR / kappa, 1) ^ alpha) *
        (pmax(LBXSCR / kappa, 1) ^ -1.200) *
        (0.9938 ^ RIDAGEYR) *
        ifelse(RIAGENDR == 2, 1.012, 1),
      NA_real_
    )
  ) %>%
  
  # --- 选择最终分析变量 ---
  select(
    SEQN, WGT4YR, SDMVPSU, SDMVSTRA, # ID & Design
    metS_status, sleep_cat, # Y & W
    PA_level, DQ_proxy_score, Age, Gender, RaceEthnicity, Education, Income, BMI, # X vars
    SmokingStatus, AlcoholCat, PHQ9_score, eGFR # X vars continued
  )


# --- 检查处理后的数据概况 ---
print("初步处理完成，数据概况：")
summary(nhanes_processed)
print(paste("处理后行数:", nrow(nhanes_processed)))

# --- 2. 处理缺失值 (多重插补 - MICE) ---
vars_to_impute <- nhanes_processed %>%
  select(-SEQN, -WGT4YR, -SDMVPSU, -SDMVSTRA) %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "na_count") %>%
  filter(na_count > 0 & na_count < nrow(nhanes_processed)) # 只显示有缺失但非完全缺失的

print("需要插补的变量（缺失比例 > 0）：")
print(vars_to_impute)

# 选择用于插补的变量 (通常是所有分析变量，除了ID和设计变量)
imputation_candidates <- nhanes_processed %>%
  select(-SEQN, -WGT4YR, -SDMVPSU, -SDMVSTRA)

# 检查常数变量或共线性问题（可能导致MICE出错）
# findCorrelation, findLinearCombos from caret package can help, but skip for now

# 执行插补 (使用较少迭代次数以加速示例)
# 注意：因子变量会自动处理。如果MICE出错，可能需要检查数据类型或共线性。
set.seed(123) # 保证可重复性
mice_obj <- mice(imputation_candidates, m = 5, method = 'pmm', maxit = 100, printFlag = TRUE)

# 提取第一个插补完整的数据集
# !! 再次强调：正式研究应汇总所有插补数据集结果 !!
nhanes_imputed <- complete(mice_obj, 1)

# 将设计变量加回
nhanes_final <- bind_cols(
  nhanes_processed %>% select(SEQN, WGT4YR, SDMVPSU, SDMVSTRA),
  nhanes_imputed
) %>% filter(!is.na(sleep_cat)) # 确保处理变量无缺失

print("--- 插补完成 (使用第一个插补集)，准备最终数据用于分析 ---")
print(paste("最终数据集维度:", dim(nhanes_final)[1], "rows,", dim(nhanes_final)[2], "columns."))
# 检查插补后的数据
summary(nhanes_final)

# 保存数据
save(nhanes_final, file = "nhanes_final.rda")
# 
# # --- 3. 准备用于 grf 包的数据 ---
# Y <- nhanes_final$metS_status
# W <- nhanes_final$sleep_cat # 因子
# X_df <- nhanes_final %>% select(PA_level, DQ_proxy_score, Age, Gender, RaceEthnicity, Education, Income, BMI, SmokingStatus, AlcoholCat, PHQ9_score, eGFR)
# 
# # 检查因子水平，确保没有问题
# # sapply(X_df, function(x) if(is.factor(x)) levels(x) else NULL)
# 
# # 创建数值矩阵 X (包含哑变量)
# # 注意：如果因子变量只有1个水平（数据问题），model.matrix会出错
# # 可以在这里添加检查或移除只有一个水平的因子
# X_df <- X_df %>% mutate(across(where(is.factor), droplevels)) # 移除未使用的因子水平
# X <- model.matrix(~ . - 1, data = X_df)
# 
# weights <- nhanes_final$WGT4YR # 使用(近似的)4年权重
# clusters <- nhanes_final$SEQN # 使用 SEQN 作为聚类ID (个体层面稳健标准误)
# 
# print("--- 数据已准备好 (Y, W, X, weights, clusters)，可以进行EDA和因果森林建模 ---")
# # 检查最终对象的维度和类型
# dim(X); class(Y); class(W); length(weights); length(clusters)
# 
# 
# 
# # --- 后续步骤 (EDA, Causal Forest, Interpretation) ---
# # 使用 Y, W, X, weights, clusters 对象进行分析
# # ... 在此添加 EDA 和 Causal Forest 建模代码 ...