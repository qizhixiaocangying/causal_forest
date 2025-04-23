# 清理环境
rm(list = ls())

# 包管理（自动安装并加载所需包）
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, survey, mice, grf, ggplot2, car, dplyr)

# 读取数据
load("nhanes_full.rda") # 请确保文件路径正确

# --- 数据预处理与变量构建 ---
nhanes_processed <- nhanes_full %>%
  # 1. 初始过滤
  filter(
    !is.na(WTINT2YR),
    !is.na(WTMEC2YR),
    RIDAGEYR >= 20,                    # 仅保留20岁及以上
    RIDEXPRG != 1 | is.na(RIDEXPRG)    # 排除孕妇
  ) %>%
  # 2. 构建权重
  mutate(WGT4YR = WTMEC2YR / 2) %>%
  filter(!is.na(WGT4YR) & WGT4YR > 0) %>%
  
  # 3. 构建代谢综合征各组分风险变量
  mutate(
    waist_risk   = as.integer((RIAGENDR == 'Male' & BMXWAIST >= 102) | (RIAGENDR == 'Female' & BMXWAIST >= 88)),
    trigly_risk  = as.integer(LBXTR >= 150),
    hdl_risk     = as.integer((RIAGENDR == 'Male' & LBDHDD < 40) | (RIAGENDR == 'Female' & LBDHDD < 50)),
    bp_risk      = as.integer(BPXSY1 >= 130 | BPXDI1 >= 85),
    glucose_risk = as.integer(LBXGLU >= 100)
  ) %>%
  
  # 4. 计算代谢综合征组分数及状态
  mutate(
    metS_components = rowSums(across(c(waist_risk, trigly_risk, hdl_risk, bp_risk, glucose_risk)), na.rm = TRUE),
    metS_status = case_when(
      metS_components >= 3 ~ 1,
      metS_components < 3 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>% 
  filter(!is.na(metS_status)) %>%
  
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
    sleep_cat = factor(case_when(
      is.na(sleep_duration) ~ NA_character_,
      sleep_duration < 6 ~ "Short sleep (<6h)",
      sleep_duration <= 9 ~ "Optimal sleep (6-9h)",
      sleep_duration > 9 ~ "Long sleep (>9h)"
    ))
  ) %>% 
  filter(!is.na(sleep_cat)) %>%
  
  # --- 构建调节变量 PA ---
  mutate(
    # 1. 处理工作相关活动（使用实际记录分钟数）
    vigorous_work_min = ifelse(PAQ605 == "Yes" & !is.na(PAD615), PAD615, 0),
    moderate_work_min = ifelse(PAQ620 == "Yes" & !is.na(PAD630), PAD630, 0),
    
    # 2. 处理休闲活动（需要结合天数和分钟数）
    vigorous_rec_min = case_when(
      PAQ650 == "Yes" & !is.na(PAQ655) & !is.na(PAD660) ~ PAQ655 * PAD660,
      TRUE ~ 0
    ),
    moderate_rec_min = case_when(
      PAQ665 == "Yes" & !is.na(PAQ670) & !is.na(PAD675) ~ PAQ670 * PAD675,
      TRUE ~ 0
    ),
    
    # 3. 计算MET-min/week（使用标准MET值）
    vigorous_MET = (vigorous_work_min + vigorous_rec_min) * 8.0,
    moderate_MET = (moderate_work_min + moderate_rec_min) * 4.0,
    total_PA_MET_min = vigorous_MET + moderate_MET,
    
    # 4. 分类体力活动水平（根据WHO标准）
    PA_level = factor(
      case_when(
        total_PA_MET_min >= 3000 ~ "High",
        total_PA_MET_min >= 600 & total_PA_MET_min < 3000 ~ "Medium",
        total_PA_MET_min < 600 ~ "Low",
        TRUE ~ NA_character_
      ),
      levels = c("Low", "Medium", "High")
    )
  ) %>%
  
  # --- 构建调节变量 DQ (极简化代理评分) ---
  mutate(
    DR1TKCAL = ifelse(is.na(DR1TKCAL) | DR1TKCAL <= 0, NA, DR1TKCAL),
    K_Na_ratio = DR1TPOTA / DR1TSODI,
    Fiber_density = DR1TFIBE / (DR1TKCAL / 1000),
    K_Na_ratio = ifelse(!is.finite(K_Na_ratio), NA, K_Na_ratio),
    Fiber_density = ifelse(!is.finite(Fiber_density), NA, Fiber_density),
    DQ_proxy_score = scale(K_Na_ratio)[,1] + scale(Fiber_density)[,1]
  ) %>%
  
  # --- 构建其他协变量 ---
  mutate(
    BMI = BMXBMI,
    Age = RIDAGEYR,
    Income = INDFMPIR,
    Gender = factor(ifelse(RIAGENDR == "Female", "Female", "Male")),
    RaceEthnicity = factor(RIDRETH3),
    Education = factor(case_when(
      DMDEDUC2 %in% c("Don't Know", "Refused") ~ NA_character_,
      TRUE ~ as.character(DMDEDUC2))
    ), 
    SmokingStatus = factor(case_when(
      SMQ020 == "No" ~ "Never Smoker",
      SMQ020 == "Yes" & SMQ040 == "Not at all" ~ "Former Smoker",
      SMQ020 == "Yes" & SMQ040 %in% c("Every day", "Some days") ~ "Current Smoker",
      TRUE ~ NA_character_
    ), levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
    AlcoholAvgDrinks = ALQ130,
    AlcoholAvgDrinks = ifelse(AlcoholAvgDrinks %in% c(777, 999), NA, AlcoholAvgDrinks),
    AlcoholCat = factor(case_when(
      AlcoholAvgDrinks > 0 & AlcoholAvgDrinks <= 1 ~ "Light",
      AlcoholAvgDrinks > 1 & AlcoholAvgDrinks <= 3 ~ "Moderate",
      AlcoholAvgDrinks > 3 ~ "Heavy",
      is.na(AlcoholAvgDrinks) ~ "Unknown/None",
      TRUE ~ NA_character_
    ), levels = c("Light", "Moderate", "Heavy","Unknown/None"))
  ) %>%
  
  # PHQ-9抑郁评分
  mutate(
    across(
      matches("^DPQ\\d{3}$"),
      ~ unname(c(
        "Not at all" = 0,
        "Several days" = 1,
        "More than half the days" = 2,
        "Nearly every day" = 3,
        "Don't know" = NA,
        "Refused" = NA
      )[.]),
      .names = "{.col}_num"
    )
  ) %>%
  mutate(
    phq9_miss_count = rowSums(is.na(select(., matches("^DPQ\\d{3}_num$")))),
    PHQ9_score = ifelse(
      phq9_miss_count <= 2,
      rowSums(select(., matches("^DPQ\\d{3}_num$")), na.rm = TRUE),
      NA
    )
  ) %>%
  
  # eGFR计算
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
    metS_status, sleep_cat, sleep_duration, waist_risk, trigly_risk, hdl_risk, bp_risk, glucose_risk, # Y & W
    bedtime_hours, PA_level, DQ_proxy_score, Age, Gender, RaceEthnicity, Education, Income, BMI, # X vars
    SmokingStatus, AlcoholCat, PHQ9_score # X vars continued
  )

# --- 缺失值处理 (多重插补) ---
set.seed(123)
mice_obj <- mice(
  nhanes_processed %>% select(-SEQN, -WGT4YR, -SDMVPSU, -SDMVSTRA, -metS_status, -sleep_duration, -bedtime_hours, -sleep_cat),
  m = 2, 
  method = 'pmm', 
  maxit = 10
)

# 提取第一个插补数据集
nhanes_imputed <- complete(mice_obj, 1)

# 合并回设计变量和处理变量
nhanes_final <- bind_cols(
  nhanes_processed %>% select(SEQN, WGT4YR, SDMVPSU, SDMVSTRA, metS_status, sleep_cat, sleep_duration, bedtime_hours),
  nhanes_imputed
) %>% 
  filter(!is.na(sleep_cat))

# --- 保存最终数据 ---
save(nhanes_final, file = "nhanes_final.rda")

# --- 检查数据 ---
print("数据处理完成，最终数据概况：")
summary(nhanes_final)
print(paste("最终样本量:", nrow(nhanes_final)))

