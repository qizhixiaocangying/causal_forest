# 数据预处理
library(tidyverse)
library(survey)
library(mice)
library(grf)
library(ggplot2)
library(car)

# 读取数据
load("nhanes_merged.rda")

# 1. 数据清洗与变量构建
nhanes_processed <- nhanes_merged %>%
  # --- 筛选人群 ---
  filter(RIDAGEYR >= 20) %>%  # 成年人
  filter(RIDEXPRG != 1 | is.na(RIDEXPRG)) %>% # 排除孕妇
  
  # --- 处理权重 ---
  # 创建近似的4年权重 WGT4YR = WTMEC2YR / 2 (因为合并了2个周期)
  # 注意：如果DEMO_J文件本身包含 WTMEC4YR，应优先使用那个变量
  mutate(WGT4YR = WTMEC2YR / 2) %>%
  filter(!is.na(WGT4YR) & WGT4YR > 0) %>%
  
  # --- 构建结局变量 Y (MetS) ---
  # 注意：药物逻辑是高度简化的示例！
  mutate(
    # 药物使用标志 (简化逻辑)
    on_bp_med = map_lgl(med_list, ~ if(is.null(.x)) FALSE else any(str_detect(.x, "pril|sartan|lol|dipine|thiazide|lisinopril|amlodipine|metoprolol"))),
    on_diab_med = map_lgl(med_list, ~ if(is.null(.x)) FALSE else any(str_detect(.x, "formin|glipizide|glyburide|pioglitazone|insulin"))),
    on_lipid_med = map_lgl(med_list, ~ if(is.null(.x)) FALSE else any(str_detect(.x, "statin|fibrate|ezetimibe|niacin|atorvastatin|simvastatin"))),
    # 处理可能因left_join产生的NA (如果没有药物记录)
    across(c(on_bp_med, on_diab_med, on_lipid_med), ~replace_na(., FALSE)),
    
    # MetS 组分判断
    waist_risk = ifelse( (RIAGENDR == 1 & BMXWAIST >= 102) | (RIAGENDR == 2 & BMXWAIST >= 88), 1, 0),
    trigly_risk = ifelse(LBXTR >= 150 | on_lipid_med, 1, 0),
    hdl_risk = ifelse( (RIAGENDR == 1 & LBDHDL < 40) | (RIAGENDR == 2 & LBDHDL < 50) | on_lipid_med, 1, 0),
    bp_risk = ifelse(BPXSY1 >= 130 | BPXDI1 >= 85 | on_bp_med, 1, 0),
    glucose_risk = ifelse(LBXGLU >= 100 | on_diab_med, 1, 0),
    glucose_risk = ifelse(!is.na(PHAFSTHR) & PHAFSTHR < 8, NA, glucose_risk),
    
    # 计算MetS状态 (处理组分判断中的NA)
    metS_components = rowSums(select(., waist_risk, trigly_risk, hdl_risk, bp_risk, glucose_risk), na.rm = FALSE), # NA if any component is NA
    metS_status = case_when(
      metS_components >= 3 ~ 1,
      metS_components < 3 ~ 0,
      TRUE ~ NA_real_ # If components have NA resulting in NA sum
    )
  ) %>%
  
  # --- 构建处理变量 W (睡眠类别) ---
  mutate(
    sleep_hours = SLD010H,
    sleep_hours = ifelse(sleep_hours %in% c(77, 99), NA, sleep_hours),
    sleep_cat = factor(case_when(
      sleep_hours < 6 ~ "Short (<6h)",
      sleep_hours >= 7 & sleep_hours <= 8 ~ "Normal (7-8h)",
      sleep_hours > 9 ~ "Long (>9h)",
      TRUE ~ NA_character_
    ), levels = c("Normal (7-8h)", "Short (<6h)", "Long (>9h)"))
  ) %>%
  filter(!is.na(sleep_cat)) %>% # 移除处理变量缺失的样本
  
  # --- 构建调节变量 PA (简化版 MET-min/周) ---
  mutate(
    across(all_of(c("PAQ605", "PAQ610", "PAQ620", "PAQ625", "PAQ650", "PAQ655", "PAQ665", "PAQ670")), ~replace_na(., 0)),
    # 确保变量存在再计算
    met_min_vig = (if_else("PAQ605" %in% names(.), PAQ605, 0) * if_else("PAQ610" %in% names(.), PAQ610, 0) +
                     if_else("PAQ650" %in% names(.), PAQ650, 0) * if_else("PAQ655" %in% names(.), PAQ655, 0)) * 8,
    met_min_mod = (if_else("PAQ620" %in% names(.), PAQ620, 0) * if_else("PAQ625" %in% names(.), PAQ625, 0) +
                     if_else("PAQ665" %in% names(.), PAQ665, 0) * if_else("PAQ670" %in% names(.), PAQ670, 0)) * 4,
    total_met_min = met_min_vig + met_min_mod,
    PA_level = factor(case_when(
      total_met_min < 500 ~ "Low",
      total_met_min >= 500 & total_met_min < 3000 ~ "Medium",
      total_met_min >= 3000 ~ "High",
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
    IsFemale = factor(ifelse(RIAGENDR == 2, "Female", "Male")),
    RaceEthnicity = factor(case_when(
      RIDRETH3 == 1 ~ "Mexican American", RIDRETH3 == 2 ~ "Other Hispanic",
      RIDRETH3 == 3 ~ "Non-Hispanic White", RIDRETH3 == 4 ~ "Non-Hispanic Black",
      RIDRETH3 == 6 ~ "Non-Hispanic Asian", RIDRETH3 == 7 ~ "Other/Mixed Race", TRUE ~ NA_character_
    )),
    Education = factor(case_when(
      DMDEDUC2 == 1 ~ "<9th grade", DMDEDUC2 == 2 ~ "9-11th grade", DMDEDUC2 == 3 ~ "High school grad/GED",
      DMDEDUC2 == 4 ~ "Some college/AA degree", DMDEDUC2 == 5 ~ "College grad or above", TRUE ~ NA_character_ # 简化：其他编码视为NA
    )),
    PIR = INDFMPIR,
    SmokingStatus = factor(case_when(
      SMQ020 == 2 ~ "Never Smoker", SMQ020 == 1 & SMQ040 == 3 ~ "Former Smoker",
      SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current Smoker", TRUE ~ NA_character_
    ), levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
    AlcoholAvgDrinks = ALQ130,
    AlcoholAvgDrinks = ifelse(AlcoholAvgDrinks %in% c(777, 999), NA, AlcoholAvgDrinks),
    AlcoholCat = factor(case_when(
      AlcoholAvgDrinks == 0 ~ "None", # 更清晰的None定义
      AlcoholAvgDrinks > 0 & AlcoholAvgDrinks <= 1 ~ "Light",
      AlcoholAvgDrinks > 1 & AlcoholAvgDrinks <= 3 ~ "Moderate", # 调整阈值示例
      AlcoholAvgDrinks > 3 ~ "Heavy",
      is.na(AlcoholAvgDrinks) ~ "Unknown/None", # 处理NA
      TRUE ~ NA_character_
    ), levels = c("None", "Light", "Moderate", "Heavy", "Unknown/None")),
    # PHQ9 (处理7, 9为NA)
    across(starts_with("DPQ"), ~ ifelse(. %in% c(7, 9), NA, .)),
    phq9_miss_count = rowSums(is.na(select(., starts_with("DPQ")))), # 计算缺失项数
    PHQ9_score = ifelse(phq9_miss_count <= 2, rowSums(select(., starts_with("DPQ")), na.rm = TRUE), NA), # 最多允许2项缺失
    # eGFR (使用 CKD-EPI 2021 - 移除非裔美国人调整)
    kappa = ifelse(RIAGENDR == 2, 0.7, 0.9),
    alpha = ifelse(RIAGENDR == 2, -0.241, -0.302),
    eGFR = ifelse(!is.na(LBXSCR) & LBXSCR > 0 & !is.na(RIDAGEYR), # 确保输入有效
                  142 * pmin(LBXSCR / kappa, 1)**alpha * pmax(LBXSCR / kappa, 1)**(-1.200) * 0.9938**RIDAGEYR * ifelse(RIAGENDR == 2, 1.012, 1),
                  NA_real_)
  ) %>%
  
  # --- 选择最终分析变量 ---
  select(
    SEQN, WGT4YR, SDMVPSU, SDMVSTRA, # ID & Design
    metS_status, sleep_cat, # Y & W
    PA_level, DQ_proxy_score, Age, IsFemale, RaceEthnicity, Education, PIR, BMI, # X vars
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
mice_obj <- mice(imputation_candidates, m = 5, method = 'pmm', maxit = 5, printFlag = TRUE)

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


# --- 3. 准备用于 grf 包的数据 ---
Y <- nhanes_final$metS_status
W <- nhanes_final$sleep_cat # 因子
X_df <- nhanes_final %>% select(PA_level, DQ_proxy_score, Age, IsFemale, RaceEthnicity, Education, PIR, BMI, SmokingStatus, AlcoholCat, PHQ9_score, eGFR)

# 检查因子水平，确保没有问题
# sapply(X_df, function(x) if(is.factor(x)) levels(x) else NULL)

# 创建数值矩阵 X (包含哑变量)
# 注意：如果因子变量只有1个水平（数据问题），model.matrix会出错
# 可以在这里添加检查或移除只有一个水平的因子
X_df <- X_df %>% mutate(across(where(is.factor), droplevels)) # 移除未使用的因子水平
X <- model.matrix(~ . - 1, data = X_df)

weights <- nhanes_final$WGT4YR # 使用(近似的)4年权重
clusters <- nhanes_final$SEQN # 使用 SEQN 作为聚类ID (个体层面稳健标准误)

print("--- 数据已准备好 (Y, W, X, weights, clusters)，可以进行EDA和因果森林建模 ---")
# 检查最终对象的维度和类型
dim(X); class(Y); class(W); length(weights); length(clusters)

# --- 后续步骤 (EDA, Causal Forest, Interpretation) ---
# 使用 Y, W, X, weights, clusters 对象进行分析
# ... 在此添加 EDA 和 Causal Forest 建模代码 ...