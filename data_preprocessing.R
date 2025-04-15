# 数据获取与预处理


# install.packages(c("nhanesA", "tidyverse", "survey", "mice", "grf", "ggplot2"))
library(nhanesA)
library(tidyverse)
library(survey)
library(mice)
library(grf)
library(ggplot2)

# --- 验证 2021-2022 ("L") 周期数据的可用性 ---
# 你可以列出某个组件在某年份可用的表格，例如:
# nhanesTables(data_group='DEMO', year=2021)
# 或者搜索特定变量:
# nhanesSearch("SLD010H") # 检查哪个周期包含这个睡眠变量

print("--- 开始下载 NHANES 2021-2022 (周期 L) 数据 ---")
print("注意：请确保周期 'L' 所需的所有文件均已在NHANES网站发布。")

# 1. 下载 NHANES 2021-2022 (周期 "L") 数据
#    将文件名中的 "_J" 替换为 "_L"。如有需要，请再次核对实际文件名。

# 人口统计学 (DEMO_L)
demo_vars_L <- nhanes("DEMO_L") %>%
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3, DMDEDUC2, INDFMPIR, RIDEXPRG, # 核心人口统计学变量 + 怀孕状态
         SDMVPSU, SDMVSTRA, WTMEC2YR) # 调查设计变量 (使用2年MEC权重)

# 睡眠 (SLQ_L)
sleep_vars_L <- nhanes("SLQ_L") %>% select(SEQN, SLD010H) # 核实 SLD010H 在 L 周期是否仍是正确的变量名

# 代谢综合征组分
# - 身体测量 (BMX_L)
bmx_vars_L <- nhanes("BMX_L") %>% select(SEQN, BMXWAIST, BMXHT, BMXWT, BMXBMI) # 腰围, 身高, 体重, BMI
# - 甘油三酯 (TRIGLY_L) & HDL (HDL_L)
trigly_vars_L <- nhanes("TRIGLY_L") %>% select(SEQN, LBXTR)
hdl_vars_L <- nhanes("HDL_L") %>% select(SEQN, LBDHDL)
# - 血压 (BPX_L)
bpx_vars_L <- nhanes("BPX_L") %>% select(SEQN, BPXSY1, BPXDI1) # 核实 BPXSY1/DI1 是否仍在使用
# - 血糖 (GLU_L)
glu_vars_L <- nhanes("GLU_L") %>% select(SEQN, LBXGLU, PHAFSTHR) # 如有需要，核实空腹小时数变量名
# -肌酐 (例如, BIOPRO_L 或类似的标准化学生化组套文件)
#   (需要确认周期 L 中包含肌酐的文件名)
#   示例占位符: 假设在 BIOPRO_L 中
biopro_vars_L <- nhanes("BIOPRO_L") %>% select(SEQN, LBXSCR) # 需要血清肌酐 LBXSCR

# 处方药 (RXQ_RX_L)
rx_vars_L <- nhanes("RXQ_RX_L") # 用于在定义MetS组分时检查相关药物使用情况

# 体力活动 (PAQ_L)
paq_vars_L <- nhanes("PAQ_L") # 包含多个体力活动相关问题

# 膳食数据 (DR1TOT_L, DR2TOT_L, DR1IFF_L, DR2IFF_L) - 仍然是最复杂的部分
diet1_tot_L <- nhanes("DR1TOT_L") # 第一天总营养素摄入
diet2_tot_L <- nhanes("DR2TOT_L") # 第二天总营养素摄入
diet1_iff_L <- nhanes("DR1IFF_L") # 第一天个体食物摄入
diet2_iff_L <- nhanes("DR2IFF_L") # 第二天个体食物摄入
# 注意: 处理膳食数据以计算饮食质量(DQ)得分仍然复杂。

# 其他协变量 (吸烟 SMQ_L, 饮酒 ALQ_L, 抑郁 DPQ_L)
smq_vars_L <- nhanes("SMQ_L") %>% select(SEQN, SMQ020, SMQ040) # 如有需要，核实变量名
alq_vars_L <- nhanes("ALQ_L") %>% select(SEQN, ALQ130) # 如有需要，核实变量名
dpq_vars_L <- nhanes("DPQ_L") %>% select(SEQN, DPQ010:DPQ090) # 如有需要，核实变量名


# 2. 合并数据
# 列出为周期 L 下载的所有数据框
data_files_L <- list(demo_vars_L, sleep_vars_L, bmx_vars_L, trigly_vars_L, hdl_vars_L,
                     bpx_vars_L, glu_vars_L, biopro_vars_L, rx_vars_L, paq_vars_L,
                     smq_vars_L, alq_vars_L, dpq_vars_L,
                     diet1_tot_L, diet2_tot_L, diet1_iff_L, diet2_iff_L /* 根据需要添加其他文件 */)

# 使用 reduce 和 full_join 按 SEQN 合并
nhanes_raw_L <- data_files_L %>% reduce(full_join, by = "SEQN")

print(paste("已合并周期 L 的原始数据。维度:", dim(nhanes_raw_L)[1], "行,", dim(nhanes_raw_L)[2], "列。"))
print("--- 后续步骤: 数据清洗, 变量构建, EDA, 模型建立 ---")

# 3. 数据清洗与变量构建 (概念性结构 - 需要详细实现)
#    这里的逻辑与之前的示例相同，但应用于 nhanes_raw_L 数据框。
#    确保使用已为周期 L 文档核实过的正确变量名。
#    示例:
nhanes_processed_L <- nhanes_raw_L %>%
  filter(RIDAGEYR >= 20 & RIDEXPRG != 1) %>% # 筛选成年人, 排除孕妇
  mutate(
    # W: 睡眠类别 (使用 _L 周期的 SLD010H)
    sleep_hours = SLD010H,
    sleep_cat = factor(case_when(
      sleep_hours < 6 ~ "Short (<6h)",
      sleep_hours >= 7 & sleep_hours <= 8 ~ "Normal (7-8h)", # 正常睡眠 (参照组)
      sleep_hours > 9 ~ "Long (>9h)", # 睡眠过长
      TRUE ~ NA_character_ # 处理 6h, 9h 或缺失值
    ), levels = c("Normal (7-8h)", "Short (<6h)", "Long (>9h)")), # 设置参照组顺序
    
    # Y: MetS 状态 (需要使用 _L 周期的变量进行复杂计算)
    # 占位符 - 替换为周期 L 的实际 MetS 计算逻辑
    metS_status = rbinom(n(), 1, 0.35),
    
    # 调节变量: PA & DQ (需要使用 _L 周期的变量进行复杂计算)
    # 占位符 - 替换为周期 L 的实际 PA/DQ 计算逻辑
    PA_level = sample(c("Low", "Medium", "High"), n(), replace = TRUE, prob = c(0.4, 0.4, 0.2)), # 低/中/高 体力活动水平 (随机示例)
    DQ_score = rnorm(n(), 60, 15), # 饮食质量得分 (随机示例)
    
    # X: 协变量 (使用 _L 周期的变量)
    BMI = BMXBMI,
    Age = RIDAGEYR,
    IsFemale = ifelse(RIAGENDR == 2, 1, 0), # 是否女性
    # ... 计算其他协变量，如PHQ9总分、吸烟状态等 ...
    
  ) %>%
  # 选择最终分析变量，包括周期 L 的调查设计变量
  select(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA, # 周期 L 的设计变量 (使用2年权重)
         metS_status, sleep_cat, # 结局 Y 和 处理 W
         PA_level, DQ_score, BMI, Age, IsFemale, # 协变量 X (包括关键调节变量)
         # ... 添加所有其他选择的 X 变量 ...
  ) %>%
  # 处理缺失值 (例如：仅保留完整案例，或进行多重插补)
  filter(complete.cases(.)) # 为简化示例，使用完整案例分析法

# --- 准备用于 grf 包的数据 (Y, W, X, weights, clusters) ---
# 这里的结构与之前的示例相同，只是使用处理过的 _L 数据
Y_L <- nhanes_processed_L$metS_status
W_L <- nhanes_processed_L$sleep_cat # 因子类型变量
X_df_L <- nhanes_processed_L %>% select(PA_level, DQ_score, BMI, Age, IsFemale, /* ... 其他 X 变量 ... */)
X_L <- model.matrix(~ . - 1, data = X_df_L) # 将因子变量转换为哑变量矩阵
weights_L <- nhanes_processed_L$WTMEC2YR # 使用2年权重
clusters_L <- nhanes_processed_L$SDMVPSU # 聚类变量

print("--- 已准备好使用 2021-2022 数据进行因果森林分析 ---")