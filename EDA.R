# 清理环境
rm(list = ls())

# 加载必要包
library(survey)
library(tableone)


load("nhanes_final.rda")

# 探索性数据分析
# 创建调查设计对象

# 1. 构建复杂抽样设计对象
nhanes_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WGT4YR,
  data = nhanes_final,
  nest = TRUE
)



# 2. 指定需要汇总的变量（除设计变量外）
vars <- c("metS_status", "sleep_duration", #"sleep_cat",
          "bedtime_hours","PA_level", "DQ_proxy_score", "Age", "Gender", "RaceEthnicity",
          "Education", "Income", "BMI", "SmokingStatus", "AlcoholCat",
          "PHQ9_score", "eGFR")

# 3. 指定分类变量（factor变量）
factorVars <- c("metS_status", "sleep_cat", "PA_level", "Gender",
                "RaceEthnicity", "Education", "SmokingStatus", "AlcoholCat")


# 4. 创建加权的TableOne对象，按metS_status分层
table1_weighted <- svyCreateTableOne(
  vars = vars,
  strata = "sleep_cat",
  data = nhanes_design,
  factorVars = factorVars,
  test = TRUE,    # 组间比较检验
  smd = TRUE      # 标准化均数差
)

# 5. 打印结果，显示分类变量所有水平，连续变量非正态时显示中位数和IQR
print(table1_weighted,
      showAllLevels = TRUE,
      #nonnormal = c("DQ_proxy_score", "PHQ9_score", "eGFR", "BMI", "Income", "Age")
)

# -------------------------------------------------------------------------------------
# 描述性统计 (加权)
# - 总体描述
svytable(~sleep_cat, design = nhanes_design)
svymean(~Age + BMI + DQ_proxy_score+sleep_duration, design = nhanes_design)
# - 按睡眠组描述
svyby(~Age + BMI + DQ_proxy_score, ~sleep_cat, design = nhanes_design, svymean)
svytable(~PA_level + sleep_cat, design = nhanes_design) # PA分布按睡眠组

# 可视化 (示例)
# - MetS 状态按睡眠时间 (加权比例)
prop_metS <- svyby(~metS_status, ~sleep_cat, design = nhanes_design, svymean)
ggplot(prop_metS, aes(x = sleep_cat, y = metS_status, fill = sleep_cat)) +
  geom_bar(stat = "identity") + labs(title = "Weighted Prevalence of MetS by Sleep Category")

# 检查重叠性 (初步)
# - 查看关键协变量（如PA, DQ, BMI, Age）在不同睡眠组的分布
ggplot(nhanes_final, aes(x = DQ_proxy_score, fill = sleep_cat)) + geom_density(alpha = 0.5)
ggplot(nhanes_final, aes(x = sleep_cat, y = BMI, fill = sleep_cat)) + geom_boxplot()




