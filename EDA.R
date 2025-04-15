# 探索性数据分析
# 创建调查设计对象
nhanes_design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR, data = nhanes_final, nest = TRUE)

# 描述性统计 (加权)
# - 总体描述
svytable(~sleep_cat, design = nhanes_design)
svymean(~Age + BMI + DQ_score, design = nhanes_design)
# - 按睡眠组描述
svyby(~Age + BMI + DQ_score, ~sleep_cat, design = nhanes_design, svymean)
svytable(~PA_level + sleep_cat, design = nhanes_design) # PA分布按睡眠组

# 可视化 (示例)
# - MetS 状态按睡眠组 (加权比例)
prop_metS <- svyby(~metS_status, ~sleep_cat, design = nhanes_design, svymean)
ggplot(prop_metS, aes(x = sleep_cat, y = metS_status, fill = sleep_cat)) +
  geom_bar(stat = "identity") + labs(title = "Weighted Prevalence of MetS by Sleep Category")

# 检查重叠性 (初步)
# - 查看关键协变量（如PA, DQ, BMI, Age）在不同睡眠组的分布
ggplot(nhanes_final, aes(x = DQ_score, fill = sleep_cat)) + geom_density(alpha = 0.5)
ggplot(nhanes_final, aes(x = sleep_cat, y = BMI, fill = sleep_cat)) + geom_boxplot()