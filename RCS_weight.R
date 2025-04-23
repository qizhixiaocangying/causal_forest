# 加载必要包
if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, survey, rms, ggplot2, broom, patchwork, ggrepel)

# 读取数据
load("nhanes_final.rda")

# 1. 构建复杂抽样设计对象
nhanes_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WGT4YR,
  data = nhanes_final,
  nest = TRUE
)

# 2. 加权描述性统计
# 计算加权患病率
svymean(~metS_status, design = nhanes_design, na.rm = TRUE)

# 按睡眠时长分组的加权患病率
svyby(~metS_status, ~cut(sleep_duration, breaks = c(4,6,8,10,12)), 
      design = nhanes_design, 
      svymean)

# 3. 准备建模环境
ddist <- datadist(nhanes_final)
options(datadist = "ddist")

# 4. 拟合加权限制性立方样条模型
weighted_model <- svyglm(
  metS_status ~ rcs(sleep_duration, 4),
  design = nhanes_design,
  family = quasibinomial()
)

# 5. 模型评估
# 模型摘要
summary(weighted_model)

# 非线性检验
anova_results <- regTermTest(weighted_model, ~rcs(sleep_duration,4), method = "Wald")
print(anova_results)

# 6. 可视化加权关系
# 生成预测数据
pred_data <- expand.grid(
  sleep_duration = seq(4, 12, length.out = 200)
)

# 获取加权预测值
pred <- predict(weighted_model, 
                newdata = pred_data, 
                type = "response",
                se.fit = TRUE)

# 正确处理预测结果
pred_data <- pred_data %>%
  mutate(
    prob = as.numeric(pred),
    se = sqrt(attributes(pred)$var) %>% as.numeric(),  # 正确提取标准误
    lower = pmax(prob - 1.96 * se, 0),  # 确保概率不小于0
    upper = pmin(prob + 1.96 * se, 1)   # 确保概率不大于1
  ) %>%
  filter(complete.cases(.))

# 找到最低风险点
min_risk <- pred_data %>% 
  filter(prob == min(prob, na.rm = TRUE)) %>%
  slice(1)  # 如果有多个取第一个

# 创建主图
main_plot <- ggplot(pred_data, aes(x = sleep_duration, y = prob)) +
  geom_line(linewidth = 1.2, color = "#2b8cbe") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#56B4E9") +
  geom_rug(data = nhanes_final, aes(x = sleep_duration, y = NULL),
           sides = "b", alpha = 0.1) +
  geom_point(data = min_risk, aes(x = sleep_duration, y = prob),
             size = 4, color = "#e34a33") +
  geom_label_repel(
    data = min_risk,
    aes(label = sprintf("最低风险:\n%.1f小时\n%.1f%%",
                        sleep_duration, prob*100)),
    nudge_y = -0.05,
    color = "#e34a33",
    size = 5
  ) +
  labs(
    title = "NHANES加权分析: 睡眠时长与代谢综合征风险",
    subtitle = "限制性立方样条模型 (95% CI)",
    x = "睡眠时长 (小时/天)",
    y = "代谢综合征患病概率"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(4, 12, 1)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 创建加权直方图
hist_data <- svyhist(~sleep_duration, design = nhanes_design)

hist_plot <- ggplot() +
  geom_bar(aes(x = hist_data$breaks[-length(hist_data$breaks)], 
               y = hist_data$counts),
           stat = "identity", fill = "#E69F00", alpha = 0.7) +
  labs(x = "睡眠时长 (小时/天)", y = "加权计数") +
  theme_minimal(base_size = 12)

# 合并图形
combined_plot <- main_plot / hist_plot + 
  plot_layout(heights = c(3, 1)) 
  # plot_annotation(tag_levels = 'A')

print(combined_plot)

# 7. 模型诊断
# 校准图
val.prob(predict(weighted_model, type = "response"), 
         as.numeric(weighted_model$y),
         smooth = TRUE,
         logistic.cal = TRUE)

# 8. 保存结果
# 保存图形
ggsave("weighted_sleep_metS_analysis.png", 
       combined_plot, 
       width = 10, height = 8, dpi = 600)

# 保存模型结果
write_csv(
  tidy(weighted_model, conf.int = TRUE) %>% 
    mutate(across(where(is.numeric), ~round(., 4))),
  "weighted_model_results.csv"
)