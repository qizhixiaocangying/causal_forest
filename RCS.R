# Load required packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, rms, ggplot2, broom, ggeffects, patchwork)

# 读取数据
# nhanes_final <- readxl::read_excel("matched_data.xlsx")
load("nhanes_final.rda")


# 1. 构建复杂抽样设计对象
nhanes_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WGT4YR,
  data = nhanes_final,
  nest = TRUE
)


# 计算加权患病率
svymean(~metS_status, design = nhanes_design, na.rm = TRUE)

# 按睡眠时长分组的加权患病率
svyby(~metS_status, ~cut(sleep_duration, breaks = c(4,6,8,10,12)), 
      design = nhanes_design, 
      svymean)


# Complete Analysis: Sleep Duration and Metabolic Syndrome
# Using Restricted Cubic Splines (RCS)


ddist <- datadist(nhanes_final)
options(datadist = "ddist")

# Fit restricted cubic spline model with 4 knots
spline_model <- lrm(metS_status ~ rcs(sleep_duration, 4), 
                    data = nhanes_final, 
                    x = TRUE, y = TRUE)


# 5. 模型评估 ----------------------------------------------------------


# 6. 可视化加权关系 ----------------------------------------------------
# 生成预测数据
pred_data <- expand.grid(
  sleep_duration = seq(4, 12, length.out = 200)
)

# 6. 可视化加权关系 (修正版) --------------------------------------------


# 3. Model Evaluation -----------------------------------------------------
# Print model summary
print(spline_model)

# Test nonlinearity
anova_results <- anova(spline_model)
print(anova_results)


# 4. Visualization --------------------------------------------------------
# Generate prediction data frame
pred_data <- Predict(spline_model, 
                     sleep_duration = seq(4, 12, length.out = 200),
                     conf.int = 0.95)

# Convert to probability scale
pred_data <- pred_data %>% 
  as_tibble() %>%
  mutate(across(c(yhat, lower, upper), plogis))

# 找到最低点索引
min_idx <- which.min(pred_data$yhat)
# 提取最低点坐标
min_point <- pred_data[min_idx, ]

# Create main plot
main_plot <- ggplot(pred_data, aes(x = sleep_duration, y = yhat)) +
  geom_line(linewidth = 1.2, color = "#0072B2") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#56B4E9") +
  geom_rug(data = nhanes_final, aes(x = sleep_duration, y = NULL),
           sides = "b", alpha = 0.1) +
  labs(title = "NHANES未加权分析: 睡眠时长与代谢综合征风险",
           subtitle = "限制性立方样条模型 (95% CI)",
           x = "睡眠时长 (小时/天)",
           y = "代谢综合征患病概率") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(4, 12, 1)) +
  geom_vline(xintercept = c(6, 9), linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = 6, y = max(pred_data$yhat) * 0.95, 
           label = "6h", color = "red", hjust = 1.1, size = 5) +
  annotate("text", x = 9, y = max(pred_data$yhat) * 0.95, 
           label = "9h", color = "red", hjust = -0.1, size = 5)+
  geom_point(data = min_point, aes(x = sleep_duration, y = yhat), 
             color = "darkred", size = 4, shape = 21, fill = "yellow") +
  annotate("text", 
           x = min_point$sleep_duration, 
           y = min_point$yhat, 
           label = paste0(round(min_point$sleep_duration, 2), "h\n", 
                          scales::percent(min_point$yhat, accuracy = 0.1)),
           color = "darkred", 
           hjust = -0.1, vjust = -1, size = 5)

# Create histogram of sleep duration
hist_plot <- ggplot(nhanes_final, aes(x = sleep_duration)) +
  geom_histogram(fill = "#E69F00", bins = 10, alpha = 0.7) +
  labs(x = "Sleep Duration (hours/day)", y = "未加权计数") +
  theme_minimal(base_size = 12)

# Combine plots
combined_plot <- main_plot / hist_plot + 
  plot_layout(heights = c(3, 1))

print(combined_plot)

# 5. Model Diagnostics ----------------------------------------------------
# Calibration plot
val.prob(predict(spline_model, type = "fitted"), 
         spline_model$y,
         smooth = TRUE,
         logistic.cal = TRUE)

# 6. Save Results ---------------------------------------------------------
# Save plot
ggsave("sleep_metabolic_syndrome_rcs.png", 
       combined_plot, 
       width = 10, height = 8, dpi = 300)


