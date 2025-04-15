# 结果解读与可视化

# 1. 平均处理效应 (ATE for Short vs Normal)
ate_short <- average_treatment_effect(cf_short_vs_normal, target.sample = "all")
print(paste("Estimated ATE (Short vs Normal sleep on MetS):", round(ate_short[1], 4),
            "SE:", round(ate_short[2], 4)))

# 2. CATE 分布 (异质性)
results_df_short <- data.frame(
  tau_hat = tau_hat_short,
  tau_se = sqrt(tau_hat_var_short)
)
# 添加协变量用于绘图
results_df_short <- bind_cols(results_df_short, as.data.frame(X_df_binary))
results_df_short$weights <- weights_binary # 添加权重用于加权绘图

# 加权直方图/密度图
ggplot(results_df_short, aes(x = tau_hat, weight = weights)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightblue", color = "black") +
  geom_density(color = "red") +
  labs(title = "Distribution of Estimated CATEs (Short vs Normal Sleep Effect)", x = "Estimated CATE (τ̂)") +
  geom_vline(xintercept = ate_short[1], linetype = "dashed", color = "blue") + # Add ATE line
  theme_minimal()

# 3. 变量重要性 (驱动异质性的因素)
var_imp <- variable_importance(cf_short_vs_normal)
# 绘制重要性排序图
var_imp_df <- data.frame(Variable = names(var_imp), Importance = var_imp) %>% arrange(desc(Importance))
ggplot(head(var_imp_df, 15), aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance for CATE Heterogeneity", x = "Variable", y = "Importance") +
  theme_minimal()
# 预期 PA_level, DQ_score, BMI, Age 等会比较重要

# 4. 可视化 CATE vs 关键调节变量 (PA, DQ)
# 需要将因子PA_level转换回来或使用哑变量
results_df_short$PA_level <- nhanes_binary_W$PA_level # Add original factor back

# CATE vs PA_level (Boxplot)
ggplot(results_df_short, aes(x = PA_level, y = tau_hat, fill = PA_level)) +
  geom_boxplot(aes(weight = weights)) + # Note: weighted boxplots might need specific packages or methods
  labs(title = "Estimated CATE (Short vs Normal) by Physical Activity Level", x = "PA Level", y = "Estimated CATE (τ̂)") +
  theme_minimal()

# CATE vs DQ_score (Scatter plot with smooth)
ggplot(results_df_short, aes(x = DQ_score, y = tau_hat)) +
  geom_point(alpha = 0.3, aes(size = weights)) + # Size points by weight (visual aid)
  geom_smooth(method = "loess", color = "red", aes(weight = weights)) + # Weighted smooth
  labs(title = "Estimated CATE (Short vs Normal) vs Diet Quality Score", x = "Diet Quality Score", y = "Estimated CATE (τ̂)") +
  theme_minimal()

# 5. 子群分析 (示例: BLP)
# - 找出哪些变量能最好地线性预测CATE
blp_results <- best_linear_projection(cf_short_vs_normal, X_binary)
print("Best Linear Projection Coefficients:")
print(blp_results) # 查看 PA, DQ, BMI 等的系数

# 重复类似分析比较 Long vs Normal sleep

# cox建模验证