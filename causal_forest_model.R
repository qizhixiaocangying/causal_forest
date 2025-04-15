# 因果森林建模

# 训练因果森林模型
# 注意: grf 的 causal_forest 对于多分类处理变量的处理方式可能需要查阅文档。
# 一种常见方法是进行一系列的二分类比较，或使用支持多分类的扩展。
# 这里以比较 Short vs Normal 为例。

# 准备二分类 W (Short vs Normal)
nhanes_binary_W <- nhanes_final %>%
  filter(sleep_cat %in% c("Normal (7-8h)", "Short (<6h)")) %>%
  mutate(W_binary = ifelse(sleep_cat == "Short (<6h)", 1, 0)) # 1=Short, 0=Normal

Y_binary <- nhanes_binary_W$metS_status
W_binary <- nhanes_binary_W$W_binary
X_df_binary <- nhanes_binary_W %>% select(PA_level, DQ_score, BMI, Age, IsFemale, /*...*/)
X_binary <- model.matrix(~ . - 1, data = X_df_binary)
weights_binary <- nhanes_binary_W$WTMEC2YR
clusters_binary <- nhanes_binary_W$SDMVPSU # Or SEQN

# 训练模型 (Short vs Normal)
# 可能需要调整参数 num.trees, min.node.size, mtry, alpha, sample.fraction
# 使用 honesty, ci.group.size 等默认设置通常是好的开始
cf_short_vs_normal <- causal_forest(
  X = X_binary,
  Y = Y_binary,
  W = W_binary,
  sample.weights = weights_binary,
  clusters = clusters_binary, # 用于稳健标准误
  num.trees = 2000, # 增加树的数量
  seed = 123
)

# 同样的方法可以比较 Long vs Normal

# 估计 CATEs
# tauhat = E[Y|W=1, X=x] - E[Y|W=0, X=x] (Short vs Normal effect)
cate_estimates_short <- predict(cf_short_vs_normal, estimate.variance = TRUE)
tau_hat_short <- cate_estimates_short$predictions
tau_hat_var_short <- cate_estimates_short$variance.estimates # 用于置信区间

# 检查模型的拟合情况 (可选，如使用 test_calibration)

