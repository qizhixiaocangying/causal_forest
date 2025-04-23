# 清理环境
rm(list = ls())

# 加载必要的包
library(grf)
library(dplyr)
library(ggplot2)
library(survey)
library(cobalt)
library(DiagrammeR)
source("best_tree.R")

# 1. 数据准备 ----------------------------------------------------------------

# 加载数据（替换为实际数据路径）
load("nhanes_final.rda")

# 数据预处理函数
preprocess_data <- function(data) {
  factorVars <- c("metS_status", "waist_risk", "trigly_risk", "hdl_risk", 
                  "glucose_risk", "bp_risk", "PA_level", "Gender", 
                  "RaceEthnicity", "Education", "SmokingStatus", "AlcoholCat")
  
  data[factorVars] <- lapply(data[factorVars], factor)
  
  data %>%
    mutate(
      SmokingStatus_num = case_when(
        SmokingStatus == "Never Smoker" ~ 0,
        SmokingStatus == "Former Smoker" ~ 1,
        SmokingStatus == "Current Smoker" ~ 2,
        TRUE ~ NA_real_
      ),
      AlcoholCat_num = case_when(
        AlcoholCat == "Unknown/None" ~ 0,
        AlcoholCat == "Light" ~ 1,
        AlcoholCat == "Moderate" ~ 2,
        AlcoholCat == "Heavy" ~ 3,
        TRUE ~ NA_real_
      ),
      PA_level_num = case_when(
        PA_level == "Low" ~ 1,
        PA_level == "Medium" ~ 2,
        PA_level == "High" ~ 3,
        TRUE ~ NA_real_ 
      ),
      Education_num = case_when(
        Education == "Less than 9th grade" ~ 1,
        Education == "9-11th grade" ~ 2,
        Education == "High school graduate" ~ 3,
        Education == "Some college" ~ 4,
        Education == "College graduate" ~ 5,
        TRUE ~ NA_real_
      ),
      Gender_num = case_when(
        Gender == "Male" ~ 1,
        Gender == "Female" ~ 2,
        TRUE ~ NA_real_
      )
    )
}

nhanes_final <- preprocess_data(nhanes_final)

# 2. 因果森林分析函数（不含双重稳健估计）----------------------------------

run_cf_analysis <- function(data, treatment_var, outcome_var, covariates, 
                            model_name = "", min_prop = 0.05) {
  
  # 准备数据
  df <- data %>%
    filter(!is.na(!!sym(treatment_var)), 
           !is.na(!!sym(outcome_var))) %>%
    mutate(W = as.numeric(!!sym(treatment_var)))
  
  Y <- as.numeric(as.character(df[[outcome_var]]))
  W <- df$W
  X <- model.matrix(as.formula(paste("~", paste(covariates, collapse = "+"))), 
                    data = df)[,-1] # 移除截距项
  weights <- df$WGT4YR
  clusters <- df$SDMVPSU
  
  # 倾向得分估计（带修剪）
  ps_model <- regression_forest(X, W, 
                                sample.weights = weights,
                                clusters = clusters,
                                min.node.size = 50)
  
  prop_scores <- predict(ps_model)$predictions
  prop_scores <- pmax(prop_scores, 0.05)  # 下限5%
  prop_scores <- pmin(prop_scores, 0.95)  # 上限95%
  
  # 计算IPW权重
  ipw_weights <- ifelse(W == 1, 1/prop_scores, 1/(1-prop_scores))
  stabilized_weights <- ifelse(W == 1, mean(W)/prop_scores, mean(1-W)/(1-prop_scores))
  
  # 诊断信息
  weight_diagnosis <- list(
    raw_weights = summary(ipw_weights),
    stabilized_weights = summary(stabilized_weights),
    effective_sample_size = list(
      original = length(W),
      weighted = sum(ipw_weights)^2 / sum(ipw_weights^2)
    )
  )
  
  # 重叠样本筛选
  overlap_idx <- which(prop_scores > min_prop & prop_scores < (1 - min_prop))
  
  # 构建三种模型
  models <- list()
  
  # 标准因果森林
  models$standard <- causal_forest(
    X = X,
    Y = Y,
    W = W,
    W.hat = prop_scores,
    sample.weights = weights,
    clusters = clusters,
    num.trees = 2000,
    min.node.size = 100,
    honesty = TRUE,
    honesty.fraction = 0.5
  )
  
  # IPW加权因果森林
  models$ipw <- causal_forest(
    X = X,
    Y = Y,
    W = W,
    W.hat = prop_scores,
    sample.weights = weights * ipw_weights,
    clusters = clusters,
    num.trees = 2000,
    min.node.size = 100,
    honesty = TRUE,
    honesty.fraction = 0.5
    
  )
  
  # 重叠样本因果森林
  models$overlap <- causal_forest(
    X = X[overlap_idx, ],
    Y = Y[overlap_idx],
    W = W[overlap_idx],
    W.hat = prop_scores[overlap_idx],
    sample.weights = weights[overlap_idx],
    clusters = clusters[overlap_idx],
    num.trees = 2000,
    min.node.size = 100,
    honesty = TRUE,
    honesty.fraction = 0.5
  )
  
  # 计算ATE（带错误处理）
  compute_ates <- function(cf) {
    list(
      all = tryCatch(average_treatment_effect(cf), 
                     error = function(e) c(estimate = NA, std.err = NA)),
      treated = tryCatch(average_treatment_effect(cf, target.sample = "treated"),
                         error = function(e) c(estimate = NA, std.err = NA))
    )
  }
  
  ate_results <- lapply(models, compute_ates)
  
  # 变量重要性（基于IPW模型）
  var_imp <- variable_importance(models$ipw)
  var_imp_df <- data.frame(
    variable = colnames(X),
    importance = var_imp[,1]
  ) %>% arrange(desc(importance))
  
  # 返回结果
  list(
    model_name = model_name,
    models = models,
    prop_scores = prop_scores,
    ipw_weights = ipw_weights,
    weight_diagnosis = weight_diagnosis,
    ate_results = ate_results,
    var_imp = var_imp_df,
    overlap_idx = overlap_idx,
    data = df,
    covariates = covariates,
    best_tree_info = find_best_tree(models$ipw, "causal")  # 添加这行
  )
}

# 3. 运行分析 --------------------------------------------------------------

# 定义协变量组
short_vars <- c("PA_level_num", "DQ_proxy_score", "Age", "Gender_num", 
                "SmokingStatus_num", "AlcoholCat_num")
long_vars <- c(short_vars, "PHQ9_score")

# 准备数据
short_optimal_data <- nhanes_final %>%
  filter(sleep_cat %in% c("Short sleep (<6h)", "Optimal sleep (6-9h)")) %>%
  mutate(short_sleep = as.numeric(sleep_cat == "Short sleep (<6h)")) 

long_optimal_data <- nhanes_final %>%
  filter(sleep_cat %in% c("Long sleep (>9h)", "Optimal sleep (6-9h)")) %>%
  mutate(long_sleep = as.numeric(sleep_cat == "Long sleep (>9h)")) 

# 执行分析
result_short <- run_cf_analysis(
  data = short_optimal_data,
  treatment_var = "short_sleep",
  outcome_var = "metS_status",
  covariates = short_vars,
  model_name = "Short vs Optimal Sleep"
)

result_long <- run_cf_analysis(
  data = long_optimal_data,
  treatment_var = "long_sleep",
  outcome_var = "metS_status",
  covariates = long_vars,
  model_name = "Long vs Optimal Sleep"
)

# 4. 结果可视化 ----------------------------------------------------------

# 倾向得分分布图
plot_prop_scores <- function(result) {
  ggplot(data.frame(prop_score = result$prop_scores, W = factor(result$data$W)), 
         aes(x = prop_score, fill = W)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = c(0.05, 0.95), linetype = "dashed", color = "red") +
    ggtitle(paste("Propensity Scores:", result$model_name)) +
    scale_fill_discrete(labels = c("Control", "Treated")) +
    theme_bw()
}

plot_prop_scores(result_short)
plot_prop_scores(result_long)

# IPW权重分布图
plot_ipw_weights <- function(result) {
  ggplot(data.frame(ipw_weight = result$ipw_weights), aes(x = ipw_weight)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    scale_x_log10() +
    ggtitle(paste("IPW Weights Distribution:", result$model_name)) +
    xlab("IPW Weight (log scale)") +
    theme_bw()
}

plot_ipw_weights(result_short)
plot_ipw_weights(result_long)

# ATE比较图
plot_ate_comparison <- function(result) {
  df <- data.frame(
    method = rep(c("Standard", "IPW", "Overlap"), each = 2),
    sample = rep(c("All", "Treated"), 3),
    estimate = c(
      result$ate_results$standard$all["estimate"],
      result$ate_results$standard$treated["estimate"],
      result$ate_results$ipw$all["estimate"],
      result$ate_results$ipw$treated["estimate"],
      result$ate_results$overlap$all["estimate"],
      result$ate_results$overlap$treated["estimate"]
    ),
    std.error = c(
      result$ate_results$standard$all["std.err"],
      result$ate_results$standard$treated["std.err"],
      result$ate_results$ipw$all["std.err"],
      result$ate_results$ipw$treated["std.err"],
      result$ate_results$overlap$all["std.err"],
      result$ate_results$overlap$treated["std.err"]
    )
  )
  
  ggplot(df, aes(x = method, y = estimate, color = sample)) +
    geom_point(position = position_dodge(width = 0.3), size = 3) +
    geom_errorbar(
      aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error),
      width = 0.2, position = position_dodge(width = 0.3)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste("ATE Comparison:", result$model_name),
         y = "ATE Estimate with 95% CI") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_ate_comparison(result_short)
plot_ate_comparison(result_long)

# 5. 结果报告 --------------------------------------------------------------

generate_results_report <- function(result) {
  cat("\n=== ", result$model_name, " Analysis Results ===\n")
  
  # 样本信息
  cat("\nSample Information:\n")
  cat("- Total sample size:", nrow(result$data), "\n")
  cat("- Treated cases:", sum(result$data$W), "\n")
  cat("- Control cases:", sum(!result$data$W), "\n")
  cat("- Overlap sample size:", length(result$overlap_idx), "\n")
  cat("- Effective sample size (IPW):", 
      round(result$weight_diagnosis$effective_sample_size$weighted), "\n")
  
  # 倾向得分诊断
  cat("\nPropensity Score Diagnostics:\n")
  cat("- Min:", round(min(result$prop_scores), 3), "\n")
  cat("- Max:", round(max(result$prop_scores), 3), "\n")
  cat("- Mean (Treated):", round(mean(result$prop_scores[result$data$W == 1]), 3), "\n")
  cat("- Mean (Control):", round(mean(result$prop_scores[result$data$W == 0]), 3), "\n")
  
  # ATE结果
  cat("\nAverage Treatment Effects:\n")
  ate_df <- data.frame(
    Method = c("Standard", "IPW", "Overlap"),
    All = sprintf("%.3f (%.3f)", 
                  c(result$ate_results$standard$all["estimate"],
                    result$ate_results$ipw$all["estimate"],
                    result$ate_results$overlap$all["estimate"]),
                  c(result$ate_results$standard$all["std.err"],
                    result$ate_results$ipw$all["std.err"],
                    result$ate_results$overlap$all["std.err"])),
    Treated = sprintf("%.3f (%.3f)", 
                      c(result$ate_results$standard$treated["estimate"],
                        result$ate_results$ipw$treated["estimate"],
                        result$ate_results$overlap$treated["estimate"]),
                      c(result$ate_results$standard$treated["std.err"],
                        result$ate_results$ipw$treated["std.err"],
                        result$ate_results$overlap$treated["std.err"]))
  )
  print(ate_df, row.names = FALSE)
  
  # 变量重要性
  cat("\nTop 5 Important Variables:\n")
  print(head(result$var_imp, 5))
}

generate_results_report(result_short)
generate_results_report(result_long)

# 6. 保存结果 --------------------------------------------------------------

save(result_short, result_long, file = "causal_forest_results.RData")

# 导出CSV结果
write.csv(
  bind_rows(
    mutate(result_short$var_imp, model = "Short vs Optimal"),
    mutate(result_long$var_imp, model = "Long vs Optimal")
  ),
  "variable_importance_results.csv",
  row.names = FALSE
)

# 导出ATE结果
ate_export <- bind_rows(
  data.frame(
    model = "Short vs Optimal",
    method = rep(c("Standard", "IPW", "Overlap"), each = 2),
    sample = rep(c("All", "Treated"), 3),
    estimate = c(
      result_short$ate_results$standard$all["estimate"],
      result_short$ate_results$standard$treated["estimate"],
      result_short$ate_results$ipw$all["estimate"],
      result_short$ate_results$ipw$treated["estimate"],
      result_short$ate_results$overlap$all["estimate"],
      result_short$ate_results$overlap$treated["estimate"]
    ),
    std.error = c(
      result_short$ate_results$standard$all["std.err"],
      result_short$ate_results$standard$treated["std.err"],
      result_short$ate_results$ipw$all["std.err"],
      result_short$ate_results$ipw$treated["std.err"],
      result_short$ate_results$overlap$all["std.err"],
      result_short$ate_results$overlap$treated["std.err"]
    )
  ),
  data.frame(
    model = "Long vs Optimal",
    method = rep(c("Standard", "IPW", "Overlap"), each = 2),
    sample = rep(c("All", "Treated"), 3),
    estimate = c(
      result_long$ate_results$standard$all["estimate"],
      result_long$ate_results$standard$treated["estimate"],
      result_long$ate_results$ipw$all["estimate"],
      result_long$ate_results$ipw$treated["estimate"],
      result_long$ate_results$overlap$all["estimate"],
      result_long$ate_results$overlap$treated["estimate"]
    ),
    std.error = c(
      result_long$ate_results$standard$all["std.err"],
      result_long$ate_results$standard$treated["std.err"],
      result_long$ate_results$ipw$all["std.err"],
      result_long$ate_results$ipw$treated["std.err"],
      result_long$ate_results$overlap$all["std.err"],
      result_long$ate_results$overlap$treated["std.err"]
    )
  )
)

write.csv(ate_export, "ate_results.csv", row.names = FALSE)






plot(tree <- get_tree(result_short$models$ipw, result_short$best_tree_info$best_tree))
plot(tree <- get_tree(result_long$models$ipw, result_long$best_tree_info$best_tree))

