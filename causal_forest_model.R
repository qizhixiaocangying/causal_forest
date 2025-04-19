# 清理环境
rm(list = ls())


# Load necessary libraries
library(grf)
library(dplyr)
library(survey) # If using survey design
library(DiagrammeR)
source("best_tree.R")

# Load your data (replace with your actual data loading)
load("nhanes_final.rda")

nhanes_final <- nhanes_final %>%
  mutate(
    SmokingStatus_num = case_when(
      SmokingStatus == "Never Smoker" ~ 0,
      SmokingStatus == "Former Smoker" ~ 1,
      SmokingStatus == "Current Smoker" ~ 2,
      TRUE ~ NA_real_
    ),
    AlcoholCat_num = case_when(
      AlcoholCat == "None" ~ 0,
      AlcoholCat == "Light" ~ 1,
      AlcoholCat == "Moderate" ~ 2,
      AlcoholCat == "Heavy" ~ 3,
      TRUE ~ NA_real_
    ),
    PA_level_num = case_when(
      PA_level == "Low" ~ 1,
      PA_level == "Medium" ~ 2,
      PA_level == "High" ~3,
      TRUE ~ NA_real_ 
  ))

# *** MODEL 1: Short Sleep vs. Optimal Sleep ***
# 1. Data Preparation for Model 1
nhanes_short_optimal <- nhanes_final %>%
  filter(sleep_cat %in% c("Short sleep (<7h)", "Optimal sleep (7-9h)")) %>%
  mutate(W_short_optimal = ifelse(sleep_cat == "Short sleep (<7h)", 1, 0))

# Define variables for Model 1
Y1 <- nhanes_short_optimal$metS_status
W1 <- nhanes_short_optimal$W_short_optimal
X1_df <- nhanes_short_optimal %>%
  select(bedtime_hours, PA_level, DQ_proxy_score, Age, Gender, RaceEthnicity,
         Education, Income, BMI, SmokingStatus, AlcoholCat, PHQ9_score, eGFR) %>%
  mutate(across(where(is.character), as.factor))



X1 <- model.matrix(~ . - 1, data = X1_df) # creates dummy variables
covars <- colnames(X1)

# Optional: Survey design elements for Model 1
weights1 <- nhanes_short_optimal$WGT4YR
clusters1 <- nhanes_short_optimal$SDMVPSU

# 2. Model Building for Model 1
cf_short_optimal <- causal_forest(
  X = X1,
  Y = Y1,
  W = W1,
  sample.weights = weights1, # if using survey design
  clusters = clusters1,       # if using survey design
  num.trees = 200,
  seed = 123,
  min.node.size = 30,
  honesty = TRUE,
  honesty.fraction = 0.5,　 #  分裂子样本的比例
  tune.parameters = c("mtry", "alpha",
                      "imbalance.penalty",
                      "honesty.prune.leaves")
)

# 3. Prediction and Interpretation for Model 1
cate_short_optimal <- predict(cf_short_optimal)
hist(cate_short_optimal$predictions, main = "CATE: Short vs Optimal Sleep", xlab = "Estimated Treatment Effect") #distribution visualization


# *** MODEL 2: Long Sleep vs. Optimal Sleep ***
# 1. Data Preparation for Model 2
nhanes_long_optimal <- nhanes_final %>%
  filter(sleep_cat %in% c("Long sleep (>9h)", "Optimal sleep (7-9h)")) %>%
  mutate(W_long_optimal = ifelse(sleep_cat == "Long sleep (>9h)", 1, 0))

# Define variables for Model 2
Y2 <- nhanes_long_optimal$metS_status
W2 <- nhanes_long_optimal$W_long_optimal
X2_df <- nhanes_long_optimal %>%
  select(bedtime_hours, PA_level_num, DQ_proxy_score, Age, Gender, RaceEthnicity,
         Education, Income, BMI, SmokingStatus_num, AlcoholCat_num, PHQ9_score, eGFR) %>%
  mutate(across(where(is.character), as.factor))

X2 <- model.matrix(~ . - 1, data = X2_df)
covars <- colnames(X2)
# Optional: Survey design elements for Model 2
weights2 <- nhanes_long_optimal$WGT4YR
clusters2 <- nhanes_long_optimal$SDMVPSU

# 2. Model Building for Model 2
cf_long_optimal <- causal_forest(
  X = X2,
  Y = Y2,
  W = W2,
  sample.weights = weights2, # if using survey design
  clusters = clusters2,       # if using survey design
  num.trees = 2000,
  seed = 123,
  sample.fraction=0.5,　#  构建每棵树时从训练集中抽样的比例
  min.node.size = 30,
  #honesty = TRUE,
  #honesty.fraction = 0.5,　 #  分裂子样本的比例
  )
  


# 3. Prediction and Interpretation for Model 2
cate_long_optimal <- predict(cf_long_optimal)
hist(cate_long_optimal$predictions,  main = "CATE: Long vs Optimal Sleep", xlab = "Estimated Treatment Effect")  #distribution visualization




var_imp <- c(variable_importance(cf_long_optimal))
names(var_imp) <- covars
sort(var_imp, decreasing=TRUE)

#-------------------------------------------------------------------------------------------------------------------------------------------------
# library(grf)


# 找最佳树
best_trees <- find_best_tree(cf_long_optimal, "causal")
plot(tree <- get_tree(cf_long_optimal, best_trees$best_tree))


best_tree_index <- best_trees$best_tree

# 提取最佳树
best_tree <- get_tree(cf_long_optimal, best_tree_index)

# 绘制最佳树
# plot(best_tree)


best_trees <- find_best_tree(cf_long_optimal, "causal")
plot(tree <- get_tree(cf_long_optimal, best_trees$best_tree))

test_calibration(cf_long_optimal)


# # 假设你已经有 cf（causal_forest对象）和 best_trees$best_tree
# tree <- get_tree(cf_long_optimal, best_trees$best_tree)
# plot(tree)        # 图形化显示分裂路径
# str(tree)         # 查看所有节点的分裂变量和阈值
# tree$split_variable  # 各节点分裂用的变量编号
# tree$split_value     # 各节点的分裂阈值
# colnames(X2)         # 变量编号和名字的对应关系

