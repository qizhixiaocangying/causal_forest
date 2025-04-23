
# 基线分析函数
perform_baseline_analysis <- function(data, columns_to_extract, output_file = "table_before_matching.csv") {
  cat("\n开始进行基线分析...\n")
  nhanes_design <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WGT4YR,
    data = data,
    nest = TRUE
  )
  table_before <- svyCreateTableOne(vars = columns_to_extract,
                                 #testExact = fisher.test,
                                 strata = "sleep_cat",
                                 data = nhanes_design,
                                 test = TRUE,
                                 addOverall = TRUE)
  
  # table_before <- CreateTableOne(vars = columns_to_extract, 
  #                                   testExact = fisher.test,
  #                                   strata = "sleep_cat", 
  #                                   data = data, 
  #                                   test = TRUE,
  #                                   addOverall = TRUE)
  
  print(table_before, smd = TRUE, showAllLevels = TRUE)
  
  # 转换为可打印格式（包括SMD）
  table_print <- print(table_before, 
                       smd = TRUE, 
                       printToggle = FALSE,
                       quote = FALSE, 
                       noSpaces = TRUE, 
                       exact = names(table_before$CatVars),
                       showAllLevels = TRUE)
  
  # 导出到Excel
  # 安装并加载readr包

  write.csv(table_print, file = output_file, fileEncoding = "UTF-8")
  cat(paste("\n基线分析表已导出到 '", output_file, "'\n", sep=""))
  
  return(table_before)
}

# 倾向性评分计算函数
calculate_propensity_scores <- function(data, ps_variables) {
  # 确保sleep_cat变量已因子化
  data$sleep_cat <- factor(data$sleep_cat)
  
  # 构建公式
  ps_formula <- as.formula(paste("sleep_cat ~", paste(ps_variables, collapse = " + ")))
  
  # 使用多项逻辑回归估计倾向性评分
  cat("\n开始进行倾向性评分估计...\n")
  multinom_model <- multinom(ps_formula, data = data, trace = FALSE)
  
  # 提取预测概率
  pred_probs <- predict(multinom_model, type = "probs")
  colnames(pred_probs) <- paste0("ps_", levels(data$sleep_cat))
  data <- cbind(data, pred_probs)
  
  # 查看结果
  cat("\n倾向性评分估计完成。前几行结果：\n")
  print(head(data[, c("sleep_cat", colnames(pred_probs))]))
  
  # 计算各组间距离
  sleep_cats <- levels(data$sleep_cat)
  if(length(sleep_cats) >= 3) {
    # 为每对处理组计算距离
    for(i in 1:(length(sleep_cats)-1)) {
      for(j in (i+1):length(sleep_cats)) {
        dist_name <- paste0("dist_", sleep_cats[i], sleep_cats[j])
        data[[dist_name]] <- calculate_gps_distance(
          data[, paste0("ps_", sleep_cats)], 
          sleep_cats[i], 
          sleep_cats[j]
        )
      }
    }
  }
  
  return(data)
}

# 广义倾向性评分距离计算函数
calculate_gps_distance <- function(ps_matrix, ref_sleep_cat, comp_sleep_cat) {
  ref_idx <- which(colnames(ps_matrix) == paste0("ps_", ref_sleep_cat))
  comp_idx <- which(colnames(ps_matrix) == paste0("ps_", comp_sleep_cat))
  
  # 计算马氏距离
  d <- sqrt((ps_matrix[, ref_idx] - ps_matrix[, comp_idx])^2 / 
              (ps_matrix[, ref_idx] * (1 - ps_matrix[, ref_idx])))
  return(d)
}

# 多组倾向性评分匹配函数
match_multi_ratio <- function(data, treatment_var = "sleep_cat", caliper = 0.1, ratio = c(2,1,1)) {
  # 获取组别
  sleep_cats <- levels(data[[treatment_var]])
  
  if(length(sleep_cats) < 3) {
    stop("需要至少3个组进行多组匹配")
  }
  
  if(length(ratio) != length(sleep_cats)) {
    stop("匹配比例数量必须与组别数量相同")
  }
  
  # 分离各处理组
  sleep_cat_data <- list()
  for(g in sleep_cats) {
    sleep_cat_data[[g]] <- data[data[[treatment_var]] == g, ]
  }
  
  # 初始化匹配结果
  matched_indices <- list()
  
  # 确定哪个组是参考组（比例小的组）
  ref_sleep_cat_idx <- 2
  ref_sleep_cat <- sleep_cats[ref_sleep_cat_idx]
  ref_data <- sleep_cat_data[[ref_sleep_cat]]
  
  # 对每个参考组的观察值找到其他组中的匹配
  for(i in 1:nrow(ref_data)) {
    # 为每个比较组找到最近的匹配
    matched_from_sleep_cats <- list()
    matched_from_sleep_cats[[ref_sleep_cat]] <- rownames(ref_data)[i]
    
    all_matched <- TRUE
    
    for(g_idx in seq_along(sleep_cats)) {
      if(g_idx == ref_sleep_cat_idx) next  # 跳过参考组
      
      comp_sleep_cat <- sleep_cats[g_idx]
      comp_ratio <- ratio[g_idx]
      
      # 计算距离名称
      if(which(sleep_cats == ref_sleep_cat) < which(sleep_cats == comp_sleep_cat)) {
        dist_name <- paste0("dist_", ref_sleep_cat, comp_sleep_cat)
      } else {
        dist_name <- paste0("dist_", comp_sleep_cat, ref_sleep_cat)
      }
      
      # 计算到比较组的距离
      dist_to_comp <- abs(ref_data[[dist_name]][i] - sleep_cat_data[[comp_sleep_cat]][[dist_name]])
      
      # 找出最小距离的索引（按照指定比例）
      if(length(dist_to_comp) >= comp_ratio) {
        # 找出comp_ratio个最小距离的索引
        min_dist_indices <- order(dist_to_comp)[1:comp_ratio]
        min_dists <- dist_to_comp[min_dist_indices]
        
        # 检查是否所有距离都在caliper内
        if(all(min_dists <= caliper)) {
          matched_from_sleep_cats[[comp_sleep_cat]] <- rownames(sleep_cat_data[[comp_sleep_cat]])[min_dist_indices]
          # 移除已匹配的观察值
          sleep_cat_data[[comp_sleep_cat]] <- sleep_cat_data[[comp_sleep_cat]][-min_dist_indices, ]
        } else {
          all_matched <- FALSE
          break
        }
      } else {
        all_matched <- FALSE
        break
      }
    }
    
    # 如果所有组都找到了匹配，添加到结果中
    if(all_matched) {
      matched_indices[[length(matched_indices) + 1]] <- matched_from_sleep_cats
    }
    
    # 如果某个组没有足够的观察值，提前退出
    insufficient_samples <- FALSE
    for(g_idx in seq_along(sleep_cats)) {
      if(g_idx == ref_sleep_cat_idx) next
      if(nrow(sleep_cat_data[[sleep_cats[g_idx]]]) < ratio[g_idx]) {
        insufficient_samples <- TRUE
        break
      }
    }
    if(insufficient_samples) break
  }
  
  # 转换为数据框
  if(length(matched_indices) > 0) {
    # 展平匹配索引
    all_indices <- unlist(matched_indices)
    
    # 提取匹配的观察值
    matched_data <- data[as.numeric(all_indices), ]
    
    return(list(
      matched_indices = matched_indices,
      matched_data = matched_data
    ))
  } else {
    return(list(
      matched_indices = NULL,
      matched_data = NULL
    ))
  }
}

# 主函数
main <- function() {
  # 清空环境
  rm(list = ls())
  
  # 加载必要的包
  library(MatchIt)      # 用于倾向性评分匹配
  library(nnet)         # 用于多项逻辑回归
  library(survival)     # 用于生存分析
  library(survminer)    # 用于生存曲线可视化
  library(tableone)     # 用于创建描述性统计表
  library(dplyr)        # 用于数据处理
  library(ggplot2)      # 用于可视化
  library(openxlsx)
  library(survey)
  
  # 读取数据
  #data <- readxl::read_excel("data_imputed.xlsx")
  load("nhanes_final.rda")
  # 定义分析所需的列（使用英文列名）
  data <-  nhanes_final
  
  
  # 2. 指定需要汇总的变量（除设计变量外）
  columns_to_extract <- c(#"metS_status", "sleep_duration", "waist_risk", "trigly_risk", "hdl_risk", "bp_risk","glucose_risk", "bedtime_hours",
               "Age", "DQ_proxy_score", "Gender","PA_level", "RaceEthnicity","Education",  "SmokingStatus", "AlcoholCat","PHQ9_score" 
               )
  
  # 3. 指定分类变量（factor变量）
  factorVars <- c(#"metS_status",  "waist_risk", "trigly_risk", "hdl_risk", "glucose_risk", "bp_risk",  
    "Gender", "PA_level", "RaceEthnicity","Education",  "SmokingStatus", "AlcoholCat")
  
  # 因子化分类变量
  data[factorVars] <- lapply(data[factorVars], factor)
  
  # 确保所有需要的列都存在
  missing_cols <- setdiff(columns_to_extract, colnames(data))
  if(length(missing_cols) > 0) {
    warning("以下列在数据中不存在：", paste(missing_cols, collapse=", "))
    # 只保留存在的列
    columns_to_extract <- intersect(columns_to_extract, colnames(data))
  }
  
  # 匹配前的基线分析
  table_before <- perform_baseline_analysis(data, columns_to_extract)
  
  # 检查是否存在sleep_cat列
  if("sleep_cat" %in% colnames(data)) {
    # 计算倾向性评分
    data_with_ps <- calculate_propensity_scores(data, columns_to_extract)
    
    # 获取组别
    sleep_cats <- levels(data_with_ps$sleep_cat)
    
    if(length(sleep_cats) >= 3) {
      # 确定匹配比例
      match_ratio <- rep(1, length(sleep_cats))
      # 第一组是要进行4:1匹配的组
      match_ratio[1] <-1
      match_ratio[2] <-1
      match_ratio[3] <-1
      
      # 进行多组匹配
      cat("\n开始进行多组倾向性评分匹配...\n")
      multi_match_result <- match_multi_ratio(data_with_ps, treatment_var = "sleep_cat", caliper = 0.1, ratio = match_ratio)
      multi_match_result$matched_indices
      # 处理匹配结果
      if(!is.null(multi_match_result$matched_data)) {
        matched_data <- multi_match_result$matched_data
        cat("多组匹配后的样本量:", nrow(matched_data), "\n")
        cat("匹配的组合数量:", length(multi_match_result$matched_indices), "\n")
        
        # 检查每个组的样本量
        sleep_cat_counts <- table(matched_data$sleep_cat)
        cat("匹配后各组的样本量:\n")
        print(sleep_cat_counts)
        write.xlsx(as.data.frame(matched_data),
                   file = "matched_data.xlsx",
                   rowNames = FALSE)
        cat(paste("\n匹配数据已保存到", "matched_data.xlsx", "\n"))
        # 匹配后的基线分析
        table_after <- perform_baseline_analysis(matched_data, columns_to_extract, "table_after_matching.csv")
      } else {
        cat("\n未找到符合条件的匹配。请尝试调整caliper值或检查数据。\n")
      }
    } else {
      warning("需要至少3个组进行多组匹配，当前只有", length(sleep_cats), "个组")
    }
  } else {
    warning("数据中不存在'sleep_cat'列，无法进行倾向性评分估计")
  }
  
  cat("\n分析完成！\n")
}

# 运行主函数
main()
