# 加载必要的包
library(nhanesA)
library(tidyverse)

print("--- 开始下载 NHANES 2015-2016 ('I') 和 2017-2018 ('J') 数据 ---")

# --- 下载周期 'I' (2015-2016) ---
cat("下载周期 I 数据...\n")
demo_I <- nhanes("DEMO_I") %>% mutate(cycle = "2015-2016")  # 不筛选变量，保留所有列
sleep_I <- nhanes("SLQ_I")
bmx_I <- nhanes("BMX_I")
trigly_I <- nhanes("TRIGLY_I")
hdl_I <- nhanes("HDL_I")
bpx_I <- nhanes("BPX_I")
glu_I <- nhanes("GLU_I")
biopro_I <- nhanes("BIOPRO_I")
paq_I <- nhanes("PAQ_I")
diet1_tot_I <- nhanes("DR1TOT_I")
smq_I <- nhanes("SMQ_I")
alq_I <- nhanes("ALQ_I")
dpq_I <- nhanes("DPQ_I")

# --- 下载周期 'J' (2017-2018) ---
cat("下载周期 J 数据...\n")
demo_J <- nhanes("DEMO_J") %>% mutate(cycle = "2017-2018")  # 不筛选变量，保留所有列
sleep_J <- nhanes("SLQ_J")
bmx_J <- nhanes("BMX_J")
trigly_J <- nhanes("TRIGLY_J")
hdl_J <- nhanes("HDL_J")
bpx_J <- nhanes("BPX_J")
glu_J <- nhanes("GLU_J")
biopro_J <- nhanes("BIOPRO_J")
paq_J <- nhanes("PAQ_J")
diet1_tot_J <- nhanes("DR1TOT_J")
smq_J <- nhanes("SMQ_J")
alq_J <- nhanes("ALQ_J")
dpq_J <- nhanes("DPQ_J")

# 2. 合并周期数据
cat("合并周期 I 和 J 的数据...\n")
demo_data <- bind_rows(demo_I, demo_J)
sleep_data <- bind_rows(sleep_I, sleep_J)
bmx_data <- bind_rows(bmx_I, bmx_J)
trigly_data <- bind_rows(trigly_I, trigly_J)
hdl_data <- bind_rows(hdl_I, hdl_J)
bpx_data <- bind_rows(bpx_I, bpx_J)
glu_data <- bind_rows(glu_I, glu_J)
biopro_data <- bind_rows(biopro_I, biopro_J)
paq_data <- bind_rows(paq_I, paq_J)
diet1_data <- bind_rows(diet1_tot_I, diet1_tot_J)
smq_data <- bind_rows(smq_I, smq_J)
alq_data <- bind_rows(alq_I, alq_J)
dpq_data <- bind_rows(dpq_I, dpq_J)

# 3. 合并所有数据到一个数据框（全连接，保留所有变量）
data_files <- list(demo_data, sleep_data, bmx_data, trigly_data, hdl_data, bpx_data,
                   glu_data, biopro_data, paq_data, diet1_data, smq_data, alq_data, dpq_data)

safe_join <- function(x, y) {
  if (!"SEQN" %in% names(y)) {
    warning("表缺少SEQN列，已跳过")
    return(x)
  }
  full_join(x, y, by = "SEQN")
}

nhanes_full <- reduce(data_files, safe_join)



print(paste("完整合并数据完成。维度:", dim(nhanes_full)[1], "rows,", dim(nhanes_full)[2], "columns."))

# 保存数据
save(nhanes_full, file = "nhanes_full.rda")

