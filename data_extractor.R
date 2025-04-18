# 数据提取
library(nhanesA)
library(tidyverse)

print("--- 开始下载 NHANES 2015-2016 ('I') 和 2017-2018 ('J') 数据 ---")

# 1. 显式下载每个周期的数据文件

# --- 定义所需变量列表 (省略周期后缀) ---
demo_vars_needed <- c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2", "INDFMPIR", "RIDEXPRG", "SDMVPSU", "SDMVSTRA", "WTMEC2YR")
sleep_vars_needed <- c("SEQN", "SLQ300", "SLQ310", "SLQ320", "SLQ330")
bmx_vars_needed <- c("SEQN", "BMXWAIST", "BMXHT", "BMXWT", "BMXBMI")
trigly_vars_needed <- c("SEQN", "LBXTR")
hdl_vars_needed <- c("SEQN", "LBDHDD")
bpx_vars_needed <- c("SEQN", "BPXSY1", "BPXDI1")
glu_vars_needed <- c("SEQN", "LBXGLU")
biopro_vars_needed <- c("SEQN", "LBXSCR")
paq_vars_needed <- c("SEQN", "PAQ605", "PAQ610", "PAQ615", "PAQ620", "PAQ625", "PAQ630", "PAQ650", "PAQ655", "PAQ660", "PAQ665", "PAQ670", "PAQ675", "PAD680")
diet_tot_needed <- c("SEQN", "DR1TKCAL", "DR1TSODI", "DR1TPOTA", "DR1TFIBE")
smq_vars_needed <- c("SEQN", "SMQ020", "SMQ040")
alq_vars_needed <- c("SEQN", "ALQ130")
dpq_vars_needed <- c("SEQN", "DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", "DPQ060", "DPQ070", "DPQ080", "DPQ090") # 直接列出

# --- 下载周期 'I' (2015-2016) ---
cat("下载周期 I 数据...\n")
demo_I <- nhanes("DEMO_I") %>% select(any_of(demo_vars_needed)) %>% mutate(cycle = "2015-2016")
sleep_I <- nhanes("SLQ_I") %>% select(any_of(sleep_vars_needed))
bmx_I <- nhanes("BMX_I") %>% select(any_of(bmx_vars_needed))
trigly_I <- nhanes("TRIGLY_I") %>% select(any_of(trigly_vars_needed))
hdl_I <- nhanes("HDL_I") %>% select(any_of(hdl_vars_needed))
bpx_I <- nhanes("BPX_I") %>% select(any_of(bpx_vars_needed))
glu_I <- nhanes("GLU_I") %>% select(any_of(glu_vars_needed))
biopro_I <- nhanes("BIOPRO_I") %>% select(any_of(biopro_vars_needed))
paq_I <- nhanes("PAQ_I") %>% select(any_of(paq_vars_needed))
diet1_tot_I <- nhanes("DR1TOT_I") %>% select(any_of(diet_tot_needed))
smq_I <- nhanes("SMQ_I") %>% select(any_of(smq_vars_needed))
alq_I <- nhanes("ALQ_I") %>% select(any_of(alq_vars_needed))
dpq_I <- nhanes("DPQ_I") %>% select(any_of(dpq_vars_needed))

# --- 下载周期 'J' (2017-2018) ---
cat("下载周期 J 数据...\n")
demo_J <- nhanes("DEMO_J") %>% select(any_of(demo_vars_needed)) %>% mutate(cycle = "2017-2018")
sleep_J <- nhanes("SLQ_J") %>% select(any_of(sleep_vars_needed))
bmx_J <- nhanes("BMX_J") %>% select(any_of(bmx_vars_needed))
trigly_J <- nhanes("TRIGLY_J") %>% select(any_of(trigly_vars_needed))
hdl_J <- nhanes("HDL_J") %>% select(any_of(hdl_vars_needed))
bpx_J <- nhanes("BPX_J") %>% select(any_of(bpx_vars_needed))
glu_J <- nhanes("GLU_J") %>% select(any_of(glu_vars_needed))
biopro_J <- nhanes("BIOPRO_J") %>% select(any_of(biopro_vars_needed))
paq_J <- nhanes("PAQ_J") %>% select(any_of(paq_vars_needed))
diet1_tot_J <- nhanes("DR1TOT_J") %>% select(any_of(diet_tot_needed))
smq_J <- nhanes("SMQ_J") %>% select(any_of(smq_vars_needed))
alq_J <- nhanes("ALQ_J") %>% select(any_of(alq_vars_needed))
dpq_J <- nhanes("DPQ_J") %>% select(any_of(dpq_vars_needed))

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


# 3. 合并所有数据到一个数据框 (与之前方法相同)
data_files <- list(demo_data, sleep_data, bmx_data, trigly_data, hdl_data, bpx_data,
                   glu_data, biopro_data, paq_data, diet1_data, smq_data, alq_data, dpq_data)
# 基础数据合并
nhanes_merged <- data_files %>% reduce(full_join, by = "SEQN")

print(paste("初步合并数据完成。维度:", dim(nhanes_merged)[1], "rows,", dim(nhanes_merged)[2], "columns."))

# 保存数据
save(nhanes_merged, file = "nhanes_merged.rda")