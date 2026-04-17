library(tidyverse)

setwd("/pastel/tools/circRNA_tools/test_data/decoupled_comparison")
setwd("/pastel/tools/circRNA_tools/test_data/decoupled_bwa_comparison")

# Compare results from joint analysis 
joint_BSJ_res = read_tsv("original/result.BSJ_Matrix")
joint_FSJ_res = read_tsv("original/result.FSJ_Matrix")

# Compare results from decoupled analysis
decoupled_BSJ_res = read_tsv("decoupled/finalize/result.BSJ_Matrix")
decoupled_FSJ_res = read_tsv("decoupled/finalize/result.FSJ_Matrix")

colnames(joint_BSJ_res) = colnames(decoupled_BSJ_res)
colnames(joint_FSJ_res) = colnames(decoupled_FSJ_res)

res_BSJ_merged = full_join(joint_BSJ_res, decoupled_BSJ_res, by = "circRNA_ID", suffix = c(".joint", ".decoupled"))
res_BSJ_merged_long = res_BSJ_merged %>%
  pivot_longer(cols = -circRNA_ID, names_to = c("Sample", "Run"), names_pattern = "(.*)\\.(.*)") %>%
  pivot_wider(names_from = Run, values_from = value)

ggplot(res_BSJ_merged_long, aes(x = joint, y = decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Sample) +
  labs(x = "Joint Analysis BSJ Count", y = "Decoupled Analysis BSJ Count") +
  theme_bw()

res_FSJ_merged = full_join(joint_FSJ_res, decoupled_FSJ_res, by = "circRNA_ID", suffix = c(".joint", ".decoupled"))
res_FSJ_merged_long = res_FSJ_merged %>%
  pivot_longer(cols = -circRNA_ID, names_to = c("Sample", "Run"), names_pattern = "(.*)\\.(.*)") %>%
  pivot_wider(names_from = Run, values_from = value)

ggplot(res_FSJ_merged_long, aes(x = joint, y = decoupled)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Sample) +
  labs(x = "Joint Analysis FSJ Count", y = "Decoupled Analysis FSJ Count") +
  theme_bw()
