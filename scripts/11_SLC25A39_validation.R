#==============================================================================
# Script 11: SLC25A39 Functional Validation
# Description: Analysis of siRNA knockdown experiments
# Corresponds to: Results 3.7 (Fig 8)
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(readr)
})

#==============================================================================
# 1. qPCR Analysis (Fig 8A)
#==============================================================================

# Load qPCR data
qpcr_data <- read_csv("data/experiments/SLC25A39_qPCR_results.csv")

# Calculate knockdown efficiency
qpcr_summary <- qpcr_data %>%
  group_by(group, replicate) %>%
  summarise(mean_ddCt = mean(ddCt), .groups = "drop")

# Normalize to control (set NC = 1)
nc_value <- qpcr_summary %>%
  filter(group == "siRNA-NC") %>%
  pull(mean_ddCt) %>%
  mean()

qpcr_summary$relative_expression <- 2^-(qpcr_summary$mean_ddCt - nc_value)

# Plot
p1 <- ggbarplot(qpcr_summary, x = "group", y = "relative_expression",
                add = c("mean_se", "jitter"),
                fill = "group") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Expression (fold change)",
       title = "SLC25A39 Knockdown Efficiency (qPCR)") +
  stat_compare_means(method = "anova")

#==============================================================================
# 2. CCK-8 Proliferation Assay (Fig 8B)
#==============================================================================

cck8_data <- read_csv("data/experiments/SLC25A39_CCK8_results.csv")

# Calculate proliferation curve
cck8_summary <- cck8_data %>%
  group_by(group, time_point) %>%
  summarise(mean_od = mean(OD450), sd_od = sd(OD450), .groups = "drop")

p2 <- ggplot(cck8_summary, aes(x = time_point, y = mean_od, group = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_od - sd_od, ymax = mean_od + sd_od), width = 0.2) +
  labs(x = "Time (hours)", y = "OD450",
       title = "Cell Proliferation (CCK-8)") +
  theme_bw()

#==============================================================================
# 3. Transwell Invasion Assay (Fig 8C)
#==============================================================================

transwell_data <- read_csv("data/experiments/SLC25A39_transwell_results.csv")

transwell_summary <- transwell_data %>%
  group_by(group, cell_line) %>%
  summarise(mean_count = mean(cell_count), sd_count = sd(cell_count), .groups = "drop")

p3 <- ggbarplot(transwell_summary, x = "group", y = "mean_count",
                facet.by = "cell_line",
                add = c("mean_se", "jitter"),
                fill = "group") +
  labs(y = "Invaded Cell Count",
       title = "Cell Invasion (Transwell)") +
  stat_compare_means(method = "t.test")

#==============================================================================
# 4. Wound Healing Assay (Fig 8D)
#==============================================================================

wound_data <- read_csv("data/experiments/SLC25A39_wound_healing_results.csv")

# Calculate migration area
wound_summary <- wound_data %>%
  mutate(migration_area = initial_area - final_area,
         migration_rate = migration_area / initial_area * 100) %>%
  group_by(group, cell_line, time_point) %>%
  summarise(mean_rate = mean(migration_rate), sd_rate = sd(migration_rate), .groups = "drop")

p4 <- ggplot(wound_summary, aes(x = time_point, y = mean_rate, group = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_rate - sd_rate, ymax = mean_rate + sd_rate), width = 2) +
  facet_wrap(~cell_line) +
  labs(x = "Time (hours)", y = "Migration Rate (%)",
       title = "Cell Migration (Wound Healing)") +
  theme_bw()

#==============================================================================
# 5. Western Blot - Pan-Kla (Fig 8E)
#==============================================================================

# Western blot quantification
wb_data <- read_csv("data/experiments/SLC25A39_WB_results.csv")

# Normalize to GAPDH and control
wb_summary <- wb_data %>%
  mutate(normalized_value = pan_kla / GAPDH) %>%
  group_by(group) %>%
  summarise(mean_norm = mean(normalized_value), sd_norm = sd(normalized_value), .groups = "drop")

# Normalize to NC
nc_norm <- wb_summary %>% filter(group == "siRNA-NC") %>% pull(mean_norm)
wb_summary$relative_level <- wb_summary$mean_norm / nc_norm

p5 <- ggbarplot(wb_summary, x = "group", y = "relative_level",
                add = c("mean_se", "jitter"),
                fill = "group") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Pan-Kla Level",
       title = "Global Lactylation Level (Western Blot)") +
  stat_compare_means(method = "anova")

#==============================================================================
# 6. Combined Summary
#==============================================================================

validation_summary <- list(
  qpcr = qpcr_summary,
  proliferation = cck8_summary,
  invasion = transwell_summary,
  migration = wound_summary,
  western_blot = wb_summary
)

saveRDS(validation_summary, "results/11_SLC25A39_validation.rds")

cat("=== Validation Summary ===\n")
cat(sprintf("Knockdown efficiency: %.1f%%\n", 
            (1 - qpcr_summary %>% filter(grepl("SLC25A39", group)) %>% 
               pull(relative_expression) %>% mean()) * 100))
cat(sprintf("Proliferation reduction: Check CCK-8 plot\n"))
cat(sprintf("Invasion reduction: Check transwell plot\n"))
cat(sprintf("Migration reduction: Check wound healing plot\n"))
cat(sprintf("Pan-Kla reduction: %.1f%%\n",
            (1 - wb_summary %>% filter(grepl("SLC25A39", group)) %>%
               pull(relative_level) %>% mean()) * 100))
