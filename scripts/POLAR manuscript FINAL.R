#===============================================================================
##                             POLAR TRIAL     
# Authors: Wungki Park, Catherine O'Connor, Walid Chatila, Marc Hilmi
#===============================================================================
# Libraries
library(readxl); library(ggplot2); library(MetBrewer); library(dplyr); library(maftools)
library(RColorBrewer); library(tidyverse); library(ggpubr); library(survival); library(survminer)
library(knitr); library(tibble); library(lubridate); library(ggsurvfit); library(gtsummary)
library(tidycmprsk); library(patchwork); library(readxl); library(patchwork)

# Run from repo root (recommended)
# If needed:
# setwd("PATH/TO/POLAR_REPO")  # user-defined

#===============================================================================
#                              Table 1 and 2
#===============================================================================

reset_gtsummary_theme()
theme_gtsummary_journal(journal = "nejm")
POLAR_clindata <- read.csv(file.path("data","POLAR_63_6_25.csv")) %>% 
  filter(Cohort %in% c("A", "B", "C")) 

A1 = POLAR_clindata %>% select('Cohort', 'Age',	'Sex', 'Race', 'ECOG', 'Initial_Stage', 'Histology', 'Previous_Surgery', 'Platinum_Type', 'CA19_9_T1', 'CEA_T1')
A2 = POLAR_clindata %>% select('Cohort', 'POLAR_BOR','TMB', 'IMPACT_HRD', 'Neoantigen', 'TIL_Density', 'NLR_T1')

#gtsummary tbl_summary for table.
tbl_summary(A1, by = Cohort) %>% add_n() %>% add_overall()
tbl_summary(A2, by = Cohort) %>% add_n() %>% add_overall()

#===============================================================================
#                              Fig 2a. WATERFALL
#===============================================================================

polar_waterfall <- read.csv(file.path("data","POLAR_RECIST.csv"))

# data is sorted by Best.change.target in descending order
ordered.waterfall <- polar_waterfall %>% arrange(desc(Best.change.target))
ordered.waterfall$Neoantigen.4Q <- as.factor(ifelse(is.na(ordered.waterfall$Neoantigen.4Q), "NA", ordered.waterfall$Neoantigen.4Q))
ordered.waterfall$IMPACT_HRD.4Q <- as.factor(ifelse(is.na(ordered.waterfall$IMPACT_HRD_4Q), "NA", ordered.waterfall$IMPACT_HRD_4Q))
# Cohort A separately 
ordered.waterfall.A <- ordered.waterfall %>% filter(Cohort == "A")
ordered.waterfall.B <- ordered.waterfall %>% filter(Cohort == "B")
ordered.waterfall.C <- ordered.waterfall %>% filter(Cohort == "C")
ordered.waterfall.BC <- ordered.waterfall %>% filter(Cohort != "A")

# Define a shared Y-axis range (e.g., from -100 to 100)
shared_y_scale <- scale_y_continuous(breaks = seq(-100, 100, by = 25), limits = c(-100, 100))

# Plot for cohort A (20 individuals)
plot_A <- ggplot(ordered.waterfall.A, aes(x = reorder(POLAR.ID, -Best.change.target), y = Best.change.target, fill = DDR_gene)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("gBRCA2m" = "#53008F", "gBRCA1m" = "#EE86BF", "gPALB2m" = "#E0A7F1", "sBRCA2m" = "#A922D3")) +
  labs(title = "Cohort A", x = "", y = "Best % Change in Target Lesions by RECIST V1.1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  shared_y_scale +
  geom_hline(yintercept = -30, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red2")

# Plot for cohort B (12 individuals)
plot_B <- ggplot(ordered.waterfall.B, aes(x = reorder(POLAR.ID, -Best.change.target), y = Best.change.target, fill = IMPACT_HRD.4Q)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("NA" = "lightgrey", "1" = "darkgrey", "2" = "#47C5FF", "3" = "#3C6BE2", "4"= "#02008F")) +
  labs(title = "Cohort B", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  shared_y_scale +
  geom_hline(yintercept = -30, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red2")

# Plot for cohort C (13 individuals)
plot_C <- ggplot(ordered.waterfall.C, aes(x = reorder(POLAR.ID, -Best.change.target), y = Best.change.target, fill = IMPACT_HRD.4Q)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("NA" = "lightgrey", "1" = "darkgrey", "2" = "#47C5FF", "3" = "#3C6BE2", "4"= "#02008F")) +
  labs(title = "Cohort C", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  shared_y_scale +
  geom_hline(yintercept = -30, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red2")

#Cohort B and C

# Arrange the plots with proportional widths using patchwork
(plot_A + plot_B + plot_C) + 
  plot_layout(widths = c(20, 12, 13))  # Proportional to number of individuals

#===============================================================================
#                       Fig.2b  POLAR: Kaplan Meier Curves PFS and OS
#===============================================================================
##survival analysis PFS
# Load required packages
library(survival); library(survminer)

POLAR <- read.csv(file.path("data","POLAR_63_6_25.csv"))
POLAR_A = POLAR %>% filter(Cohort == "A")
POLAR_2 <- read.csv(file.path("data","POLAR_63_6_25.csv"))
POLAR_BC = POLAR  %>% filter(Cohort != "A")

# PFS Fit survival curves using the Kaplan-Meier method
fit <- survfit(Surv(PFS_months, PFS_event) ~ Cohort, data = POLAR)

#===============================================================================
#                 Fig. 2b PFS Plot the Kaplan-Meier curves
#===============================================================================
p <- ggsurvplot(fit, data = POLAR, 
                censor = TRUE, censor.shape = "|", censor.size = 10,
                pval = FALSE,
                fun = "pct",
                xlab = "Months since the first cycle", 
                ylab = "Probability of Progression-free Survival (%)",
                font.x = 12, font.y = 16,
                break.x.by = 3, xlim = c(0, 52),
                size = 2,  
                legend = "none",
                risk.table = TRUE,
                risk.table.height = 0.15,
                legend.labs = c('A', 'B', 'C'),
                tables.theme = theme_cleantable(),
                palette = c("A" = "#A922D3", "B" = "#22E3CD", "C" = "#02008F"))

# Add post-hoc padj values
p$plot <- p$plot +
  annotate("text", x = 30, y = 85, 
           label = expression("A vs B:" ~ p[adj] ~ "=" ~ 0.147), size = 6, hjust = 0) +
  annotate("text", x = 30, y = 75, 
           label = expression("A vs C:" ~ p[adj] ~ "≤ 0.001"), size = 6, hjust = 0) +
  theme(axis.title.y = element_text(face = "bold"))

# Print the updated plot
print(p)
#===============================================================================
#                 Fig. 2c,d OS Plot the Kaplan-Meier curves
#===============================================================================

# OS Fit survival curves using the Kaplan-Meier method
fit_OS_A <- survfit(Surv(OS, OS_Event) ~ Cohort, data = POLAR_A)
fit_OS_BC <- survfit(Surv(OS, OS_Event) ~ Cohort, data = POLAR_BC)

# Extract survival probabilities
surv_summary <- summary(fit_OS_A, times = c(24, 36))
os_2yr <- round(surv_summary$surv[1] * 100, 1)
os_3yr <- round(surv_summary$surv[2] * 100, 1)

# Generate the survival plot
km_plot <- ggsurvplot(
  fit_OS_A, data = POLAR_A, 
  censor = TRUE, censor.shape = "|", censor.size = 10,
  pval = FALSE,
  fun = "pct",
  xlab = "Months since the first cycle", 
  ylab = "Probability of Overall Survival (%)",
  font.x = 12, font.y = 16,
  break.x.by = 3, xlim = c(0, 52),
  size = 2,
  legend = "none",
  risk.table = TRUE,
  risk.table.height = 0.15,
  legend.labs = c('A'),
  tables.theme = theme_cleantable(),
  palette = c("A" = "#A922D3")
)

# Text to annotate
annotation_text <- paste0(
  "Median follow-up: 26.0 months (range 1.4–52.5)\n",
  "2Y-OS rate: 56% (95% CI: 41–76%)\n",
  "3Y-OS rate: 44% (95% CI: 28–69%)"
)

# Add annotation, vertical lines, and bold y-axis label
km_plot$plot <- km_plot$plot +
  annotate("text", x = 0, y = 20, 
           label = annotation_text, 
           hjust = 0, vjust = 1,
           size = 6, color = "black") +
  geom_vline(xintercept = 24, linetype = "dashed", color = "black", size = 0.4) +
  geom_vline(xintercept = 36, linetype = "dashed", color = "black", size = 0.4) +
  theme(axis.title.y = element_text(face = "bold"))

# Print the plot
print(km_plot)

# OS for POLAR BC Plot the Kaplan-Meier curves
p_BC <- ggsurvplot(
  fit_OS_BC, data = POLAR_BC, 
  censor = TRUE, censor.shape = "|", censor.size = 12,
  pval = TRUE, pval.coord = c(0, 10), pval.size = 8,
  fun = "pct",
  xlab = "Months since the first cycle", 
  ylab = "Probability of Overall Survival (%)",
  font.x = 12, font.y = 16,
  break.x.by = 3, xlim = c(0, 43),
  size = 2,
  legend = "none",
  risk.table = TRUE,
  risk.table.height = 0.15,
  legend.labs = c('B', 'C'),
  tables.theme = theme_cleantable(),
  palette = c("B" = "#22E3CD", "C" = "#02008F")
)

# Add bold styling to y-axis label
p_BC$plot <- p_BC$plot + theme(axis.title.y = element_text(face = "bold"))

# Print the styled plot
print(p_BC)
#===============================================================================
library(survival); library(survminer)

# Fit survival curves for all three cohorts
fit_OS_all <- survfit(Surv(OS, OS_Event) ~ Cohort, data = POLAR)

# Extract survival for Cohort A
fit_OS_A <- survfit(Surv(OS, OS_Event) ~ 1, data = subset(POLAR, Cohort == "A"))
surv_summary <- summary(fit_OS_A, times = c(24, 36))
os_2yr <- round(surv_summary$surv[1] * 100, 1)
os_3yr <- round(surv_summary$surv[2] * 100, 1)

# Plot KM curves for all cohorts
km_plot_all <- ggsurvplot(fit_OS_all, data = POLAR, 
                          censor = TRUE, censor.shape = "|", censor.size = 9,
                          pval = FALSE, 
                          fun = "pct",
                          xlab = "Months since the first cycle", 
                          ylab = "Probability of Overall Survival (%)",
                          font.x = 12, font.y = 16,
                          break.x.by = 3, xlim = c(0, 52),
                          size = 2,
                          legend = "none",
                          risk.table = TRUE,
                          risk.table.height = 0.15,
                          legend.labs = c("A", "B", "C"),
                          tables.theme = theme_cleantable(),
                          palette = c("A" = "#A922D3", "B" = "#22E3CD", "C" = "#02008F")
)

# Add manual annotations
km_plot_all$plot <- km_plot_all$plot +
  annotate("text", x = 20, y = 100, label = "Median follow-up: 26.0 months (range 1.4–52.5)", hjust = 0, size = 5.5) +
  annotate("text", x = 20, y = 94, label = paste0("2Y-OS rate (Cohort A): ", os_2yr, "% (95% CI: 41–76%)"), hjust = 0, size = 5.5) +
  annotate("text", x = 20, y = 88, label = paste0("3Y-OS rate (Cohort A): ", os_3yr, "% (95% CI: 28–69%)"), hjust = 0, size = 5.5) +
  annotate("text", x = 38, y = 76, label = expression("OS log-rank" ~ p == 0.005), hjust = 0, size = 5.5) +
  annotate("text", x = 40, y = 70, label = expression("A vs B:" ~ p[adj] == 0.709), hjust = 0, size = 5.5) +
  annotate("text", x = 40, y = 64, label = expression("A vs C:" ~ p[adj] == 0.004), hjust = 0, size = 5.5) +
  annotate("segment", x = 24, xend = 24, y = -5, yend = 56, 
           linetype = "dashed", color = "black", size = 0.4) +
  annotate("segment", x = 36, xend = 36, y = -5, yend = 44, 
           linetype = "dashed", color = "black", size = 0.4) +
  theme(axis.title.y = element_text(face = "bold"))

# Print the combined plot
print(km_plot_all)

#===============================================================================
#                 Fig 3a. Swimmer's plot by Ramzi Homsi, Walid Chatila
#===============================================================================
# Load data ----------------------------------------------------------------
# Edit path as needed to point to your local data file
swimmer_data <- read.csv(file.path("data","swimmer_data_POLAR.csv"))

# Data preparation ---------------------------------------------------------

# Create median split groups for CA19-9 and NLR at T1 (baseline)
swimmer_data$CA199_T1_median <- ntile(swimmer_data$CA19_9_T1, n = 2)
swimmer_data$NLR_T1_median <- ntile(swimmer_data$NLR_T1, n = 2)

# Sort patients by cohort (descending) and PFS duration (ascending)
# This ordering determines the y-axis arrangement in the swimmer plot
swimmer_data <- swimmer_data %>%
  arrange(desc(Cohort), PFS_months) %>%
  mutate(POLAR_ID = factor(POLAR_ID, levels = POLAR_ID))

# Create OS variable for swimmer plot display
# For patients with ongoing PFS (PFS_event = 0), censor OS to match PFS
swimmer_data$OS_swimmer <- swimmer_data$OS
swimmer_data[swimmer_data$PFS_event == 0, "OS_swimmer"] <-
  swimmer_data[swimmer_data$PFS_event == 0, "PFS_months"]

# Generate individual annotation plots -------------------------------------

# Gene panel plot (left strip #7)
# Color palette for different gene categories
genes_plot <- ggplot(
  swimmer_data,
  aes(x = POLAR_ID, y = 0.5, fill = as.factor(gene))
) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "#717171", size = 0.1, width = 0.95) +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "#c2c1c0", "#6cbdea", "#b5ddbe", "#c13d8f",
      "#815ca6", "#94cb7e", "#eedfa2", "#f9c01f", "#dd797d"
    ),
    na.value = "grey30"
  ) +
  ggtitle("Genes")

# CA19-9 at baseline (T1) - median split (left strip #3)
# Low = light orange, High = dark orange
CA199_T1_plot <- ggplot(
  swimmer_data,
  aes(x = POLAR_ID, y = 0.5, fill = as.factor(CA199_T1_median))
) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "#717171", size = 0.1, width = 0.95) +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("#ffd48f", "#da7127"),
    na.value = "#c2c1c0"
  ) +
  ggtitle("CA 19-9")

# Neutrophil-to-Lymphocyte Ratio at baseline (T1) - median split (left strip #2)
# Low = light orange, High = dark orange
NLR_T1_plot <- ggplot(
  swimmer_data,
  aes(x = POLAR_ID, y = 0.5, fill = as.factor(NLR_T1_median))
) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "#717171", size = 0.1, width = 0.95) +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("#ffd48f", "#da7127"),
    na.value = "#c2c1c0"
  ) +
  ggtitle("NLR")

# Immune-related adverse events (irAE) plot (left strip #1)
# White = no irAE, Black = irAE present
irAE_plot <- ggplot(
  swimmer_data,
  aes(x = POLAR_ID, y = 0.5, fill = as.factor(irAE))
) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "#717171", size = 0.1, width = 0.95) +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("white", "black"),
    na.value = "#c2c1c0"
  ) +
  ggtitle("irAE")

# KRAS mutation status (left strip #5)
# Gradient of greens representing different KRAS mutation categories
KRAS_plot <- ggplot(
  swimmer_data,
  aes(x = POLAR_ID, y = 0.5, fill = as.factor(KRAS_Mutation_Final_Grouping))
) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "#717171", size = 0.1, width = 0.95) +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "#0d4420", "#086d38", "#40ab5d",
      "#e6f2e0", "#a2d39a", "white"
    ),
    na.value = "#c2c1c0"
  ) +
  ggtitle("KRAS")

# CDKN2A homozygous deletion status (left strip #4)
# White = wild-type, Dark blue = homozygous deletion
CDKN2A_plot <- ggplot(
  swimmer_data,
  aes(x = POLAR_ID, y = 0.5, fill = as.factor(CDKN2A_HomeDel_Final))
) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "#717171", size = 0.1, width = 0.95) +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("white", "#0c2c84"),
    na.value = "#c2c1c0"
  ) +
  ggtitle("CDKN2A")

# PFS event indicator (right strip #9)
# Black = event occurred, White = censored
p_event_plot <- ggplot(
  swimmer_data,
  aes(x = POLAR_ID, y = 0.5, fill = as.factor(PFS_event))
) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "#717171", size = 0.1, width = 0.95) +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("black", "white"),
    na.value = "#c2c1c0"
  ) +
  ggtitle("Still on\ntreatment")

# Zygosity status (left strip #6)
# Color palette for biallelic, monoallelic, indeterminate, and NA
zyg_plot <- ggplot(
  swimmer_data,
  aes(
    x = POLAR_ID,
    y = 0.5,
    fill = factor(Zygosity, levels = c("Biallelic", "Monoallelic", "Indeterminate", "NA"))
  )
) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "#717171", size = 0.1, width = 0.95) +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Biallelic"     = "#2e1a58",  # Dark navy
      "Monoallelic"   = "#63499e",  # Lighter navy
      "Indeterminate" = "#af9ec5",  # Light purple
      "NA"            = "#c2c1c0"   # Light grey
    ),
    drop = FALSE,
    na.value = "#cfcfcf"
  ) +
  ggtitle("Zygosity")

# Main swimmer plot (central panel #8) ------------------------------------
# Shows PFS (solid bars) and OS (semi-transparent overlay)
# Different colors represent different treatment cohorts

swim_combo_plot <- ggplot(swimmer_data, aes(x = POLAR_ID)) +
  # PFS base layer (solid color)
  geom_col(
    aes(y = PFS_months, fill = Cohort),
    show.legend = FALSE,
    color = "#717171",
    size = 0.1,
    width = 0.95
  ) +
  # OS overlay layer (lighter alpha to show both PFS and OS)
  geom_col(
    aes(y = OS_swimmer, fill = Cohort),
    alpha = 0.35,
    show.legend = FALSE,
    color = "#717171",
    size = 0.1,
    width = 0.95
  ) +
  # Add grey dashed line at 6 months (use hline after coord_flip)
  geom_hline(yintercept = 6, linetype = "dashed", color = "grey50", size = 0.5) +
  # Cohort color scheme
  scale_fill_manual(values = c("#8c5ca6", "#6ac7ba", "#556587")) +
  coord_flip() +
  theme_pubr() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(hjust = 0.5, size = 11)
  ) +
  # X-axis shows time in months (0-60 months, breaks every 6 months)
  # Remove expansion to eliminate gap
  scale_y_continuous(
    limits = c(0, 60),
    breaks = seq(0, 60, 6),
    expand = c(0, 0)
  ) +
  labs(y = "Progression Free Survival and Overall Survival")




##################################
# Create individual legend components using dummy plots

# 1. Cohort PFS legend
legend_cohort_pfs <- ggplot(
  data.frame(
    x = rep(1, 3),
    y = rep(1, 3),
    cohort = factor(
      c("Cohort A", "Cohort B", "Cohort C"),
      levels = c("Cohort A", "Cohort B", "Cohort C")
    )
  ),
  aes(x = x, y = y, fill = cohort)
) +
  geom_tile() +
  scale_fill_manual(
    name = "Cohort\nPFS",
    values = c(
      "Cohort A" = "#8c5ca6",
      "Cohort B" = "#6ac7ba",
      "Cohort C" = "#556587"
    ),
    labels = c("Cohort A: HRD", "Cohort B: ncHRD", "Cohort C: HRP"),
    drop = FALSE
  ) +
  guides(fill = guide_legend(title = "Cohort\nPFS", title.position = "top", title.hjust = 0)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

# 2. Cohort OS legend (lighter version)
legend_cohort_os <- ggplot(
  data.frame(
    x = rep(1, 3),
    y = rep(1, 3),
    cohort = factor(
      c("Cohort A", "Cohort B", "Cohort C"),
      levels = c("Cohort A", "Cohort B", "Cohort C")
    )
  ),
  aes(x = x, y = y, fill = cohort)
) +
  geom_tile(alpha = 0.35) +
  scale_fill_manual(
    name = "OS",
    values = c(
      "Cohort A" = "#8c5ca6",
      "Cohort B" = "#6ac7ba",
      "Cohort C" = "#556587"
    ),
    labels = c("Cohort A: HRD", "Cohort B: ncHRD", "Cohort C: HRP"),
    drop = FALSE
  ) +
  guides(fill = guide_legend(
    title = "OS",
    title.position = "top",
    title.hjust = 0,
    override.aes = list(alpha = 0.35)
  )) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )
# 3. DDR Gene legend - Create manual data to ensure correct order
legend_genes <- ggplot(
  data.frame(
    x = rep(1, 8),
    y = rep(1, 8),
    gene = factor(
      c("BRCA1", "BRCA2", "PALB2", "ATM", "BLM", "CHEK2", "FANCC", "MUTYH"),
      levels = c("BRCA1", "BRCA2", "PALB2", "ATM", "BLM", "CHEK2", "FANCC", "MUTYH")
    )
  ),
  aes(x = x, y = y, fill = gene)
) +
  geom_tile() +
  scale_fill_manual(
    name = "DDR Gene",
    values = c(
      "BRCA1" = "#c13d8f",
      "BRCA2" = "#815ca6",
      "PALB2" = "#dd797d",
      "ATM" = "#6cbdea",
      "BLM" = "#b5ddbe",
      "CHEK2" = "#94cb7e",
      "FANCC" = "#eedfa2",
      "MUTYH" = "#f9c01f"
    ),
    drop = FALSE
  ) +
  guides(fill = guide_legend(title = "DDR Gene", title.position = "top", title.hjust = 0)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(face = "italic", size = 9)
  )

# 4. Zygosity legend
legend_zygosity <- ggplot(swimmer_data, aes(x = 1, y = 1, fill = factor(Zygosity))) +
  geom_tile() +
  scale_fill_manual(
    name = "Zygosity",
    values = c(
      "Biallelic" = "#2e1a58",
      "Monoallelic" = "#63499e",
      "Indeterminate" = "#af9ec5"
    ),
    labels = c("Biallelic", "Monoallelic", "Indeterminate"),
    na.translate = FALSE
  ) +
  guides(fill = guide_legend(title = "Zygosity", title.position = "top", title.hjust = 0)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

# 5. CA19-9/NLR legend
legend_ca199_nlr <- ggplot(swimmer_data, aes(x = 1, y = 1, fill = as.factor(CA199_T1_median))) +
  geom_tile() +
  scale_fill_manual(
    name = "CA19-9/NLR",
    values = c("2" = "#da7127", "1" = "#ffd48f"),
    labels = c("High", "Low"),
    na.translate = FALSE
  ) +
  guides(fill = guide_legend(title = "CA19-9/NLR", title.position = "top", title.hjust = 0)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

# 6. irAE legend
legend_irae <- ggplot(swimmer_data, aes(x = 1, y = 1, fill = as.factor(irAE))) +
  geom_tile() +
  scale_fill_manual(
    name = "irAE",
    values = c("1" = "black"),
    labels = c("Yes"),
    na.translate = FALSE
  ) +
  guides(fill = guide_legend(title = "irAE", title.position = "top", title.hjust = 0)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

# 7. KRAS legend - Fixed Q61 and NOS color assignments
legend_kras <- ggplot(
  data.frame(
    x = rep(1, 6),
    y = rep(1, 6),
    kras = factor(
      c("G12D", "G12V", "G12R", "Q61", "NOS", "WT"),
      levels = c("G12D", "G12V", "G12R", "Q61", "NOS", "WT")
    )
  ),
  aes(x = x, y = y, fill = kras)
) +
  geom_tile() +
  scale_fill_manual(
    name = "KRAS",
    values = c(
      "G12D" = "#0d4420",   # Darkest green
      "G12V" = "#086d38",   # Dark green
      "G12R" = "#40ab5d",   # Medium green
      "Q61" = "#a2d39a",    # Light green (darker than NOS)
      "NOS" = "#e6f2e0",    # Very light green (lightest)
      "WT" = "white"
    ),
    drop = FALSE
  ) +
  guides(fill = guide_legend(
    title = "KRAS",
    title.position = "top",
    title.hjust = 0,
    override.aes = list(color = "black", size = 0.5)
  )) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold.italic", size = 10),
    legend.text = element_text(size = 9)
  )

# 8. CDKN2A legend
legend_cdkn2a <- ggplot(data.frame(
  x = rep(1, 2),
  y = rep(1, 2),
  cdkn2a = factor(c("0", "1"), levels = c("1", "0"))
), aes(x = x, y = y, fill = cdkn2a)) +
  geom_tile() +
  scale_fill_manual(
    name = "CDKN2A",
    values = c("1" = "#0c2c84", "0" = "white"),
    labels = c("Homozygous Deletion", "WT"),
    drop = FALSE
  ) +
  guides(fill = guide_legend(
    title = "CDKN2A",
    title.position = "top",
    title.hjust = 0,
    override.aes = list(color = "black", size = 0.5)
  )) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold.italic", size = 10),
    legend.text = element_text(size = 9)
  )


# 9. Not Available legend (create manually)
legend_na <- ggplot(data.frame(x = 1, y = 1, fill = "NA"), aes(x, y, fill = fill)) +
  geom_tile() +
  scale_fill_manual(
    name = "",
    values = c("NA" = "#c2c1c0"),
    labels = c("Not Available")
  ) +
  guides(fill = guide_legend(title = NULL)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.text = element_text(size = 9)
  )

# 10. Still on treatment legend (using point shape for arrow)
legend_treatment <- ggplot(data.frame(x = 1, y = 1, shape = "treatment"), aes(x, y, shape = shape)) +
  geom_point(size = 3) +
  scale_shape_manual(
    name = "",
    values = c("treatment" = 16),  # Use 16 for filled circle, or could use arrow
    labels = c("Still on treatment")
  ) +
  guides(shape = guide_legend(title = NULL)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.text = element_text(size = 9)
  )

# Extract all legends
leg_cohort_pfs <- cowplot::get_legend(legend_cohort_pfs)
leg_cohort_os <- cowplot::get_legend(legend_cohort_os)
leg_genes <- cowplot::get_legend(legend_genes)
leg_zygosity <- cowplot::get_legend(legend_zygosity)
leg_ca199_nlr <- cowplot::get_legend(legend_ca199_nlr)
leg_irae <- cowplot::get_legend(legend_irae)
leg_kras <- cowplot::get_legend(legend_kras)
leg_cdkn2a <- cowplot::get_legend(legend_cdkn2a)
leg_na <- cowplot::get_legend(legend_na)
leg_treatment <- cowplot::get_legend(legend_treatment)

# Combine all legends vertically
combined_legend <- cowplot::plot_grid(
  leg_cohort_pfs,
  leg_cohort_os,
  leg_genes,
  leg_zygosity,
  leg_ca199_nlr,
  leg_irae,
  leg_kras,
  leg_cdkn2a,
  leg_na,
  leg_treatment,
  ncol = 1,
  align = "v",
  rel_heights = c(1.2, 1.2, 2.5, 1.2, 0.8, 0.3, 2.2, 1.0, 0.5, 0.5)
)

# Combine main plot with legend
final_plot_with_legend <- cowplot::plot_grid(
  cowplot::plot_grid(
    irAE_plot,
    NLR_T1_plot,
    CA199_T1_plot,
    CDKN2A_plot,
    KRAS_plot,
    zyg_plot,
    genes_plot,
    swim_combo_plot,
    p_event_plot,
    ncol = 9,
    nrow = 1,
    rel_widths = c(0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.64, 0.08),
    align = "hv",
    axis = "l"
  ),
  combined_legend,
  ncol = 2,
  rel_widths = c(0.80, 0.20)
)

# Display the final plot
print(final_plot_with_legend)

#===============================================================================
#                      Fig 3b. POLAR WES Oncoprint
#===============================================================================

POLAR_Baseline_WES <- read.csv(
  file.path("data", "POLAR_MS_Baseline_WES_ClinData.csv"),
  header = TRUE,
  check.names = FALSE
)

POLAR_Baseline <- read.maf(maf = file.path("data","POLAR_Baseline_WES.maf"),
                           clinicalData = POLAR_Baseline_WES)

#Shows sample summry.
getSampleSummary(POLAR_Baseline)

#Shows gene summary.
getGeneSummary(POLAR_Baseline)

#shows clinical data associated with samples
getClinicalData(POLAR_Baseline)

#Shows all fields in MAF
getFields(POLAR_Baseline)

#Writes maf summary to an output file with basename iPC.
write.mafSummary(maf = POLAR_Baseline, basename = 'POLAR_Baseline')

plotmafSummary(maf = POLAR_Baseline, rmOutlier = TRUE, addStat = 'median', 
               dashboard = TRUE, titvRaw = FALSE)
#oncoplot for top ten mutated genes.
oncoplot(maf = POLAR_Baseline, draw_titv = TRUE)

#Shows sample summry.
getSampleSummary(POLAR_Baseline)

#Shows gene summary.
getGeneSummary(POLAR_Baseline)

#shows clinical data associated with samples
getClinicalData(POLAR_Baseline)

#Shows all fields in MAF
getFields(POLAR_Baseline)

#Writes maf summary to an output file with basename iPC.
write.mafSummary(maf = POLAR_Baseline, basename = 'POLAR_Baseline')

plotmafSummary(maf = POLAR_Baseline, rmOutlier = TRUE, addStat = 'median', 
               dashboard = TRUE, titvRaw = FALSE)
#oncoplot for top ten mutated genes.
oncoplot(maf = POLAR_Baseline, draw_titv = TRUE)

#Color map for annotations
iPC_cust_cols = c("Missense_Mutation" = "#9d9cd5", "Nonsense_Mutation" = "#b3e1bf", "Multi_Hit" = "#299d3e", "Frame_Shift_Ins" = "#ff2a67", "Frame_Shift_Del" = "#011a52",
                  "In_Frame_Ins" = "#A6CEE3", "In_Frame_Del" = "deepskyblue", "Splice_Site" = "#355828")
cohort_palette = c("A" = "#A922D3", "B" = "#22E3CD", "C" = "#02008F")

kRAS = c("G12D" = "#17154f", "G12V" = "#98d1d9", "G12R" = "#9d9cd5", "Q61H" = "#a5506d", "G12C" = "#edf181", "Q61R" = "#ada43b", "G13D" = "#ff2a67", "I171Nfs*14" = "#2f357c")

dDRm = c("gATMm" = "#299d3e", "gATMm_gMUTYHm" =  "#ada43b",  "gBLMm" =  "#e69b00", "gBRCA1m" = "#011a52", "gBRCA2m" =  "#ff2a67", "gBRCA2m_gFLCNm"= "#ff3100", "gCHEK2m" = "#9d9cd5",
         "gFANCCm" = "#355828", "gPALB2m" = "#a5506d", "gMUTYHm" = "#f6b3b0", "sBRCA2m" = "#ff2a67")

wGD = c("TRUE" = "#ff2a67", "FALSE" = "#011a52")
quartile = c("1" = "#2f357c", "2" = "#A6CEE3",  "3" = "#a5506d", "4"= "#ff2a67")
titv = c("C>T" = "#9d9cd5", "C>G" = "#ff2a67", "C>A" = "#011a52", "T>A" =  "#318c97", "T>C" = "#98d1d9", "T>G" = "#fba702")
titvcolors = list(titv_cols = titv)
iPC_genes = c("KRAS", "TP53", "CDKN2A", "SMAD4", "BRCA2", "BRCA1", "PALB2", "ATM")

#Definte annocolors before making maf oncoplot
annocolors = list(Cohort = cohort_palette, HRD_gene = dDRm, KRAS = kRAS, WGD = wGD, Time_point= time_point, Pre_Post = pre_Post, IMPACT_HRD_Quartile = quartile, Neoantigen_Quartile = quartile)

#Variant allele frequencies (Left bar plot)
POLAR_Baseline_genes_vaf = subsetMaf(maf = POLAR_Baseline, genes = iPC_genes, fields = "t_var_freq", mafObj = FALSE)[,mean(t_var_freq, na.rm = TRUE), Hugo_Symbol]

dir.create("outputs", showWarnings = FALSE)
dir.create(file.path("outputs","figures"), recursive = TRUE, showWarnings = FALSE)

outfile <- file.path("outputs","figures","Figure_3b.pdf")  # safer filename

pdf(outfile,
    width = 14,
    height = 10,
    useDingbats = FALSE,
    family = "Helvetica",
    pointsize = 10)

oncoplot(
  maf = POLAR_Baseline,
  removeNonMutated = FALSE,
  genes = iPC_genes,
  keepGeneOrder = TRUE,
  logColBar = TRUE,
  sepwd_samples = 1,
  fontSize = 0.9,
  SampleNamefontSize = 1.6,
  titleFontSize = 1.6,
  legendFontSize = 1.4,
  drawColBar = TRUE,
  annotationFontSize = 1.4,
  anno_height = 3.5,
  colors = iPC_cust_cols,
  leftBarData = POLAR_Baseline_genes_vaf,
  leftBarLims = c(0, 0.5),
  clinicalFeatures = c("Cohort", "KRAS", "HRD_gene", "WGD", "IMPACT_HRD_Quartile", "Neoantigen_Quartile"),
  sortByAnnotation = TRUE,
  annotationOrder = c("A", "B", "C"),
  annotationColor = annocolors,
  annoBorderCol = "white",
  titv_col = titv,
  additionalFeatureCol = "gray70",
  bgCol = "#CCCCCC",
  borderCol = "white",
  draw_titv = TRUE
)

dev.off()

#===============================================================================
# Fig 3C. Spaghetti plot: Tumor response change with diff TIL
#===============================================================================

create_spider_plot <- function(tumorresponse, cohort_label,
                               y_limits = c(-100, 110),
                               y_breaks = 10,
                               months_max = 36,
                               month_breaks = 2,
                               weeks_per_month = 4.345) {
  
  cohort_data <- tumorresponse %>% filter(Cohort == cohort_label)
  
  # Time columns can be "0","9",... OR "X0","X9",...
  time_cols <- names(cohort_data)[str_detect(names(cohort_data), "^X?\\d+$")]
  
  if (length(time_cols) == 0) {
    stop("No timepoint columns detected. Expected columns like 0/9/18... or X0/X9/X18...")
  }
  
  data_long <- cohort_data %>%
    pivot_longer(
      cols = all_of(time_cols),
      names_to = "Time_raw",
      values_to = "Change"
    ) %>%
    mutate(
      Time_weeks  = parse_number(Time_raw),          # works for "9" or "X9"
      Time_months = Time_weeks / weeks_per_month,    # convert to months
      Change      = as.numeric(Change)
    ) %>%
    filter(!is.na(Time_months), !is.na(Change))
  
  ggplot(data_long, aes(x = Time_months, y = Change, group = POLAR_ID)) +
    geom_line(aes(color = TIL), linewidth = 0.7, alpha = 0.9) +
    geom_point(aes(color = TIL), size = 1.8, alpha = 0.9) +
    scale_y_continuous(
      limits = y_limits,
      breaks = seq(y_limits[1], y_limits[2], by = y_breaks)
    ) +
    scale_x_continuous(
      limits = c(0, months_max),
      breaks = seq(0, months_max, by = month_breaks)
    ) +
    scale_color_gradient(low = "blue", high = "red", name = "TIL") +
    scale_fill_gradient(low = "blue", high = "red", name = "TIL") +  # kept to match your prior block
    geom_hline(yintercept = 20,  linetype = "dashed", color = "red2") +
    geom_hline(yintercept = -30, linetype = "dashed", color = "green4") +
    labs(
      title = paste0("Cohort ", cohort_label, " % Change of Target Lesions"),
      x = "Months",
      y = "% Change"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
      axis.title.y = element_text(vjust = 2),
      axis.text.y = element_text(vjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# usage
tumorresponse <- read.csv(file.path("data","POLAR_RECIST_CohortA_Fig3C.csv"))
pA <- create_spider_plot(tumorresponse, "A")
print(pA)

# ---- Save as vector PDF (Nature-safe)
outfile <- file.path("outputs","figures","Figure_3c.pdf")

pdf(outfile,
    width = 7.5, height = 5.5,          # adjust if you need exact panel sizing
    useDingbats = FALSE,
    family = "Helvetica",
    pointsize = 10)

print(pA)
dev.off()

#===============================================================================
#          Figure 3D - Immunofluorescent  Imaging by Marc Hilmi
#===============================================================================
df0 <- read_excel(file.path("data","POLAR_IF_final_NMV2.xlsx")) %>%
  filter(!is.na(CD3_CD8))

# Helpers
add_bracket <- function(p, x1, x2, y, y_text, label, tick,
                        lw = 0.8, text_size = 4.2) {
  p +
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = lw) +
    annotate("segment", x = x1, xend = x1, y = y - tick, yend = y, linewidth = lw) +
    annotate("segment", x = x2, xend = x2, y = y - tick, yend = y, linewidth = lw) +
    annotate("text", x = (x1 + x2) / 2, y = y_text, label = label, size = text_size)
}

# returns a y upper limit that leaves room for annotations
ylim_with_headroom <- function(y, top_frac = 0.28) {
  y_max <- max(y, na.rm = TRUE)
  c(0, y_max * (1 + top_frac))
}

base_theme <- theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(size = 13),
    axis.title.y = element_text(size = 14),
    axis.line    = element_line(linewidth = 0.8),
    axis.ticks   = element_line(linewidth = 0.8),
    plot.margin  = margin(12, 14, 8, 12)  # extra top margin helps
  )

# Plot P (PFS group)
df_pfs <- df0 %>%
  mutate(
    PFS_group = ifelse(PFS <= 4, "PFS ≤ 4 months", "PFS > 4 months"),
    PFS_group = factor(PFS_group, levels = c("PFS ≤ 4 months", "PFS > 4 months"))
  ) %>%
  filter(!is.na(PFS_group))

wilcox_pfs <- wilcox.test(CD3_CD8 ~ PFS_group, data = df_pfs, exact = FALSE)
p_label <- paste0("Wilcoxon p = ", formatC(wilcox_pfs$p.value, format = "e", digits = 2))

yl_p <- ylim_with_headroom(df_pfs$CD3_CD8, top_frac = 0.22)
y_text_p <- yl_p[2] * 0.96  # place text within reserved headroom

group_colors <- c("PFS ≤ 4 months" = "#4E79A7", "PFS > 4 months" = "#E15759")

p <- ggplot(df_pfs, aes(x = PFS_group, y = CD3_CD8, color = PFS_group)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, fill = "white", linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 2.3, alpha = 0.85) +
  annotate("text", x = 1.5, y = y_text_p, label = p_label, size = 4.4) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "CD3+CD8+ infiltrate",
    x = NULL,
    y = "Stromal immune infiltrate (% positive cells)"
  ) +
  base_theme +
  coord_cartesian(ylim = yl_p, clip = "off")

print(p)

# Plot Q (Cohorts A/B/C with brackets vs C)
df_cohort <- df0 %>%
  mutate(Cohort = factor(Cohort, levels = c("A","B","C"))) %>%
  filter(!is.na(Cohort))

p_CA <- wilcox.test(CD3_CD8 ~ Cohort, data = df_cohort %>% filter(Cohort %in% c("A","C")), exact = FALSE)$p.value
p_CB <- wilcox.test(CD3_CD8 ~ Cohort, data = df_cohort %>% filter(Cohort %in% c("B","C")), exact = FALSE)$p.value

label_CA <- paste0("p = ", formatC(p_CA, format = "e", digits = 2))
label_CB <- paste0("p = ", formatC(p_CB, format = "e", digits = 2))

yl_q <- ylim_with_headroom(df_cohort$CD3_CD8, top_frac = 0.35)
y_max_q <- max(df_cohort$CD3_CD8, na.rm = TRUE)

# bracket geometry anchored to y_max (stable) but stays within yl_q
tick <- y_max_q * 0.03
y1 <- y_max_q * 1.06
y1t <- y_max_q * 1.11
y2 <- y_max_q * 1.20
y2t <- y_max_q * 1.25

cohort_colors <- c("A"="#A922D3","B"="#22E3CD","C"="#02008F")

q <- ggplot(df_cohort, aes(Cohort, CD3_CD8, color = Cohort)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 2.3, alpha = 0.85) +
  scale_color_manual(values = cohort_colors) +
  labs(
    title = "CD3+CD8+ infiltrate",
    x = NULL,
    y = "Stromal immune infiltrate (% positive cells)"
  ) +
  base_theme +
  coord_cartesian(ylim = yl_q, clip = "off")

q <- add_bracket(q, x1 = 1, x2 = 3, y = y1, y_text = y1t, label = label_CA, tick = tick)
q <- add_bracket(q, x1 = 2, x2 = 3, y = y2, y_text = y2t, label = label_CB, tick = tick)

print(q)

# Export both as vector PDFs
dir.create("outputs", showWarnings = FALSE)
dir.create(file.path("outputs","figures"), recursive = TRUE, showWarnings = FALSE)

pdf(file.path("outputs","figures","Figure_3d_p.pdf"),
    width = 6.5, height = 5.5, useDingbats = FALSE, family = "Helvetica", pointsize = 10)
print(p)
dev.off()

pdf(file.path("outputs","figures","Figure_3d_q.pdf"),
    width = 5.5, height = 5.5, useDingbats = FALSE, family = "Helvetica", pointsize = 10)
print(q)
dev.off()

#===============================================================================
#                   Supplementary Fig 4. cfDNA Figures
#===============================================================================
library(tidyverse); library(ggplot2); library(dplyr)

# Load and filter data
cfDNA_df <- read.csv(file.path("data","POLAR_MSK_ACCESS.csv"), stringsAsFactors = FALSE)

# Prepare long format for T1 and T2 values
cfDNA_long <- cfDNA_df %>%
  filter(!is.na(cfDNA_T1) & !is.na(cfDNA_T2)) %>%
  select(cfDNA_T1, cfDNA_T2) %>%
  pivot_longer(cols = everything(), names_to = "Timepoint", values_to = "cfDNA")

# Clean up labels
cfDNA_long <- cfDNA_long %>%
  mutate(Timepoint = recode(Timepoint, cfDNA_T1 = "T1", cfDNA_T2 = "T2"))

# Plot
ggplot(cfDNA_long, aes(y = cfDNA, fill = Timepoint, color = Timepoint)) +
  geom_density(alpha = 0.5, size = 1.2, orientation = "y") +
  scale_fill_manual(values = c("T1" = "#1f78b4", "T2" = "#33a02c")) +
  scale_color_manual(values = c("T1" = "#1f78b4", "T2" = "#33a02c")) +
  labs(
    title = "Density of Mean VAF at T1 and T2",
    x = "Density",
    y = "Mean Variant Allele Frequency (%)",
    fill = "Timepoint",
    color = "Timepoint"
  ) +
  theme_classic(base_size = 14)

# Filter the data (remove undefined RECIST and ensure cfDNA_T1_T2_dynamic is not NA)
plot_df <- cfDNA_df %>%
  filter(Best_RECIST != "Undefined" & !is.na(cfDNA_T1_T2_dynamic))

# Ensure Best_RECIST is numeric
plot_df$Best_RECIST <- as.numeric(plot_df$Best_RECIST)

# Scatterplot: cfDNA VAF Change vs. Best RECIST, colored by Cohort
ggplot(plot_df, aes(x = cfDNA_T1_T2_dynamic, y = Best_RECIST, color = Cohort)) +
  geom_point(size = 5, alpha = 0.6) +
  scale_color_manual(values = c("A" = "#A922D3", "B" = "#22E3CD", "C" = "#02008F")) +
  labs(
    x = "cfDNA mean VAF Change (T2 - T1)",
    y = "Best RECIST (% Tumor Change)",
    color = "Cohort"
  ) +
  theme_classic(base_size = 14) + 
  geom_vline(xintercept = 0, color = "red", size = 1)

# PLOT ALL PATIENTS
cfDNA_df <- cfDNA_df %>%
  mutate(
    PFS_group = if_else(
      PFS_months > 6,
      "PFS > 6 months",
      "PFS ≤ 6 months"
    )
  )

delta_df <- cfDNA_df %>%
  filter(!is.na(cfDNA_T1) & !is.na(cfDNA_T2)) %>%
  mutate(delta_cfDNA = cfDNA_T2 - cfDNA_T1)

p_val <- wilcox.test(delta_cfDNA ~ PFS_group, data = delta_df)$p.value
p_label <- paste0("p = ", signif(p_val, 2))

plotData <- cfDNA_df %>%
  pivot_longer(cols = c(cfDNA_T1, cfDNA_T2),
               names_to = "Timepoint", values_to = "cfDNA") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("cfDNA_T1", "cfDNA_T2"),
                            labels = c("T1", "T2")))

pfs_colors <- c(
  "PFS > 6 months" = "#ff2a67",
  "PFS ≤ 6 months" = "#9d9cd5"
)

ggplot(plotData, aes(x = Timepoint, y = cfDNA,
                     group = cfDNA_CMO_PT_ID,
                     color = PFS_group)) +
  geom_line(lineend = "round", size = 1.5) +
  geom_point(size = 3, shape = 16) +
  scale_color_manual(values = pfs_colors, name = "") +
  theme_classic() +
  labs(
    title = "cfDNA mean VAF Change\n 30 PFS outlier pairs from all cohorts",
    x = "6 weeks",
    y = "Mean Variant Allele Frequency %"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 14),
    legend.text  = element_text(size = 14),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim = c(0, max(plotData$cfDNA, na.rm = TRUE) * 1.1)) +
  annotate("text", x = 1.5, y = max(plotData$cfDNA, na.rm = TRUE) * 1.07,
           label = p_label, size = 7)

#===============================================================================
#         Suppl Fig 6. POLAR HRD. neoantigen, immune context correlation 
#===============================================================================
library(ggplot2); library(ggpubr); library(cowplot); library(dplyr)

# Read POLAR MAF with clinical data
POLAR_Baseline_WES <- read.csv(
  file.path("data", "POLAR_MS_Baseline_WES_ClinData.csv"),
  header = TRUE, check.names = FALSE)

POLAR_WES_df_NA <- read.csv(
  file.path("data", "POLAR_WES_df_NA.csv"),
  header = TRUE,
  check.names = FALSE
)

POLAR_WES_df_NA$Cohort <- factor(POLAR_WES_df_NA$Cohort, levels = c("A", "B", "C"))

cohort_palette <- c("A" = "#A922D3", "B" = "#22E3CD", "C" = "#02008F")
comparisons <- list(c("A", "B"), c("A", "C"))

make_box <- function(df, var, ylab) {
  mx <- max(df[[var]], na.rm = TRUE)
  ggplot(df, aes(x = Cohort, y = .data[[var]], fill = Cohort)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.9) +
    stat_compare_means(
      comparisons = comparisons, method = "wilcox.test",
      method.args = list(exact = FALSE), label = "p.signif",
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns")),
      tip.length = 0.02, label.y = mx * c(1.05, 1.12), size = 5
    ) +
    scale_fill_manual(values = cohort_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)),
                       breaks = scales::pretty_breaks(n = 5)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = var, title = ylab) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 6)),
      axis.title.y = element_text(size = 12, margin = margin(r = 5))
    )
}
#===============================================================================
#             Suppl Fig 6a. Suppl b. Purity, TMB, Ploidy, FGA
#===============================================================================
p_purity <- make_box(POLAR_WES_df_NA, "Purity", "Purity")
p_tmb    <- make_box(POLAR_WES_df_NA, "WES_TMB", "WES_TMB")
p_ploidy <- make_box(POLAR_WES_df_NA, "Ploidy", "Ploidy")
p_fga    <- make_box(POLAR_WES_df_NA, "FGA", "FGA")

top_title_1 <- ggdraw() + draw_label(
  "POLAR WES (N=35): Baseline Genomic Characteristics",
  fontface = "bold", size = 14, x = 0.5, hjust = 0.5)

grid_1 <- plot_grid(p_purity, p_tmb, p_ploidy, p_fga,
                    ncol = 4, align = "hv",
                    labels = c("i", "ii", "iii", "iv"),
                    label_fontface = "bold", label_size = 12,
                    label_x = 0, label_y = 1)

final_plot_1 <- plot_grid(top_title_1, grid_1, ncol = 1, rel_heights = c(0.1, 1))
print(final_plot_1)

#===============================================================================
#             Suppl Fig 6b. nsSNVs, Indel, fsindel, Neoantigen
#===============================================================================
# Create plots
p1 <- make_box(POLAR_WES_df_NA, "nsSNVs", "nsSNVs")
p2 <- make_box(POLAR_WES_df_NA, "Total_Indel_Burden", "Indels")
p3 <- make_box(POLAR_WES_df_NA, "Frameshift_Indel_Burden", "Frameshift Indels")
p4 <- make_box(POLAR_WES_df_NA, "Neoantigen_Burden", "Neoantigens")

# Assemble title + grid
top_title2 <- ggdraw() +
  draw_label("POLAR WES (N=35): Mutation Types and Neoantigen Burden",
             fontface = "bold", size = 14, x = 0.5, hjust = 0.5)

grid2 <- plot_grid(p1, p2, p3, p4, ncol = 4, align = "hv",
                   labels = c("i", "ii", "iii", "iv"),
                   label_fontface = "bold", label_size = 12,
                   label_x = 0, label_y = 1)

final_plot2 <- plot_grid(top_title2, grid2, ncol = 1, rel_heights = c(0.1, 1))
print(final_plot2)

#==================================================================================
#           Suppplementary Fig 5c. IMPACT_TMB (N=42), IMPACT_HRD (N=34)
#                   Neoantigen (N=35), TIL (N=38) by H&E 
#==================================================================================

# Load POLAR baseline WES clinical data
POLAR_WES_df <- read.csv(file.path("data","POLAR_63_6_25.csv"), check.names = FALSE)

# Define Cohort order
POLAR_WES_df$Cohort <- factor(POLAR_WES_df$Cohort, levels = c("A","B","C"))

# Define a quick plot function

make_box <- function(df, var, ylab) {
  mx <- max(df[[var]], na.rm = TRUE)
  ggplot(df, aes(x = Cohort, y = .data[[var]], fill = Cohort)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.9) +
    stat_compare_means(
      comparisons = comparisons, method = "wilcox.test",
      method.args = list(exact = FALSE), label = "p.signif",
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns")),
      tip.length = 0.02, label.y = mx * c(1.05, 1.12), size = 5
    ) +
    scale_fill_manual(values = cohort_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)),
                       breaks = scales::pretty_breaks(n = 5)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = var, title = ylab) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 6)),
      axis.title.y = element_text(size = 12, margin = margin(r = 5))
    )
}

# Individual plots
p_hrd <- make_box(POLAR_WES_df, "IMPACT_HRD", "IMPACT-HRD Score")
p_tmb <- make_box(POLAR_WES_df, "TMB", "TMB")
p_neo <- make_box(POLAR_WES_df, "Neoantigen", "Neoantigen")
p_til <- make_box(POLAR_WES_df, "TILs", "TIL")
p_nlr <- make_box(POLAR_WES_df, "NLR_T1", "Baseline NLR")

# Arrange plots with title
plot_combined <- ggarrange(p_hrd, p_tmb, p_neo, p_til, p_nlr, 
                           ncol = 5, align = "hv")
annotate_figure(plot_combined, top = text_grob("POLAR Baseline Genomic and Immune Features\nWES (N=35), IMPACT (N=42), H&E (N=38), Baseline blood (N=63)", face = "bold", size = 16))

#===============================================================================
#          Suppplementary Fig 6 - Immunofluorescent Imaging by Marc Hilmi
#===============================================================================
# 1) Load data

POLAR <- read_excel(file.path("data","POLAR_IF_final_NMV2.xlsx")) %>%
  transmute(
    POLAR_ID,
    Cohort = as.factor(Cohort),
    IMPACT_HRD,
    Neoantigen_Burden,
    CD3,
    CD3_CD8,
    CD68
  ) %>%
  mutate(
    CD4 = CD3 - CD3_CD8
  ) %>%
  filter(!is.na(Cohort))

cohort_palette <- c(
  "A" = "#A922D3",
  "B" = "#22E3CD",
  "C" = "#02008F"
)

# Generic scatter function
make_scatter <- function(df, x, y, xlabel, ylabel, title) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = Cohort)) +
    geom_smooth(
      method = "lm", se = TRUE,
      linetype = "dashed", color = "black"
    ) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = cohort_palette) +
    stat_cor(
      aes(group = 1),
      method = "pearson",
      label.x.npc = "left",
      label.y.npc = "top",
      size = 4,
      color = "black"
    ) +
    labs(
      title = title,
      x = xlabel,
      y = ylabel,
      color = "Cohort"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5)
    )
}

# a) Neoantigen vs IMPACT_HRD
p1 <- make_scatter(
  POLAR,
  x = "IMPACT_HRD",
  y = "Neoantigen_Burden",
  xlabel = "IMPACT_HRD",
  ylabel = "Neoantigen Burden",
  title  = "Neoantigen vs HRD"
)

# b) CD3 vs Neoantigen (Neoantigen on X)
p2 <- make_scatter(
  POLAR,
  x = "Neoantigen_Burden",
  y = "CD3",
  xlabel = "Neoantigen Burden",
  ylabel = "CD3+ T cells",
  title  = "CD3+T vs Neoantigen"
)

# c) Macrophages vs IMPACT_HRD
p3 <- make_scatter(
  POLAR,
  x = "IMPACT_HRD",
  y = "CD68",
  xlabel = "IMPACT_HRD",
  ylabel = "CD68+ Macrophages",
  title  = "Macrophages vs HRD"
)

# d) CD4 vs Neoantigen (Neoantigen on X)
p4 <- make_scatter(
  POLAR,
  x = "Neoantigen_Burden",
  y = "CD4",
  xlabel = "Neoantigen Burden",
  ylabel = "CD4+ T cells",
  title  = "CD4+T vs Neoantigen"
)

# 5) Combine panels (a–d)
(p1 + p2) / (p3 + p4) +
  plot_annotation(tag_levels = "a")