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

POLAR_AE <- read_excel(file.path("data","POLAR_AEs.xlsx"))

A1 = POLAR_clindata %>% select('Cohort', 'Age',	'Sex', 'Race', 'ECOG', 'Initial_Stage', 'Histology', 'Previous_Surgery', 'Platinum_Type', 'CA19_9_T1', 'CEA_T1')
A2 = POLAR_clindata %>% select('Cohort', 'FHx_of_PDAC','POLAR_BOR','TMB', 'IMPACT_HRD', 'Neoantigen', 'TIL_Density', 'NLR_T1', 'Mean_VAF_T1')
# POLAR_AE %>% select('Study_ID',	'Event_Name',	'Grade',	'Attribution',	'Intervention')

#gtsummary tbl_summary for table.
tbl_summary(A1, by = Cohort) %>% add_n() %>% add_overall()
tbl_summary(A2, by = Cohort) %>% add_n() %>% add_overall()

#===============================================================================
#                              Fig 2a. WATERFALL
#===============================================================================

polar_waterfall <- read.csv(file.path("data","POLAR.RECIST.csv"))

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

library(survival)
cox1 = coxph(Surv())

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

# Set the theme to remove background gridlines
theme_set(theme_bw() + theme(panel.grid = element_blank()))

tumorresponse <- read.csv(file.path("data","POLAR.RECIST.csv"))

create_spider_plot <- function(tumorresponse, cohort_label, y_limits, y_breaks) {
  # Subset the data for the specified cohort
  cohort_data <- subset(tumorresponse, Cohort == cohort_label)
  
  # Ensure Time is properly formatted and remove NA values
  data_long <- cohort_data %>%
    gather(key = "Time", value = "Change", -POLAR.ID, -Cohort, -TIL) %>%
    filter(!is.na(Change)) %>%  # Filter out rows where Change is NA
    mutate(Time = as.numeric(as.character(Time))) %>%  # Convert Time to numeric
    filter(!is.na(Time))  # Remove rows where Time couldn't be converted to numeric
  
  # Ensure Change column is numeric
  data_long$Change <- as.numeric(as.character(data_long$Change))
  
  # Filter out any rows where Change could not be converted to numeric
  data_long <- data_long %>% filter(!is.na(Change))
  
  # Create the plot
  ggplot(data_long, aes(x = (Time), y = Change, group = POLAR.ID)) +
    geom_line(aes(color = TIL)) + 
    geom_point(aes(color = TIL)) +
    scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], by = y_breaks)) +
    scale_x_continuous(breaks = seq(0, 140, by = 12)) +  
    labs(title = paste("Cohort", cohort_label, "% Change of Target Lesions"), x = "Weeks", y = "% Change") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
      axis.title.y = element_text(vjust = 2),
      axis.text.y = element_text(vjust = 0.5),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank()   
    ) +
    scale_color_gradient(low = "grey", high = "red", name = "TIL") +
    scale_fill_gradient(low = "grey", high ="red", name = "TIL") +
    geom_hline(yintercept = 20, linetype = "dashed", color = "red2") +
    geom_hline(yintercept = -30, linetype = "dashed", color = "blue")
}

# Create plots for each cohort
print(create_spider_plot(tumorresponse, "A", c(-100, 110), 10))
create_spider_plot(tumorresponse, "B", c(-60, 60), 10)
create_spider_plot(tumorresponse, "C", c(-100, 70), 10)

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
#                 Fig 3a. Swimmer's plot by Walid Chatila
#===============================================================================
#===============================================================================
#                      Fig 3b. POLAR WES Oncoprint
#===============================================================================
POLAR_Baseline_WES <- read.csv(file.path("data","POLAR_MS_Baseline_WES_ClinData.csv"), check.names = FALSE)
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
#primary_Metastasis = c("Primary" = "#9d9cd5", "Metastasis" = "#299d3e")
#sample_Origin = c("Liver" = "#b3e1bf", "Lung" =  "#011a52", "Pancreas" = "#9d9cd5", "Ovary" = "#916c36", "Lymph_node" =  "#a31300", "Peritoneum" = "#f6b3b0")
kRAS = c("G12D" = "#17154f", "G12V" = "#98d1d9", "G12R" = "#9d9cd5", "Q61H" = "#a5506d", "G12C" = "#edf181", "Q61R" = "#ada43b", "G13D" = "#ff2a67", "I171Nfs*14" = "#2f357c")

dDRm = c("gATMm" = "#299d3e", "gATMm_gMUTYHm" =  "#ada43b",  "gBLMm" =  "#e69b00", "gBRCA1m" = "#011a52", "gBRCA2m" =  "#ff2a67", "gBRCA2m_gFLCNm"= "#ff3100", "gCHEK2m" = "#9d9cd5",
         "gFANCCm" = "#355828", "gPALB2m" = "#a5506d", "gMUTYHm" = "#f6b3b0", "sBRCA2m" = "#ff2a67")

#zygosity = c("Monoallelic"= "#17154f", "Biallelic"= "#ff2a67")
#plasticity = c("Classical" ="#9d9cd5", "Basal-like" =  "#fba702", "N/A" =  "#98d1d9")
#hRDLOH = c("LOH" = "#ff3100", "HOMDEL" = "#a31300", "LOH_LOH" = "#ff2a67", "None" = "#9d9cd5")
wGD = c("TRUE" = "#ff2a67", "FALSE" = "#011a52")
quartile = c("1" = "#2f357c", "2" = "#A6CEE3",  "3" = "#a5506d", "4"= "#ff2a67")
titv = c("C>T" = "#9d9cd5", "C>G" = "#ff2a67", "C>A" = "#011a52", "T>A" =  "#318c97", "T>C" = "#98d1d9", "T>G" = "#fba702")
titvcolors = list(titv_cols = titv)
iPC_genes = c("KRAS", "TP53", "CDKN2A", "SMAD4", "BRCA2", "BRCA1", "PALB2", "ATM")
time_point = c("T0" = "deepskyblue", "T1" = "#A6CEE3", "T2" = "#011a52", "T3" = "purple")
pre_Post = c("PRE" ="#ff2a67", "POST" = "#011a52")

#Definte annocolors before making maf oncoplot
annocolors = list(Cohort = cohort_palette, HRD_gene = dDRm, KRAS = kRAS, WGD = wGD, Time_point= time_point, Pre_Post = pre_Post, IMPACT_HRD_Quartile = quartile, Neoantigen_Quartile = quartile)

#Variant allele frequencies (Left bar plot)
POLAR_Baseline_genes_vaf = subsetMaf(maf = POLAR_Baseline, genes = iPC_genes, fields = "t_var_freq", mafObj = FALSE)[,mean(t_var_freq, na.rm = TRUE), Hugo_Symbol]

#Oncoplot
dir.create("outputs", showWarnings = FALSE)
dir.create(file.path("outputs","figures"), recursive = TRUE, showWarnings = FALSE)

oncoplot(maf = POLAR_Baseline, 
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
         numericAnnoCol = numericcol,
         clinicalFeatures = c('Cohort', 'KRAS', 'HRD_gene', 'WGD', 'IMPACT_HRD_Quartile', 'Neoantigen_Quartile'),
         sortByAnnotation = TRUE,
         annotationOrder = c('A', 'B', 'C'),
         annotationColor = annocolors,
         annoBorderCol = "white",
         titv_col = titv,
         additionalFeatureCol = "gray70",
         bgCol = "#CCCCCC",
         borderCol = "white",
         draw_titv = TRUE)

pdf(file.path("outputs","figures","POLAR_Oncoplot_Final.pdf"), width=14, height=10, useDingbats=FALSE)
out_pdf <- file.path("outputs","figures","3dii_CD8T_boxplot_by_Cohort.pdf")

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
POLAR_Baseline_WES <- read.csv("/XX/POLAR_MS_Baseline_WES_ClinData.csv", header = TRUE, check.names = FALSE)
POLAR_WES_df_NA <- read.csv("/XX/POLAR_WES_df_NA.csv", header = TRUE)
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

final_plot2 <- plot_grid(top_title, grid2, ncol = 1, rel_heights = c(0.1, 1))
print(final_plot2)
#==================================================================================
#           Suppplementary Fig 5c. IMPACT_TMB (N=42), IMPACT_HRD (N=34)
#                   Neoantigen (N=35), TIL (N=38) by H&E 
#==================================================================================

# Load POLAR baseline WES clinical data
POLAR_WES_df <- read.csv("/XX/POLAR_63_6_25.csv", header = TRUE, check.names = FALSE)

# Define Cohort order
POLAR_WES_df$Cohort <- factor(POLAR_WES_df$Cohort, levels = c("A", "B", "C"))

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
#          Figure 3D - Immunofluorescent  Imaging by Marc Hilmi
#===============================================================================
# -------------------------
# Import & preparation - Survival groups
# -------------------------
df <- read_excel('/XX/POLAR_IF_final_NMV2.xlsx') 

df <- df %>%
  mutate(
    PFS_group = ifelse(PFS <= 4, "PFS ≤ 4 months", "PFS > 4 months"),
    PFS_group = factor(PFS_group, levels = c("PFS ≤ 4 months", "PFS > 4 months"))
  ) %>%
  filter(!is.na(CD3_CD8))

# -------------------------
# Statistics
# -------------------------
wilcox_res <- wilcox.test(CD3_CD8 ~ PFS_group, data = df, exact = FALSE)
p_label <- paste0("Wilcoxon p = ", formatC(wilcox_res$p.value, format = "e", digits = 2))

y_max <- max(df$CD3_CD8, na.rm = TRUE)

# -------------------------
# Colors 
# -------------------------
group_colors <- c(
  "PFS ≤ 4 months" = "#4E79A7",   # muted blue
  "PFS > 4 months" = "#E15759"    # muted red
)

# -------------------------
# Plot
# -------------------------
p <- ggplot(df, aes(x = PFS_group, y = CD3_CD8, color = PFS_group)) +
  
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.8
  ) +
  
  geom_jitter(
    width = 0.15,
    size = 2.3,
    alpha = 0.85
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = y_max * 1.05,
    label = p_label,
    size = 4.5
  ) +
  
  scale_color_manual(values = group_colors) +
  
  labs(
    title = "CD3+CD8+ infiltrate",
    x = NULL,
    y = "Stromal immune infiltrate (% positive cells)"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(size = 13),
    axis.title.y = element_text(size = 14),
    axis.line    = element_line(linewidth = 0.8),
    axis.ticks   = element_line(linewidth = 0.8),
    plot.margin  = margin(10, 10, 10, 10)
  ) +
  
  coord_cartesian(ylim = c(0, y_max * 1.15))

p

out_pdf <- "/xx/3di_CD3_CD8_boxplot_by_PFS4M.pdf"

ggsave(
  filename = out_pdf,
  plot     = p,
  width    = 6.5,
  height   = 5.2,
  units    = "in"
)

# -------------------------
# Import & preparation - Cohorts
# -------------------------

df <- df %>%
  mutate(Cohort = factor(Cohort, levels = c("A","B","C"))) %>%
  filter(!is.na(CD3_CD8), !is.na(Cohort))

# p values
p_CA <- wilcox.test(CD3_CD8 ~ Cohort, data = df %>% filter(Cohort %in% c("A","C")), exact = FALSE)$p.value
p_CB <- wilcox.test(CD3_CD8 ~ Cohort, data = df %>% filter(Cohort %in% c("B","C")), exact = FALSE)$p.value

label_CA <- paste0("p = ", formatC(p_CA, format = "e", digits = 2))
label_CB <- paste0("p = ", formatC(p_CB, format = "e", digits = 2))

y_max <- max(df$CD3_CD8, na.rm = TRUE)

# heights
y_CA <- y_max * 1.05
y_CB <- y_max * 1.20
y_CA_text <- y_max * 1.10
y_CB_text <- y_max * 1.25
tick <- y_max * 0.02

cohort_colors <- c("A"="#A922D3","B"="#22E3CD","C"="#02008F")

# -------------------------
# Helper: add bracket + label
# -------------------------
add_bracket <- function(p, x1, x2, y, y_text, label, tick, lw = 0.8, text_size = 4.3) {
  p +
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = lw) +
    annotate("segment", x = x1, xend = x1, y = y - tick, yend = y, linewidth = lw) +
    annotate("segment", x = x2, xend = x2, y = y - tick, yend = y, linewidth = lw) +
    annotate("text", x = (x1 + x2) / 2, y = y_text, label = label, size = text_size)
}

# -------------------------
# Plot
# -------------------------
q <- ggplot(df, aes(Cohort, CD3_CD8, color = Cohort)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 2.3, alpha = 0.85) +
  scale_color_manual(values = cohort_colors) +
  labs(title = "CD3+CD8+ infiltrate", x = NULL, y = "Stromal immune infiltrate (% positive cells)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(size = 13),
    axis.title.y = element_text(size = 14),
    axis.line    = element_line(linewidth = 0.8),
    axis.ticks   = element_line(linewidth = 0.8),
    plot.margin  = margin(10, 10, 10, 10)
  ) +
  coord_cartesian(ylim = c(0, y_max * 1.32))

# add the two comparisons
q <- q %>%
  add_bracket(x1 = 1, x2 = 3, y = y_CA, y_text = y_CA_text, label = label_CA, tick = tick) %>%
  add_bracket(x1 = 2, x2 = 3, y = y_CB, y_text = y_CB_text, label = label_CB, tick = tick)

q

out_pdf <- "/xx/3dii_CD8T_boxplot_by_Cohort.pdf"

ggsave(
  filename = out_pdf,
  plot     = q,
  width    = 6.5,
  height   = 5.2,
  units    = "in"
)


# -------------------------
# Import & preparation
# -------------------------
df <- read_excel(file.path("data","POLAR_IF_final_NMV2.xlsx"))

df <- df %>%
  mutate(
    Cohort = factor(Cohort, levels = c("A", "B", "C"))
  ) %>%
  filter(!is.na(CD3_CD8), !is.na(Cohort))

# -------------------------
# Statistics
# -------------------------
# Global test
kw_res <- kruskal.test(CD3_CD8 ~ Cohort, data = df)

# Pairwise tests of interest
p_CA <- wilcox.test(CD3_CD8 ~ Cohort, data = df %>% filter(Cohort %in% c("A", "C")), exact = FALSE)$p.value
p_CB <- wilcox.test(CD3_CD8 ~ Cohort, data = df %>% filter(Cohort %in% c("B", "C")), exact = FALSE)$p.value

label_CA <- paste0("p = ", formatC(p_CA, format = "e", digits = 2))
label_CB <- paste0("p = ", formatC(p_CB, format = "e", digits = 2))

y_max <- max(df$CD3_CD8, na.rm = TRUE)

# -------------------------
# Colors 
# -------------------------
cohort_colors <- c(
  "A" = "#A922D3",
  "B" = "#22E3CD",
  "C" = "#02008F"
)

# -------------------------
# Plot
# -------------------------
p <- ggplot(df, aes(x = Cohort, y = CD3_CD8, color = Cohort)) +
  
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.8
  ) +
  
  geom_jitter(
    width = 0.15,
    size = 2.3,
    alpha = 0.85
  ) +
  
  # --- C vs A ---
  annotate("segment",
           x = 1, xend = 3,
           y = y_max * 1.05, yend = y_max * 1.05,
           linewidth = 0.8) +
  annotate("segment",
           x = 1, xend = 1,
           y = y_max * 1.03, yend = y_max * 1.05,
           linewidth = 0.8) +
  annotate("segment",
           x = 3, xend = 3,
           y = y_max * 1.03, yend = y_max * 1.05,
           linewidth = 0.8) +
  annotate("text",
           x = 2,
           y = y_max * 1.07,
           label = label_CA,
           size = 4.3) +
  
  # --- C vs B ---
  annotate("segment",
           x = 2, xend = 3,
           y = y_max * 1.13, yend = y_max * 1.13,
           linewidth = 0.8) +
  annotate("segment",
           x = 2, xend = 2,
           y = y_max * 1.11, yend = y_max * 1.13,
           linewidth = 0.8) +
  annotate("segment",
           x = 3, xend = 3,
           y = y_max * 1.11, yend = y_max * 1.13,
           linewidth = 0.8) +
  annotate("text",
           x = 2.5,
           y = y_max * 1.15,
           label = label_CB,
           size = 4.3) +
  
  scale_color_manual(values = cohort_colors) +
  
  labs(
    title = "CD3+CD8+ infiltrate",
    x = NULL,
    y = "Stromal immune infiltrate (% positive cells)"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(size = 13),
    axis.title.y = element_text(size = 14),
    axis.line    = element_line(linewidth = 0.8),
    axis.ticks   = element_line(linewidth = 0.8),
    plot.margin  = margin(10, 10, 10, 10)
  ) +
  
  coord_cartesian(ylim = c(0, y_max * 1.22))

p

#===============================================================================
#          Suppplementary Fig 6 - Immunofluorescent Imaging by Marc Hilmi
#===============================================================================
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# -----------------------------
# 1) Load data
# -----------------------------
POLAR <- read_excel('/xx/POLAR_IF_final_NMV2.xlsx') %>%
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

# -----------------------------
# 2) Cohort colors
# -----------------------------
cohort_palette <- c(
  "A" = "#A922D3",
  "B" = "#22E3CD",
  "C" = "#02008F"
)

# -----------------------------
# 3) Generic scatter function
# -----------------------------
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

# -----------------------------
# 4) Plots (axes EXACTLY as requested)
# -----------------------------

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

# -----------------------------
# 5) Combine panels (a–d)
# -----------------------------
(p1 + p2) / (p3 + p4) +
  plot_annotation(tag_levels = "a")


