# =====================================================================
# THC / cortisol DDI — statistical analysis
### Summary: geometric mean (90% geometric CI)
### Paired t-test on log transformed data 
### Note that this script does run sequentially. 
# =====================================================================

# ---- Clean session ---------------------------------------------------
rm(list = ls(all.names = TRUE)); gc()
if (!is.null(dev.list())) graphics.off()
cat("\014")

# ---- Libraries -------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(outliers)
library(writexl)

# ---- Working directory -----------------------------------------------
#setwd() #set to working directory with code and dataframe

# ---- Load NCA output -------------------------------------------------
# One sheet per analyte. Each row is a participant (SID).

THC  <- read_excel("all_data_wide_HC.xlsx")                  # first sheet = parent THC
OH   <- read_excel("all_data_wide_HC.xlsx", sheet = "OH")    # 11-OH-THC
COOH <- read_excel("all_data_wide_HC.xlsx", sheet = "COOH")  # 11-COOH-THC
Gluc <- read_excel("all_data_wide_HC.xlsx", sheet = "Gluc")  # 11-COOH-THC-glucuronide
Misc <- read_excel("all_data_wide_HC.xlsx", sheet = "Misc")  # AUC ratios, COOH+Gluc

sheets <- list(THC = THC, OH = OH, COOH = COOH, Gluc = Gluc, Misc = Misc)

# =====================================================================
# Table 2. Geometric mean and 90% geometric confidence interval
# =====================================================================

# ---- Settings --------------------------------------------------------
CONF_LEVEL <- 0.90
OUT_FILE   <- "THC_HC_DDI_summary_table_GM_90CI.xlsx"

# ---- Exclusions ------------------------------------------------------
# One row per exclusion. param = NA drops that SID from the whole sheet.
# Otherwise only <param>_Control/_Treated/_TC are set to NA for that SID.
exclusions <- tribble(
  ~sheet, ~SID, ~param
  # "OH",   2,   "t_half",
  # "THC",  8,   NA
)

apply_exclusions <- function(df, sheet_name, excl) {
  ex <- filter(excl, sheet == sheet_name)
  for (i in seq_len(nrow(ex))) {
    rows <- df$SID == ex$SID[i]
    if (is.na(ex$param[i])) {
      df <- df[!rows, ]
    } else {
      cols <- paste0(ex$param[i], c("_Control", "_Treated", "_TC"))
      df[rows, intersect(cols, names(df))] <- NA
    }
  }
  df
}

sheets <- imap(sheets, ~ apply_exclusions(.x, .y, exclusions))

# ---- Helpers ---------------------------------------------------------

gm_ci <- function(x, conf = CONF_LEVEL) {
  x  <- x[is.finite(x) & x > 0]
  n  <- length(x)
  lx <- log(x)
  se <- sd(lx) / sqrt(n)
  tq <- qt(1 - (1 - conf) / 2, df = n - 1)
  tibble(
    gm = exp(mean(lx)),
    lo = exp(mean(lx) - tq * se),
    hi = exp(mean(lx) + tq * se),
    n  = n,
    p  = if (n > 2) t.test(lx)$p.value else NA_real_
  )
}

# 2 significant figures, keep trailing zeros (0.070, 1.0), integers at >= 100
fmt <- function(x) {
  out <- formatC(signif(x, 2), digits = 2, format = "fg", flag = "#")
  out <- sub("\\.$", "", out)
  big <- !is.na(x) & abs(x) >= 100
  out[big] <- format(round(x[big]), trim = TRUE)
  out
}

# ---- Table layout ----------------------------------------------------
# scale converts units. CLr is exported in L/h, so scale = 1000 gives mL/h.
table_spec <- tribble(
  ~Section,                    ~sheet, ~param,            ~label,                     ~scale,
  "THC",                       "THC",  "AUCinf",          "AUCinf (nM×h)",            1,
  "THC",                       "THC",  "t_half",          "t1/2 (h)",                 1,
  "THC",                       "THC",  "ke",              "ke (1/h)",                 1,
  "THC",                       "THC",  "CL_F",            "CL/F (L/h)",               1,
  "THC",                       "THC",  "Vb_F",            "Vb/F (L)",                 1,
  "11-OH-THC",                 "OH",   "AUCinf",          "AUCinf (nM×h)",            1,
  "11-OH-THC",                 "OH",   "t_half",          "t1/2 (h)",                 1,
  "11-OH-THC",                 "OH",   "ke",              "ke (1/h)",                 1,
  "11-OH-THC",                 "Misc", "OH_THC_AUCr",     "AUC11-OH-THC/AUCTHC",      1,
  "11-COOH-THC",               "COOH", "AUCinf",          "AUCinf (nM×h)",            1,
  "11-COOH-THC",               "COOH", "t_half",          "t1/2 (h)",                 1,
  "11-COOH-THC",               "COOH", "ke",              "ke (1/h)",                 1,
  "11-COOH-THC",               "Misc", "COOH_THC_AUCr",   "AUC11-COOH-THC/AUCTHC",    1,
  "11-COOH-THC",               "COOH", "CLr_COOH_24",     "CLr (mL/h)",               1000,
  "11-COOH-THC-glucuronide",   "Gluc", "AUCinf",          "AUCinf (nM×h)",            1,
  "11-COOH-THC-glucuronide",   "Gluc", "t_half",          "t1/2 (h)",                 1,
  "11-COOH-THC-glucuronide",   "Gluc", "ke",              "ke (1/h)",                 1,
  "11-COOH-THC-glucuronide",   "Misc", "Gluc_THC_AUCr",   "AUCGluc/AUCTHC",           1,
  "11-COOH-THC-glucuronide",   "Gluc", "CLr_Gluc_24",     "CLr (mL/h)",               1000
)

# ---- Summarise -------------------------------------------------------
summarise_param <- function(sheet, param, scale) {
  df <- sheets[[sheet]]
  map_dfr(c("Control", "Treated", "TC"), function(cond) {
    col <- paste0(param, "_", cond)
    if (!col %in% names(df)) stop("Missing column ", col, " in sheet ", sheet)
    x <- as.numeric(df[[col]])
    if (cond != "TC") x <- x * scale    # ratios are unitless
    gm_ci(x) %>% mutate(cond = cond, .before = 1)
  })
}

results_long <- table_spec %>%
  mutate(stats = pmap(list(sheet, param, scale), summarise_param)) %>%
  unnest(stats)

# ---- Format and pivot ------------------------------------------------
ci_lab <- paste0(CONF_LEVEL * 100, "% CI")

results_table <- results_long %>%
  mutate(cell = paste0(fmt(gm), " (", fmt(lo), ", ", fmt(hi), ")")) %>%
  select(Section, label, cond, cell) %>%
  pivot_wider(names_from = cond, values_from = cell) %>%
  rename(
    Parameter                                    = label,
    !!paste0("Control GM (", ci_lab, ")")         := Control,
    !!paste0("Treated GM (", ci_lab, ")")         := Treated,
    !!paste0("Treated/Control GMR (", ci_lab, ")") := TC
  )

print(results_table, n = Inf)


# ---- Save to Excel ---------------------------------------------------
write_xlsx(
  list(Table = results_table),
  OUT_FILE
)
cat("Saved:", file.path(getwd(), OUT_FILE), "\n")



# =====================================================================
# Paired t-test on log transformed data 
# =====================================================================

# ---- Diagnostics -----------------------------------------------------
# Text in value columns that becomes NA (e.g. #NUM!, #DIV/0!, embedded headers)
text_cells <- imap_dfr(sheets, function(d, sheet) {
  d %>%
    select(matches("_(Control|Treated|TC)$")) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(everything(), names_to = "Column", values_to = "Raw") %>%
    filter(!is.na(Raw), is.na(suppressWarnings(as.numeric(Raw)))) %>%
    distinct(Column, Raw) %>%
    mutate(Sheet = sheet, .before = 1)
})
print(text_cells, n = Inf)   # should be empty or error codes only

# Rows missing SID or Analyte (dropped below)
unlabeled <- imap_dfr(sheets, function(d, sheet) {
  d %>%
    mutate(Excel_row = row_number() + 1, .before = 1) %>%
    filter(is.na(SID) | is.na(Analyte)) %>%
    mutate(across(everything(), as.character), Sheet = sheet, .before = 1)
})
if (nrow(unlabeled) > 0) warning(nrow(unlabeled), " unlabeled row(s) dropped. See 'unlabeled'.")

# ---- Reshape to long -------------------------------------------------
# Split on the final suffix only, so CL_F and CLr_COOH_24 stay intact
df_long <- map_dfr(sheets, function(d) {
  d %>%
    filter(!is.na(SID), !is.na(Analyte)) %>%
    pivot_longer(
      cols          = matches("_(Control|Treated|TC)$"),
      names_to      = c("Parameter", "Condition"),
      names_pattern = "^(.*)_(Control|Treated|TC)$",
      values_to     = "Value",
      values_transform = list(Value = ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    mutate(SID = as.numeric(SID)) %>%
    select(Analyte, Study_Arm, SID, Parameter, Condition, Value)
}) %>%
  filter(Condition != "TC", !is.na(Value))   # TC recomputed from paired logs

# ---- Exclusions: profiles < 10 h -------------------------------------
flags <- tribble(
  ~Study_Arm, ~SID, ~Flag_Analyte,
  "HC",       10,   "OH",
  "HC",       17,   "COOH"
)

direct_params <- c("AUClast", "AUCinf", "ke", "p_extrapolated", "t_half")

derived_params <- list(
  THC  = c("CL_F", "Vb_F", "OH_THC_AUCr", "COOH_THC_AUCr", "Gluc_THC_AUCr"),
  OH   = c("OH_THC_AUCr", "COOH_OH_AUCr", "CLf_COOH_24"),
  COOH = c("COOH_OH_AUCr", "COOH_THC_AUCr", "Gluc_COOH_AUCr", "CLr_COOH_24", "CLf_Gluc_24")
)

direct_drop <- flags %>%
  crossing(Parameter = direct_params) %>%
  rename(Analyte = Flag_Analyte)

derived_drop <- flags %>%
  mutate(Parameter = derived_params[Flag_Analyte]) %>%
  unnest(Parameter) %>%
  select(Study_Arm, SID, Parameter)

excluded <- bind_rows(
  semi_join(df_long, direct_drop,  by = c("Study_Arm", "SID", "Analyte", "Parameter")),
  semi_join(df_long, derived_drop, by = c("Study_Arm", "SID", "Parameter"))
) %>%
  distinct() %>%
  arrange(SID, Analyte, Parameter, Condition)
print(excluded, n = Inf)

df_long <- df_long %>%
  anti_join(direct_drop,  by = c("Study_Arm", "SID", "Analyte", "Parameter")) %>%
  anti_join(derived_drop, by = c("Study_Arm", "SID", "Parameter"))

# ---- Assign AUC ratios to the metabolite; log transform --------------
df_long <- df_long %>%
  mutate(
    Analyte = case_when(
      Parameter == "OH_THC_AUCr"                          ~ "OH",
      Parameter %in% c("COOH_OH_AUCr", "COOH_THC_AUCr")   ~ "COOH",
      Parameter %in% c("Gluc_COOH_AUCr", "Gluc_THC_AUCr") ~ "Gluc",
      TRUE ~ Analyte
    ),
    logValue = log(Value)
  )

dups <- count(df_long, Analyte, SID, Parameter, Condition) %>% filter(n > 1)
if (nrow(dups) > 0) { print(dups, n = Inf); stop("Duplicate rows found. See above.") }

# ---- Pre-specified parameters ----------------------------------------
# FDR families per Methods. ke omitted: same test as t_half.
df_ttest <- df_long %>%
  filter(
    (Parameter == "t_half" & Analyte %in% c("THC", "OH", "COOH", "Gluc")) |
      (Analyte == "THC"  & Parameter %in% c("Vb_F", "CL_F")) |
      (Analyte == "OH"   & Parameter %in% c("OH_THC_AUCr")) |
      (Analyte == "COOH" & Parameter %in% c("CLr_COOH_24", "COOH_OH_AUCr", "COOH_THC_AUCr")) |
      (Analyte == "Gluc" & Parameter %in% c("CLr_Gluc_24", "Gluc_COOH_AUCr", "Gluc_THC_AUCr"))
  )

# ---- Paired t-test on log values (Treated vs Control) ----------------
# GMR with 90% CI. Grubbs on paired differences is diagnostic only.
paired_test <- function(d) {
  w <- d %>%
    select(SID, Condition, logValue) %>%
    pivot_wider(names_from = Condition, values_from = logValue)
  if (!all(c("Control", "Treated") %in% names(w))) return(tibble())
  
  w <- filter(w, is.finite(Control), is.finite(Treated))
  if (nrow(w) < 3) return(tibble())
  
  tt   <- t.test(w$Treated, w$Control, paired = TRUE, conf.level = 0.90)
  diff <- w$Treated - w$Control
  g    <- tryCatch(grubbs.test(diff), error = function(e) NULL)
  
  tibble(
    N             = nrow(w),
    SIDs          = paste(sort(w$SID), collapse = ", "),
    GMR           = exp(unname(tt$estimate)),
    GMR_90CI_low  = exp(tt$conf.int[1]),
    GMR_90CI_high = exp(tt$conf.int[2]),
    t_stat        = unname(tt$statistic),
    p_raw         = tt$p.value,
    G_diff        = if (is.null(g)) NA_real_ else unname(g$statistic[["G"]]),
    p_G_diff      = if (is.null(g)) NA_real_ else g$p.value
  )
}

paired_HC <- df_ttest %>%
  group_by(Study_Arm, Analyte, Parameter) %>%
  group_modify(~ paired_test(.x)) %>%
  group_by(Study_Arm, Analyte) %>%
  mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>%   # BH within analyte
  ungroup()

# ---- Clean p-value table (full precision) ----------------------------
compound_labels <- c(
  THC  = "THC",
  OH   = "11-OH-THC",
  COOH = "11-COOH-THC",
  Gluc = "11-COOH-THC-glucuronide"
)

param_labels <- c(
  t_half         = "t1/2",
  CL_F           = "CL/F",
  Vb_F           = "Vb/F",
  OH_THC_AUCr    = "AUC11-OH-THC/AUCTHC",
  CLr_COOH_24    = "CLr",
  COOH_OH_AUCr   = "AUC11-COOH-THC/AUC11-OH-THC",
  COOH_THC_AUCr  = "AUC11-COOH-THC/AUCTHC",
  CLr_Gluc_24    = "CLr",
  Gluc_COOH_AUCr = "AUCGluc/AUC11-COOH-THC",
  Gluc_THC_AUCr  = "AUCGluc/AUCTHC"
)

p_table <- paired_HC %>%
  mutate(
    Compound  = factor(Analyte, levels = names(compound_labels), labels = compound_labels),
    order     = match(Parameter, names(param_labels)),
    Parameter = unname(param_labels[Parameter])
  ) %>%
  arrange(Compound, order) %>%
  select(Compound, Parameter, p_raw, p_FDR)

print(as.data.frame(p_table), digits = 6)

# ---- Export ----------------------------------------------------------
write_xlsx(
  list(P_values = p_table, Full_results = paired_HC),
  "paired_ttests_HC.xlsx"
)

# =====================================================================
# Supplemental Table 4 AUClast and percent extrapolated
# =====================================================================

# All participants included (no exclusions applied)
analyte_labels <- c(THC  = "THC",
                    OH   = "11-OH-THC",
                    COOH = "11-COOH-THC",
                    Gluc = "11-COOH-THC-glucuronide")

param_labels <- c(AUClast        = "AUClast (nMxh)",
                  p_extrapolated = "Percent extrapolated (%)")

# 2 significant figures, whole numbers at 100 and above
fmt <- function(x) {
  ifelse(x >= 100,
         formatC(round(x), format = "d"),
         sub("\\.$", "", formatC(signif(x, 2), digits = 2, format = "fg", flag = "#")))
}
fmt_ci <- function(est, lo, hi) paste0(fmt(est), " (", fmt(lo), ", ", fmt(hi), ")")

# Rebuild from the sheets so the exclusions above are not applied
df_s4 <- map_dfr(sheets, to_long) %>%
  filter(Condition != "TC",
         !is.na(Value),
         Analyte %in% names(analyte_labels),
         Parameter %in% names(param_labels)) %>%
  mutate(logValue = log(Value)) %>%
  filter(is.finite(logValue))

# Geometric mean (90% CI) per study day, all available participants
arm_tbl <- df_s4 %>%
  group_by(Analyte, Parameter, Condition) %>%
  summarise(n  = n(),
            m  = mean(logValue),
            hw = qt(0.95, n - 1) * sd(logValue) / sqrt(n),
            .groups = "drop") %>%
  mutate(cell = fmt_ci(exp(m), exp(m - hw), exp(m + hw))) %>%
  select(Analyte, Parameter, Condition, n, cell) %>%
  pivot_wider(names_from = Condition, values_from = c(cell, n))

# Treated/Control GMR (90% CI), complete pairs only
pair_tbl <- df_s4 %>%
  select(Analyte, Parameter, SID, Condition, logValue) %>%
  pivot_wider(names_from = Condition, values_from = logValue) %>%
  filter(is.finite(Control), is.finite(Treated)) %>%
  mutate(d = Treated - Control) %>%
  group_by(Analyte, Parameter) %>%
  summarise(n_pairs = n(),
            m  = mean(d),
            hw = qt(0.95, n_pairs - 1) * sd(d) / sqrt(n_pairs),
            .groups = "drop") %>%
  mutate(TC = fmt_ci(exp(m), exp(m - hw), exp(m + hw))) %>%
  select(Analyte, Parameter, n_pairs, TC)

supp4 <- arm_tbl %>%
  left_join(pair_tbl, by = c("Analyte", "Parameter")) %>%
  arrange(match(Analyte, names(analyte_labels)),
          match(Parameter, names(param_labels)))

# N used per cell
print(supp4 %>% select(Analyte, Parameter, n_Control, n_Treated, n_pairs), n = Inf)

# Formatted table with analyte header rows
supp4_out <- map_dfr(names(analyte_labels), function(a) {
  s <- supp4 %>% filter(Analyte == a)
  bind_rows(
    tibble(Row = analyte_labels[[a]], Control = "", Treated = "", `Treated/Control` = ""),
    tibble(Row               = unname(param_labels[s$Parameter]),
           Control           = s$cell_Control,
           Treated           = s$cell_Treated,
           `Treated/Control` = s$TC)
  )
})

print(supp4_out, n = Inf)
write_xlsx(supp4_out, "Supp_Table4_AUClast_pext.xlsx")
