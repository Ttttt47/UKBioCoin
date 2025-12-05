
source("functions.R")

### Evaluation 
setwd("simulation")
library(data.table)
library(ggplot2)

B_start = 1
B_end = 100

# simulation settings
all_traits = c("continuous")
all_sig_mafs = c(0.1)
all_h2s = c(0.3)
all_mcars = c(0.1)
all_ns = c(1000, 5000, 10000, 100000)
B_start = 1
B_end = 100
B = B_end - B_start + 1
prefix = "beta_compare"


## --- make root dir for evaluation -----------------------------------------
eval_root <- file.path("evaluation", prefix)
dir.create(eval_root, recursive = TRUE, showWarnings = FALSE)
simu_res_root <- file.path("simulation_results", prefix)

## --- iterate over the factorial design ---------------------------------
p_thres = 0.05/5000
eval_results = list()

for (trait   in all_traits)
for (sig_maf in all_sig_mafs)
for (h2      in all_h2s)
for (mcar    in all_mcars)
for (n       in all_ns)
 {

  one_errors = matrix(0, nrow=B, ncol=3)
  signed_errors = matrix(0, nrow=B, ncol=3)
  for (b        in B_start:B_end) {
    out_pref_b = paste0(trait, "_", sig_maf, "_", h2, "_", mcar, "_", b)
    res_dir = paste0(simu_res_root,'/',paste0(trait, "_", sig_maf, "_", h2, "_", mcar, "_n", n),'/',out_pref_b)
    # load simulation result
    sim <- readRDS(paste0(res_dir,".rds"))
    
    # load PLINK result
    if (trait == "continuous") {
        plink_res_correct <- fread(paste0(res_dir,"_correct.PHENO.glm.linear"), header=T)
    } else {
        # binary trait
        plink_res_correct <- fread(paste0(res_dir,"_correct.PHENO.glm.logistic.hybrid"), header=T)
    }    # load UKC result
    ukc_res <- fread(paste0(res_dir,"_UKC_results.table"), header=T)
    
    # calculate error of beta
    true_beta = sim$beta/sqrt(2*plink_res_correct$A1_FREQ*(1-plink_res_correct$A1_FREQ))[1:10]
    plink_beta_correct = as.numeric(plink_res_correct$BETA)[1:length(true_beta)]
    ukc_beta = as.numeric(ukc_res$BETA)[1:length(true_beta)]
    one_errors[b-B_start+1,1] = mean(abs(true_beta - plink_beta_correct))
    signed_errors[b-B_start+1,1] = mean(true_beta - plink_beta_correct)
    one_errors[b-B_start+1,2] = mean(abs(plink_beta_correct - ukc_beta))
    signed_errors[b-B_start+1,2] = mean(plink_beta_correct - ukc_beta)
    one_errors[b-B_start+1,3] = mean(abs(true_beta - ukc_beta))
    signed_errors[b-B_start+1,3] = mean(true_beta - ukc_beta)


  }
  # calculate TPR and FPR
  eval_results[[paste0(trait, "_", sig_maf, "_", h2, "_", mcar, "_n", n)]] = list(
    abs_sample_error = mean(one_errors[,1]),
    abs_approx_error = mean(one_errors[,2]),
    abs_total_error = mean(one_errors[,3]),
    signed_sample_error = mean(signed_errors[,1]),
    signed_approx_error = mean(signed_errors[,2]),
    signed_total_error = mean(signed_errors[,3]),
    sd_signed_sample_error = sd(signed_errors[,1]),
    sd_signed_approx_error = sd(signed_errors[,2]),
    sd_signed_total_error = sd(signed_errors[,3])
  ) 

  print(paste0(trait, "_", sig_maf, "_", h2, "_", mcar, "_n", n))
}

# save results
saveRDS(eval_results, file=paste0("evaluation/eval_results_",prefix,".rds"))

# load results
eval_results = readRDS(paste0("evaluation/eval_results_",prefix,".rds"))

# plot results
library(ggplot2)
library(cowplot)

# plot error accros different n
error_df = do.call(rbind, lapply(names(eval_results), function(x) {
  n = strsplit(x, "_")[[1]][5]
  n = substr(n, 2, nchar(n))  # remove 'n' prefix
  n = as.numeric(n)
  data.frame(
    trait = strsplit(x, "_")[[1]][1],
    sig_maf = as.numeric(strsplit(x, "_")[[1]][2]),
    h2 = as.numeric(strsplit(x, "_")[[1]][3]),
    mcar = as.numeric(strsplit(x, "_")[[1]][4]),
    n = n,
    abs_sample_error = eval_results[[x]]$abs_sample_error,
    abs_approx_error = eval_results[[x]]$abs_approx_error,
    abs_total_error = eval_results[[x]]$abs_total_error
  )
}))

error_df_long = reshape2::melt(error_df, id.vars = c("trait", "sig_maf", "h2", "mcar", "n"),
                               measure.vars = c("abs_sample_error", "abs_approx_error", "abs_total_error"),
                               variable.name = "error_type", value.name = "error_value")

pdf(file = paste0("evaluation/error_plot_", prefix, ".pdf"), width = 7, height = 5)
p1 = ggplot(error_df_long, aes(x = n, y = error_value, color = error_type)) +
  geom_line() +
  geom_point() +
  scale_x_log10(
    breaks = c(1000, 5000, 10000, 100000),
    labels = function(x) format(x, scientific = FALSE)
  ) +
  labs(title = "Absolute Error Comparison across Sample Sizes",
       x = "Sample Size (n)",
       y = "Absolute Error") +
  theme_bw() +
  theme(legend.position = "bottom")
p1
dev.off()

#############################
# signed error and sd.
#############################

signed_error_df <- do.call(rbind, lapply(names(eval_results), function(x) {
  
  parts <- strsplit(x, "_")[[1]]
  n_val <- as.numeric(sub("n", "", parts[5]))

  data.frame(
    trait  = parts[1],
    sig_maf = as.numeric(parts[2]),
    h2      = as.numeric(parts[3]),
    mcar    = as.numeric(parts[4]),
    n       = n_val,

    # Signed errors
    signed_sample_error = eval_results[[x]]$signed_sample_error,
    signed_approx_error = eval_results[[x]]$signed_approx_error,
    signed_total_error  = eval_results[[x]]$signed_total_error,

    # SDs
    sd_sample_error = eval_results[[x]]$sd_sample_error,
    sd_approx_error = eval_results[[x]]$sd_approx_error,
    sd_total_error  = eval_results[[x]]$sd_total_error
  )
}))

#########################################################
# 2. Convert to long form for plotting with ribbons
#########################################################
# Long data for signed errors
signed_long <- melt(
  signed_error_df,
  id.vars = c("trait","sig_maf","h2","mcar","n"),
  measure.vars = c("signed_sample_error","signed_approx_error","signed_total_error"),
  variable.name = "error_type",
  value.name = "signed_error"
)

# Long data for corresponding SDs
sd_long <- melt(
  signed_error_df,
  id.vars = c("trait","sig_maf","h2","mcar","n"),
  measure.vars = c("sd_sample_error","sd_approx_error","sd_total_error"),
  variable.name = "sd_type",
  value.name = "sd"
)

# Make the labels match so we can merge
sd_long$sd_type <- gsub("sd_", "signed_", sd_long$sd_type)

# Merge signed error with SD
signed_plot_df <- left_join(signed_long, sd_long,
                            by = c("trait","sig_maf","h2","mcar","n",
                                   "error_type" = "sd_type"))

# Compute ±2 SD bands
signed_plot_df <- signed_plot_df %>%
  mutate(
    upper = signed_error + 1.96 * sd,
    lower = signed_error - 1.96 * sd
  )

#########################################################
# 3. Colors consistent with absolute error plot
#########################################################
signed_colors <- c(
  "signed_sample_error" = "#F8766D",  # red
  "signed_approx_error" = "#00BA38",  # green
  "signed_total_error"  = "#619CFF"   # blue
)

#########################################################
# 4. Plot with ±2 SD ribbons
#########################################################
pdf(file = paste0("evaluation/signed_error_plot_with_sd_", prefix, ".pdf"),
    width = 7, height = 5)

p2 = ggplot(signed_plot_df,
       aes(x = n, y = signed_error, color = error_type)) +

  # ±2 SD shaded region
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = error_type),
              alpha = 0.20, color = NA) +

  geom_line() +
  geom_point() +

  scale_color_manual(values = signed_colors) +
  scale_fill_manual(values = signed_colors) +

  scale_x_log10(
    breaks = c(1000, 5000, 10000, 100000),
    labels = function(x) format(x, scientific = FALSE)
  ) +

  labs(
    title = "Error with 95% CI Bands Across Sample Sizes",
    x = "Sample Size (n)",
    y = "Error",
    color = "error_type",
    fill  = "error_type"
  ) +

  theme_bw() +
  theme(legend.position = "bottom")
p2
dev.off()

# Combine both plots
pdf(file = paste0("evaluation/combined_error_plots_", prefix, ".pdf"),
    width = 10, height = 4)
combined_plot <- plot_grid(p1, p2, ncol = 2)
print(combined_plot)
dev.off()