library(glue)

setwd('simulation')

ncores = 1

# simulation settings
all_traits = c("continuous")
all_sig_mafs = c(0.1)
all_h2s = c(0.3)
all_mcars = c(0.1)
all_ns = c(100000)
B_start = 1
B_end = 100
prefix = "beta_compare"

## --- make root dir for scripts -----------------------------------------
qsub_root <- file.path("qsub_scripts", prefix)
dir.create(qsub_root, recursive = TRUE, showWarnings = FALSE)

## --- make root dir for results -----------------------------------------
simu_res_root <- file.path("simulation_results", prefix)
dir.create(simu_res_root, recursive = TRUE, showWarnings = FALSE)


## --- iterate over the factorial design ---------------------------------
for (trait   in all_traits)
for (maf_sig in all_sig_mafs)
for (h2      in all_h2s)
for (mcar    in all_mcars)
for (n       in all_ns)
 {

    task_name <- glue(
    "{trait}_maf{maf_sig}_h2{h2}_mcar{mcar}_n{n}"
    )

    ## filenames
    qsub_file  <- glue("{task_name}.pbs")
    log_file   <- glue("{task_name}.log")
    err_file   <- glue("{task_name}.err")
    rout_file  <- glue("{task_name}.Rout")
    qsub_path  <- file.path(qsub_root, qsub_file)

    ## PBS script -----------------------------------------------------------
    qsub_script <- glue(
    '#!/bin/bash
    #PBS -N {task_name}
    #PBS -o {qsub_root}/{log_file}
    #PBS -e {qsub_root}/{err_file}
    #PBS -l nodes=1:ppn={ncores}
    #PBS -l walltime=1000:00:00
    #PBS -m abe

    source ~/.bashrc
    cd {getwd()}/simulation_results
    mamba activate r4.3

    Rscript ../beta_sub_task.R {trait} {maf_sig} {h2} {mcar} {B_start} {B_end} {n} {prefix} > ../{qsub_root}/{rout_file} 2>&1
    ')

    writeLines(qsub_script, qsub_path)
      system(glue("qsub {qsub_path}"))
}