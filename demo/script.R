#!/usr/bin/env Rscript

options(warn = 0)
required_packages <- c("data.table", "getopt", "jsonlite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". Create the environment from environment.yml.")
}
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(getopt))

cat_header <- function() {
  cat("###################################################\n")
  cat("## Generating NSS for UKBioCoin (UKC)\n")
  cat("## Version: 2.0\n")
  cat("## Written by: Jing-cheng He, Guo-An Qi, Zhejiang University\n")
  cat("## Bug report: jc_he@zju.edu.cn\n")
  cat("###################################################\n\n")
}

command <- matrix(c(
  "pfile", "pf", 2, "character", "PLINK2 genotype prefix (use this or --bfile)",
  "bfile", "bf", 2, "character", "PLINK1 binary genotype prefix (use this or --pfile)",
  "pheno", "p", 2, "character", "Phenotype table with FID, IID and phenotype columns",
  "novisualize", "v", 0, "logical", "Skip comparison plots",
  "threads", "t", 2, "integer", "PLINK/data.table threads (default: 4)",
  "memory", "m", 2, "integer", "PLINK memory limit in MB (default: 8000)",
  "prefix", "pr", 2, "character", "NSS prefix (default: UKC)",
  "ukc", "ue", 2, "character", "Path to the UKBioCoin executable",
  "out", "o", 2, "character", "Output directory (default: current directory)",
  "nss-format", "nf", 2, "character", "NSS output format: v2 (default) or legacy",
  "force", "f", 0, "logical", "Replace an existing NSS output",
  "help", "h", 0, "logical", "Show this help"
), byrow = TRUE, ncol = 5)
args <- getopt(command)

if (!is.null(args$help) || (is.null(args$pfile) && is.null(args$bfile)) ||
    is.null(args$pheno)) {
  cat_header()
  cat(getopt(command, usage = TRUE), "\n")
  quit(status = ifelse(is.null(args$help), 1, 0))
}
if (!is.null(args$pfile) && !is.null(args$bfile)) {
  stop("Specify only one of --pfile and --bfile")
}

threads <- if (is.null(args$threads)) 4L else as.integer(args$threads)
memory <- if (is.null(args$memory)) 8000L else as.integer(args$memory)
prefix <- if (is.null(args$prefix)) "UKC" else args$prefix
ukc <- if (is.null(args$ukc)) "UKBioCoin" else args$ukc
out <- if (is.null(args$out)) getwd() else normalizePath(args$out, mustWork = FALSE)
nss_format <- if (is.null(args$`nss-format`)) "v2" else tolower(args$`nss-format`)
novisualize <- !is.null(args$novisualize)
force <- !is.null(args$force)
if (!nss_format %in% c("v2", "legacy")) stop("--nss-format must be v2 or legacy")
if (threads < 1) stop("--threads must be a positive integer")
if (memory < 1) stop("--memory must be a positive integer")

pheno_input <- normalizePath(args$pheno)
if (is.null(args$pfile)) {
  genotype_prefix <- sub("\\.fam$", "", normalizePath(paste0(args$bfile, ".fam")))
  plink_input <- c("--bfile", genotype_prefix)
  sample_path <- paste0(genotype_prefix, ".fam")
  sample_file <- fread(sample_path, header = FALSE, data.table = FALSE,
                       nThread = threads)
} else {
  genotype_prefix <- sub("\\.psam$", "", normalizePath(paste0(args$pfile, ".psam")))
  plink_input <- c("--pfile", genotype_prefix)
  sample_path <- paste0(genotype_prefix, ".psam")
  sample_file <- fread(sample_path, data.table = FALSE, nThread = threads)
}
sample_count <- nrow(sample_file)
if (sample_count < 1) stop("No samples found in ", sample_path)

cat_header()
cat("Options:\n",
    "  genotype: ", genotype_prefix, "\n",
    "  phenotype: ", pheno_input, "\n",
    "  sample_count: ", sample_count, "\n",
    "  nss-format: ", nss_format, "\n",
    "  threads: ", threads, "\n",
    "  memory: ", memory, " MB\n",
    "  prefix: ", prefix, "\n",
    "  output: ", out, "\n\n", sep = "")

dir.create(out, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(out)
dir.create("0.basic", showWarnings = FALSE)
dir.create("1.plink_temp", showWarnings = FALSE)

run_command <- function(arguments) {
  display <- paste(shQuote(arguments), collapse = " ")
  cat(display, "\n")
  status <- system2(arguments[1], arguments[-1])
  if (status != 0) stop("Command failed with status ", status, ": ", display)
}

total_minutes <- 0
stage_start <- proc.time()[["elapsed"]]
run_command(c("plink2", plink_input, "--silent", "--geno-counts",
              "--out", "0.basic/0.count", "--threads", threads,
              "--memory", memory))
run_command(c("plink2", plink_input, "--silent", "--freq",
              "--out", "0.basic/0.freq", "--threads", threads,
              "--memory", memory))
minutes <- (proc.time()[["elapsed"]] - stage_start) / 60
total_minutes <- total_minutes + minutes
cat(sprintf("Generating basic statistics done. Elapsed time: %.3f minute.\n\n", minutes))

phenotype_table <- fread(pheno_input, data.table = FALSE, nThread = threads)
if (ncol(phenotype_table) < 3 || !all(c("FID", "IID") %in% names(phenotype_table))) {
  stop("Phenotype table must contain FID, IID and at least one phenotype")
}
sample_iids <- as.character(sample_file[[2]])
phenotype_table <- phenotype_table[as.character(phenotype_table$IID) %in% sample_iids, , drop = FALSE]
if (nrow(phenotype_table) == 0) stop("No phenotype rows match genotype sample IIDs")
phenotypes <- setdiff(names(phenotype_table), c("FID", "IID"))
keep <- vapply(phenotype_table[phenotypes], function(value) mean(is.na(value)) < 0.8,
               logical(1))
if (any(!keep)) {
  cat("Dropping phenotypes with missing rate >=80%: ",
      paste(phenotypes[!keep], collapse = ", "), "\n", sep = "")
}
phenotypes <- phenotypes[keep]
if (length(phenotypes) == 0) stop("No phenotype passed the missingness filter")
phenotype_table <- phenotype_table[c("FID", "IID", phenotypes)]
for (name in phenotypes) phenotype_table[[name]] <- as.numeric(scale(phenotype_table[[name]]))
scaled_path <- file.path(out, "scaled_pheno.table")
fwrite(phenotype_table, scaled_path, sep = "\t", quote = FALSE, na = "NA")

stage_start <- proc.time()[["elapsed"]]
run_command(c("plink2", plink_input, "--glm", "allow-no-covars",
              "skip-invalid-pheno", "--no-input-missing-phenotype",
              "--out", "1.plink_temp/single_reg", "--pheno", scaled_path,
              "--pheno-name", paste(phenotypes, collapse = ","),
              "--threads", threads, "--memory", memory,
              "--read-freq", "0.basic/0.freq.afreq"))
minutes <- (proc.time()[["elapsed"]] - stage_start) / 60
total_minutes <- total_minutes + minutes
cat(sprintf("Generating PLINK statistics done. Elapsed time: %.3f minute.\n\n", minutes))

write_npy <- function(path, values, dtype, shape) {
  if (!dtype %in% c("<f4", "<f8")) stop("Unsupported NPY dtype: ", dtype)
  if (prod(shape) != length(values)) stop("NPY shape/value length mismatch for ", path)
  shape_text <- if (length(shape) == 1) {
    paste0("(", shape[[1]], ",)")
  } else {
    paste0("(", paste(shape, collapse = ", "), ")")
  }
  dictionary <- paste0("{'descr': '", dtype,
                       "', 'fortran_order': False, 'shape': ", shape_text, ", }")
  padding <- (16 - ((10 + nchar(dictionary, type = "bytes") + 1) %% 16)) %% 16
  header <- paste0(dictionary, paste(rep(" ", padding), collapse = ""), "\n")
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(as.raw(c(0x93, charToRaw("NUMPY"), 1, 0)), connection)
  writeBin(as.integer(nchar(header, type = "bytes")), connection,
           size = 2, endian = "little")
  writeBin(charToRaw(header), connection)
  if (dtype == "<f4") {
    writeBin(as.numeric(values), connection, size = 4, endian = "little")
  } else {
    writeBin(as.double(values), connection, size = 8, endian = "little")
  }
}

validate_npy <- function(path, dtype, shape) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  magic <- readBin(connection, "raw", n = 6)
  if (!identical(magic, as.raw(c(0x93, charToRaw("NUMPY"))))) {
    stop("Invalid NPY magic in ", path)
  }
  version <- readBin(connection, "raw", n = 2)
  if (length(version) != 2 || version[[1]] != as.raw(1)) {
    stop("Unsupported NPY version in ", path)
  }
  header_length <- readBin(connection, integer(), n = 1, size = 2,
                           endian = "little", signed = FALSE)
  header_raw <- readBin(connection, "raw", n = header_length)
  if (length(header_raw) != header_length) stop("Truncated NPY header in ", path)
  header <- rawToChar(header_raw)
  expected_shape <- if (length(shape) == 1) {
    paste0("(", shape[[1]], ",)")
  } else {
    paste0("(", paste(shape, collapse = ", "), ")")
  }
  if (!grepl(paste0("'descr': '", dtype, "'"), header, fixed = TRUE) ||
      !grepl("'fortran_order': False", header, fixed = TRUE) ||
      !grepl(paste0("'shape': ", expected_shape), header, fixed = TRUE)) {
    stop("NPY header does not match expected dtype/shape in ", path)
  }
  item_size <- ifelse(dtype == "<f4", 4, 8)
  expected_bytes <- 10 + header_length + prod(shape) * item_size
  actual_bytes <- file.info(path)$size
  if (actual_bytes != expected_bytes) {
    stop("NPY payload size mismatch in ", path, ": expected ", expected_bytes,
         ", found ", actual_bytes)
  }
}

npy_descriptor <- function(path, dtype, shape) {
  list(path = path, dtype = dtype, shape = as.list(as.integer(shape)), order = "C")
}

indexed_cov_path <- function(index) sprintf("cov_xy/%08d.npy", index - 1L)

stage_start <- proc.time()[["elapsed"]]
gcount <- fread("0.basic/0.count.gcount", nThread = threads)
required_count <- c("ID", "HOM_REF_CT", "HET_REF_ALT_CTS", "TWO_ALT_GENO_CTS",
                    "MISSING_CT")
if (!all(required_count %in% names(gcount))) {
  stop("PLINK genotype-count output lacks required columns: ",
       paste(setdiff(required_count, names(gcount)), collapse = ", "))
}
called_count <- rowSums(gcount[, c("HOM_REF_CT", "HET_REF_ALT_CTS",
                                   "TWO_ALT_GENO_CTS")])
var_x <- (4 * gcount$HOM_REF_CT + gcount$HET_REF_ALT_CTS) / called_count -
  (2 * gcount$HOM_REF_CT / called_count + gcount$HET_REF_ALT_CTS / called_count)^2
variant_ids <- as.character(gcount$ID)
variant_count <- length(var_x)
if (variant_count == 0 || anyDuplicated(variant_ids)) stop("Invalid/duplicate variant IDs in gcount")

afreq <- fread("0.basic/0.freq.afreq", nThread = threads)
if (!identical(as.character(afreq$ID), variant_ids)) {
  stop("Variant ID/order mismatch between frequency and genotype-count outputs")
}
if (is.null(args$bfile)) {
  pvar <- fread(paste0(genotype_prefix, ".pvar"), nThread = threads)
} else {
  pvar <- fread(paste0(genotype_prefix, ".bim"), header = FALSE, nThread = threads)
  setnames(pvar, c("#CHROM", "ID", "GENETICPOS", "POS", "ALT", "REF"))
}
if (!identical(as.character(pvar$ID), variant_ids)) {
  stop("Variant ID/order mismatch between pvar/bim and genotype-count outputs")
}

scaled_phenotypes <- fread(scaled_path, data.table = FALSE,
                           nThread = threads)[phenotypes]
cov_yy <- cov(scaled_phenotypes, use = "pairwise")
y_missing <- vapply(scaled_phenotypes, function(value) mean(is.na(value)), numeric(1))
x_missing <- gcount$MISSING_CT / (gcount$MISSING_CT + called_count)

if (!all(c("REF", "ALT", "ALT_FREQS") %in% names(afreq))) {
  stop("PLINK frequency output lacks REF, ALT, or ALT_FREQS")
}
afreq$REF_FREQ <- 1 - afreq$ALT_FREQS
afreq$POS <- pvar$POS
afreq$Effect_Allele <- ifelse(afreq$REF_FREQ > 0.5, afreq$ALT, afreq$REF)
metadata <- afreq[, .(`#CHROM`, ID, POS, Effect_Allele,
                      REF_Allele = REF, ALT_Allele = ALT, REF_FREQ)]
if (!identical(as.character(metadata$ID), variant_ids)) stop("Internal metadata order mismatch")

if (nss_format == "v2") {
  dir.create("2.nss", showWarnings = FALSE)
  target <- file.path("2.nss", paste0(prefix, ".nss"))
  temporary <- paste0(target, ".tmp.", Sys.getpid())
  if (dir.exists(target) && !force) stop(target, " already exists; use --force to replace it")
  if (dir.exists(temporary)) unlink(temporary, recursive = TRUE)
  dir.create(file.path(temporary, "cov_xy"), recursive = TRUE)
  completed <- FALSE
  on.exit(if (!completed && dir.exists(temporary)) unlink(temporary, recursive = TRUE), add = TRUE)

  for (i in seq_along(phenotypes)) {
    name <- phenotypes[[i]]
    glm_path <- file.path("1.plink_temp", paste0("single_reg.", name, ".glm.linear"))
    if (!file.exists(glm_path)) {
      warning("Missing PLINK GLM for ", name, "; writing an all-NaN COVxy vector")
      cov_xy <- rep(NaN, variant_count)
    } else {
      glm <- fread(glm_path, nThread = threads)
      if ("TEST" %in% names(glm)) glm <- glm[TEST == "ADD"]
      if (!identical(as.character(glm$ID), variant_ids)) {
        stop("Variant ID/order mismatch in ", glm_path)
      }
      if (!"BETA" %in% names(glm)) stop(glm_path, " lacks BETA")
      cov_xy <- var_x * glm$BETA
    }
    write_npy(file.path(temporary, indexed_cov_path(i)), cov_xy, "<f4", variant_count)
    cat(sprintf("  COVxy vectors: %d/%d (%d%%)\n", i, length(phenotypes),
                floor(100 * i / length(phenotypes))))
  }
  fwrite(metadata, file.path(temporary, "meta.tsv"), sep = "\t",
         quote = FALSE, na = "NA")
  write_npy(file.path(temporary, "cov_yy.npy"), as.vector(t(cov_yy)), "<f8",
            c(length(phenotypes), length(phenotypes)))
  write_npy(file.path(temporary, "var_x.npy"), var_x, "<f8", variant_count)
  write_npy(file.path(temporary, "x_missing.npy"), x_missing, "<f8", variant_count)
  write_npy(file.path(temporary, "y_missing.npy"), y_missing, "<f8", length(phenotypes))

  manifest <- list(
    schema = "org.ukbiocoin.nss",
    schema_version = 2L,
    ukc_version = "2.0",
    byte_order = "little",
    dimensions = list(sample_count = sample_count,
                      variant_count = variant_count,
                      phenotype_count = length(phenotypes)),
    metadata = list(path = "meta.tsv", format = "tsv", rows = variant_count),
    arrays = list(
      cov_yy = npy_descriptor("cov_yy.npy", "<f8",
                              c(length(phenotypes), length(phenotypes))),
      var_x = npy_descriptor("var_x.npy", "<f8", variant_count),
      x_missing = npy_descriptor("x_missing.npy", "<f8", variant_count),
      y_missing = npy_descriptor("y_missing.npy", "<f8", length(phenotypes))
    ),
    phenotypes = lapply(seq_along(phenotypes), function(i) {
      list(index = i - 1L, name = phenotypes[[i]],
           cov_xy = npy_descriptor(indexed_cov_path(i), "<f4", variant_count))
    })
  )
  jsonlite::write_json(manifest, file.path(temporary, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA,
                       null = "null")
  validate_npy(file.path(temporary, "cov_yy.npy"), "<f8",
               c(length(phenotypes), length(phenotypes)))
  validate_npy(file.path(temporary, "var_x.npy"), "<f8", variant_count)
  validate_npy(file.path(temporary, "x_missing.npy"), "<f8", variant_count)
  validate_npy(file.path(temporary, "y_missing.npy"), "<f8", length(phenotypes))
  for (i in seq_along(phenotypes)) {
    validate_npy(file.path(temporary, indexed_cov_path(i)), "<f4", variant_count)
  }
  checked_manifest <- jsonlite::read_json(file.path(temporary, "manifest.json"),
                                          simplifyVector = TRUE)
  checked_dimensions <- unlist(checked_manifest$dimensions, use.names = TRUE)
  if (!identical(as.integer(checked_dimensions[c("sample_count", "variant_count",
                                                  "phenotype_count")]),
                 as.integer(c(sample_count, variant_count, length(phenotypes))))) {
    stop("Generated manifest dimensions failed self-validation")
  }
  if (nrow(fread(file.path(temporary, "meta.tsv"), nThread = threads)) != variant_count) {
    stop("Generated metadata row count failed self-validation")
  }
  backup <- paste0(target, ".backup.", Sys.getpid())
  moved_existing <- FALSE
  if (dir.exists(backup)) unlink(backup, recursive = TRUE)
  if (dir.exists(target)) {
    if (!file.rename(target, backup)) stop("Unable to move existing output to backup")
    moved_existing <- TRUE
  }
  if (!file.rename(temporary, target)) {
    if (moved_existing) file.rename(backup, target)
    stop("Unable to atomically rename ", temporary, " to ", target)
  }
  if (moved_existing) unlink(backup, recursive = TRUE)
  completed <- TRUE
  nss_input <- target
} else {
  dir.create("2.matrix", showWarnings = FALSE)
  target_prefix <- file.path("2.matrix", prefix)
  target_files <- paste0(target_prefix, c("_cov_xy.table", "_cov_yy.table",
                                          "_var_x.table", "_meta.table",
                                          "_x_missing.table", "_y_missing.table"))
  if (any(file.exists(target_files)) && !force) {
    stop("Legacy NSS output exists; use --force to replace it")
  }
  cov_xy_columns <- vector("list", length(phenotypes))
  for (i in seq_along(phenotypes)) {
    name <- phenotypes[[i]]
    glm_path <- file.path("1.plink_temp", paste0("single_reg.", name, ".glm.linear"))
    if (!file.exists(glm_path)) {
      warning("Missing PLINK GLM for ", name, "; writing an all-NA COVxy column")
      cov_xy_columns[[i]] <- rep(NA_real_, variant_count)
    } else {
      glm <- fread(glm_path, nThread = threads)
      if ("TEST" %in% names(glm)) glm <- glm[TEST == "ADD"]
      if (!identical(as.character(glm$ID), variant_ids)) stop("Variant ID/order mismatch in ", glm_path)
      cov_xy_columns[[i]] <- var_x * glm$BETA
    }
  }
  names(cov_xy_columns) <- phenotypes
  cov_path <- paste0(target_prefix, "_cov_xy.table")
  writeLines(paste(sprintf('"%s"', phenotypes), collapse = " "), cov_path)
  cov_table <- as.data.table(cov_xy_columns)
  cov_table <- cbind(row_label = sprintf('"%d"', seq_len(variant_count)), cov_table)
  fwrite(cov_table, cov_path, sep = " ", quote = FALSE, na = "NA",
         col.names = FALSE, append = TRUE)
  write.table(cov_yy, paste0(target_prefix, "_cov_yy.table"), quote = TRUE,
              row.names = TRUE, col.names = TRUE, sep = " ")
  write.table(var_x, paste0(target_prefix, "_var_x.table"), row.names = TRUE,
              col.names = "var_x")
  fwrite(metadata, paste0(target_prefix, "_meta.table"), sep = " ",
         quote = FALSE, na = "NA")
  write.table(x_missing, paste0(target_prefix, "_x_missing.table"),
              row.names = TRUE, col.names = "x_missing")
  write.table(y_missing, paste0(target_prefix, "_y_missing.table"),
              row.names = TRUE, col.names = "missing_rate")
  nss_input <- target_prefix
}

minutes <- (proc.time()[["elapsed"]] - stage_start) / 60
total_minutes <- total_minutes + minutes
cat(sprintf("Generating NSS %s done. Elapsed time: %.3f minute.\n\n",
            nss_format, minutes))

if (!novisualize) {
  dir.create("3.analysis", showWarnings = FALSE)
  stage_start <- proc.time()[["elapsed"]]
  result_prefix <- file.path("3.analysis", paste0("test.", phenotypes[[1]]))
  run_command(c(ukc, "--file", nss_input, "--phe", phenotypes[[1]],
                "--use-missing-rate-estimate", "--threads", threads,
                "--out", result_prefix,
                if (nss_format == "legacy") c("--totalsize", sample_count) else NULL))
  result <- fread(paste0(result_prefix, "_results.table"), nThread = threads)
  glm <- fread(file.path("1.plink_temp", paste0("single_reg.", phenotypes[[1]],
                                                ".glm.linear")), nThread = threads)
  if ("TEST" %in% names(glm)) glm <- glm[TEST == "ADD"]
  result$BETA <- ifelse(result$Effect_Allele == glm$A1, result$BETA, -result$BETA)
  comparison <- merge(
    result[, .(ID, UKC_BETA = BETA, UKC_SE = SE, UKC_TSTAT = `T-STAT`,
               UKC_P = 10^(-`-log10_P`))],
    glm[, .(ID, UKB_BETA = BETA, UKB_SE = SE, UKB_TSTAT = T_STAT, UKB_P = P)],
    by = "ID")
  comparison <- na.omit(comparison)
  if (nrow(comparison) > 10000) comparison <- comparison[sample(.N, 10000)]
  png(file.path("3.analysis", paste0("Validation.", phenotypes[[1]], ".png")),
      width = 4800, height = 4800, res = 300)
  par(mfrow = c(2, 2))
  plot(comparison$UKB_BETA, comparison$UKC_BETA, xlab = "PLINK BETA", ylab = "UKC BETA"); abline(0, 1, col = "red")
  plot(comparison$UKB_SE, comparison$UKC_SE, xlab = "PLINK SE", ylab = "UKC SE"); abline(0, 1, col = "red")
  plot(comparison$UKB_TSTAT, comparison$UKC_TSTAT, xlab = "PLINK T", ylab = "UKC T"); abline(0, 1, col = "red")
  plot(-log10(comparison$UKB_P), -log10(comparison$UKC_P), xlab = "PLINK -log10(P)", ylab = "UKC -log10(P)"); abline(0, 1, col = "red")
  dev.off()
  minutes <- (proc.time()[["elapsed"]] - stage_start) / 60
  total_minutes <- total_minutes + minutes
  cat(sprintf("Validation plot done. Elapsed time: %.3f minute.\n\n", minutes))
}

cat(sprintf("ALL UKC NSS generation done. Elapsed time: %.3f minute.\n", total_minutes))
