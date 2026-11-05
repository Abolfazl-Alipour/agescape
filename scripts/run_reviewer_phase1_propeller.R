library(speckle)
library(limma)

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
outfile <- args[2]
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
results <- list()

for (i in seq_len(nrow(manifest))) {
  key <- manifest$key[i]
  counts <- read.csv(file.path(root, key, "counts.csv"), row.names = 1, check.names = FALSE)
  metadata <- read.csv(file.path(root, key, "metadata.csv"), row.names = 1)
  counts <- as.matrix(counts)
  if (ncol(counts) == 0 || nrow(counts) < 5) {
    cat("Skipping", manifest$species[i], manifest$organ[i], "because no adequately covered cell types remain\n")
    next
  }
  metadata <- metadata[rownames(counts), , drop = FALSE]
  props.pseudo <- (counts + 0.5) / rowSums(counts + 0.5)
  transformed <- t(log(props.pseudo / (1 - props.pseudo)))
  design <- data.frame(age = metadata$age, log_total = log(rowSums(counts)))
  if (length(unique(metadata$sex_male)) > 1) {
    design$sex_male <- metadata$sex_male
  }
  if ("platform_nextseq" %in% colnames(metadata) && length(unique(metadata$platform_nextseq)) > 1) {
    design$platform_nextseq <- metadata$platform_nextseq
  }
  design <- model.matrix(~ ., data = design)
  fit <- lmFit(transformed, design)
  fit <- eBayes(fit[, "age"], robust = TRUE, trend = FALSE)
  df <- fit$df.residual + fit$df.prior
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, 1]
  estimate <- fit$coefficients[, 1]
  critical <- qt(0.975, df)
  result <- data.frame(
    species = manifest$species[i], organ = manifest$organ[i], cohort = manifest$cohort[i],
    cell_type = rownames(fit$coefficients), effect_logit = estimate, se = se,
    ci95_low = estimate - critical * se, ci95_high = estimate + critical * se,
    pvalue = fit$p.value[, 1], q_bh = p.adjust(fit$p.value[, 1], method = "BH"),
    df_total = df, row.names = NULL
  )
  results[[key]] <- result
  cat(manifest$species[i], manifest$organ[i], manifest$cohort[i], "\n")
}

output <- do.call(rbind, results)
write.csv(output, outfile, row.names = FALSE)
