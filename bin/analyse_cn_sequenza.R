#!/usr/bin/env Rscript

if (!require(sequenza)) stop("Package 'sequenza' missing\n.")

args <- commandArgs(TRUE)
input <- args[1]
output_dir <- args[2]
output_prefix <- args[3]
gender <- args[4]
ploidy <- as.integer(args[5])
ccf <- as.numeric(args[6])
gam <- as.integer(args[7])

print(paste("full input:", input))
print(paste("output dir:", output_dir))
print(paste("output_prefix:", output_prefix))
print(paste("gender:", gender))
print(paste("ploidy:", ploidy))
print(paste("ccf:", ccf))
print(paste("gamma:", gam))

# Compute ploidy range
if (ploidy == 7) {
  low_p <- 1
  up_p <- 7
} else {
  low_p <- ploidy - 1
  up_p <- ploidy + 1
}

# Compute CCF range
if (ccf == 100) {
  high_ccf <- 1.0
  low_ccf <- 0.1
} else {
  high_ccf <- ccf / 100 + 0.1
  low_ccf <- ccf / 100 - 0.1
}

params_list <- list(
  input = input,
  output_prefix = output_prefix,
  output_dir = output_dir,
  gender = gender,
  gamma = gam,
  low_p = low_p,
  up_p = up_p,
  low_ccf = low_ccf,
  high_ccf = high_ccf
)

sequenzaAnalysis <- function(input,
                             output_prefix,
                             output_dir,
                             gender,
                             gamma,
                             low_p,
                             up_p,
                             low_ccf,
                             high_ccf,
                             window = 1e5,
                             overlap = 1,
                             kmin = 300,
                             min_reads = 40,
                             min_reads_normal = 10,
                             min_reads_baf = 1,
                             max_mut_types = 1,
                             breaks = NULL,
                             assembly = "hg38",
                             weighted_mean = TRUE,
                             normalization_method = "mean",
                             is_female = TRUE,
                             segment_filter = 3e6,
                             ratio_priority = FALSE,
                             method = "baf",
                             CNt_max = 20) {

    cat(sprintf("- Starting analysis for %s\n", input))
    cat("- Calculating gc-stats\n")
    gc_stats <- gc.sample.stats(input)
    chr_list <- gc_stats$file.metrics$chr
    cat("- Loading data\n")
    modDat <- sequenza.extract(input,
        window = window,
        overlap = overlap,
        gamma = gamma,
        kmin = kmin,
        min.reads = min_reads,
        min.reads.normal = min_reads_normal,
        min.reads.baf = min_reads_baf,
        max.mut.types = max_mut_types,
        chromosome.list = chr_list,
        breaks = breaks,
        assembly = assembly,
        weighted.mean = weighted_mean,
        normalization.method = normalization_method,
        gc.stats = gc_stats
    )

    cat("- Fitting the model\n")
    cells <- seq(low_ccf, high_ccf, 0.01)
    plo <- seq(low_p, up_p, 0.1)
    fit <- sequenza.fit(modDat,
        female = is_female,
        segment.filter = segment_filter,
        cellularity = cells,
        ploidy = plo,
        XY = c(X = "chrX", Y = "chrY"),
        ratio.priority = ratio_priority,
        method = method
    )

    cat("- Exporting results\n")
    sequenza.results(modDat,
        fit,
        output_prefix,
        out.dir = output_dir,
        female = is_female,
        CNt.max = CNt_max,
        XY = c(X = "chrX", Y = "chrY"),
        ratio.priority = ratio_priority
    )
}

do.call(sequenzaAnalysis, params_list)
