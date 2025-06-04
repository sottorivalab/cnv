#!/usr/bin/env Rscript
# Load libs:
if (!require(sequenza)) stop("Package 'sequenza' missing\n.")

args <- commandArgs(TRUE)
print(args)
input <- args[1]
output_dir <- args[2]
output_prefix <- args[3]
gender <- args[4]
ploidy <- as.integer(args[5])
ccf_val <- as.numeric(args[6])
gam <- as.integer(args[7])

if (ploidy == 7) {
  low_p <- 1
  up_p <- 7
} else  {
  low_p <- ploidy - 1
  up_p <- ploidy + 1
}

if (ccf == 100) {
  high_ccf <- 1.0
  low_ccf <- 0.1
} else {
  high_ccf <- ccf/100 + 0.1
  low_ccf <- ccf/100 - 0.1
}
params_list <- list("input" = input, "output_prefix" = output_prefix)
# Function:
sequenzaAnalysis <- function(input,
                             output_prefix,
                             window = 1e5,
                             overlap = 1,
                             gamma = gam,
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
                             low_cell = low_ccf,
                             up_cell = high_ccf,
                             low_ploidy = low_p,
                             up_ploidy = up_p,
                             CNt_max = 20) {

    #  Define chromosomes to analyse (note these will subset to those that
    # are available for sequenza:
    chr_list <- c(1:22, "X")
    if (gender != "XX") {
        chr_list <- c(chr_list, "Y")
    }

    chr_list <- c(chr_list, paste0("chr", chr_list))

    # Extract sequenza data for model:
    cat(sprintf("- Starting analysis for %s\n", input))
    cat("- Calculating gc-stats\n")
    gc_stats <- gc.sample.stats(input)

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
        parallel = 8,
        gc.stats = gc_stats
    )

    # Fit the model:
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

    # Export the data:
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

