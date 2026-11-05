#!/usr/bin/env Rscript

render_qcplots <- function(qc_dir, nClusters, GMKey, has_qc_report){
    workd <- getwd()

    output.path <- file.path(workd, snakemake@output[["dashboard"]])
    script.in <- file.path(snakemake@scriptdir, "qc_dashboard_interactive.Rmd")

    # Render straight to the Snakemake output, with knit intermediates in a
    # private directory. Rendering beside the Rmd wrote into the workflow's
    # own scripts directory, so projects run from one clone collided and a
    # read-only install failed.
    intermediates <- tempfile("qc_dashboard_")
    dir.create(intermediates)
    on.exit(unlink(intermediates, recursive = TRUE), add = TRUE)

    rmarkdown::render(script.in,
                    params = list(qc_dir = qc_dir,
                                  nClusters = nClusters,
                                  GMKey = GMKey,
                                  has_qc_report = has_qc_report),
                    knit_root_dir = workd,
                    output_dir = dirname(output.path),
                    output_file = basename(output.path),
                    intermediates_dir = intermediates)
}

render_qcplots(snakemake@params[["qc_dir"]],
               snakemake@params[["clusters"]],
               snakemake@params[["google_api_key"]],
               snakemake@params[["has_qc_report"]])
