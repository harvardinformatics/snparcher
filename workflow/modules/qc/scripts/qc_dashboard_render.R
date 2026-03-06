#!/usr/bin/env Rscript

render_qcplots <- function(qc_dir, nClusters, GMKey, has_qc_report){
    workd <- getwd()

    # Output path for the dashboard HTML
    output.path <- snakemake@output[["dashboard"]]

    script.in <- paste0(snakemake@scriptdir, "/qc_dashboard_interactive.Rmd")
    script.out <- gsub(".Rmd", ".html", script.in)

    rmarkdown::render(script.in,
                    params = list(qc_dir = qc_dir,
                                  nClusters = nClusters,
                                  GMKey = GMKey,
                                  has_qc_report = has_qc_report),
                    knit_root_dir = workd)

    copy_successful <- file.copy(script.out, output.path)

    if (copy_successful) {
        file.remove(script.out)
    } else {
        cat("snpArcher: Failed to move the qc dashboard html.\n")
    }
}

render_qcplots(snakemake@params[["qc_dir"]],
               snakemake@params[["clusters"]],
               snakemake@params[["google_api_key"]],
               snakemake@params[["has_qc_report"]])
