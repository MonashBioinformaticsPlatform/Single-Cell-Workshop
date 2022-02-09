
build_html <- function() {
    rmarkdown::render("vignettes/pbmc3k_tutorial.Rmd", output_dir="docs")
}

if (!interactive()) {
    build_html()
    print(warnings())
    cat("Done.\n")
}

invisible()
