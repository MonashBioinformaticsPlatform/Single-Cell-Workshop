
library(stringr)

build_html <- function() {
    rmarkdown::render("vignettes/installation.Rmd", output_dir="docs")
    rmarkdown::render("vignettes/solutions.Rmd", output_dir="docs")
    rmarkdown::render("vignettes/pbmc3k_tutorial.Rmd", output_dir="docs")
}

build_r <- function() {
    lines <- readLines("vignettes/pbmc3k_tutorial.Rmd")
    lines <- str_trim(lines, "right")

    state <- "start1"
    need_blank <- FALSE
    output <- c()
    show <- function(...) output[length(output)+1] <<- paste0(...)
    for(line in lines) {
        if (state == "start1" && str_detect(line,"^---")) {
            state <- "start2"
        } else if (state == "start2" && str_detect(line,"^---")) {
            state <- "main"
            need_blank <- FALSE
        } else if (state == "main" && str_detect(line,"^```.*include=FALSE")) {
            state <- "hidden_code"
        } else if (state == "main" && str_detect(line,"^```")) {
            state <- "code"
            show("")
        } else if (state %in% c("code","hidden_code") && str_detect(line,"^```")) {
            state <- "main"
            need_blank <- TRUE
        } else if (state == "main") {
            clean_line <- str_replace_all(line,"`|<details>|</details>|<summary>|</summary>|\\{\\..*\\}","")
            clean_line <- str_trim(clean_line, "right")
            if (clean_line=="" || clean_line=="\\" || clean_line=="***") {
                need_blank <- TRUE
            } else if (str_detect(clean_line,"^#+\\s+\\w+")) {
                show("")
                show("")
                show(clean_line, " --------")
                need_blank <- TRUE
            } else {
                if (need_blank)
                    show("")
                show("# ", clean_line)
                need_blank <- FALSE
            }
        } else if (state == "code") {
            show(line)
        }
    }

    writeLines(output, "vignettes/pbmc3k_tutorial.R")
}

if (!interactive()) {
    build_r()
    build_html()
    cat("Done.\n")
}

invisible()
