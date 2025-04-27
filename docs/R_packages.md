# R Packages overview code

This R script installs and loads a set of specified R packages, retrieves their citation information, and generates a tidytable of package versions and citations formatted in Markdown, which can be directly used in a GitHub README.

## Install Required Packages

The script starts by installing and loading the required R packages. If any package is missing, it will be automatically installed. 

```r
packages <- c(
    "broom", "broom.mixed", "cowplot", "dbplyr", "doSNOW", "extrafont", 
    "flextable", "foreach", "ftExtra", "ggfortify", "ggtext", "grid", 
    "kableExtra", "lmerTest", "MuMIn", "officer", "parallel", "patchwork", 
    "purrr", "purrrlyr", "ramboseli", "rptR", "RPostgreSQL", "scales", 
    "showtext", "survminer", "survMisc", "survival", "systemfonts", 
    "tidyverse", "zoo"
)

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
}

devtools::install_github("amboseli/ramboseli")

lapply(packages, require, character.only = TRUE, warn.conflicts = TRUE, quietly = FALSE)
```

```r
get_citation_text <- function(library) {
    library %>%
        citation() %>%
        capture.output() %>%
        as_tibble(.name_repair = "unique") %>%
        filter(!str_starts(value, "To cite"),
               !str_starts(value, "A BibTeX entry for LaTeX users is"),
               !str_starts(value, "To see these entries in BibTeX format"),
               !str_starts(value, "bibtex=TRUE")
               ) %>%
        pull(value) %>%
        str_c(collapse = " ") %>%
        str_remove("@Manual\\{.*") %>%
        str_remove("@Article\\{.*") %>%
        str_remove("R package version [^,]+,\\s*") %>%
        str_remove("https?://CRAN\\.R-project\\.org/package=[^\\s]+") %>%
        str_remove_all("‘.*options.*’") %>%
        str_remove_all("options\\(citation\\.bibtex\\.max=999\\)") %>%
        str_replace_all("_", "") %>%
        str_squish() %>%
        str_remove("<")  # remove any lingering
}
```


```r
ciation_data <- tibble(library = packages) %>%
    mutate(version = map(.x = library, .f = packageVersion)) %>%
    mutate(reference = map_chr(.x = library, .f = get_citation_text)) %>%
    unnest(cols = c(version)) %>%
    mutate(version = as.character(version))
    
    # Convert the tibble to a Markdown table
markdown_table <- citation_data %>%
    mutate(Markdown = str_glue("| {library} | {version} | {reference} |")) %>%
    pull(Markdown) %>%
    str_c(collapse = "\n")

# Output the Markdown table for use in README
cat(markdown_table)
    # flextable() %>%
    # flextable::width(j = 3, 8.5) %>%
    # ftExtra::colformat_md()
```
