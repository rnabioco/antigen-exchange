renv::use(lockfile = "../../renv.lock")

library(here)
library(rmarkdown)
library(tidyverse)

slides <- c(
  "1" = "output-XETG00230__0022624__[A-Z]__",
  "2" = "output-XETG00230__0022841__[A-Z]__"
)

slides %>%
  iwalk(~ {
    output <- str_c("xenium-slide-", .y, ".html")

    render(
      input       = "xenium-analysis.Rmd",
      output_file = output,
      params      = list(
        slide_id    = .y,
        slide_regex = .x 
      )
    )
  })

