# flowFun

R package containing functions for analysis of flow cytometry data. This
package can be installed with vignettes from the Rstudio console as
follows:

\`\`\`{r, eval = FALSE} install.packages(“BiocManager”)
install.packages(“devtools”)

# Install with vignettes

# The package ‘flowFunData’ must first be installed

devtools::install_github(“00berst33/flowFunData”) \# Then install
flowFun devtools::install_github(“00berst33/flowFun”, dependencies =
TRUE, build_vignettes = TRUE)

    Or to install without vignettes:

    ```{r, eval = FALSE}
    install.packages("devtools")

    # Install without vignettes
    #   No installation of 'flowFunData', set argument `build_vignettes` to FALSE
    devtools::install_github("00berst33/flowFun", dependencies = TRUE, build_vignettes = FALSE)
