
Templates with ![](https://img.shields.io/badge/status-stable-green) revision indicates that the components or processes have undergone comprehensive parameterization and testing.

Templates with ![](https://img.shields.io/badge/status-alpha-yellow) revision indicates that the components or processes are currently being tested. There is some test data available, but there are parameters that need to be set up manually within the code.

Templates with ![](https://img.shields.io/badge/status-draft-grey) revision indicates that the components or processes are not fully tested. There is no test data available, parameters need to be set up manually within the code, and specific code changes are required based on the data used.

# Guidelines for downstream analysis

- Set the working directory to the directory containing this README. We recommend using a [Project](https://support.posit.co/hc/en-us/articles/200526207-Using-RStudio-Projects) in Rstudio.
- Use [install_dependencies.R](install_dependencies.R) to install all packages used in these reports.

## nf-core pipelines

These templates assume that the raw data has been processed by 

- [nf-core/chipseq](https://nf-co.re/chipseq/2.1.0/docs/usage)
- EPIC methylation array

## Quick Start

### With Rstudio

```
source(install_depedencies.R)
rmarkdown::render("QC.Rmd", quiet = TRUE
```

## Downstream analysis

Before using any template:
1. **Modify** [information.R](information.R) with the right information. You can use this file with any template to include the project/analysis information.
2. **Modify** *params.R* with the locations of select files/folders from the output of [nf-core/chipseq](https://nf-co.re/chipseq/2.1.0/docs/output). These nf-core outputs will become inputs to various templates.
3. **Modify** the `YAML` header of the `Rmd` files to choose the right parameters for that report.

Additional useful info:
- `params*example.R` are files containing parameters pointing to a small, simple dataset that can be used to test the report code and see how the fully rendered report looks.

### Quality assessment

![](https://img.shields.io/badge/status-draft-grey) [chipseq/01_quality_assessment/QC.Rmd](chipseq/01_quality_assessment/QC.Rmd)

### Differential expression

![](https://img.shields.io/badge/status-stable-green) [chipseq/02_diffbind/diffbind.Rmd](chipseq/02_diffbind/diffbind.Rmd) 

On the `YAML` header file of the `Rmd` you can specify some parameters or just set them up in the second chunk of code of the template. 

