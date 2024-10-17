# acd7

## Refactoring of the Acadian Variant of FVS

### Introduction

The Acadian Variant of FVS ($ACD_0$) officially exists as R[^readme-1] code and is maintained by [Aaron Weiskittel](mailto:aaron.weiskittel@umaine.edu) and [Ben Rice](mailto:midgard.natural.resources@gmail.com) (among others).

[^readme-1]: R Core Team (2023). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>

While $ACD_0$ performs well and is better than FVSne[^readme-2], there are several things of concern:

[^readme-2]: FVS Staff. 2008 (revised April 16, 2024). Northeast (NE) Variant Overview – Forest Vegetation Simulator. Internal Rep. Fort Collins, CO: U. S. Department of Agriculture, Forest Service, Forest Management Service Center. 56p

-   Use of a plot-level basal area increment equation to limit basal area growth emerging from tree-level diameter growth estimates;
-   Over-estimation of survival for over-stocked stands;
-   Use of a cap on height growth to regulate unreasonably high height growth rates;
-   Some problematic formulations of component models;
-   A number of bugs and inconsistencies in the R code; and
-   Slow execution times.

### Methodology: Developing $ACD_7$

-   Re-coded the model in C++ (to the C++ 23 standard) and built an R package interface;
-   Refactoring:
    -   Abandoned mixed-model approach and fit species independently
    -   Fit equations for:
        -   [diameter growth](./documentation/Diameter%20Growth%20Equation.pdf);
        -   [individual tree survival](./documentation/Species%20Specific%20Survival%20Equation.pdf);
        -   [height growth](./documentation/Height%20Growth%20Equation.pdf);
        -   [height to crown base change](./documentation/Height%20to%20Crown%20Base%20Change%20Equation.pdf);
        -   [height prediction](./documentation/Height-Imputation-for-the-Acadian-Variant-of-FVS--ACD-.pdf);
-   Developed an [algorithm](./documentation/Tree-List-Expansion-and-Contraction-Algorithm.pdf) to solve high expansion factor bias;
-   Worked with Ben Rice's FIA plot benchmarking data set to [validate](./documentation/Acadian-Benchmarking-v7.pdf) and [test](./documentation/Acadian-Benchmarking-100-Year-Projections-v7.pdf) the model.

### R Package

$ACD_7$ is accessible through an R package `acdR` available in this repository as both a [Windows binary](./acdR_0.76.zip) and
[Source tarball](./acdR_0.76.tar.gz). 
The package depends on the following libraries:

* Rcpp (>= 1.0.9)
* RSQLite
* dplyr
 
### Repo Contents

**acd7**\
|\
|-- **acdR**: R package source\
|\
|-- **cl_acd**: command line interface source \
|\
|-- **db_acd**: SQLite interface source and example files\
|\
|-- **documentation**: pdfs of fitting documentation for component equations\

Note: all C++ source code tested with g++ 13.2.0.

Examples for running the R package and the command line interfaces can be found in `README.txt`.

### Repo Development Path

1. Insure cross-platform build integrity
2. Build in runtime optimizations
3. Chasing model residuals


