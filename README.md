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
        -   diameter growth;
        -   individual tree survival;
        -   height growth;
        -   height to crown base change;
-   Developed algorithm to solve high expansion factor bias;
-   Worked with Ben Rice's FIA plot benchmarking data set to validate and test the model.

### R Package

$ACD_7$ is accessible through an R package (`acdR`) available in this repository as both a Windows binary and Source tarball. 
The package depends on the following libraries:

* Rcpp (>= 1.0.9)
* RSQLite
* dplyr
 
### Repo Development Path

1. Include C++ and R source directories
2. Insure cross-platform build integrity

