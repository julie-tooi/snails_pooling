---
output: html_document
---

# Effect of the proteomic sample pooling on statistical power and type I error rate
## Project description

Pulling - mixing biological material of samples, for example, when combining several individual organisms in one sample for analysis. 

As a result, the expression level of most proteins in the pool will be equal to the average expression level of the corresponding proteins in individual samples, the total amount of biological material in the sample increases, the quality of comparison of protein spots on different gels increases. Finally, the level of biological variability between individual pools is reduced, and this affects the power of statistical tests.

In this case, questions arise:

(1) How much the biological variability of the pool is equal to the real biological variability
(2) Are the average values of individual samples really equal to the average value of the pool
(3) When combining samples possibility of data loss

This problem was considered by the authors from a mathematical point of view [1], experimental work was also carried out on a small sample [2]. We tried to answer questions using modeling with randomization. We studied the behavior of the moderated t-test, which is often used in proteomic studies.

## Goal and tasks

Goal: Analyze how pool size can affect the power of a moderated t-test

Tasks:

- Simulation for test power (Koyenov effect)
- Different scenarios
- Simulation of differential expression data. Technical and biological variability
- Sample pooling
- Raw data (with technical replications)
- Pooled data (different pool sizes)
- Simulation of data with different sample sizes

## Methods

R script performs simulation with different parameters (user-defined):

- Amount of protein
- Level of technical variability
- Level of biological variability
- Log expression level in group 1
- Difference between the groups (the parameter through which the log-expression is calculated in group 2)
- Number of generation repeats for each scenario
- Number of repetitions of generations

Mathematical modeling using the limma package, occurs according to the principle:
Model is being built for one protein with different levels of expression in groups; from the model, we extract the moderated t-test
then the script generations are repeated a specified number of times, confidence intervals are determined.
Graphing is done using the ggplot2 package
To speed up calculations, parallelization of calculations across different kernels using the dbframe package is used.

## Requirements
R 3.6
Rstudio 1.2

Tested on Ubuntu 19.04, R 3.6.1, Rstudio 1.2.5019.

Libraries: plyr [3], dplyr [4], limma [5], ggplot2 [6], tidyr [7], dbframe [8]

## Running script
The script get started in the Rstudio environment.

Variables must be specified:

- Number of Generations <mark> (expression set) </mark>
- Number of proteins <mark> (number of lines in expression set) </mark>
- Expression level of one of the groups <mark> (beta 0) </mark>
- Level of differences in the expression of groups <mark> (beta 1) </mark>
- Degree of biological variability <mark> (sig B) </mark>
- Degree of technical variability <mark> (sig e) </mark>

## Results
![](https://s113sas.storage.yandex.net/rdisk/5d962a5a69df8da17fbb07c5591514a3e4f10d5e89200e691090b49cda5c63ad/5e3df18d/fKqInKw3d7bLFOeFnMGnhOb3edIyamFhWf4Uwg5ySpDJSEmrILl_KY9P2tbiGk3ZgCL4254OKF73jPp2DqhqreQJXnv3OeZc6ErgGSxXO6ar8npumZHI4midPdWhecNq?uid=1130000038919147&filename=fin_r_plot.png&disposition=inline&hash=&limit=0&content_type=image%2Fpng&owner_uid=1130000038919147&fsize=132636&hid=846f9b03b0452173074931151a4d9c5b&media_type=image&tknv=v2&etag=c8fc833b07c9484b5d166351a5ed16f3&rtoken=ZExucjrFXtWE&force_default=yes&ycrid=na-8e4cdf4e8cb823e078ac6261dc1bec7a-downloader23h&ts=59e04af05bd40&s=263f935ff8425893cfe02ba6b256e5f62fa92855926b0a15326ceb448a0745a6&pb=U2FsdGVkX1_BPAvv_xUU5Q4ZfvkBe676aqG7_14xgmYJM2d24yGKmh1XAB8jc7G4e7oM2Y9JICRDuwu4wG2GEYs30gaAce-K_d7qVVCeKvwSRunr-hfHZfMZbdUFmeZK)

## References 
[1] Shu-Dong Zhang, Timothy W. Gant, Effect of pooling samples on the efficiency of comparative studies using microarrays, Bioinformatics, Volume 21, Issue 24, 15 December 2005, Pages 4378–4383, https://doi.org/10.1093/bioinformatics/bti717.

[2] Diz, A. P., Truebano, M., & Skibinski, D. O. (2009). The consequences of sample pooling in proteomics: an empirical study. Electrophoresis, 30(17), 2967–2975.

[3] Wickham H (2011). “The Split-Apply-Combine Strategy for Data Analysis.” Journal of Statistical Software, 40(1), 1–29. http://www.jstatsoft.org/v40/i01/.

[4] Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2018). dplyr: A Grammar of Data Manipulation. R package version
  0.7.6. https://CRAN.R-project.org/package=dplyr.
  
[5] Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47. doi: 10.1093/nar/gkv007.

[6] H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

[7] Hadley Wickham and Lionel Henry (2019). tidyr: Tidy Messy Data. R package version 1.0.0. https://CRAN.R-project.org/package=tidyr.

[8] Gray Calhoun (2010). dbframe: An R to SQL interface. R package version 0.4.0.
