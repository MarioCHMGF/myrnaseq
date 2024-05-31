# A function to calculate DEGS from raw count matrix. It creates all additional figures and normalizes the data using TMM and VOOM.
## Example:
```
countdata # dataframe of raw counts
sampleinfo # sampleinfo file (needs columns: sample and group aka. outcome)
# please double check and perform sanity check that table(names(countdata) == sampleinfo$sample) is all TRUE!
samplenum # number of all samples together
groupnum # minimal samples per group
my.object <- myrnaseq(countdata, sampleinfo, samplenum, groupnum)
my.object$degs - returns calculated degs
my.object$dgeObj - returns dgeObj from edgeR for additional analyses e.g. adjusting for covariates
```
### Example of sampleinfo.txt
```
sample  group
x1  0
x2  0
x3  0
x4  1
x5  1
x6  1
# if not specified in factor levels, it is always 1 relative to 0. If specified with levels, it is always level2 relative to level1; levels=c(level1, level2) 
```
# Additonally it can calculate degs adjusted up to four covariates
## Example:
```
dgeObj # dgeObj from your analysis (can be from my.object$dgeObj from myrnaseq function)
outcome # outcome variable from sampleinfo file (factor, ordered factor or numeric e.g. sampleinfo$group)
cov1 # first covariate to adjust for e.g. sampleinfo$cov1
cov2 # second covariate to adjust for
cov3 # third covariate to adjust for
cov4 # fourth covariate to adjust for
my.adjusted.results <- calcdeg(dgeObj, outcome, cov1, cov2, cov3, cov4)
# defualt for all covs is NULL
```
# Draw a volcano plot
## Example:
```
myresults # annotated results from calcdeg
Q # if TRUE take Q values, default is FALSE
top # number of top genes to annotate, default is 5
xl # lower x axis value, must be negative, default is NULL
xu # upper x axis value, must be positive, default is NULL
plot <-  volcano(myresults, Q=FALSE, top=5, xl=NULL, xu=NULL)
```
