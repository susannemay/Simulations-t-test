* call poissonshell multiple times (which calls poissonpower) and requires 5 input parameters
*
*       1 = number of observations total (to be divided into 2 groups)
*       2 = exp(beta0) for the SOC group
*       3 = exp(beta1) for the tx group
*       4 = percent of observations with zero count/outcome value in SOC group (1)
*       5 = percent of observations with zero count/outcome value in tx group (2)
*       6 = random number seed
*       7 = number of replications 
*
*       . do poissonshell 1280 3.1 0.31 0.69 0.69 539083 10000 
*

capture log close
log using runlots.log, replace

do poissonshell 1280 2.82 0.335 0.69 0.69 105637 10000

*do poissonshell 1280 2.82 0.27 0.69 0.67 427433 1000

log close
