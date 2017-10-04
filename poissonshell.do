**************************************************************************
*
* Program to try out power and sample size for COS project
*
* this is the "shell" program for poissonpower 
*
*
* Author: Susanne May
* Date: 09/16/2017 
*
*
* this program requires 5 input parameters, NOTE input parameter do not seem to work, directly code the values instead
*
*       1 = number of observations total (to be divided into 2 groups)
*       2 = exp(beta0) for the SOC group 
*       3 = exp(beta1) for the tx group 
*       4 = percent of observations with zero count/outcome value in SOC group (1)
*       5 = percent of observations with zero count/outcome value in tx group (2)
*       6 = random number seed
*       7 = number of replications 
*       8 = number of log output file
*
*       . do poissonshell 1280 3.1 0.31 0.69 0.69 39083 500
*
**************************************************************************

    display "Simulation started  $S_DATE   $S_TIME"
    
    clear
    use one
    replace nummer=1 in 1
    save,replace
    clear

    simulate obs=r(obs) meanSOC1=r(meanSOC1) sdSOC1=r(sdSOC1) n1=r(n1) n1_0s=r(n1_0s) perc0sSOC1=r(perc0sSOC1) ///
        meantx2=r(meantx2) sdtx2=r(sdtx2) n2=r(n2) n2_0s=r(n2_0s) perc0stx2=r(perc0stx2) ///
        meandiff=r(meandiff) realdiff=r(realdiff) ///
        pvallrtest=r(pvallrtest) chi2lrtest=r(chi2lrtest) reject05lrtest=r(reject05lrtest) ///
        pvalttestu=r(pvalttestu) tttestu=r(tttestu) reject05tttestu=r(reject05tttestu) incittu=r(incittu) ///
        pvaltteste=r(pvaltteste) ttteste=r(ttteste) reject05ttteste=r(reject05ttteste) incitte=r(incitte) ///
        b1wo=r(b1wo) b1sewo=r(b1sewo) b1pvalwo=r(b1pvalwo) reject05regrwo=r(reject05regrwo) b1uciwo=r(b1uciwo) b1lciwo=r(b1lciwo) incilrwo=r(incilrwo) ///
        b1w=r(b1w) b1sew=r(b1sew) b1pvalw=r(b1pvalw) reject05regrw=r(reject05regrw) b1uciw=r(b1uciw) b1lciw=r(b1lciw) incilrw=r(incilrw) ///
	    , reps(`7'): simula, nn(`1') beta0(`2') beta1(`3') perc0s1(`4') perc0s2(`5') seed(`6')

    display "Program was run as: do poissonshell `*'"
    display "with arguments: n beta0 beta1 perc0s1 perc0s2 seed reps"
    sum obs-realdiff
    sum pvallrtest-reject05lrtest
    sum pvaltteste-incitte
    sum pvalttestu-incittu
    sum b1wo-incilrwo
    sum b1w-incilrw
    
    capture save results.dta, replace

    display "Simulation ended  $S_DATE   $S_TIME"



