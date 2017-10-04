**************************************************************************
*
* Program to compare power and size for t-test vs (zero inflated) poisson regression  
* for Steroid trial (COS) to respond to reviewer's comments
*
*
* Author: Susanne May
* Date: 09/15/2017
*
*
* this program requires input parameters, which need to be supplied via the 
* poissonshell do file
*
* These input parameters are: 
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
* run this file once as 
*  
*    do poissonpower
*
* before running dothedofile.do
*
**************************************************************************
/*
*** for singlerun
  local nn=1000
  local beta0=log(5.7)
  local beta1=-0.3
  local perc0s1=0.6
  local perc0s2=0.6
  
  local seed=389201
  local reps=1
  local no=01
  
  set more off
  set seed `seed'

*** end for singlerun
*/
clear

capture log close
*** for single run use next line
* capture log using sims_onerun.log, replace


capture program drop simula
program define simula, rclass

        version 15.0
        
        *----------------
        syntax [, nn(integer 1) beta0(real 1) beta1(real 1) ///
             perc0s1(real 1) perc0s2(real 1) seed(integer 1) ]
        *----------------
        drop _all
        
        use one
        local nummerin=nummer[1]
        local nummerout=`nummerin'+1
        replace nummer=`nummerout' in 1
        save,replace
        clear
                
        set obs `nn'
        
        gen x=group(2)
        tab x, gen(x)
        
        *** generate linear predictor 
        gen lp=`beta0'+`beta1'*x2
        gen mu=exp(lp)
        
        *** generate percent of zeros for each group separately from Poisson distribution
        *** to achieve zero inflated Poisson
        
        gen uni=uniform()
        gen     xp=0 if uni<`perc0s1' & x==1
        replace xp=0 if uni<`perc0s2' & x==2
        
        replace xp=rpoisson(mu) if xp!=0

        sum xp if x==1
        local meanSOC1=r(mean)
        local sdSOC1=r(sd)
        count if x==1
        local n1=r(N)
        count if x==1 & xp==0
        local n1_0s=r(N)
        local perc0sSOC1=`n1_0s'/`n1'

        sum xp if x==2
        local meantx2=r(mean)
        local sdtx2=r(sd)
        count if x==2
        local n2=r(N)
        count if x==2 & xp==0
        local n2_0s=r(N)
        local perc0stx2=`n2_0s'/`n2'
        
        local meandiff=`meantx2'-`meanSOC1'
        
        *** poisson regression to test for differences        
        poisson xp x2
        estimates store model1
        poisson xp
        lrtest model1

        local pvallrtest=r(p)
        local chi2lrtest=r(chi2)
        local reject05lrtest=(r(p)<=0.05)

        
        *** t-test to test for differences, unequal variances
        ttest xp, by(x) unequal
 
        *** calculate upper and lower boundaries for CI for difference
        local sdp=sqrt(r(sd_1)^2/r(N_1)+r(sd_2)^2/r(N_2))
        local dfu=r(df_t)
        local cittu_low=-(r(mu_1)-r(mu_2)+invt(`dfu',0.975)*`sdp')
        local cittu_upp=-(r(mu_1)-r(mu_2)-invt(`dfu',0.975)*`sdp')
        local onem1=1-`perc0s1'
        local onem2=1-`perc0s2'
        local txgrp=`beta1'+`beta0'
        
        local realdiff=exp(`txgrp')*`onem2'-exp(`beta0')*(`onem1')
        local incittu=(`realdiff'<=`cittu_upp' & `realdiff'>=`cittu_low')

        local pvalttestu=r(p)
        local tttestu=r(t)
        local reject05tttestu=(r(p)<=0.05)
        

        *** t-test to test for differences, equal variances
        ttest xp, by(x)

        local sdp=sqrt(((r(N_1)-1)*r(sd_1)^2+(r(N_2)-1)*r(sd_2)^2)/(r(N_1)+r(N_2)-2))
        local df=r(N_1)+r(N_2)-2
        local citte_low=-(r(mu_1)-r(mu_2)+invt(`df',0.975)*`sdp'*sqrt(1/r(N_1)+1/r(N_2)))
        local citte_upp=-(r(mu_1)-r(mu_2)-invt(`df',0.975)*`sdp'*sqrt(1/r(N_1)+1/r(N_2)))
        *** onem1, onem2, txgrp and realdiff already calculated above
        local incitte=(`realdiff'<=`citte_upp' & `realdiff'>=`citte_low')

        local pvaltteste=r(p)
        local ttteste=r(t)
        local reject05ttteste=(r(p)<=0.05)

        
        *** linear regression withOUT robust standard errors 
        regress xp x
        matrix regrHFSworobust=r(table)
        matrix list regrHFSworobust
        local b1wo=regrHFSworobust[1,1]
        local b1sewo=regrHFSworobust[2,1]
        local b1pvalwo=regrHFSworobust[4,1]
        local reject05regrwo=(`b1pvalwo'<=0.05)
        local b1uciwo=regrHFSworobust[5,1]
        local b1lciwo=regrHFSworobust[6,1]
*        display "b1wo=" `b1wo' " b1sewo=" `b1sewo' " b1uciwo=" `b1uciwo' " b1lciwo=" `b1lciwo'
        local incilrwo=(`realdiff'<=`b1lciwo' & `realdiff'>=`b1uciwo')



        *** linear regression with robust standard errors 
        regress xp x, vce(robust)
        matrix regrHFSwrobust=r(table)
        matrix list regrHFSwrobust
        local b1w=regrHFSwrobust[1,1]
        local b1sew=regrHFSwrobust[2,1]
        local b1pvalw=regrHFSwrobust[4,1]
        local reject05regrw=(`b1pvalw'<=0.05)
        local b1uciw=regrHFSwrobust[5,1]
        local b1lciw=regrHFSwrobust[6,1]
*        display "b1w=" `b1w' " b1sew=" `b1sew' " b1uciw=" `b1uciw' " b1lciw=" `b1lciw'
        local incilrw=(`realdiff'<=`b1lciw' & `realdiff'>=`b1uciw')




        count 
        local obss=r(N)

*       *** for testing
*        if `nummerin'==7 {
*         capture save data`nummerin', replace
*        }


        return scalar obs=`obss'
        return scalar meanSOC1=`meanSOC1'
        return scalar sdSOC1=`sdSOC1'
        return scalar n1=`n1'
        return scalar n1_0s=`n1_0s'
        return scalar perc0sSOC1=`perc0sSOC1'
        return scalar meantx2=`meantx2'
        return scalar sdtx2=`sdtx2'
        return scalar n2=`n2'
        return scalar n2_0s=`n2_0s'
        return scalar perc0stx2=`perc0stx2'
        return scalar meandiff=`meandiff'
        return scalar realdiff=`realdiff'

        return scalar pvallrtest=`pvallrtest'
        return scalar chi2lrtest=`chi2lrtest'
        return scalar reject05lrtest=`reject05lrtest'
        
        return scalar pvalttestu=`pvalttestu'
        return scalar tttestu=`tttestu'
        return scalar reject05tttestu=`reject05tttestu'
        return scalar incittu=`incittu'
        
        return scalar pvaltteste=`pvaltteste'
        return scalar ttteste=`ttteste'
        return scalar reject05ttteste=`reject05ttteste'
        return scalar incitte=`incitte'
        
        return scalar b1wo=`b1wo'
        return scalar b1sewo=`b1sewo'
        return scalar b1pvalwo=`b1pvalwo'
        return scalar reject05regrwo=`reject05regrwo'
        return scalar b1uciwo=`b1uciwo'
        return scalar b1lciwo=`b1lciwo'
        return scalar incilrwo=`incilrwo'
        
        return scalar b1w=`b1w'
        return scalar b1sew=`b1sew'
        return scalar b1pvalw=`b1pvalw'
        return scalar reject05regrw=`reject05regrw'
        return scalar b1uciw=`b1uciw'
        return scalar b1lciw=`b1lciw'
        return scalar incilrw=`incilrw'
        
        
    end        

* see poissonshell for how this program is called
