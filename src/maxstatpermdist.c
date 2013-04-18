/*

  $Id: maxstatpermdist.c 372 2003-08-30 14:29:45Z hothorn $
  
  maxstatpermdist : Simulate Distribution of Maximally Selected Rank Statistics
  Copyright (C) 2003  Torsten Hothorn 
                      <Torsten.Hothorn@rzmail.uni-erlangen.de>
    
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or (at
  your option) any later version.
            
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
                   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
                          
  SYNOPSIS
 
    SEXP maxstatpermdist(SEXP scores, SEXP msample, SEXP expect, SEXP variance,
                         SEXP Nsim, SEXP pvalonly, SEXP ostat)
                                             
  DESCRIPTION
  
    maxstatpermdist:	permutation distribution of maximally selected 
    			rank statistics
                                                                  
*/

                                                                   
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>

/* this will move to exactRankTests */

SEXP sim2distr(SEXP distr) {

  int ns;		/* number of simulation runs: Nsim */
  
  SEXP counts;		/* vector of length ns for counting duplicates */
  SEXP dlist;		/* a list of two elements: */
  SEXP T, prob;		/* the vector of possible realizations T and
  			   the vector of corresponding probabilities */

  double stat;

  /*
  	little helpers
  */
  
  int i, lasti, countus = 0, this;

  ns = LENGTH(distr);
              
  PROTECT(counts = allocVector(INTSXP, ns));

  for (i = 0; i < ns; i++) INTEGER(counts)[i] = 0;
  /* we need to compute the distribution based on the simulation */
  
  /* 
     <FIXME>
     all of this is the overkill for simulating p-values only
     but currently I'll have the whole distribution at hand 
     (for conf.int's for example)
     </FIXME>
  */
    
  R_rsort(REAL(distr), ns);

  /* determine duplicates and its frequency */
  
  stat = REAL(distr)[0];
  lasti = 0;
  for (i = 0; i < ns; i++) {
    if (stat == REAL(distr)[i]) {
      INTEGER(counts)[lasti]++;
    } else {
      INTEGER(counts)[i]++;
      lasti = i;
    }
    if (INTEGER(counts)[i] == 0) countus++;
    stat = REAL(distr)[i];
  }
  
  /* 
     allocate memory for the return values: we do have ns-countus 
     unique statistics
  */
  
  countus = ns - countus;
  PROTECT(dlist = allocVector(VECSXP, 2));
  PROTECT(T = allocVector(REALSXP, countus));  
  PROTECT(prob = allocVector(REALSXP, countus));  

  /* we are ready to compute the (estimated) density now */
  
  this = 0;
  for (i = 0; i < ns; i++) {
    if (INTEGER(counts)[i] != 0) {
      REAL(T)[this] = REAL(distr)[i];
      REAL(prob)[this] = (double) INTEGER(counts)[i]/ns;
      this++;
    }
  }
  
  /* store the distribution in a lists with two elements */
  
  SET_VECTOR_ELT(dlist, 0, T);
  SET_VECTOR_ELT(dlist, 1, prob);
  
  UNPROTECT(4);
  return(dlist);
}


SEXP maxstatpermdist(SEXP scores, SEXP msample, SEXP expect, SEXP variance, 
                     SEXP Nsim, SEXP pvalonly, SEXP ostat) {

  /*
  
    Permutation distribution of maximally selected rank statistics
    
    scores:	N-vector of scores of the response variable
    msample:	sample sizes associated with all possible cutpoints
    expect:	vector of expecations of the statistics for each possible cutpoint
    variance:	vector of variances of the statistics for each possible cutpoint
    Nsim:	number of Monte-Carlo replications to be performed
    pvalonly:	logical whether a p-value is required only
    ostat:	in case pvalonly is TRUE, the observed statistic
    
  */

  SEXP ans;		/* SEXP for the return value: either p-value of distribution */  
  SEXP distr;		/* placeholder for the simulated statistics */
      
  /* mothers little helpers */

  SEXP myscores, mymsample, myexpect, myvariance;
  double *y, thisstat = 0.0;
  int *m, N, mN, ns;
  double *e, *v, *urand, *stat, dummy = 0.0;
  int i, j, *perm, k;

  PROTECT(myscores = coerceVector(scores, REALSXP));
  PROTECT(mymsample = coerceVector(msample, INTSXP));
  PROTECT(myexpect = coerceVector(expect, REALSXP));
  PROTECT(myvariance = coerceVector(variance, REALSXP));

  /* for readability reasons only */
  
  y = REAL(myscores);
  m = INTEGER(mymsample);
  e = REAL(myexpect);
  v = REAL(myvariance);
  ns = INTEGER(Nsim)[0];

  N = LENGTH(scores);
  mN = LENGTH(msample);
  
  /* this function is for internal use only, so no excessive checking here */
  
  if (mN != LENGTH(expect))
    error("expect not of length msample");

  if (mN != LENGTH(variance))
    error("variance not of length msample");
    
  urand = (double *) R_alloc(N, sizeof(double));
  stat = (double *) R_alloc(mN, sizeof(double));
  perm = (int *) R_alloc(N, sizeof(int));
  
  if (LOGICAL(pvalonly)[0] == TRUE) {
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = 0.0;
  } else {    
    PROTECT(ans = allocVector(VECSXP, 2));
  }
  PROTECT(distr = allocVector(REALSXP, ns));

  for (k = 0; k < mN; k++) v[k] = sqrt(v[k]);
                    
  for (j = 0; j < N; j++) perm[j] = j;

  GetRNGstate();
      
  for (i = 0; i < ns; i++) {
  
    /* a new permutation of the response */
    
    for (j = 0; j < N; j++)
      urand[j] = unif_rand();
    rsort_with_index(urand, perm, N);
    
    dummy = 0.0;
    k = 0;
    for (j = 0; j < N; j++) {
    
      if (k == mN) break;
    
      /* the statistic at cutpoint k is the standardized sum of the scores */
    
      dummy = dummy + y[perm[j]];
      if ((j + 1) == m[k]) {
        stat[k] = fabs((dummy - e[k])/v[k]);

        /*  | stat | > ostat is sufficient to count when a 
            p-value is to be computed only
        */

        if (LOGICAL(pvalonly)[0] == TRUE) {
          if (stat[k] > REAL(ostat)[0]) {
            REAL(ans)[0] = REAL(ans)[0] + 1.0;
            break;
          }
        }                         
        k++;
      }
    }
    
    /* get the maximum if the distribution is of special interest */

    if (LOGICAL(pvalonly)[0] == FALSE) {
      thisstat = 0.0;
      for (k = 0; k < mN; k++)
        if (stat[k] > thisstat) thisstat = stat[k];
      REAL(distr)[i] = thisstat;
    }
  }
  
  PutRNGstate();

  /* either return the p-value of the density */
  if (LOGICAL(pvalonly)[0] == TRUE) {
    REAL(ans)[0] = REAL(ans)[0] / ns;
  } else {
     ans = sim2distr(distr);
  }
  UNPROTECT(6);
  return(ans);
}

