/* $Id: newcorr.c 372 2003-08-30 14:29:45Z hothorn $ 
**
** Correlation matrix for maximally selected rank statistics.
**
** Source:
**
** Lausen, Hothorn, Bretz, Schumacher (2002): Assessment of optimal 
** selected prognostic factors, submitted.
**
** Input variables:
** 	ilist: 	a list of ranks of observations of each prognostic factor
**		under test (computed using irank, package exactRankTests)
**
**	prop:	the proportions of the minimum sample size in the first group
**		usually (0.1, 0.9).
**
** Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>
**
*/

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

SEXP newcorr (SEXP ilist, SEXP prop) {
  int p, nobs, kx, ky, i, j, k, l, mn, m, n;
  int nobsmm, nobsmn, index;
  int *sortx, *sorty, help;
  int minprop, maxprop;
  double fact;
  int *dupx, *dupy;
  int getoutahere=0;
  SEXP cm, outlist, rowdel, coldel;

  /*
  ** check if ilist is a list and get the number of prognostic factors (p)
  ** and observations (nobs)
  */
  if(!isNewList(ilist)) error("`ilist' must be a list");
  p = LENGTH(ilist);
  nobs = LENGTH(VECTOR_ELT(ilist, 0));
  
  /*
  ** calculate the minimal and maximal sample size allowed for the first group
  */
  minprop = floor(nobs*REAL(prop)[0]);
  if (minprop == 0) minprop=1;
  maxprop = floor(nobs*REAL(prop)[1]);
  
  /* 
  ** helpers for determination of ties
  */
  dupx = (int *) R_alloc(nobs, sizeof(int));
  dupy = (int *) R_alloc(nobs, sizeof(int));
  for (l = 0; l < nobs; l++) {
    dupx[l] = 0;
    dupy[l] = 0;
  }
  
  /*
  ** cm is the correlation matrix and rowdel/coldel are 0/1 vectors indicating
  ** rows/cols to be deleted in R (because of ties or sample 
  ** sizes > < min/maxprop)
  */
  PROTECT(cm = allocMatrix(REALSXP, nobs*p, nobs*p));
  PROTECT(rowdel = allocVector(INTSXP, nobs*p));
  PROTECT(coldel = allocVector(INTSXP, nobs*p));
  PROTECT(outlist = allocVector(VECSXP, 3));
  for (i = 0; i < nobs*p; i++) {
    INTEGER(rowdel)[i] = 0;
    INTEGER(coldel)[i] = 0;
    for (j = 0; j < nobs*p; j++) {
      if (i == j) {
        /* diagonal elements */
        REAL(cm)[i + nobs*p*j] = 1.0; 
      } else {
        REAL(cm)[i + nobs*p*j] = 0.0;
      }
    }
  }

  for (kx = 0; kx < p; kx++) {
  
    /* extract the first prognostic factor: "x" */
    sortx = INTEGER(VECTOR_ELT(ilist, kx));
    for (ky = 0; ky < p; ky++) {

      /* and compute the correlations for the second one: "y"*/
      sorty = INTEGER(VECTOR_ELT(ilist, ky));

      /* reset helpers */
      for (l=0; l < nobs; l++) {
        dupx[l] = 0;
        dupy[l] = 0;
      }

      /* for each cutpoint in x */
      for (i=0; i < nobs; i++) {
        m = sortx[i];
        
        /* 
        ** check if is tied and do not consider this correlation
        */
        if (dupx[m-1] == 1) {
          INTEGER(rowdel)[kx*nobs + i] = 1;
          continue;
        } else {
          dupx[m-1] = 1;
        }
        
        /* reset helper for y */
        for (l=0; l < nobs; l++) {
          dupy[l] = 0;
        }
        
        /* for each cutpoint in y */
        for (j=0; j < nobs; j++) {
          n = sorty[j];

          /* check for ties */
          if (dupy[n-1] == 1) {
            INTEGER(coldel)[ky*nobs + j] = 1;
            continue;
          } else {
            dupy[n-1] = 1;
          }                   
          
          /*
          ** those splits are not allowed, therefore no correlation.
          */
          if (m < minprop || m > maxprop) {
            INTEGER(rowdel)[kx*nobs + i] = 1;
            getoutahere = 1;
          }
          if (n < minprop || n > maxprop) {
            INTEGER(coldel)[ky*nobs + j] = 1;
            getoutahere = 1;
          }
          if (getoutahere == 1) {
            getoutahere = 0;
            continue;
          }

          /* 
          ** determine the correct entry in cm (no array indexing here)
          */                 
          index =  (kx*nobs + i)*p*nobs + ky*nobs+j;

          /*
          ** count the elements both in the first group
          */
          mn = 0;
          for (k=0; k < nobs; k++) {
            if (sortx[k] <= m && sorty[k] <= n) mn++;
          }
          
          /*
          ** and calculate the correlation
          */
          nobsmm = nobs - m;
          nobsmn = nobs - n;
          fact = nobs*sqrt(m)*sqrt(nobsmm)*sqrt(n)*sqrt(nobsmn);
          fact = 1/(double)fact;
          help = mn*nobsmm*nobsmn;
          help += -(m -mn)*nobsmm*n;
          help += -(n-mn)*nobsmn*m;
          help += (nobsmm - (n- mn))*m*n;
          
          /*
          ** ok: place yourself here
          */
          REAL(cm)[index] = fact*help;
        }
      }
    }
  }
  
  /*
  ** we need to return an array and two vectors, use a list
  */
  SET_VECTOR_ELT(outlist, 0, cm);
  SET_VECTOR_ELT(outlist, 1, rowdel);
  SET_VECTOR_ELT(outlist, 2, coldel);
  UNPROTECT(4);
  return(outlist);
}
