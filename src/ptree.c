#include <R.h>
#include <Rinternals.h>
#include <R_ext/Mathlib.h>

int cardunion(int *x, int *y, int m, int n)
{
  /*
    cardinality of the union of x and y
  */
  
  int count = 0;
  int i, j;
  
  for (i=0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (x[i] == y[j]) { count++; break; }
      if (x[i] < y[j]) break; /* x, y are ordered */
    }
  }
  return(count);
}

void setminus(int *index, int *x, int *y, int m, int n)
{

  /*
    x \ y
  */
  
  int i, j, count = 0;
  
  for (i=0; i < m; i++) index[i] = 1;
  
  for (i=0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (x[i] == y[j]) { index[i] = 0; count++; break; }
      if (x[i] < y[j]) break; /* x, y are ordered */
    }
  }
}

double corrgauss(int *x, int *y, int m, int n, int N)
{
  int i, j, cxy, cnxy;
  int *Nindex, *notx, *noty, *notxindx, *notyindx;
  int mnotx = 0, nnoty = 0;
  double fact, rhs;
  
  Nindex = (int *) calloc(N, sizeof(int));
  notxindx = (int *) calloc(N, sizeof(int));
  notyindx = (int *) calloc(N, sizeof(int));


  for (i = 0; i < N; i++) Nindex[i] = i + 1;

  setminus(notxindx, Nindex, x, N, m);
  setminus(notyindx, Nindex, y, N, n);
   
  for (i = 0; i < N; i++) { mnotx += notxindx[i]; nnoty += notyindx[i]; }
  
  notx = (int *) calloc(mnotx, sizeof(int));
  noty = (int *) calloc(nnoty, sizeof(int));
  
  j = 0;
  for (i = 0; i < N; i++) {
    if (notxindx[i] == 1) {
       notx[j] = Nindex[i];
       j++;
    }
  }
  
  j = 0;
  for (i = 0; i < N; i++) {
    if (notyindx[i] == 1) {
       noty[j] = Nindex[i];
       j++;
    }
  }
  

  fact = (sqrt(m)*sqrt(N - m)*sqrt(n)*sqrt(N - n)) / N;
  cxy = cardunion(x, y, m, n);
  cnxy = n - cxy; 
 
  rhs = 1 / (double)(m * n) * cxy;
  rhs = rhs - 1/ (double) (m * (N-n)) * (m - cxy); 
  rhs = rhs - 1/ (double) ((N - m)*n) * cnxy;
  rhs = rhs + 1/ (double) ((N - m)*(N - n)) * (mnotx - cnxy); 
  free((void *) Nindex);
  free((void *) notxindx);
  free((void *) notyindx);
  free((void *) notx);
  free((void *) noty);
  return(fact*rhs);
}

SEXP corr(SEXP ilist, SEXP N) {
  int n, i, j, m, mm;
  SEXP cm;
  int *x;
  int *y;
  
  n = length(ilist);
  
  PROTECT(cm = allocMatrix(REALSXP, n, n));
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      REAL(cm)[i + n*j] = 0;
    }
  }
  
  for (i = 0; i < n; i++) {
    m = length(VECTOR_ELT(ilist, i));
    x = INTEGER(VECTOR_ELT(ilist, i));
    for (j = i; j < n; j++) {
      mm = length(VECTOR_ELT(ilist, j));
      y = INTEGER(VECTOR_ELT(ilist, j));
      REAL(cm)[i + n*j] = corrgauss(x, y, m, mm, INTEGER(N)[0]);
    }
  }
  UNPROTECT(1);
  return(cm);
}
