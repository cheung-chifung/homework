/* sample program: for average-case analysis of MAX-3LIN solvers
   file:     sample3gwalk.c
   update:   2010/7/7 revised 7/16
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt.h"

#define K         3       /* # of variables in each eq */
#define MAXvar    200     /* max # of variables */
#define MAXeq     10000   /* max # of equations */
#define MAXv2e    150     /* max # of equations containing each var. */

#define Nex       5       /* # of example instances per each param. setting */

int main(int argc, char *argv[])
{
  /* random seeds */
  int rand_seed;

  /* parameters */
  int nvar, ncombi;
  int p100;  /* clause density parameter p = p100/(100*n*ln(n)) */
  int q100;  /* noise prob. q = q100/100 */
  double p, q;

  /* for generating formulas */
  int neq, nerreq, maxv2e;
  int E[MAXeq][3], R[MAXeq];
  int V2E[MAXvar][MAXv2e], nV2E[MAXvar];
  int v1, v2, v3, rval;
  double r;

  /* for heuristics (greedy random walk) */
  int A[MAXvar], Vadv[MAXvar], Esat[MAXeq];
  int tt, maxtt;
  int minadv, nminadv, Vlist[MAXvar], vflip;
  int dist, nunsat;

  /* work var.s */
  int iex;
  int i, j, v, eq, val;

  /****************/
  /* Program Body */
  /****************/

  rand_seed = 1234;
  sgenrand(rand_seed);

  printf("# n = 100, 110, ..., 200\n");
  printf("# P = 0.2/nlog n, 0.3/nlog n, ..., 10/nlog n, Q = 0\n");
  printf("# random walk steps = 1.5 * n\n");

  for(nvar = 100; nvar <= 200; nvar = nvar + 10) {
    printf("> n = %d\n", nvar);

    for(p100 = 200; p100 <= 300; p100 = p100 + 100) {
      q100 = 0;

      p = (double)p100 / ( 100.0 * (double)nvar * log((double)nvar) );
      q = (double)q100 / 100.0;
      ncombi = nvar * (nvar - 1) * (nvar - 2) / 6;

      printf("> P, Q = %2.4f, %2.4f:\n", p, q);
      /*
      printf("> average  n.eqs = %d, n.err.eqs = %d, n.eqs.per.var = %d\n",
	     (int)((float)ncombi * p),
	     (int)((float)ncombi * p * q),
	     (int)((float)(nvar * (nvar - 1)) / 2.0 * p));
      */

      for(iex = 0; iex < Nex; iex++) {

	/****************/
	/* generate eqs */
	/****************/
	neq = nerreq = maxv2e = 0;
	for(v = 1; v <= nvar; v++) nV2E[v] = 0;
	for(v1 = 1; v1 <= nvar; v1++)
	  for(v2 = v1 + 1; v2 <= nvar; v2++)
	    for(v3 = v2 + 1; v3 <= nvar; v3++) {
	      r = genrand();
	      if(r < p) {
		E[neq][0] = v1; E[neq][1] = v2; E[neq][2] = v3;
		V2E[v1][nV2E[v1]] = neq;
		nV2E[v1]++;  if(nV2E[v1] > maxv2e) maxv2e = nV2E[v1];
		V2E[v2][nV2E[v2]] = neq;
		nV2E[v2]++;  if(nV2E[v2] > maxv2e) maxv2e = nV2E[v2];
		V2E[v3][nV2E[v3]] = neq;
		nV2E[v3]++;  if(nV2E[v3] > maxv2e) maxv2e = nV2E[v3];
		r = genrand();
		if(r < q) {
		  R[neq] = -1;
		  nerreq++;
		}
		else
		  R[neq] = +1;
		neq++;
	      }
	    }
	if(neq > MAXeq || maxv2e > MAXv2e) {
	  printf("error too many: %d, %d\n", neq, maxv2e);
	  exit(-1);
	}
	/*
	printf("> obtained n.eqs = %d, n.err.eqs = %d, n.eqs.per.var = %d\n",
	       neq, nerreq, maxv2e);
	*/

	/*********************************/
	/* heuristics solver simultation */
	/* greedy random walk            */
	/*********************************/
	maxtt = 1.5 * (float)nvar;

	/* initial assignment */
	dist = 0;
	for(v = 1; v <= nvar; v++) {
	  if(genrand() > 0.5) {
	    A[v] = +1;
	  } 
	  else {
	    A[v] = -1;
	    dist++;
	  }
	  Vadv[v] = 0;
	}
	nunsat = 0;
	for(eq = 1; eq <= neq; eq++) {
	  val = +1;
	  for(i = 0; i < 3; i++) val = val * A[E[eq][i]];
	  if(val == R[eq]) {
	    Esat[eq] = +1;
	    for(i = 0; i < 3; i++) Vadv[E[eq][i]]++;
	  }
	  else {
	    Esat[eq] = -1;
	    for(i = 0; i < 3; i++) Vadv[E[eq][i]]--;
	    nunsat++;
	  }
	}

	/* one execution */
	printf("(%d, %d) ==> ", dist, nunsat);
	for(tt = 1; tt <= maxtt; tt++) {
      
	  /* var.s with max. disadvantage */
	  minadv = nvar;
	  nminadv = 0;
	  for(v = 1; v <= nvar; v++) {
	    if(Vadv[v] < minadv) {
	      minadv = Vadv[v];
	      nminadv = 0;
	      Vlist[nminadv] = v;
	      nminadv++;
	    }
	    else if(Vadv[v] == minadv) {
	      Vlist[nminadv] = v;
	      nminadv++;
	    }
	  }
	  /* choose one flipped var. */
	  i = genrand() * (float)nminadv;
	  vflip = Vlist[i];

	  /* update vflip itself */
	  A[vflip] = -A[vflip];
	  if(A[vflip] == +1) dist--; else dist++;

	  /* update the other related var.s */
	  for(i = 0; i < nV2E[vflip]; i++) {
	    eq = V2E[vflip][i];
	    if(Esat[eq] < 0) {
	      Esat[eq] = +1;
	      nunsat--;
	      for(j = 0; j < 3; j++) Vadv[E[eq][j]] = Vadv[E[eq][j]] + 2;
	    }
	    else { /* Esat[eq] > 0 */
	      Esat[eq] = -1;
	      nunsat++;
	      for(j = 0; j < 3; j++) Vadv[E[eq][j]] = Vadv[E[eq][j]] - 2;
	    }
	  }
	  if(dist < 2 || nunsat < 10) break;
	}
	printf("(%d, %d), t = %d\n", dist, nunsat, tt);
      }
    }
  }
}
