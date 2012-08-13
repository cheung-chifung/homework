/* sample program: for average-case analysis
                   of greedy random walk solver of MAX-2LIN
   file:     sample2gwalk.c
   update:   2010/7/7 revised 7/16
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt.h"

#define K         2       /* # of variables in each eq */
#define MAXvar    400     /* max # of variables */
#define MAXeq     1000000   /* max # of equations */
#define MAXv2e    300     /* max # of equations containing each var. */

#define Nex       3       /* # of example instances per each param. setting */

int main(int argc, char *argv[])
{
  /* random seeds */
  int rand_seed;

  /* parameters */
  int nvar, ncombi;
  int p10;  /* clause density parameter p = p10/(10*n) */
  int q10;  /* noise prob. q = q10/10 */
  double p, q;

  /* for generating formulas */
  int neq, nerreq, maxv2e;
  int E[MAXeq][2], R[MAXeq];
  int V2E[MAXvar][MAXv2e], nV2E[MAXvar];
  int v1, v2, v3, rval;
  double r;
  long mt;

  /* for heuristics (greedy random walk) */
  int A[MAXvar], Vadv[MAXvar], Esat[MAXeq];
  int tt, maxtt;
  int minadv, nminadv, Vlist[MAXvar], vflip;
  int dist, nunsat;

  /* for heuristics (message passing algo.) */

  /* work var.s */
  int iex;
  int i, j, v, eq, val;

  /****************/
  /* Program Body */
  /****************/

  rand_seed = 1234;
  sgenrand(rand_seed);

  printf("# n = 100, 110, ..., 200\n");
  printf("# P = 0.2/n, 0.3/n, ..., 10/n, Q = 0\n");
  printf("# random walk steps = 1.5 * n\n");

  for(nvar = 100; nvar <= 200; nvar = nvar + 10) {
    printf("> n = %d\n", nvar);

    for(p10 = 20; p10 <= 20; p10 = p10 + 10) {
      q10 = 0;
      p = (double)p10 / ( 10.0 * (double)nvar );
      q = (double)q10 / 10.0;
      ncombi = nvar * (nvar - 1) / 2;

      printf("> P, Q = %2.3f, %2.3f:\n", p, q);
      /*
      printf("> average  n.eqs = %d, n.err.eqs = %d, n.eqs.per.var = %d\n",
	     (int)((float)ncombi * p),
	     (int)((float)ncombi * p * q),
	     (int)((float)nvar * p));
      */

      for(iex = 0; iex < Nex; iex++) {

	/****************/
	/* generate eqs */
	/****************/
	neq = nerreq = maxv2e = 0;
	for(v = 1; v <= nvar; v++) nV2E[v] = 0;
	for(v1 = 1; v1 <= nvar; v1++)
	  for(v2 = v1 + 1; v2 <= nvar; v2++) {
	      r = genrand();
	      if(r < p) {
		E[neq][0] = v1; E[neq][1] = v2;
		V2E[v1][nV2E[v1]] = neq;
		nV2E[v1]++;  if(nV2E[v1] > maxv2e) maxv2e = nV2E[v1];
		V2E[v2][nV2E[v2]] = neq;
		nV2E[v2]++;  if(nV2E[v2] > maxv2e) maxv2e = nV2E[v2];
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
  // A[]を初期化
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
  //すべての組み合せを計算する
	for(eq = 1; eq <= neq; eq++) {
	  val = +1;
	  for(i = 0; i < 2; i++) val = val * A[E[eq][i]]; //論理式を計算
	  if(val == R[eq]) {
	    Esat[eq] = +1;
	    for(i = 0; i < 2; i++) Vadv[E[eq][i]]++; //優位度を計算
	  }
	  else {
	    Esat[eq] = -1;
	    for(i = 0; i < 2; i++) Vadv[E[eq][i]]--;
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
	      for(j = 0; j < 2; j++) Vadv[E[eq][j]] = Vadv[E[eq][j]] + 2;
	    }
	    else { /* Esat[eq] > 0 */
	      Esat[eq] = -1;
	      nunsat++;
	      for(j = 0; j < 2; j++) Vadv[E[eq][j]] = Vadv[E[eq][j]] - 2;
	    }
	  }
	  if(dist < 2 || nunsat < 10) break;
	}
	printf("(%d, %d), t = %d\n", dist, nunsat, tt);

      } /* end of  for(iex = 0; iex < Nex; iex++) */
    } /* end of for(p10 = 20; p10 <= 100; p10 = p10 + 10) */
  } /* end of   for(nvar = 100; nvar <= 200; nvar = nvar + 10) */
}
