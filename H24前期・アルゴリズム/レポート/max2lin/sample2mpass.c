/* sample program: for average-case analysis
                   of message passing solver of MAX-2LIN
   file:     sample2mpass.c
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

#define Nex       5       /* # of example instances per each param. setting */

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

  /* for heuristics (message passing algo.) */
  int A[MAXvar], Atmp[MAXvar];
  int tt, maxtt;
  int message, tmessage;
  int vtmp;
  int dist, nunsat;

  /* work var.s */
  int iex;
  int i, j, v, eq, val;

  /****************/
  /* Program Body */
  /****************/

  rand_seed = 1234;
  mt=0;
  sgenrand(rand_seed);

  printf("# n = 100, 110, ..., 200\n");
  printf("# P = 0.2/n, 0.3/n, ..., 10/n, Q = 0\n");
  printf("# assignment update steps = 1.5 * n\n");

  for(nvar = 400; nvar <= 400; nvar = nvar + 10) {
    printf("> n = %d\n", nvar);

    for(p10 = 100; p10 <= 100; p10 = p10 + 10) {
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
	/* simple message passing        */
	/*********************************/
	maxtt = 1.5 * (float)nvar;
  // A[]を初期化
	for(v = 2; v <= nvar; v++) A[v] = 0;
	A[1] = +1;

	for(tt = 1; tt <= maxtt; tt++) {

	  /* compute new assignments */
	  for(v = 2; v <= nvar; v++) {

	    /* collect messages from all related eqs */
	    tmessage = 0;
	    for(i = 0; i < nV2E[v]; i++) {
	      eq = V2E[v][i];
	      for(j = 0; j < 2; j++) {
		vtmp = E[eq][j];
		if(vtmp != v) {
		  message = A[vtmp] * R[eq];
		  break;
		}
	      }
	      tmessage = tmessage + message;
	    }

	    /* determine the new assignment for v */
	    if(tmessage == 0) Atmp[v] = 0;
	    else {
	      if(tmessage > 0) Atmp[v] = +1;
	      else             Atmp[v] = -1;
	    }
	  }

	  /* update assignments */
	  dist = 0;
	  for(v = 2; v <= nvar; v++) {
	    A[v] = Atmp[v];
	    if(Atmp[v] <= 0) dist++;
	  }
	  //printf("dist:(%d) t:(%d)\n", dist, tt);
	  if(dist < 2) break;
	}
	printf("t = %d, dist = %d\n", tt, dist);
  mt=mt+tt;
  printf("mt=%ld\n",mt);
      } /* end of  for(iex = 0; iex < Nex; iex++) */
    } /* end of for(p10 = 20; p10 <= 100; p10 = p10 + 10) */
  } /* end of   for(nvar = 100; nvar <= 200; nvar = nvar + 10) */
}
