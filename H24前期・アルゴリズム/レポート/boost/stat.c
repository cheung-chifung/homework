/* program:  Preprocess the mushroom dataset as an example of AdaBoost
   file:     preprocess.c
   update:   2010/6/25
*/

#include <stdio.h>
#include <stdlib.h>
#include "mt.h"

#define AttrN      21      /* # of attributes */
#define MaxNumData 10000   /* max # of data items */
#define MaxNumTest 3000    /* max # of test data items */

int main(int argc, char *argv[])
{
  /* for data */
  char data[MaxNumData][AttrN+1];
  double f[MaxNumData];
  char test[MaxNumTest][AttrN+1];
  double w[MaxNumTest];
  int m;

  /* weak hyptotheses */
  int Ha[119]
  = {  1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
       3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
       4, 4, 5, 5, 5, 5, 5, 5, 5, 5,
       5, 6, 6, 6, 6, 7, 7, 7, 8, 8,
       9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
       9, 9,10,10,11,11,11,11,12,12,
      12,12,13,13,13,13,13,13,13,13,
      13,14,14,14,14,14,14,14,14,14,
      15,15,16,16,16,16,17,17,17,18,
      18,18,18,18,18,18,18,19,19,19,
      19,19,19,19,19,19,20,20,20,20,
      20,20,21,21,21,21,21,21,21 };
  char Hv[119]
  = {  'b','c','x','f','k','s','f','g','y','s',
       'n','b','c','g','r','p','u','e','w','y',
       't','f','a','l','c','y','f','m','n','p',
       's','a','d','f','n','c','w','d','b','n',
       'k','n','b','h','g','r','o','p','u','e',
       'w','y','e','t','f','y','k','s','f','y',
       'k','s','n','b','c','g','o','p','e','w',
       'y','n','b','c','g','o','p','e','w','y',
       'p','u','n','o','w','y','n','o','t','c',
       'e','f','l','n','p','s','z','k','n','b',
       'h','r','o','u','w','y','a','c','n','s',
       'v','y','g','l','m','p','u','w','d' };

  /* stat data */
  int num_poison;
  int numA4, numA8, numA10, numA15;
  int numH20, numH38, numH52, numH80;
  double prob_poison;
  double errA4, errA8, errA10, errA15;
  double advH20, advH20N, advH38, advH52, advH80;

  /* random seeds */
  int rand_seed;

  /* data files */
  char fnameI[256];
  char instr[600];
  FILE *fpI;
  int num;

  /* work */
  int i, j;
  char inputchar;

  /****************/
  /* Program Body */
  /****************/

  /* read data from the file */
  sprintf(fnameI, "data.txt");
  fpI = fopen(fnameI, "r");
  if(!fpI){
    printf("File open error! %s\n", fnameI);
    exit -1;
  }
  num = 0;
  while(1) {
    fscanf(fpI, "%c", &inputchar);
    if(inputchar == 'x') break; /* end to read */
    data[num][0] = inputchar;
    if(inputchar == 'p') f[num] = 1.0; else f[num] = -1.0;
    /*
     *  p -> posionous
     *  otherwises -> edible
     */
    for(j = 1; j <= AttrN; j++) {
      fscanf(fpI, "%c", &inputchar);  /* ignore one blank */
      fscanf(fpI, "%c", &inputchar);
      data[num][j] = inputchar;
    }
    fgets(instr, 600+1, fpI);         /* skip to the next line */
    num++;
  }
  printf("num of data items = %d\n", num);

  /* prepare test data and their weights */
  m = 1000;     /* m <= MaxNumTest */
  for(i = 0; i < m; i++) {
    for(j = 0; j <= AttrN; j++)
      test[i][j] = data[i][j];
    w[i] = 1 / (double)m;
  }

  /* stat review */
  num_poison = 0;
  numA4 = numA8 = numA15 = 0;
  numH20 = numH38 = numH80 = 0;
  prob_poison = 0;
  errA4 = errA8 = errA15 = 0;
  advH20 = advH20N = advH38 = advH80 = 0;
  for(i = 0; i < m; i++) {
    if(test[i][0] == 'p')
      { num_poison++; prob_poison = prob_poison + w[i]; }
    if(test[i][4] == 't') {
      numA4++;
      if(test[i][0] == 'e') errA4 = errA4 + w[i];
    }
    else {
      if(test[i][0] == 'p') errA4 = errA4 + w[i];
    }
    if(test[i][Ha[20]] == Hv[20]) {
      numH20++;
      advH20 = advH20 + f[i] * w[i];
      advH20N = advH20N - f[i] * w[i];
    }
    else {
      advH20 = advH20 - f[i] * w[i];
      advH20N = advH20N + f[i] * w[i];
    }
    if(test[i][8] == 'b') {
      numA8++;
      if(test[i][0] == 'e') errA8 = errA8 + w[i];
    }
    else {
      if(test[i][0] == 'p') errA8 = errA8 + w[i];
    }
    if(test[i][Ha[38]] == Hv[38]) {
      numH38++;
      advH38 = advH38 + f[i] * w[i];
    }
    else {
      advH38 = advH38 - f[i] * w[i];
    }
    if(test[i][15] == 'p') {
      numA15++;
      if(test[i][0] == 'e') errA15 = errA15 + w[i];
    }
    else {
      if(test[i][0] == 'p') errA15 = errA15 + w[i];
    }
    if(test[i][Ha[80]] == Hv[80]) {
      numH80++;
      advH80 = advH80 + f[i] * w[i];
    }
    else {
      advH80 = advH80 - f[i] * w[i];
    }
  }

  printf("stat review in the test data (%d items):\n", m);
  printf("   # of poisonous mushrooms: number = %d, prob. = %f\n",
	 num_poison, prob_poison);
  /* for debugging
  printf("   #4 attr.  = t: number = %d, error prob. = %f\n",
	 numA4, errA4);
  */
  printf("   #4 attr.  = t: number = %d, advantage   = %f\n",
	 numH20, advH20 / 2.0);
  printf("   #4 attr. != t: number = %d, advantage   = %f\n",
	 m - numH20, advH20N / 2.0);
  /* for debugging
  printf("   #8 attr.  = b: number = %d, error prob. = %f\n",
	 numA8, errA8);
  */
  printf("   #8 attr.  = b: number = %d, advantage   = %f\n",
	 numH38, advH38 / 2.0);
  /* for debugging
  printf("   #15 attr. = p: number = %d, error prob. = %f\n",
	 numA15, errA15);
  */
  printf("   #15 attr. = p: number = %d, advantage   = %f\n",
	 numH80, advH80 / 2.0);

  return 0;
}
