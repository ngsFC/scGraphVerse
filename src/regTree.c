/*******************************************************************
   This file is a modified version of file regTree.c contained in the R package 
   randomForest.
   

   Copyright (C) 2001-7 Leo Breiman, Adele Cutler and Merck & Co., Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*******************************************************************/

/******************************************************************
 * buildtree and findbestsplit routines translated from Leo's
 * original Fortran code.
 *
 *      copyright 1999 by leo Breiman
 *      this is free software and can be used for any purpose.
 *      It comes with no guarantee.
 *
 ******************************************************************/
#include <Rmath.h>
#include <R.h>
#include "rf.h"

void regTree(double *x, double *y, int mdim, int *sampsize,int nsample, int *lDaughter,
             int *rDaughter, double *upper, double *avnode, int *nodestatus, int nrnodes,
             int *treeSize, int nthsize, int mtry, int *mbest, int *cat,
            double *tgini, int *varUsed, int nclasses,double *weight) {
              
             
   int i, j, k, m, ncur, *jdex, *nodestart, *nodepop;
   int *ndstart, *ndend, *ndendl, *nodecnt, jstat, msplit, ind, s;
   double *d, *ss, *av, *decsplit, *ubest, *sumnode;


    sumnode = (double *) S_alloc(nclasses, sizeof(double));
    d = (double *) S_alloc(nclasses, sizeof(double));
    ubest = (double *) S_alloc(nclasses, sizeof(double));
    ndendl = (int *) Calloc(nclasses, int);
    nodecnt = (int *) Calloc(nclasses, int);
    decsplit = (double *) S_alloc(nclasses, sizeof(double));
    
    
    nodestart = (int *) Calloc(nclasses * nrnodes, int);
    nodepop   = (int *) Calloc(nclasses * nrnodes, int);
    av         = (double *) S_alloc(nclasses, sizeof(double)); /* average for each class */
    ss         = (double *) S_alloc(nclasses, sizeof(double)); /* standard deviation for each class */
    avnode     = (double *) S_alloc(nclasses * nrnodes, sizeof(double)); /* matrix average node x classes */
    ndstart = (int *) Calloc(nclasses, int);
    ndend = (int *) Calloc(nclasses, int);
    
    /* initialize some arrays for the tree */
    zeroInt(nodestatus, nrnodes * nclasses);
    zeroInt(nodestart, nrnodes * nclasses);
    zeroInt(nodepop, nrnodes * nclasses);
    
   /* zeroDouble(avnode, nrnodes); */

    jdex = (int *) Calloc(nclasses * nsample, int);
    
    ncur = 0;
    for (s = 0; s < nclasses; ++s) {
      for (i = 1; i <= nsample; ++i) { 
        jdex[(i-1) * nclasses + s] = i;
    }
    
    nodepop[0 + s * nrnodes] = sampsize[s]; /* number of sample in each node */
    nodestart[0 + s * nrnodes] = 0;
    nodestatus[0 + s * nrnodes] = NODE_TOSPLIT;
   
    av[s] = 0.0;
    ss[s] = 0.0; 
    
/*    printf("E22 = %lf <------\n",weight[0]);*/
  for (i = 0; i < sampsize[s]; ++i) {
    d[s] = y[(jdex[i * nclasses + s] - 1) * nclasses + s];
    ss[s] += i * (av[s] - d[s]) * (av[s] - d[s]) / (i + 1);
    av[s] = (i * av[s] + d[s]) / (i + 1);
  }

  avnode[nrnodes * s] = av[s];   

} 

/** START MAIN loop over nodes **/
  /* for (k = 0; k < nrnodes - 2; ++k) { */
  for (k = 0; k < nrnodes - 2; ++k) {
                                             
                                             if (k > ncur || ncur >= nrnodes - 2) break;
                                             
                                             ind = 0;
                                             for (s = 0; s < nclasses; ++s) {
                                               
                                               if (nodestatus[s * nrnodes + k] != NODE_TOSPLIT) {
                                                 ind = ind + 1; 
                                               }
                                             }
                                             
                                             if (ind == nclasses) continue;
                                             
                                             /*   Rprintf("nodi = %d\n",k);*/
                                               /* initialize for next call to findbestsplit */
                                               for (s = 0; s < nclasses; ++s) {
                                                 if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) {
                                                   ndstart[s] = nodestart[k + s * nrnodes];
                                                   ndend[s] = ndstart[s] + nodepop[k + s * nrnodes] - 1;
                                                   nodecnt[s] = nodepop[k + s * nrnodes];
                                                   sumnode[s] = nodecnt[s] * avnode[k + s * nrnodes];
                                                   decsplit[s] = 0.0;
                                                 }
                                               }
                                             
                                             jstat = 0;
                                             
                                             findBestSplit(x, jdex, y, mdim, nsample, ndstart, ndend, &msplit,
                                                           decsplit, ubest, ndendl, &jstat, mtry, sumnode,
                                                           nodecnt, cat, nclasses, nodestatus, nrnodes,k,weight);
                                             
    /*  Rprintf("jstat = %d\n",jstat);*/
                        
      if (jstat == 1) {
  		for (s = 0; s < nclasses; ++s) {nodestatus[k + s * nrnodes] = NODE_TERMINAL;}
       continue;
      }
       
         
        varUsed[msplit - 1] = 1;
		    mbest[k] = msplit;
          
        for (s = 0; s < nclasses; ++s) { 
            if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) { 
      
          upper[k * nclasses + s] = ubest[s];
		      tgini[(msplit - 1) + s * mdim] += decsplit[s];
          
         
          nodepop[s * nrnodes + ncur + 1] = ndendl[s] - ndstart[s] + 1;
  	      nodepop[s * nrnodes + ncur + 2] = ndend[s] - ndendl[s];
		      nodestart[s * nrnodes + ncur + 1] = ndstart[s];
		      nodestart[s * nrnodes + ncur + 2] = ndendl[s] + 1;
  	      nodestatus[s * nrnodes + k] = NODE_INTERIOR;
              
/*        Rprintf("pop_left = %d\n",nodepop[s * nrnodes + ncur + 1]);
        Rprintf("pop_right = %d\n",nodepop[s * nrnodes + ncur + 2]);*/

            }             
        }            

		
    for (s = 0; s < nclasses; ++s) { 
      
      if (nodestatus[s * nrnodes + k] == NODE_INTERIOR) {
      
    av[s] = 0.0;
		ss[s] = 0.0;
		for (j = ndstart[s]; j <= ndendl[s]; ++j) {
			d[s] = y[(jdex[j * nclasses + s] - 1) * nclasses + s];
  
			m = j - ndstart[s];
			ss[s] += m * (av[s] - d[s]) * (av[s] - d[s]) / (m + 1);
			av[s] = (m * av[s] + d[s]) / (m+1);
		}
		avnode[nrnodes * s + ncur + 1] = av[s];
    
		nodestatus[nrnodes * s + ncur + 1] = NODE_TOSPLIT;
		if (nodepop[nrnodes * s + ncur + 1] <= nthsize) {			nodestatus[nrnodes * s + ncur + 1] = NODE_TERMINAL;}

		av[s] = 0.0;
		ss[s] = 0.0;
		for (j = ndendl[s] + 1; j <= ndend[s]; ++j) {
     
     d[s] = y[(jdex[j * nclasses + s] - 1) * nclasses + s];
 			m = j - (ndendl[s] + 1);
			ss[s] += m * (av[s] - d[s]) * (av[s] - d[s]) / (m + 1);
			av[s] = (m * av[s] + d[s]) / (m + 1);
		}
		avnode[nrnodes * s + ncur + 2] = av[s];
		nodestatus[nrnodes * s + ncur + 2] = NODE_TOSPLIT;
		if (nodepop[nrnodes * s + ncur + 2] <= nthsize) {		nodestatus[nrnodes * s + ncur + 2] = NODE_TERMINAL;}
    }   
}
		lDaughter[k] = ncur + 1 + 1;
		rDaughter[k] = ncur + 2 + 1;
	 
		ncur += 2;
    }
   
   
    *treeSize = nrnodes;
    for (k = nrnodes - 1; k >= 0; --k) {
        ind = 0;
        for (s = 0; s < nclasses; ++s) { 
        if (nodestatus[nrnodes * s + k] == 0) { ind++; }
        
         if (nodestatus[nrnodes * s + k] == NODE_TOSPLIT) { nodestatus[nrnodes * s + k] = NODE_TERMINAL;}
        }
        
        if (ind == nclasses) (*treeSize)--;
        

    }




    
    Free(nodestart);
    Free(jdex);
    Free(nodepop);
    Free(ndendl);
    Free(nodecnt);
    Free(nodestart);
    Free(ndend);
      
}



/*--------------------------------------------------------------*/



void findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample,
		   int *ndstart, int *ndend, int *msplit, double *decsplit,
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double *sumnode, int *nodecnt, int *cat, int nclasses, int *nodestatus,  int nrnodes, int k,
       double *weight) {
         


    int last, ncat[32], icat[32], lc, nl, nr, npopl, npopr;
    int i, j, kv, l, *mind, *ncase, s;
    double *xt, *ut, *v, *yl, sumcat[32], avcat[32], tavcat[32], *ubestt;
    double crit, *critmax, *critvar, suml, sumr, d, critParent, sumcritvar, sumcritmax;

    
    critvar = (double *) Calloc(nclasses, double);
    critmax = (double *) Calloc(nclasses, double);
    ubestt = (double *) Calloc(nclasses, double);
  
    ut = (double *) Calloc(nclasses * nsample, double);
    xt = (double *) Calloc(nclasses * nsample, double);
    v  = (double *) Calloc(nsample, double);
    yl = (double *) Calloc(nsample * nclasses, double);
    mind  = (int *) Calloc(mdim, int);
    ncase = (int *) Calloc(nsample, int);
    zeroDouble(avcat, 32);
    zeroDouble(tavcat, 32);

    /* START BIG LOOP */
    *msplit = -1;
    
    
    last = mdim - 1;
    sumcritmax = 0.0;
    for (s=0; s < nclasses; ++s) { /* initialize */
      ubestt[s] = 0.0;
      critmax[s] = 0.0;
    }

    for (i=0; i < mdim; ++i) mind[i] = i;
    
    /** START MAIN loop over variables to consider for the split **/
    for (i = 0; i < mtry; ++i) { 


    for (s=0; s < nclasses; ++s) critvar[s] = 0.0; /* initialize */
    
    
    /* select variable */
  	j = (int) (unif_rand() * (last+1));
		kv = mind[j];
    swapInt(mind[j], mind[last]);
		last--;
		
		lc = cat[kv];
	
			
  		for (s = 0; s < nclasses; ++s) { /* xt: value of selected variable [kv] for each class each sample  */
       
       if (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
			for (j = ndstart[s]; j <= ndend[s]; ++j) {
				xt[s + j * nclasses ] = x[ (jdex[j * nclasses + s] - 1) * mdim * nclasses + s * mdim + kv ];
				yl[s + j * nclasses] = y[nclasses * (jdex[j * nclasses + s] - 1) + s ]; 
			}
       }
  		}
	  
              
    for (s = 0; s < nclasses; ++s) {     /** START loop over classes: for each variable kv & for each class find best threshold */          

     if (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
       
		 for (j = ndstart[s]; j <= ndend[s]; ++j) {
			 v[j] = xt[s  + j * nclasses ];  
       
		 }
     
    
		for (j = 1; j <= nsample; ++j) ncase[j - 1] = j;
		
    R_qsort_I(v, ncase, ndstart[s] + 1, ndend[s] + 1); /* sort v and return index ncase */
		
    
    if (v[ndstart[s]] >= v[ndend[s]]) continue;
    
		/* ncase(n)=case number of v nth from bottom */
		/* Start from the right and search to the left. */
		critParent = sumnode[s] * sumnode[s] / nodecnt[s];
		suml = 0.0;
		sumr = sumnode[s];
		npopl = 0;
		npopr = nodecnt[s];
		crit = 0.0;
		
		for (j = ndstart[s]; j <= ndend[s] - 1; ++j) { /* Search through the "gaps" in the x-variable. FIND THE TRESHOLD */
			d = yl[(ncase[j] - 1) * nclasses + s];	
    
    suml += d;
			sumr -= d;
			npopl++;
			npopr--;
			if (v[j] < v[j+1]) {
				crit = (suml * suml / npopl) + (sumr * sumr / npopr) -
					critParent;
				if (crit > critvar[s]) {
					ubestt[s] = (v[j] + v[j+1]) / 2.0;
					critvar[s] = crit;
				}
			}
		}
     }
    } /** END loop over classes */
    
		/* Find the best variables */
    sumcritvar=0.0;
    for (s = 0; s < nclasses; ++s) {
      if (nodestatus[s * nrnodes + k] == NODE_TOSPLIT) { sumcritvar = sumcritvar + weight[s]*critvar[s];
      } /* weight decrease in node impurity based on sample size */
    
    }
    
    sumcritvar=sumcritvar/nclasses;
    
    if (sumcritvar > sumcritmax) {
			
      for (s = 0; s < nclasses; ++s) {

       if (nodestatus[s * nrnodes + k]  == NODE_TOSPLIT) { 
        ubest[s] = ubestt[s];  /* ubest is the threshold */
        critmax[s] = critvar[s];
        
        for (j = ndstart[s]; j <= ndend[s]; ++j) { 	ut[s  + j * nclasses ] = xt[s  + j * nclasses ]; /* ut stores the values of variables used for the best split */
	     }
       }
      }
			*msplit = kv + 1; /* variable used for the best split */		
       sumcritmax = sumcritvar;       
		}
    } /** END MAIN loop over variables **/
    
    
  
/*     Rprintf("best split = %d\n",*msplit);*/
  
          
    if (*msplit != -1) { /* best variable has been found */
    
    /* divide samples into groups based on best splitting variable */
      for (s = 0; s < nclasses; ++s) {
        if (nodestatus[s * nrnodes + k]  == NODE_TOSPLIT) { 
        nl = ndstart[s];
      
        decsplit[s] = critmax[s];

       for (j =ndstart[s]; j <= ndend[s]; ++j) {
            if (ut[s + j * nclasses] <= ubest[s]) {
                nl++;
                ncase[nl-1] = jdex[s + j * nclasses]; 
            }
        }
        ndendl[s] = imax2(nl - 1, ndstart[s]);
        nr = ndendl[s] + 1;
        for (j = ndstart[s]; j <= ndend[s]; ++j) {
            if (ut[s + j * nclasses] > ubest[s]) {
                if (nr >= nsample) break;
                nr++;
                ncase[nr - 1] = jdex[s + j * nclasses];
            }
	    }
        if (ndendl[s] >= ndend[s]) ndendl[s] = ndend[s] - 1;
  
        for (j = ndstart[s]; j <= ndend[s]; ++j) jdex[s + j * nclasses] = ncase[j]; /* update jdex; left leave obs first */
		   
        }
      }
    } else *jstat = 1;      /* If best split can not be found, set to terminal node and return. */

  
	
    Free(ncase);
    Free(mind);
    Free(v);
    Free(yl);
    Free(xt);
    Free(ut);
    Free(critvar);
    Free(critmax);  
    Free(ubestt);
}

void zeroInt(int *x, int length) {
    memset(x, 0, length * sizeof(int));
}

void zeroDouble(double *x, int length) {
    memset(x, 0, length * sizeof(double));
}
/*==============================================================================================================================*/

void predictRegTree(double *x, int nsample, int mdim,
		    int *lDaughter, int *rDaughter, int *nodestatus,
                    double *ypred, double *split, double *nodepred,
                    int *splitVar, int treeSize, int *cat, int maxcat,
                    int *nodex, int nclasses, int nrnodes) {
    int i, j, k, m, *cbestsplit, s;
	unsigned int npack;

    /* decode the categorical splits */
    
     for (s = 0; s < nclasses; ++s) { /* loop over classes */
       
/*    if (maxcat > 1) {
        cbestsplit = (int *) Calloc(maxcat * treeSize, int);
        zeroInt(cbestsplit, maxcat * treeSize);
       
        for (i = 0; i < nrnodes; ++i) {
          
            if (nodestatus[s * nrnodes + i] != NODE_TERMINAL && cat[splitVar[i] - 1] > 1) {
                npack = (unsigned int) split[i * nclasses + s];
                
                for (j = 0; npack; npack >>= 1, ++j) {
                    cbestsplit[j + i*maxcat] = npack & 1;
                }
            }
        }
        }*/
    
   for (i = 0; i < nsample; ++i) {
	k = 0;
	while (nodestatus[s * nrnodes + k] != NODE_TERMINAL) { 
	    m = splitVar[k] - 1;
	    k = (x[m + i * mdim * nclasses + s * mdim] <= split[k * nclasses + s]) ?
		    lDaughter[k] - 1 : rDaughter[k] - 1;
	    } 

 
	ypred[i * nclasses + s] = nodepred[k + s * nrnodes];
 
	nodex[i * nclasses + s] = k + 1;
    } 
     }
    if (maxcat > 1) Free(cbestsplit);
}


void permuteOOB(int m, double *x, int *in, int nsample, int mdim, int s, int nclasses) {
/* Permute the OOB part of a variable in x.
 * Argument:
 *   m: the variable to be permuted
 *   x: the data matrix (variables in rows)
 *   in: vector indicating which case is OOB
 *   nsample: number of cases in the data
 *   mdim: number of variables in the data
 */
    double *tp, tmp;
    int i, last, k, nOOB = 0;

    tp = (double *) Calloc(nsample, double);

    for (i = 0; i < nsample; ++i) {
  	/* make a copy of the OOB part of the data into tp (for permuting) */
		if (in[i] == 0) {
            tp[nOOB] = x[m + i * mdim * nclasses + mdim * s];
            nOOB++;
        }
    }
    /* Permute tp */
    last = nOOB;
    for (i = 0; i < nOOB; ++i) {
		k = (int) last * unif_rand();
		tmp = tp[last - 1];
		tp[last - 1] = tp[k];
		tp[k] = tmp;
		last--;
    }

    /* Copy the permuted OOB data back into x. */
    nOOB = 0;
    for (i = 0; i < nsample; ++i) {
		if (in[i] == 0) {
            x[m + i * mdim * nclasses + mdim * s] = tp[nOOB];
            nOOB++;
		}
    }
    Free(tp);
}

