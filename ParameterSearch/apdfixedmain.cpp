#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

//#define TT
#ifdef TT
#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
#include "Cells/TTCellIto.cpp"
#elif TTSK
#define TYPECELL TTCellItoSK
#define TYPECELLSTRING "TTSK"
#include "Cells/TTCellItoSK.cpp"
#elif LR1
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"
#elif LR1SK
#define TYPECELL LR1CellItoSK
#define TYPECELLSTRING "LR1SK"
#include "Cells/LR1CellItoSK.cpp"
#elif LR1_taud
#define TYPECELL LR1CellIto_taud
#define TYPECELLSTRING "LR1_taud"
#include "Cells/LR1CellIto_taud.cpp"
#else
#define TYPECELL FentonCell
#define TYPECELLSTRING "Fenton"
#include "Cells/FentonCell.cpp"
#endif

#ifndef numvar1
#define numvar1 20
#endif

#ifndef numpcl
#define numpcl 1
#endif

#define NCELLS numvar1*numpcl

#ifndef BEATS
#define BEATS 134
#endif

#ifndef SAVEBEATS
#define SAVEBEATS 50
#endif

#ifndef stimt
#define stimt 100.0L
#endif

#ifndef curfac1
#define curfac1 ibarcafac
#endif

#ifndef TYPEVAR1
#define TYPEVAR1 "ICa"
#endif

#ifndef curfac2
#define curfac2 ikrfac
#endif

#ifndef TYPEVAR2
#define TYPEVAR2 "IKr"
#endif

#ifndef TOL
#define TOL 0.1
#endif

#include "APDBifurcation.cpp"

int main(int argc, char *argv[])
{
    const double minvar1 = atof(argv[1]);
    const double maxvar1 = atof(argv[2]);
    const double minpcl = atof(argv[3]);
    const double maxpcl = atof(argv[4]);
    const long double dt = atof(argv[5]);
    
    double byvar1 = (numvar1 > 1) ? (maxvar1 - minvar1)/(numvar1 - 1) : 1;
    double bypcl = (numpcl > 1) ? (maxpcl - minpcl)/(numpcl - 1) : 1;
    
    APDBifurcation<TYPECELL, 1, BEATS>* h_cells;
    
    h_cells = new APDBifurcation<TYPECELL, 1, BEATS>();
    
    FILE *paramsearch;
    char fileSpec1[100];
#ifdef TTSK
    snprintf(fileSpec1, 100, "%sPSearch_PCL_%g_%s_%g_%s.txt", TYPECELLSTRING, minpcl, TYPEVAR1, minvar1, TYPEVAR2);
#elif LR1SK
    snprintf(fileSpec1, 100, "%sPSearch_PCL_%g_%s_%g_%s.txt", TYPECELLSTRING, minpcl, TYPEVAR1, minvar1, TYPEVAR2);
#elif LR1_taud
    snprintf(fileSpec1, 100, "%sPSearch_PCL_%g_%s_%g_%s.txt", TYPECELLSTRING, minpcl, TYPEVAR1, minvar1, TYPEVAR2);
#else
    snprintf(fileSpec1, 100, "%sPSearch_PCL_%g_%s_%g_%s.txt", TYPECELLSTRING, minpcl, TYPEVAR1, minvar1, TYPEVAR2);
#endif
    paramsearch = fopen(fileSpec1, "w");
    
    FILE *paramsearch_apd;
    char fileSpec2[100];
#ifdef TTSK
    snprintf(fileSpec2, 100, "%sPSearch_APD_PCL_%g_%s_%g_%s.txt", TYPECELLSTRING, minpcl, TYPEVAR1, minvar1, TYPEVAR2);
#elif LR1SK
    snprintf(fileSpec2, 100, "%sPSearch_APD_PCL_%g_%s_%g_%s.txt", TYPECELLSTRING, minpcl, TYPEVAR1, minvar1, TYPEVAR2);
#elif LR1_taud
    snprintf(fileSpec2, 100, "%sPSearch_APD_PCL_%g_%s_%g_%s.txt", TYPECELLSTRING, minpcl, TYPEVAR1, minvar1, TYPEVAR2);
#else
    snprintf(fileSpec2, 100, "%sPSearch_APD_PCL_%g_%s_%g_%s.txt", TYPECELLSTRING, minpcl, TYPEVAR1, minvar1, TYPEVAR2);
#endif
    paramsearch_apd = fopen(fileSpec2, "w");
#ifdef TT
    FILE *paramsearch_cai;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%sPSearch_cai_PCL_%g_%s_%g_%s.txt", TYPECELLSTRING, minpcl, TYPEVAR1, minvar1, TYPEVAR2);
    paramsearch_cai = fopen(fileSpec3, "w");
#endif

    //printf("Byte Size of Cells on Device: %lu\n", sizeof(APDBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    double minapd;
    double maxapd;
    double saveapd;
    
    double curfac2min;
    double curfac2max;
    double curfac2cur;
    int counter;
    
    int j_sim;
    long double t;
    
    for (int i = 0; i < numpcl; i++) {
        h_cells = new APDBifurcation<TYPECELL, 1, BEATS>();
        h_cells->Cells.itofac[0] = 0.0;
        h_cells->pcls[0] = bypcl*i + minpcl;
        t = -stimt;
        h_cells->dobif(dt,t);
        minapd = h_cells->apds[BEATS - SAVEBEATS];
        maxapd = h_cells->apds[BEATS - SAVEBEATS];
        
        for (j_sim = BEATS - SAVEBEATS + 1; j_sim < BEATS; j_sim++) {
            if (h_cells->apds[j_sim] < minapd) {
                minapd = h_cells->apds[j_sim];
            }
            if (h_cells->apds[j_sim] > maxapd) {
                maxapd = h_cells->apds[j_sim];
            }
        }
        
        if (maxapd - minapd > 1) {
            printf("WARNING: PCL = %g LACKS STABILITY, MINAPD = %g, MAXAPD = %g",bypcl*i + minpcl,minapd,maxapd);
        }
        
        saveapd = maxapd;
        
        printf("PCL = %g, SAVEAPD = %g\n",bypcl*i + minpcl,saveapd);
        
        for (int j = 0; j < numvar1; j++) {
            h_cells = new APDBifurcation<TYPECELL, 1, BEATS>();
            h_cells->Cells.itofac[0] = 0.0;
            t = -stimt;
#ifdef LR1_taud
            h_cells->Cells.tauXfac[0] = 5.0;
            //h_cells->Cells.itofac[0] = 1.05;
#endif
            h_cells->Cells.curfac1[0] = byvar1*j + minvar1;
            h_cells->pcls[0] = bypcl*i + minpcl;
            
            h_cells->Cells.curfac2[0] = 1.0;
            curfac2min = 1.0;
            curfac2max = 1.0;
            curfac2cur = 1.0;
            
            h_cells->dobif(dt, t);            
 
            minapd = h_cells->apds[BEATS - SAVEBEATS];
            maxapd = h_cells->apds[BEATS - SAVEBEATS];
            
            for (j_sim = BEATS - SAVEBEATS + 1; j_sim < BEATS; j_sim++) {
                if (h_cells->apds[j_sim] < minapd) {
                    minapd = h_cells->apds[j_sim];
                }
                if (h_cells->apds[j_sim] > maxapd) {
                    maxapd = h_cells->apds[j_sim];
                }
            }
            
            if (maxapd - minapd > 1) {
                printf("WARNING: PCL = %g LACKS STABILITY, MINAPD = %g, MAXAPD = %g\n",bypcl*i + minpcl,minapd,maxapd);
            }
            
            printf("PCL = %g, CURFAC1 = %g, CURFAC2 = %g, SAVEAPD = %g, THISAPD = %g\n",bypcl*i + minpcl,byvar1*j + minvar1,curfac2cur,saveapd,maxapd);
            
            counter = 0;
            while (fabs(maxapd - saveapd) > TOL) {
                if (counter == 0) {
                    if (maxapd > saveapd) {
                        while (maxapd > saveapd) {
                            curfac2min = curfac2cur;
                            curfac2cur = curfac2cur*2.0;
                            curfac2max = curfac2cur;
                            h_cells = new APDBifurcation<TYPECELL, 1, BEATS>();
                            h_cells->Cells.itofac[0] = 0.0;
#ifdef LR1_taud
                            h_cells->Cells.tauXfac[0] = 5.0;
                            //h_cells->Cells.itofac[0] = 1.05;
#endif
                            h_cells->Cells.curfac1[0] = byvar1*j + minvar1;
                            h_cells->pcls[0] = bypcl*i + minpcl;
                        
                            h_cells->Cells.curfac2[0] = curfac2cur;
                        
                            h_cells->dobif(dt, t);
                        
                            minapd = h_cells->apds[BEATS - SAVEBEATS];
                            maxapd = h_cells->apds[BEATS - SAVEBEATS];
                        
                            for (j_sim = BEATS - SAVEBEATS + 1; j_sim < BEATS; j_sim++) {
                                if (h_cells->apds[j_sim] < minapd) {
                                    minapd = h_cells->apds[j_sim];
                                }
                                if (h_cells->apds[j_sim] > maxapd) {
                                    maxapd = h_cells->apds[j_sim];
                                }
                            }
                        
                            if (maxapd - minapd > 1) {
                                printf("WARNING: PCL = %g LACKS STABILITY, MINAPD = %g, MAXAPD = %g\n",bypcl*i + minpcl,minapd,maxapd);
                            }
                            printf("PCL = %g, CURFAC1 = %g, CURFAC2 = %g, SAVEAPD = %g, THISAPD = %g\n",bypcl*i + minpcl,byvar1*j + minvar1,curfac2cur,saveapd,maxapd);
                        }
                    }
                    else {
                        while (maxapd < saveapd) {
                            curfac2max = curfac2cur;
                            curfac2cur = curfac2cur/2.0;
                            curfac2min = curfac2cur;
                            h_cells = new APDBifurcation<TYPECELL, 1, BEATS>();
                            h_cells->Cells.itofac[0] = 0.0;
#ifdef LR1_taud
                            h_cells->Cells.tauXfac[0] = 5.0;
                            //h_cells->Cells.itofac[0] = 1.05;
#endif
                            h_cells->Cells.curfac1[0] = byvar1*j + minvar1;
                            h_cells->pcls[0] = bypcl*i + minpcl;
                            
                            h_cells->Cells.curfac2[0] = curfac2cur;
                            
                            h_cells->dobif(dt, t);
                            
                            minapd = h_cells->apds[BEATS - SAVEBEATS];
                            maxapd = h_cells->apds[BEATS - SAVEBEATS];
                            
                            for (j_sim = BEATS - SAVEBEATS + 1; j_sim < BEATS; j_sim++) {
                                if (h_cells->apds[j_sim] < minapd) {
                                    minapd = h_cells->apds[j_sim];
                                }
                                if (h_cells->apds[j_sim] > maxapd) {
                                    maxapd = h_cells->apds[j_sim];
                                }
                            }
                            
                            if (maxapd - minapd > 1) {
                                printf("WARNING: PCL = %g LACKS STABILITY, MINAPD = %g, MAXAPD = %g\n",bypcl*i + minpcl,minapd,maxapd);
                            }
                            printf("PCL = %g, CURFAC1 = %g, CURFAC2 = %g, SAVEAPD = %g, THISAPD = %g\n",bypcl*i + minpcl,byvar1*j + minvar1,curfac2cur,saveapd,maxapd);
                        }
                    }
                    counter = counter + 1;
                    printf("COUNTER: %d\n",counter);
                }
                else {
                    curfac2cur = (curfac2max + curfac2min)/2;
                    h_cells = new APDBifurcation<TYPECELL, 1, BEATS>();
                    h_cells->Cells.itofac[0] = 0.0;
#ifdef LR1_taud
                    h_cells->Cells.tauXfac[0] = 5.0;
                    //h_cells->Cells.itofac[0] = 1.05;
#endif
                    h_cells->Cells.curfac1[0] = byvar1*j + minvar1;
                    h_cells->pcls[0] = bypcl*i + minpcl;
                    
                    h_cells->Cells.curfac2[0] = curfac2cur;
                    
                    h_cells->dobif(dt, t);
                    
                    minapd = h_cells->apds[BEATS - SAVEBEATS];
                    maxapd = h_cells->apds[BEATS - SAVEBEATS];
                    
                    for (j_sim = BEATS - SAVEBEATS + 1; j_sim < BEATS; j_sim++) {
                        if (h_cells->apds[j_sim] < minapd) {
                            minapd = h_cells->apds[j_sim];
                        }
                        if (h_cells->apds[j_sim] > maxapd) {
                            maxapd = h_cells->apds[j_sim];
                        }
                    }
                    
                    if (maxapd - minapd > 1) {
                        printf("WARNING: PCL = %g LACKS STABILITY, MINAPD = %g, MAXAPD = %g\n",bypcl*i + minpcl,minapd,maxapd);
                    }
                    printf("PCL = %g, CURFAC1 = %g, CURFAC2 = %g, SAVEAPD = %g, THISAPD = %g\n",bypcl*i + minpcl,byvar1*j + minvar1,curfac2cur,saveapd,maxapd);
                    
                    if (maxapd > saveapd) {
                        curfac2min = curfac2cur;
                    }
                    else {
                        curfac2max = curfac2cur;
                    }
                    
                    counter = counter + 1;
                    printf("COUNTER: %d\n",counter);
                }
            }
            
            fprintf(paramsearch,"%g\t%g\t%g\t%g\t%g\n",h_cells->pcls[0],h_cells->Cells.curfac1[0],h_cells->Cells.curfac2[0],saveapd,maxapd);
	    printf("%g\t%g\t%g\t%g\n",h_cells->apds[BEATS - SAVEBEATS],h_cells->apds[BEATS - SAVEBEATS+1],h_cells->apds[BEATS - SAVEBEATS+2],h_cells->apds[BEATS - SAVEBEATS+3]);
            
            for (double pclnow = 4000; pclnow < 5001; pclnow = pclnow + 10) {
            
                h_cells = new APDBifurcation<TYPECELL, 1, BEATS>();
                h_cells->Cells.itofac[0] = 0.9;
#ifdef LR1_taud
                h_cells->Cells.tauXfac[0] = 5.0;
                //h_cells->Cells.itofac[0] = 1.05;
#endif
                h_cells->Cells.curfac1[0] = byvar1*j + minvar1;
                h_cells->pcls[0] = pclnow; //bypcl*i + minpcl;
            
                h_cells->Cells.curfac2[0] = curfac2cur;
            
                h_cells->dobif(dt, t);
            
                fprintf(paramsearch_apd,"%g\t%g\t%g\t%g",h_cells->pcls[0],h_cells->Cells.curfac1[0],h_cells->Cells.curfac2[0],h_cells->Cells.itofac[0]);
#ifdef TT
		fprintf(paramsearch_cai,"%g\t%g\t%g\t%g",h_cells->pcls[0],h_cells->Cells.curfac1[0],h_cells->Cells.curfac2[0],h_cells->Cells.itofac[0]);
#endif
                for (int i = 0; i < SAVEBEATS; i++) {
                    fprintf(paramsearch_apd,"\t%g",h_cells->apds[BEATS - SAVEBEATS + i]);
#ifdef TT
		    fprintf(paramsearch_cai,"\t%g\t%g",h_cells->cais[2*(BEATS - SAVEBEATS + i)],h_cells->cais[2*(BEATS - SAVEBEATS + i) + 1]);
#endif 
                }
                fprintf(paramsearch_apd,"\n");
#ifdef TT
		fprintf(paramsearch_cai,"\n");
#endif
            }
        }
    }
    
    delete h_cells;
    fclose(paramsearch);
    fclose(paramsearch_apd);
    fclose(paramsearch_cai);
    return 0;
}
