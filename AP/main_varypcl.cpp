//
//  apdbifurcationmain.cpp
//
//  Main file for generating APD bifurcation diagrams.
//  .txt output files are organized as:
//  1st column:         pacing cycle length (PCL)
//  2nd column:         multiplication factor of Ito (called "itofac")
//  Remaining columns:  beat-to-beat series of action potential durations (APDs)
//  The number of remaining columns depends on the macros BEATS and REMOVEBEATS
//  where #remaining columns = (BEATS - REMOVEBEATS)
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

#ifndef PRECTYPE
#define PRECTYPE double
#endif

//#define TT
#ifdef OHara
#define TYPECELL OHaraCell
#define TYPECELLSTRING "OHara"
#include "Cells/OHaraCell.cpp"

#elif UCLA
#define TYPECELL UCLACell
#define TYPECELLSTRING "UCLA"
#include "Cells/UCLACell.cpp"

#elif TT
#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
#define gnaped 0.01
#include "Cells/TTCellIto.cpp"

#elif LR1
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"

#elif LR1SK
#define TYPECELL LR1CellItoSK
#define TYPECELLSTRING "LR1SK"
#include "Cells/LR1CellItoSK.cpp"

#elif LR2
#define TYPECELL LR2CellIto
#ifdef SK
#define TYPECELLSTRING "LR2SK"
#else
#define TYPECELLSTRING "LR2"
#endif
#include "Cells/LR2CellIto.cpp"

#elif Fenton
#define TYPECELL FentonCell
#define TYPECELLSTRING "Fenton"
#include "Cells/FentonCell.cpp"

#elif LR1CellMod
#define TYPECELL LR1CellMod
#define TYPECELLSTRING "LR1Mod"
#include "Cells/LR1CellItoMod.cpp"

#else
#define TYPECELL TTCellMod
#define TYPECELLSTRING "TTMod"
#include "Cells/TTCellMod.cpp"
#endif

#ifndef numitofac
#define numitofac 1
#endif

#ifndef numpcl
#define numpcl 1
#endif

#define NCELLS numitofac*numpcl

#ifndef BEATS
#define BEATS 20
#endif

#ifndef REMOVEBEATS
#define REMOVEBEATS 10
#endif

#ifndef stimt
#define stimt 100.0L
#endif


#ifndef threshold
#define threshold -75.0
#endif

#ifndef curfac
#define curfac itofac
#endif

int main(int argc, char *argv[])
{
    const double itofac = atof(argv[1]);
    const double pcl1 = atof(argv[2]);
    const double pcl2 = atof(argv[3]);
    const double dt = atof(argv[4]);
    
    TYPECELL<1>* h_cell;
    h_cell = new TYPECELL<1>();
    
    FILE *ap;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sap_PCL_%g_%g_ITO_%g.txt", TYPECELLSTRING, pcl1, pcl2, itofac);
    ap = fopen(fileSpec1, "w");
    
    FILE *saveconditions;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sVARS_PCL_%g_%g_ITO_%g.txt", TYPECELLSTRING, pcl1, pcl2, itofac);
    saveconditions = fopen(fileSpec2, "w");
#ifndef Fenton
    h_cell->curfac[0] = itofac;
#endif
#ifdef TT
    h_cell->itofac[0] = 0.0;
    h_cell->iskfac[0] = 0.0;
    h_cell->skh[0] = 0.0015;
#ifdef clampki
    h_cell->ki[0] = 138.0;
#endif
#ifdef clampnai
    h_cell->nai[0] = 12.0;
#endif
    //h_cell->ibarcafac[0] = 6e-4/(gcal);
    //h_cell->ikrfac[0] = 0.01/(gkr);
    //h_cell->iksfac[0] = 0.036/(gks);
    //h_cell->nacafac[0] = 2.0;
#endif
#ifdef LR2
    h_cell->skh[0] = 0.002;
    h_cell->itofac[0] = 0.0;
    h_cell->iupfac[0] = 3.0;
#endif
#ifdef UCLA
    h_cell->itofac[0] = 2.0;
    h_cell->iskfac[0] = 0.0;
    h_cell->skh[0] = 2.0;
#endif
#ifdef OHara
    h_cell->itofac[0] = 0.0;
    h_cell->iskfac[0] = 0.0;
    h_cell->inafac[0] = 2.0;
    h_cell->skh[0] = 0.01;
#endif
#ifdef TTMod
    h_cell->ibarcafac[0] = 5.0;
#endif
    h_cell->curfac[0] = itofac;
    //printf("Byte Size of Cells on Device: %lu\n", sizeof(APDBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    double t_savestart = 0;
    for (int i = 0; i < REMOVEBEATS; i++) {
       t_savestart += pcl1 + i*(pcl2 - pcl1)/(BEATS - 1);
    }
//#define t_savestart (pcl1*REMOVEBEATS)
    double t_save = t_savestart - 100; //200; //pcl*200 - 200;
    double other_tsave = 0;
#ifdef LR1
    double dzdv, dydv;
#endif
#ifdef LR1Modified
    double dzdv, dydv;
#endif
    int curbeat = -1;
    double vold, vnew;
    double stimtime = 0;
    while (curbeat < BEATS - 1) {
        vold = h_cell->v[0];
        if (t > stimtime - dt/4.0 && t < stimtime + stimduration - dt/4.0) {
            h_cell->stepdt(0,dt,stimulus);
#ifdef TTMod
            h_cell->ki[0] = 140.0;
            h_cell->nai[0] = 10.0;
#endif
        }
        else {
#ifdef TTMod
            h_cell->ki[0] = 140.0;
            h_cell->nai[0] = 10.0;
#endif
            h_cell->stepdt(0,dt,0);
#ifdef TTMod
            h_cell->ki[0] = 140.0;
            h_cell->nai[0] = 10.0;
#endif
        }
        t = t + dt;
        vnew = h_cell->v[0];
        if (vold >= threshold && vnew < threshold) {
           curbeat += 1;  
           stimtime += pcl1 + curbeat*(pcl2 - pcl1)/(BEATS - 1); 
           while (stimtime < t - dt/4.0) {
               stimtime += pcl1 + curbeat*(pcl2 - pcl1)/(BEATS - 1);
           }
           printf("END AP, TIME = %g, CURBEAT = %d, %g\n",(double)t,curbeat,pcl1 + curbeat*(pcl2 - pcl1)/(BEATS - 1)); 
        }
        if (t > t_save - dt/4) {
#ifdef TT
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0],h_cell->casr[0],h_cell->ki[0],h_cell->nai[0]);
#endif
#ifdef LR2
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0],h_cell->ki[0],h_cell->nai[0]);
#endif
#ifdef TTSK
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0],h_cell->casr[0],h_cell->ki[0],h_cell->nai[0]);
#endif
#ifdef TTMod
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0],h_cell->casr[0],h_cell->ki[0],h_cell->nai[0]);
#endif
#ifdef LR1
            fprintf(ap,"%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->comp_ito(0,dt,dzdv,dydv));
#endif
#ifdef LR1Modified
            fprintf(ap,"%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->z[0],h_cell->y[0]);
#endif
#ifdef LR1SK
            fprintf(ap,"%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0]);
#endif
#ifdef Fenton
            fprintf(ap,"%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->h[0]);
#endif
#ifdef UCLA
            h_cell->saveconditions(ap,0,0,(double)t - t_savestart);
	    //fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart, h_cell->v[0], h_cell->ci[0], h_cell->cs[0], h_cell->cp[0], h_cell->nai[0]);
#endif
#ifdef OHara
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart, h_cell->v[0], h_cell->cai[0], h_cell->cass[0], h_cell->nai[0]);
#endif
            t_save = t_save + 1;
        }
        if (t > stimtime - dt/4.0 & t < stimtime + dt/4.0) {
            printf("%g\n",(double)t);
        }
    }
    
    //h_cell->saveconditions(saveconditions,0,true,t - pcl*1900);

    delete h_cell;
    fclose(ap);
    fclose(saveconditions);
    return 0;
}


