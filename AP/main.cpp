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
//#define gnaped 0.01
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

#ifndef curfac
#define curfac itofac
#endif

int main(int argc, char *argv[])
{
    const double itofac = atof(argv[1]);
    const double pcl = atof(argv[2]);
    const double dt = atof(argv[3]);
    
    TYPECELL<1>* h_cell;
    h_cell = new TYPECELL<1>();
    
    FILE *ap;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sap_PCL_%g_ITO_%g.txt", TYPECELLSTRING, pcl, itofac);
    ap = fopen(fileSpec1, "w");
    
    FILE *saveconditions;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sVARS_PCL_%g_ITO_%g.txt", TYPECELLSTRING, pcl, itofac);
    saveconditions = fopen(fileSpec2, "w");
#ifndef Fenton
    h_cell->curfac[0] = itofac;
#endif
#ifdef LR1
    h_cell->itofac[0] = 0;
    h_cell->itoslowfac[0] = 5.0;
    h_cell->icalfac[0] = 1.0;
    h_cell->ikfac[0] = 1.0;
    h_cell->ikifac[0] = 1.0;
    h_cell->yshift[0] = 0.0; //-8.0;
    h_cell->tauXfac[0] = 1;
#endif
#ifdef TT
    h_cell->itofac[0] = 0.0;
    h_cell->iskfac[0] = 0.0;
    h_cell->skh[0] = 0.0015;
    h_cell->ki[0] = 138.0;
    h_cell->nai[0] = 12.0;
    //h_cell->ibarcafac[0] = 6e-4/(gcal);
    //h_cell->ikrfac[0] = 0.01/(gkr);
    //h_cell->iksfac[0] = 0.036/(gks);
    h_cell->nacafac[0] = 4.0;
#endif
#ifdef LR2
    h_cell->skh[0] = 0.0008; //0.0025
    h_cell->iskfac[0] = 0.0;
    h_cell->itofac[0] = 0.0;
    h_cell->iupfac[0] = 3.0;
#endif
#ifdef UCLA
    h_cell->itofac[0] = 1.0;
    h_cell->iskfac[0] = 0.0;
    h_cell->skh[0] = 0.6;
    h_cell->nacafac[0] = 0.3;
    h_cell->iksfac[0] = 0.5;
    //h_cell->nacafac[0] = 0.15;
    //h_cell->iksfac[0] = 0.9;
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
#define t_savestart (pcl*REMOVEBEATS)
    double t_save = t_savestart - 100; //200; //pcl*200 - 200;
    double other_tsave = 0;

    double ito, dxtos, dytos, dxtof, dytof, isk;
    double stimtime = 0;
    while (t < pcl*BEATS + pcl - dt/2) {
        if (t > -dt/4.0 && t > stimtime - dt/4.0 && t < stimtime + stimduration - dt/4.0) {
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
	if (t > stimtime + pcl/2.0) {
	    stimtime = stimtime + pcl;
        }	
        if (t > t_save - dt/4) {
#ifndef LR1
#ifdef UCLA
	    ito = h_cell->comp_ito(0,dt,dxtos,dytos,dxtof,dytof);
	    isk = h_cell->comp_isk(0);
#elif TT
	    ito = h_cell->comp_ito(0,dt,dxtof,dytof);
	    isk = h_cell->comp_isk(0);
#else
	    h_cell->comp_ito(0,dt,dxtos,dytos,dxtof,dytof,ito);
	    h_cell->comp_isk(0,isk);
#endif
#endif
#ifdef TT
	    fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0],ito,isk);
#endif
#ifdef LR2
	    double dxs1, dxs2, iks;
            h_cell->comp_iks(0,dt,dxs1,dxs2,iks);
	    fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0],ito,isk,iks);
#endif
#ifdef TTSK
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0],ito,isk);
#endif
#ifdef TTMod
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0],ito,isk);
#endif
#ifdef LR1
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->comp_ito(0,dt,dxtos,dytos,dxtof,dytof),h_cell->xr[0],h_cell->ytos[0]);
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
			fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart, h_cell->v[0], h_cell->ci[0], h_cell->cs[0], h_cell->nai[0], h_cell->ytos[0], isk);
#endif
#ifdef OHara
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart, h_cell->v[0], h_cell->cai[0], h_cell->cass[0], ito, isk);
#endif
            t_save = t_save + 1;
        }
        if (t > other_tsave - dt/4) {
            printf("%g\t%g\n",(double)t,h_cell->v[0]);
            other_tsave = other_tsave + pcl;
        }
    }

    delete h_cell;
    fclose(ap);
    fclose(saveconditions);
    return 0;
}


