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

#elif TP06
#define TYPECELL TP06Cell
#define TYPECELLSTRING "TP06"
#include "Cells/TP06Cell.cpp"

#elif TT
#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
//#define gnaped 0.01
#include "Cells/TTCellIto.cpp"

#elif TTHopen
#define TYPECELL TTCellHopen
#define TYPECELLSTRING "TTHopen"
#include "Cells/TTCellHopen.cpp"

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

#ifndef VARIABLE1
#define VARIABLE1 itofac
#endif

#ifndef VARIABLE2
#define VARIABLE2 iskfac
#endif

int main(int argc, char *argv[])
{
    const double var1 = atof(argv[1]);
    const double var2 = atof(argv[2]);
    const double pcl = atof(argv[3]);
    const double dt = atof(argv[4]);
    
    TYPECELL<1>* h_cell;
    h_cell = new TYPECELL<1>();
    
    FILE *ap;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sap_PCL_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl, var1, var2);
    ap = fopen(fileSpec1, "w");
    
    FILE *saveconditions;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sVARS_PCL_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl, var1, var2);
    saveconditions = fopen(fileSpec2, "w");
#ifndef Fenton
    h_cell->VARIABLE1[0] = var1;
    h_cell->VARIABLE2[0] = var2;
#endif
#ifdef LR1
    h_cell->itofac[0] = 0;
    h_cell->itoslowfac[0] = 0.0;
    h_cell->icalfac[0] = 1.15;
    h_cell->ikfac[0] = 1.0;
    h_cell->ikifac[0] = 2.2; //1.0;
    h_cell->yshift[0] = 8.0; //0.0; //-8.0;
    h_cell->tauXfac[0] = 5; //1;
#endif
#ifdef TT
    h_cell->itofac[0] = 0.0;
    h_cell->iskfac[0] = 0.0;
    h_cell->skh[0] = 0.0015;
    h_cell->ki[0] = 138.0;
    h_cell->nai[0] = 12.0;
    h_cell->ikrfac[0] = 0.0;
    h_cell->iksfac[0] = 0.5;
    h_cell->ibarcafac[0] = 1.0;
    //h_cell->ibarcafac[0] = 6e-4/(gcal);
    //h_cell->ikrfac[0] = 0.01/(gkr);
    //h_cell->iksfac[0] = 0.036/(gks);
    h_cell->nacafac[0] = 1.0; //4.0;
#endif

#ifdef TTHopen
    h_cell->itofac[0] = 0.0;
    h_cell->iskfac[0] = 0.0;
    h_cell->skh[0] = 0.0015;
    h_cell->ki[0] = 138.0;
    h_cell->nai[0] = 12.0;

#endif
#ifdef LR2
    h_cell->skh[0] = 0.002; //0.0025; //0.0025; //0.0008; //0.0025
    h_cell->iskfac[0] = 0.0;
    h_cell->itofac[0] = 0.0;
    h_cell->iupfac[0] = 1.0;
#endif
#ifdef UCLA
    h_cell->itofac[0] = 1.0;
    h_cell->iskfac[0] = 0.0;
    h_cell->skh[0] = 0.6;
    h_cell->nacafac[0] = 1.0; //0.3;
    h_cell->iksfac[0] = 1.0; //0.5;
    //h_cell->nacafac[0] = 0.15;
    //h_cell->iksfac[0] = 0.9;
#endif
#ifdef OHara
    h_cell->itofac[0] = 0.0;// 1.9;
    h_cell->iskfac[0] = 0.0; //20.0;
    h_cell->inafac[0] = 2.0;
    h_cell->skh[0] = 0.001; // 0.0005; //0.006; //0.01;
    h_cell->nai[0] = 7.96; //9.0;
    h_cell->nass[0] = h_cell->nai[0];
    h_cell->ki[0] = 146.0;
    h_cell->kss[0] = h_cell->ki[0];
    h_cell->yshift[0] = -8.0;
    h_cell->ikrfac[0] = 1.0;
    h_cell->iksfac[0] = 1.0; //1.0;
    h_cell->vssfac[0] = 1.0; //5.0;
    h_cell->alphask[0] = 0.0;
#endif
#ifdef TTMod
    h_cell->ibarcafac[0] = 5.0;
#endif
    h_cell->VARIABLE1[0] = var1;
    h_cell->VARIABLE2[0] = var2;
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
#elif OHara
	    ito = h_cell->comp_ito(0,dt,dxtos,dytos,dxtof,dytof,dxtof,dytof);
	    isk = h_cell->comp_isk(0,dt,dytof);
#else
	    //h_cell->comp_ito(0,dt,dxtos,dytos,dxtof,dytof,ito);
	    //h_cell->comp_isk(0,isk);
#endif
#endif
#ifdef TT
	    fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0], h_cell->nai[0], ito,isk);
#endif
#ifdef TP06
	    fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n", (double)t - t_savestart,h_cell->v[0],h_cell->cai[0],h_cell->cass[0], h_cell->nai[0], h_cell->ki[0]);
#endif
#ifdef LR2
	    double dxs1, dxs2, iks;
	    double dxsk, isk;
            h_cell->comp_iks(0,dt,dxs1,dxs2,iks);
	    h_cell->comp_isk(0,dt,dxsk,isk);
	    fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart,h_cell->v[0],h_cell->cai[0],ito,isk,h_cell->ki[0]);
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
	    double dum, ICaL, ICaNa, ICaK, INaCa_i, INaCa_ss, INa, INaL, Ikr, Iks;
	    h_cell->comp_ical(0,dt,ICaL, ICaNa, ICaK, dum, dum, dum, dum, dum, dum, dum, dum, dum);
	    h_cell->comp_inaca(0,INaCa_i, INaCa_ss);
        h_cell->comp_ina(0,dt,INa,INaL,dum,dum,dum,dum,dum,dum,dum,dum,dum);
        Ikr = h_cell->comp_ikr(0,dt,dum,dum);
        Iks = h_cell->comp_iks(0,dt,dum,dum);
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - t_savestart, h_cell->v[0], h_cell->cai[0], h_cell->cass[0], h_cell->CaMKt[0], h_cell->nai[0], h_cell->nass[0], h_cell->ki[0], h_cell->kss[0], h_cell->m[0],h_cell->hf[0], h_cell->hs[0], h_cell->j[0], h_cell->hsp[0], h_cell->jp[0], h_cell->mL[0], h_cell->hL[0], h_cell->hLp[0], h_cell->nca[0], h_cell->jca[0], h_cell->d[0], h_cell->ff[0], h_cell->fs[0], h_cell->fcaf[0], h_cell->fcas[0], h_cell->ffp[0], h_cell->fcafp[0], h_cell->fcas[0], h_cell->xrf[0], h_cell->xrs[0], h_cell->xs1[0], h_cell->xs2[0], h_cell->xk1[0], ICaL, ICaNa, ICaK, INaCa_i, INaCa_ss, INa, INaL, Ikr, Iks, ito, isk);
#endif
            //t_save = t_save + 1;
	    t_save = t_save + 1.0; //dt; // dt; //1;
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


