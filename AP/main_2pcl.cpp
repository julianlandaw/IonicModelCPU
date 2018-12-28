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

//#define TT
#ifdef TT
#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
#include "Cells/TTCellIto.cpp"
#elif LR1
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"
#elif LR1SK
#define TYPECELL LR1CellItoSK
#define TYPECELLSTRING "LR1SK"
#include "Cells/LR1CellItoSK.cpp"
#else
#define TYPECELL FentonCell
#define TYPECELLSTRING "Fenton"
#include "Cells/FentonCell.cpp"
#endif

#ifndef numitofac
#define numitofac 1
#endif

#ifndef numpcl
#define numpcl 1
#endif

#define NCELLS numitofac*numpcl

#ifndef BEATS
#define BEATS 50
#endif

#ifndef REMOVEBEATS
#define REMOVEBEATS 1000
#endif

#ifndef stimt
#define stimt 100.0L
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
    snprintf(fileSpec1, 100, "%sap_PCLs_%g_%g_ITO_%g.txt", TYPECELLSTRING, pcl1, pcl2, itofac);
    ap = fopen(fileSpec1, "w");
    
    FILE *saveconditions;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sVARS_PCLs_%g_%g_ITO_%g.txt", TYPECELLSTRING, pcl1, pcl2, itofac);
    saveconditions = fopen(fileSpec2, "w");
#ifndef Fenton
    h_cell->itofac[0] = itofac;
#endif
    //printf("Byte Size of Cells on Device: %lu\n", sizeof(APDBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    double t_save = 0; //pcl*200 - 200;
    double other_tsave = 0;
#ifdef LR1
    double dzdv, dydv;
#endif
    while (t < pcl1*(BEATS + REMOVEBEATS) - dt/2) {
#ifdef TT
        h_cell->ki[0] = 138;
        h_cell->nai[0] = 12;
#endif
        if (t > -dt/4.0 && fmod((t + dt/4.0), pcl1) < (stimduration - dt/4.0)) {
            h_cell->stepdt(0,dt,stimulus);
        }
        else {
            h_cell->stepdt(0,dt,0);
        }
        t = t + dt;
        if (t > t_save - dt/4) {
            t_save = t_save + 1;
            if (t > pcl1*REMOVEBEATS - 10.0) {
#ifdef TT
                fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - pcl1*REMOVEBEATS,h_cell->v[0],h_cell->cai[0],h_cell->casr[0],h_cell->ki[0],h_cell->nai[0]);
#endif
#ifdef LR1
                fprintf(ap,"%g\t%g\t%g\n",(double)t - pcl*200,h_cell->v[0],h_cell->comp_ito(0,dt,dzdv,dydv));
#endif
#ifdef LR1SK
                fprintf(ap,"%g\t%g\t%g\n",(double)t - pcl*200,h_cell->v[0],h_cell->cai[0]);
#endif
#ifdef Fenton
                fprintf(ap,"%g\t%g\t%g\n",(double)t - pcl*200,h_cell->v[0],h_cell->h[0]);
#endif
            }
        }
        if (t > other_tsave - dt/4) {
            printf("%g\n",(double)t);
            other_tsave = other_tsave + 10;
        }
    }
    long double lasttime = t;
    
    while (t - lasttime < pcl2*BEATS - dt/2) {
#ifdef TT
        h_cell->ki[0] = 138;
        h_cell->nai[0] = 12;
#endif
        if (fmod((t - lasttime + dt/4.0), pcl2) < (stimduration - dt/4.0)) {
            h_cell->stepdt(0,dt,stimulus);
        }
        else {
            h_cell->stepdt(0,dt,0);
        }
        t = t + dt;
        if (t > t_save - dt/4) {
#ifdef TT
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t - pcl1*REMOVEBEATS, h_cell->v[0],h_cell->cai[0],h_cell->casr[0],h_cell->ki[0],h_cell->nai[0]);
#endif
#ifdef LR1
            fprintf(ap,"%g\t%g\t%g\n",(double)t - pcl*200,h_cell->v[0],h_cell->comp_ito(0,dt,dzdv,dydv));
#endif
#ifdef LR1SK
            fprintf(ap,"%g\t%g\t%g\n",(double)t - pcl*200,h_cell->v[0],h_cell->cai[0]);
#endif
#ifdef Fenton
            fprintf(ap,"%g\t%g\t%g\n",(double)t - pcl*200,h_cell->v[0],h_cell->h[0]);
#endif
            t_save = t_save + 1;
        }
        if (t > other_tsave - dt/4) {
            printf("%g\n",(double)t);
            other_tsave = other_tsave + 10;
        }
    }
    
    //h_cell->saveconditions(saveconditions,0,true,t - pcl*1900);

    delete h_cell;
    fclose(ap);
    fclose(saveconditions);
    return 0;
}


