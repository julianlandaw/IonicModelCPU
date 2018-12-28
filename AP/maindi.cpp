//
//  maindi.cpp
//
//  Simulate APs with a fixed DI after every beat.
//
//  Created by Julian Landaw on 11/08/17.
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

#ifndef numdi
#define numdi 1
#endif

#define NCELLS numitofac*numdi

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

#ifndef threshold
#define threshold -75.0
#endif

int main(int argc, char *argv[])
{
    const double itofac = atof(argv[1]);
    const double di = atof(argv[2]);
    const double dt = atof(argv[3]);
    
    TYPECELL<1>* h_cell;
    h_cell = new TYPECELL<1>();
    
    FILE *ap;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sap_DI_%g_ITO_%g.txt", TYPECELLSTRING, di, itofac);
    ap = fopen(fileSpec1, "w");
/*
    FILE *saveconditions;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sVARS_DI_%g_ITO_%g.txt", TYPECELLSTRING, di, itofac);
    saveconditions = fopen(fileSpec2, "w");
*/
    FILE *apds;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%sAPD_DI_%g_ITO_%g.txt", TYPECELLSTRING, di, itofac);
    apds = fopen(fileSpec3, "w");
    fprintf(apds,"%g\t%g",di,itofac);
#ifndef Fenton
    h_cell->curfac[0] = itofac;
#endif
#ifdef TT
    h_cell->ibarcafac[0] = 6e-4/(gcal);
    h_cell->ikrfac[0] = 0.01/(gkr);
    h_cell->iksfac[0] = 0.036/(gks);
    //h_cell->nacafac[0] = 2.0;
#endif
#ifdef TTMod
    h_cell->ibarcafac[0] = 5.0;
#endif
    //printf("Byte Size of Cells on Device: %lu\n", sizeof(APDBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    double t_save = 0;
    double other_tsave = 0;
    double t_endap = di;
    double vold;
    double vnew;
    //int flag = 0;
#ifdef LR1
    double dzdv, dydv;
#endif
#ifdef LR1Modified
    double dzdv, dydv;
#endif
    while (t < (di + 600)*400 - dt/2) {
#ifdef TT
        h_cell->ki[0] = 138.0;
        h_cell->nai[0] = 12.0;
#endif
#ifdef TTSK
        h_cell->ki[0] = 138.0;
        h_cell->nai[0] = 12.0;
#endif
       // h_cell->stepdt(0,dt,stimulus);
#ifdef TT
        h_cell->ki[0] = 138.0;
        h_cell->nai[0] = 12.0;
#endif
#ifdef TTSK
        h_cell->ki[0] = 138.0;
        h_cell->nai[0] = 12.0;
#endif
#ifdef TTMod
        h_cell->ki[0] = 140.0;
        h_cell->nai[0] = 10.0;
#endif
        vold = h_cell->v[0];
        if (t > -dt/4.0 && (t > (t_endap + di - dt/4.0)) && (t < (t_endap + di + stimduration - dt/4.0))) {
            h_cell->stepdt(0,dt,stimulus);
        }
        else {
            h_cell->stepdt(0,dt,0);
        }
        t = t + dt;
        vnew = h_cell->v[0];
        if (vnew < threshold && vold > threshold) {
            fprintf(apds,"\t%g",(double)t - t_endap - di);
            t_endap = t;
        }
        if (t > t_save - dt/4) {
#ifdef TT
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t,h_cell->v[0],h_cell->cai[0],h_cell->casr[0],h_cell->ki[0],h_cell->nai[0]);
#endif
#ifdef TTSK
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t,h_cell->v[0],h_cell->cai[0],h_cell->casr[0],h_cell->ki[0],h_cell->nai[0]);
#endif
#ifdef TTMod
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\t%g\n",(double)t,h_cell->v[0],h_cell->cai[0],h_cell->casr[0],h_cell->ki[0],h_cell->nai[0]);
#endif
#ifdef LR1
            fprintf(ap,"%g\t%g\t%g\n",(double)t,h_cell->v[0],h_cell->comp_ito(0,dt,dzdv,dydv));
#endif
#ifdef LR1Modified
            fprintf(ap,"%g\t%g\t%g\t%g\n",(double)t,h_cell->v[0],h_cell->z[0],h_cell->y[0]);
#endif
#ifdef LR1SK
            fprintf(ap,"%g\t%g\t%g\n",(double)t,h_cell->v[0],h_cell->cai[0]);
#endif
#ifdef Fenton
            fprintf(ap,"%g\t%g\t%g\n",(double)t,h_cell->v[0],h_cell->h[0]);
#endif
            t_save = t_save + 1;
        }
        if (t > other_tsave - dt/4) {
            printf("%g\t%g\n",(double)t, t_endap);
            other_tsave = other_tsave + 10;
        }
    }
    
    //h_cell->saveconditions(saveconditions,0,true,t - pcl*1900);

    delete h_cell;
    fclose(ap);
    //fclose(saveconditions);
    fclose(apds);
    return 0;
}


