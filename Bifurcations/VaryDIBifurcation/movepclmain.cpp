//
//  dibifurcationmain.cpp
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
//#define gnaped 0.01
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

#ifndef numvar
#define numvar 1
#endif

#ifndef numdi
#define numdi 1
#endif

#ifndef numalpha
#define numalpha 1
#endif

#define NCELLS numvar*numdi

#ifndef BEATS
#define BEATS 20
#endif

#ifndef REMOVEBEATS
#define REMOVEBEATS 10
#endif

#ifndef stimt
#define stimt 100.0L
#endif

#ifndef VARIABLE
#define VARIABLE itofac
#endif

#include "VaryDIBifurcation.cpp"

int main(int argc, char *argv[])
{
    const double itofac = atof(argv[1]);
    const double di0 = atof(argv[2]);
    const double alphanew = atof(argv[3]);
    const int numapdsave = atoi(argv[4]);
    const long double dt = atof(argv[5]);
    
    VaryDIBifurcation<TYPECELL, NCELLS, BEATS>* h_cells;
    VaryDIBifurcation<TYPECELL, NCELLS, BEATS>* h_cellsnew;
    VaryDIBifurcation<TYPECELL, NCELLS, BEATS>* h_cellsnew2;
    
    h_cells = new VaryDIBifurcation<TYPECELL, NCELLS, BEATS>();
    h_cellsnew = new VaryDIBifurcation<TYPECELL, NCELLS, BEATS>();
    h_cellsnew2 = new VaryDIBifurcation<TYPECELL, NCELLS, BEATS>();

    FILE *allbifs;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sbifsDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    allbifs = fopen(fileSpec1, "w");
    
    
#ifdef LR1
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    xrbifs = fopen(fileSpec2, "w");
#endif
    
#ifdef LR1_taud
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrDI_%g_taud_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    xrbifs = fopen(fileSpec2, "w");
#endif
    
#ifdef TT
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    naibifs = fopen(fileSpec5, "w");
#endif
    
#ifdef TTMod
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiDI_%g_ITO_%g_ALPHA_%g_APDs_%d.txt", TYPECELLSTRING, di0, itofac, alphanew, numapdsave);
    naibifs = fopen(fileSpec5, "w");
#endif
    
    int index = 0;
#ifdef LR1
#ifndef icafac
#define icafac 1.0
#endif
#ifndef taufac
#define taufac 1.0
#endif                
    h_cells->Cells.ibarcafac[index] = icafac;
    h_cells->Cells.tauXfac[index] = taufac;
#endif
            
#ifdef LR1_taud
    h_cells->Cells.tauXfac[index] = 5.0;
    h_cells->Cells.itofac[index] = 1.05;
#endif
#ifdef TTMod
    h_cells->Cells.ibarcafac[index] = 5.0;
#endif

#ifdef TT
#ifdef clampnai
    h_cells->Cells.nai[index] = 12.0;
#endif
#ifdef clampki
    h_cells->Cells.ki[index] = 138.0;
#endif
#endif
    h_cells->Cells.VARIABLE[index] = itofac;
    h_cells->d0[index] = di0;
    h_cells->apdfac[index] = 0;
    h_cells->numapdsave = numapdsave;
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(VaryDIBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    h_cells->dobif(dt, t);
    
    index = 0;           
    fprintf(allbifs, "%g\t%g\t%g\t%d", h_cells->d0[0], h_cells->apdfac[0], itofac, h_cells->numapdsave);
    for (int l = REMOVEBEATS; l < BEATS; l++) {
         fprintf(allbifs, "\t%g", h_cells->apds[BEATS*(index) + l]);
    }
    fprintf(allbifs, "\n");

    double apd_min = h_cells->apds[REMOVEBEATS]; 
    double apd_max = h_cells->apds[REMOVEBEATS];
    double apd_temp;
    for (int i = REMOVEBEATS + 1; i < BEATS; i++) {
        apd_temp = h_cells->apds[i];
        if (apd_temp > apd_max) {
            apd_max = apd_temp;
        }
        if (apd_temp < apd_min) {
            apd_min = apd_temp;
        }
    }
    if (apd_max - apd_min < 5) {
        printf("STABLE!\n");
    }
    double apd_star = h_cells->apds[BEATS - 1];    
    
    h_cellsnew->Cells.setcell(0,&(h_cells->Cells));

    if (numapdsave == 0) {    
        h_cellsnew->d0[index] = di0 + apd_star*h_cells->apdfac[index] - alphanew*apd_star;
    }
    else {
        h_cellsnew->d0[index] = di0;
    }
    h_cellsnew->apdfac[index] = alphanew;
    h_cellsnew->numapdsave = numapdsave;
    t = -stimt;

    h_cellsnew->dobif(dt, t);
    
    fprintf(allbifs, "%g\t%g\t%g\t%d", h_cellsnew->d0[0], h_cellsnew->apdfac[0], itofac, h_cellsnew->numapdsave);
    for (int l = REMOVEBEATS; l < BEATS; l++) {
         fprintf(allbifs, "\t%g", h_cellsnew->apds[BEATS*(index) + l]);
    }
    fprintf(allbifs, "\n");
    
    for (int i = 1; i < 200; i = i + 1) {
        printf("\n\n%d\n\n",i);
        h_cellsnew->d0[index] = h_cellsnew->d0[index] - 20;
        h_cellsnew->startapds[index] = 0.0;
        h_cellsnew->st[index] = 0.0;
        h_cellsnew->curbeat[index] = -1;
        h_cellsnew->stimtime[index] = 0.0;
        h_cellsnew->donebifurcation = 0;
        
        t = -stimt;
        h_cellsnew->dobif(dt,t);

        fprintf(allbifs, "%g\t%g\t%g\t%d", h_cellsnew->d0[0], h_cellsnew->apdfac[0], itofac, h_cellsnew->numapdsave);
        for (int l = REMOVEBEATS; l < BEATS; l++) {
             fprintf(allbifs, "\t%g", h_cellsnew->apds[BEATS*(index) + l]);
        }
        fprintf(allbifs, "\n");
    }
            
            
#ifdef LR1
                fprintf(xrbifs, "%g\t%g\t%g\t%d", di0, alphanew, itofac, numapdsave);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(index) + l]);
                }
                fprintf(xrbifs, "\n");
#endif
            
#ifdef LR1_taud
                fprintf(xrbifs, "%g\t%g\t%g", di0, alphanew, itofac, numapdsave);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(index) + l]);
                }
                fprintf(xrbifs, "\n");
#endif
#ifdef TT
                fprintf(caibifs, "%g\t%g\t%g\t%d", di0, alphanew, itofac, numapdsave);
                fprintf(casrbifs, "%g\t%g\t%g\t%d", di0, alphanew, itofac, numapdsave);
                fprintf(kibifs, "%g\t%g\t%g\t%d", di0, alphanew, itofac, numapdsave);
                fprintf(naibifs, "%g\t%g\t%g\t%d", di0, alphanew, itofac, numapdsave);
            
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(caibifs, "\t%g", h_cells->cais[2*BEATS*(index) + l]);
                    fprintf(casrbifs, "\t%g", h_cells->casrs[2*BEATS*(index) + l]);
                    fprintf(kibifs, "\t%g", h_cells->kis[2*BEATS*(index) + l]);
                    fprintf(naibifs, "\t%g", h_cells->nais[2*BEATS*(index) + l]);
                }
                fprintf(caibifs, "\n");
                fprintf(casrbifs, "\n");
                fprintf(kibifs, "\n");
                fprintf(naibifs, "\n");
#endif
#ifdef TTMod
                fprintf(caibifs, "%g\t%g\t%g", di0, alphanew, itofac, numapdsave);
                fprintf(casrbifs, "%g\t%g\t%g", di0, alphanew, itofac, numapdsave);
                fprintf(kibifs, "%g\t%g\t%g", di0, alphanew, itofac, numapdsave);
                fprintf(naibifs, "%g\t%g\t%g", di0, alphanew, itofac, numapdsave);
            
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(caibifs, "\t%g", h_cells->cais[2*BEATS*(index) + l]);
                    fprintf(casrbifs, "\t%g", h_cells->casrs[2*BEATS*(index) + l]);
                    fprintf(kibifs, "\t%g", h_cells->kis[2*BEATS*(index) + l]);
                    fprintf(naibifs, "\t%g", h_cells->nais[2*BEATS*(index) + l]);
                }
                fprintf(caibifs, "\n");
                fprintf(casrbifs, "\n");
                fprintf(kibifs, "\n");
                fprintf(naibifs, "\n");
#endif

    delete h_cells;
    delete h_cellsnew;
    fclose(allbifs);
#ifdef LR1
    fclose(xrbifs);
#endif
#ifdef LR1_taud
    fclose(xrbifs);
#endif
#ifdef TT
    fclose(caibifs);
    fclose(casrbifs);
    fclose(kibifs);
    fclose(naibifs);
#endif
#ifdef TTMod
    fclose(caibifs);
    fclose(casrbifs);
    fclose(kibifs);
    fclose(naibifs);
#endif
    return 0;
}


