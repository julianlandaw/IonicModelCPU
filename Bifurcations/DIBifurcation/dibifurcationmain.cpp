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

#include "DIBifurcation.cpp"

int main(int argc, char *argv[])
{
    const double minitofac = atof(argv[1]);
    const double maxitofac = atof(argv[2]);
    const double mindi = atof(argv[3]);
    const double maxdi = atof(argv[4]);
    const long double dt = atof(argv[5]);
    
    double byitofac = (numvar > 1) ? (maxitofac - minitofac)/(numvar - 1) : 1;
    double bydi = (numdi > 1) ? (maxdi - mindi)/(numdi - 1) : 1;
    
    DIBifurcation<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new DIBifurcation<TYPECELL, NCELLS, BEATS>();
    
    FILE *allbifs;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sbifsDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    allbifs = fopen(fileSpec1, "w");
    
    
#ifdef LR1
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    xrbifs = fopen(fileSpec2, "w");
#endif
    
#ifdef LR1_taud
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrDI_%g_taud_%g.txt", TYPECELLSTRING, mindi, minitofac);
    xrbifs = fopen(fileSpec2, "w");
#endif
    
#ifdef TT
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    naibifs = fopen(fileSpec5, "w");
#endif
    
#ifdef TTMod
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiDI_%g_ITO_%g.txt", TYPECELLSTRING, mindi, minitofac);
    naibifs = fopen(fileSpec5, "w");
#endif
    
    for (int i = 0; i < numdi; i++) {
        for (int j = 0; j < numvar; j++) {
            h_cells->Cells.VARIABLE[numvar*i + j] = byitofac*j + minitofac;
            h_cells->dis[numvar*i + j] = bydi*i + mindi;
#ifdef LR1
            h_cells->Cells.tauXfac[numvar*i + j] = 10.0;
#endif
            
#ifdef LR1_taud
            h_cells->Cells.tauXfac[numvar*i + j] = 5.0;
            h_cells->Cells.itofac[numvar*i + j] = 1.05;
#endif
#ifdef TTMod
            h_cells->Cells.ibarcafac[numvar*i + j] = 5.0;
#endif
//#ifdef TT
//	    h_cells->Cells.ibarcafac[numvar*i + j] = 6e-4/(gcal);
//	    h_cells->Cells.ikrfac[numvar*i + j] = 0.01/(gkr);
//	    h_cells->Cells.iksfac[numvar*i + j] = 0.036/(gks);
//#endif
        }
    }
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(DIBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    h_cells->dobif(dt, t);
    
    for (int i = 0; i < numdi; i++) {
        for (int j = 0; j < numvar; j++) {
            fprintf(allbifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            for (int k = REMOVEBEATS; k < BEATS; k++) {
                fprintf(allbifs, "\t%g", h_cells->apds[BEATS*(numvar*i + j) + k]);
            }
            fprintf(allbifs, "\n");
            
            
#ifdef LR1
            fprintf(xrbifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            for (int k = 2*REMOVEBEATS; k < 2*BEATS; k++) {
                fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(numvar*i + j) + k]);
            }
            fprintf(xrbifs, "\n");
#endif
            
#ifdef LR1_taud
            fprintf(xrbifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            for (int k = 2*REMOVEBEATS; k < 2*BEATS; k++) {
                fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(numvar*i + j) + k]);
            }
            fprintf(xrbifs, "\n");
#endif
#ifdef TT
            fprintf(caibifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            fprintf(casrbifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            fprintf(kibifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            fprintf(naibifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            
            for (int k = 2*REMOVEBEATS; k < 2*BEATS; k++) {
                fprintf(caibifs, "\t%g", h_cells->cais[2*BEATS*(numvar*i + j) + k]);
                fprintf(casrbifs, "\t%g", h_cells->casrs[2*BEATS*(numvar*i + j) + k]);
                fprintf(kibifs, "\t%g", h_cells->kis[2*BEATS*(numvar*i + j) + k]);
                fprintf(naibifs, "\t%g", h_cells->nais[2*BEATS*(numvar*i + j) + k]);
            }
            fprintf(caibifs, "\n");
            fprintf(casrbifs, "\n");
            fprintf(kibifs, "\n");
            fprintf(naibifs, "\n");
#endif
#ifdef TTMod
            fprintf(caibifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            fprintf(casrbifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            fprintf(kibifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            fprintf(naibifs, "%g\t%g", bydi*i + mindi, byitofac*j + minitofac);
            
            for (int k = 2*REMOVEBEATS; k < 2*BEATS; k++) {
                fprintf(caibifs, "\t%g", h_cells->cais[2*BEATS*(numvar*i + j) + k]);
                fprintf(casrbifs, "\t%g", h_cells->casrs[2*BEATS*(numvar*i + j) + k]);
                fprintf(kibifs, "\t%g", h_cells->kis[2*BEATS*(numvar*i + j) + k]);
                fprintf(naibifs, "\t%g", h_cells->nais[2*BEATS*(numvar*i + j) + k]);
            }
            fprintf(caibifs, "\n");
            fprintf(casrbifs, "\n");
            fprintf(kibifs, "\n");
            fprintf(naibifs, "\n");
#endif
        }
    }

    delete h_cells;
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


