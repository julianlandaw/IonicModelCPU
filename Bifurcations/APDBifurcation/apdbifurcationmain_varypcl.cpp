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

#ifdef OHara
#define TYPECELL OHaraCell
#define TYPECELLSTRING "OHara"
#ifndef VARIABLE1
#define VARIABLE1 iskfac
#endif
#ifndef VARIABLE2
#define VARIABLE2 skh
#endif
#include "Cells/OHaraCell.cpp"

#elif UCLA
#define TYPECELL UCLACell
#define TYPECELLSTRING "UCLA"
#ifndef VARIABLE1
#define VARIABLE1 iskfac
#endif
#ifndef VARIABLE2
#define VARIABLE2 skh
#endif
#include "Cells/UCLACell.cpp"

#elif TT
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

#elif LR2
#ifdef SK
#define TYPECELL LR2CellIto
#define TYPECELLSTRING "LR2SK"
#ifndef VARIABLE1
#define VARIABLE1 iskfac
#endif
#ifndef VARIABLE2
#define VARIABLE2 skh
#endif
#include "Cells/LR2CellIto.cpp"
#else

#define TYPECELL LR2CellIto
#define TYPECELLSTRING "LR2"
#include "Cells/LR2CellIto.cpp"
#endif

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

#ifndef numvar1
#define numvar1 1
#endif

#ifndef numvar2
#define numvar2 1
#endif

//#ifndef numpcl
//#define numpcl 1
//#endif

#define NCELLS (numvar1*numvar2)

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
#define VARIABLE2 tauXfac
#endif

#include "APDBifurcation_varypcl.cpp"

int main(int argc, char *argv[])
{
    const double minvar1 = atof(argv[1]);
    const double maxvar1 = atof(argv[2]);
    const double minvar2 = atof(argv[3]);
    const double maxvar2 = atof(argv[4]);
    const double pcl1 = atof(argv[5]);
    const double pcl2 = atof(argv[6]);
    const long double dt = atof(argv[7]);
    
    double byvar1 = (numvar1 > 1) ? (maxvar1 - minvar1)/(numvar1 - 1) : 1;
    double byvar2 = (numvar2 > 1) ? (maxvar2 - minvar2)/(numvar1 - 1) : 1;
    //double bypcl = (numpcl > 1) ? (maxpcl - minpcl)/(numpcl - 1) : 1;
    
    APDBifurcation_varypcl<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new APDBifurcation_varypcl<TYPECELL, NCELLS, BEATS>();
    
    
    FILE *allbifs;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sbifsPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    allbifs = fopen(fileSpec1, "w");
    
    
#ifdef LR1
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    xrbifs = fopen(fileSpec2, "w");
#endif
    
#ifdef LR1_taud
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    xrbifs = fopen(fileSpec2, "w");
#endif
    
#ifdef TT
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    naibifs = fopen(fileSpec5, "w");
#endif
    
#ifdef UCLA
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *cassbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scassPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    cassbifs = fopen(fileSpec3, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    naibifs = fopen(fileSpec5, "w");
#endif    
    
#ifdef LR2
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    naibifs = fopen(fileSpec5, "w");
#endif
    
#ifdef TTMod
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, pcl1, pcl2, minvar1,minvar2);
    naibifs = fopen(fileSpec5, "w");
#endif
    int index;
    //for (int i = 0; i < numpcl; i++) {
        for (int j = 0; j < numvar1; j++) {
            for (int k = 0; k < numvar2; k++) {
                index = (numvar2*j) + k; //(numvar1*numvar2*i) + (numvar2*j) + k; //numvar*i + j
#ifdef LR2
                h_cells->Cells.itofac[index] = 0;
                h_cells->Cells.iskfac[index] = 0;
                h_cells->Cells.iupfac[index] = 3.0;
#endif
#ifdef TT
                h_cells->Cells.itofac[index] = 0.0;
                h_cells->Cells.iskfac[index] = 0.0;
#ifdef clampnai
                h_cells->Cells.nai[index] = 12.0;
#endif
#ifdef clampki
                h_cells->Cells.ki[index] = 138.0;
#endif
#endif

#ifdef OHara
                h_cells->Cells.itofac[index] = 0.0;
                h_cells->Cells.iskfac[index] = 0.0;
                h_cells->Cells.inafac[index] = 2.0;
#endif
#ifdef UCLA
                h_cells->Cells.itofac[index] = 0.0;
                h_cells->Cells.iskfac[index] = 0.0;
                h_cells->Cells.skh[index] = 2.0;
#endif
#ifdef LR1
                h_cells->Cells.tauXfac[index] = 1.0;
#endif
            
#ifdef LR1_taud
                h_cells->Cells.tauXfac[index] = 5.0;
                h_cells->Cells.itofac[index] = 1.05;
#endif
#ifdef TTMod
                h_cells->Cells.ibarcafac[index] = 1.0;
#endif
                h_cells->Cells.VARIABLE1[index] = byvar1*j + minvar1;
                h_cells->Cells.VARIABLE2[index] = byvar2*k + minvar2;
                h_cells->pcls1[index] = pcl1;
                h_cells->pcls2[index] = pcl2;
            }
        }
    //}
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(APDBifurcation_varypcl<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    h_cells->dobif(dt, t);
    
    //for (int i = 0; i < numpcl; i++) {
        for (int j = 0; j < numvar1; j++) {
            for (int k = 0; k < numvar2; k++) {
                index = (numvar2*j) + k; //(numvar1*numvar2*i) + (numvar2*j) + k;
                fprintf(allbifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                for (int l = REMOVEBEATS; l < BEATS; l++) {
                    fprintf(allbifs, "\t%g", h_cells->apds[BEATS*(index) + l]);
                }
                fprintf(allbifs, "\n");
            
            
#ifdef LR1
                fprintf(xrbifs, "%g\t%g\t%g", pcl1, pcl2, byitofac*j + minitofac);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(index) + l]);
                }
                fprintf(xrbifs, "\n");
#endif
            
#ifdef LR1_taud
                fprintf(xrbifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(index) + l]);
                }
                fprintf(xrbifs, "\n");
#endif
#ifdef TT
                fprintf(caibifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(casrbifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(kibifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(naibifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
            
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
                
#ifdef UCLA
                fprintf(caibifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(naibifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(cassbifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(caibifs, "\t%g", h_cells->cais[2*BEATS*(index) + l]);
                    fprintf(naibifs, "\t%g", h_cells->nais[2*BEATS*(index) + l]);
                    fprintf(cassbifs, "\t%g", h_cells->cass[2*BEATS*(index) + l]);
                }
                fprintf(caibifs, "\n");
                fprintf(naibifs, "\n");
                fprintf(cassbifs, "\n");
#endif                
                
#ifdef LR2
                fprintf(caibifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(kibifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(naibifs, "%g\t%g\t%g\t%g", pcl1, pcl2, byvar1*j + minvar1, byvar2*k + minvar2);
            
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(caibifs, "\t%g", h_cells->cais[2*BEATS*(index) + l]);
                    fprintf(kibifs, "\t%g", h_cells->kis[2*BEATS*(index) + l]);
                    fprintf(naibifs, "\t%g", h_cells->nais[2*BEATS*(index) + l]);
                }
                fprintf(caibifs, "\n");
                fprintf(kibifs, "\n");
                fprintf(naibifs, "\n");
#endif            
#ifdef TTMod
                fprintf(caibifs, "%g\t%g\t%g", byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(casrbifs, "%g\t%g\t%g", byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(kibifs, "%g\t%g\t%g", byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(naibifs, "%g\t%g\t%g", byvar1*j + minvar1, byvar2*k + minvar2);
            
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
            }
        }
    //}

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
#ifdef UCLA
    fclose(caibifs);
    fclose(naibifs);
    fclose(cassbifs);
#endif
#ifdef LR2
    fclose(caibifs);
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


