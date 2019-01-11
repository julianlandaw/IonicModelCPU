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

#ifndef numpcl
#define numpcl 1
#endif

#define NCELLS (numvar1*numvar2*numpcl)

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

#include "NaBifurcation.cpp"

int main(int argc, char *argv[])
{
    const double startna = atof(argv[1]);
    const double endna = atof(argv[2]);
    const double minvar1 = atof(argv[3]);
    const double maxvar1 = atof(argv[4]);
    const double minvar2 = atof(argv[5]);
    const double maxvar2 = atof(argv[6]);
    const double minpcl = atof(argv[7]);
    const double maxpcl = atof(argv[8]);
    const long double dt = atof(argv[9]);
    
    	#ifdef LR1
    double _itofac = 0;
    double _icalfac = 1.0; //1.15; //1.0;
    double _ikfac = 1.0;
    double _ikifac = 1.0; //2.2; //1.0;
    double _tauXfac = 1.0; //5.0;
    double _yshift = 0.0; //8.0; //-8.0;
        #endif 
        #ifdef LR2
    double _itofac = 0.0;
    double _iskfac = 0.0;
    double _iupfac = 1.0;
    double _skh = 0.0025;
        #endif
        #ifdef TT
    double _itofac = 0.0;
    double _iskfac = 0.0;
    double _nai = 12.0;
    double _ki = 138.0;
    double _nacafac = 5.75;

    double _ikrfac = 1.0;
    double _iksfac = 1.0;
    double _icalfac = 1.0;
    //double _ikrfac = 0.01/(gkr);  //1.0;
    //double _iksfac = 0.036/(gks); //1.0;
    //double _icalfac = 0.0006/(gcal); //1.0;
        #endif
        #ifdef OHara
    double _itofac = 0.0;
    double _iskfac = 0.0;
    double _inafac = 2.0; 
        #endif
        #ifdef UCLA
    double _itofac = 1.0;
    double _iskfac = 0.0;
    double _skh = 2.0;
    double _icalfac = 1.0;
    double _ikrfac = 1.0;
    double _iksfac = 1.0;
    double _nai = 8.0;
    double _nacafac = 0.3;
        #endif
    
    double byvar1 = (numvar1 > 1) ? (maxvar1 - minvar1)/(numvar1 - 1) : 1;
    double byvar2 = (numvar2 > 1) ? (maxvar2 - minvar2)/(numvar1 - 1) : 1;
    double bypcl = (numpcl > 1) ? (maxpcl - minpcl)/(numpcl - 1) : 1;
    
    NaBifurcation<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new NaBifurcation<TYPECELL, NCELLS, BEATS>();
    
    
    FILE *allbifs;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sbifsPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    allbifs = fopen(fileSpec1, "w");
    
    
#ifdef LR1
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    xrbifs = fopen(fileSpec2, "w");

    FILE *ytosbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%sYtosPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1, minvar2,startna);
    ytosbifs = fopen(fileSpec3, "w");
#endif
    
#ifdef LR1_taud
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    xrbifs = fopen(fileSpec2, "w");
#endif
    
#ifdef TT
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    naibifs = fopen(fileSpec5, "w");
#endif
    
#ifdef UCLA
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *cassbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scassPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    cassbifs = fopen(fileSpec3, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    naibifs = fopen(fileSpec5, "w");
#endif    
    
#ifdef LR2
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    naibifs = fopen(fileSpec5, "w");
#endif
    
#ifdef TTMod
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrPCL_%g_VAR1_%g_VAR2_%g_Na_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g.txt", TYPECELLSTRING, minpcl, minvar1,minvar2,startna);
    naibifs = fopen(fileSpec5, "w");
#endif
    int index;
    for (int i = 0; i < numpcl; i++) {
        for (int j = 0; j < numvar1; j++) {
            for (int k = 0; k < numvar2; k++) {
                index = (numvar1*numvar2*i) + (numvar2*j) + k; //numvar*i + j
		    #ifdef LR1
                h_cells->Cells.itofac[index] = _itofac;
                h_cells->Cells.icalfac[index] = _icalfac;
                h_cells->Cells.ikfac[index] = _ikfac;
                h_cells->Cells.ikifac[index] = _ikifac;
                h_cells->Cells.tauXfac[index] = _tauXfac;
                h_cells->Cells.yshift[index] = _yshift;
		    #endif
		    #ifdef LR2
                h_cells->Cells.itofac[index] = _itofac;
                h_cells->Cells.iskfac[index] = _iskfac;
                h_cells->Cells.iupfac[index] = _iupfac;
                h_cells->Cells.skh[index] = _skh;
            #endif
            #ifdef TT
                h_cells->Cells.itofac[index] = _itofac;
                h_cells->Cells.iskfac[index] = _iskfac;
                h_cells->Cells.nai[index] = _nai;
                h_cells->Cells.ki[index] = _ki;
                h_cells->Cells.nacafac[index] = _nacafac;

                h_cells->Cells.ikrfac[index] = _ikrfac;
                h_cells->Cells.iksfac[index] = _iksfac;
                h_cells->Cells.ibarcafac[index] = _icalfac;
                    #endif

#ifdef OHara
                h_cells->Cells.itofac[index] = _itofac;//0.0;
                h_cells->Cells.iskfac[index] = _iskfac;//0.0;
                h_cells->Cells.inafac[index] = _inafac; //2.0;
#endif
#ifdef UCLA
                h_cells->Cells.itofac[index] = _itofac; //0.0;
                h_cells->Cells.iskfac[index] = _iskfac; //0.0;
                h_cells->Cells.skh[index] = _skh; //2.0;
                h_cells->Cells.icalfac[index] = _icalfac;
                h_cells->Cells.ikrfac[index] = _ikrfac;
                h_cells->Cells.iksfac[index] = _iksfac;
                h_cells->Cells.nai[index] = _nai;
		h_cells->Cells.nacafac[index] = _nacafac;
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
                h_cells->pcls[index] = bypcl*i + minpcl;
                h_cells->Cells.nai[index] = startna;
                h_cells->startna[index] = startna;
                h_cells->endna[index] = endna;
                h_cells->Cells.naiclamped[index] = true;
            }
        }
    }
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(NaBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    h_cells->dobif(dt, t);
    
    for (int i = 0; i < numpcl; i++) {
        for (int j = 0; j < numvar1; j++) {
            for (int k = 0; k < numvar2; k++) {
                index = (numvar1*numvar2*i) + (numvar2*j) + k;
                fprintf(allbifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                for (int l = REMOVEBEATS; l < BEATS; l++) {
                    fprintf(allbifs, "\t%g", h_cells->apds[BEATS*(index) + l]);
                }
                fprintf(allbifs, "\n");
            
            
#ifdef LR1
                fprintf(xrbifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(ytosbifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(index) + l]);
                fprintf(ytosbifs, "\t%g", h_cells->ytoss[2*BEATS*(index) + l]);
                }
                fprintf(xrbifs, "\n");
                fprintf(ytosbifs, "\n");
#endif
            
#ifdef LR1_taud
                fprintf(xrbifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(index) + l]);
                }
                fprintf(xrbifs, "\n");
#endif
#ifdef TT
                fprintf(caibifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(casrbifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(kibifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(naibifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
            
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
                fprintf(caibifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(naibifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(cassbifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
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
                fprintf(caibifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(kibifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
                fprintf(naibifs, "%g\t%g\t%g", bypcl*i + minpcl, byvar1*j + minvar1, byvar2*k + minvar2);
            
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
    }

    delete h_cells;
    fclose(allbifs);
#ifdef LR1
    fclose(xrbifs);
    fclose(ytosbifs);
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


