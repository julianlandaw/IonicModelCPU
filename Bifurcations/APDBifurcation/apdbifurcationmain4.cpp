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

#elif TP06
#define TYPECELL TP06Cell
#define TYPECELLSTRING "TP06"
#include "Cells/TP06Cell.cpp"

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

#elif Hund
#define TYPECELL HundCell
#define TYPECELLSTRING "Hund"
#include "Cells/HundCell.cpp"

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
#define VARIABLE1 icalfac
#endif

#ifndef VARIABLE2
#define VARIABLE2 ikrfac
#endif

#ifndef VARIABLE3
#define VARIABLE3 iksfac
#endif

#ifndef VARIABLE4
#define VARIABLE4 ik1fac
#endif

#include "APDBifurcation.cpp"

int main(int argc, char *argv[])
{
    const double var1 = atof(argv[1]);
    const double var2 = atof(argv[2]);
    const double var3 = atof(argv[3]);
    const double var4 = atof(argv[4]);
    const double pcl = atof(argv[5]);
    const long double dt = atof(argv[6]);
    
    	#ifdef LR1
    double _itofac = 0.0; //1.0;
    double _icalfac = 1.0; //1.15; //1.0; //1.15; //1.0;
    double _ikfac = 1.0;
    double _ikifac = 1.0; //2.2; //1.0; //2.2; //1.0;
    double _tauXfac = 1.0; //5.0; //5.0;
    double _yshift = -8.0; //8.0; //-8.0; //8.0; //-8.0;
        #endif 
        #ifdef LR2
    double _itofac = 0.0;
    double _iskfac = 0.0;
    double _iupfac = 1.0; //3.0; //3.0;
    double _skh = 0.002; //0.0025
    double _tauyfac = 1.0;
    //double _nai = 20.0;
    //double _ki = 134.2;
        #endif
        #ifdef TT
    double _itofac = 0.0;
    double _iskfac = 0.0;
    //double _nai = 12.0;
    //double _ki = 138.0;
    double _nacafac = 1.0; //5.75;

    double _ikrfac = 0.0; //1.0;
    double _iksfac = 1;
    double _icalfac = 1; //1.0;
    //double _ikrfac = 0.01/(gkr);  //1.0;
    //double _iksfac = 0.036/(gks); //1.0;
    //double _icalfac = 0.0006/(gcal); //1.0;
        #endif
    #ifdef Hund
    double _itofac = 1.0;
    double _inafac = 1.0;
    double _ikrfac = 1.0;
    double _iksfac = 1.0;
    double _icalfac = 1.0;
    #endif
        #ifdef OHara
    double _itofac = 0.0; //0.0; //2.0;
    double _iskfac = 0.0;
    double _skh = 0.001; //0.0005; //0.006;
    double _inafac = 2.0; //2.0; //2.0;
    double _nai = 7.92; //8.73; 
    double _ki = 146.0;
    double _vssfac = 1.0;
    double _alphask = 0.1; //0.1; //0.1;
    double _ikrfac = 1.0;
    double _iksfac = 1.0;
        #endif
        #ifdef UCLA
    double _itofac = 0.0;
    double _iskfac = 0.0;
    double _skh = 2.0;
    double _icalfac = 1.0;
    double _ikrfac = 1.0;
    double _iksfac = 1.0;
    double _nai = 12.0;
    double _nacafac = 1.0;//1.0;
        #endif
    
    APDBifurcation<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new APDBifurcation<TYPECELL, NCELLS, BEATS>();
    
    
    FILE *allbifs;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sbifsPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    allbifs = fopen(fileSpec1, "w");
    
    
#ifdef LR1
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    xrbifs = fopen(fileSpec2, "w");

    FILE *ytosbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%sYtosPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    ytosbifs = fopen(fileSpec3, "w");
#endif
    
#ifdef LR1_taud
    FILE *xrbifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sXrPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    xrbifs = fopen(fileSpec2, "w");
#endif
    
#ifdef TT
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    naibifs = fopen(fileSpec5, "w");
    
    FILE *cpeakbifs;
    char fileSpec6[100];
    snprintf(fileSpec6, 100, "%scpeakPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    cpeakbifs = fopen(fileSpec6, "w");
#endif
    
#ifdef UCLA
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *cassbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scassPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    cassbifs = fopen(fileSpec3, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    naibifs = fopen(fileSpec5, "w");
#endif  

#ifdef OHara
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *cassbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scassPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    cassbifs = fopen(fileSpec3, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    naibifs = fopen(fileSpec5, "w");
#endif  
    
#ifdef LR2
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    naibifs = fopen(fileSpec5, "w");
#endif
    
#ifdef TTMod
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    naibifs = fopen(fileSpec5, "w");
#endif
    int index = 0;
    //for (int i = 0; i < numpcl; i++) {
 //       for (int j = 0; j < numvar1; j++) {
   //         for (int k = 0; k < numvar2; k++) {
    //            index = (numvar1*numvar2*i) + (numvar2*j) + k; //numvar*i + j
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
		h_cells->Cells.tauyfac[index] = _tauyfac;
		//h_cells->Cells.nai[index] = _nai;
		//h_cells->Cells.ki[index] = _ki;
            #endif
		#ifdef Hund
		h_cells->Cells.itofac[index] = _itofac;
		h_cells->Cells.inafac[index] = _inafac;
		h_cells->Cells.icalfac[index] = _icalfac;
		h_cells->Cells.ikrfac[index] = _ikrfac;
		h_cells->Cells.iksfac[index] = _iksfac;
		#endif
            #ifdef TT
                h_cells->Cells.itofac[index] = _itofac;
                h_cells->Cells.iskfac[index] = _iskfac;
                //h_cells->Cells.nai[index] = _nai;
                //h_cells->Cells.ki[index] = _ki;
		h_cells->Cells.nacafac[index] = _nacafac;

                h_cells->Cells.ikrfac[index] = _ikrfac;
                h_cells->Cells.iksfac[index] = _iksfac;
                h_cells->Cells.ibarcafac[index] = _icalfac;
                    #endif

#ifdef OHara
                h_cells->Cells.itofac[index] = _itofac;//0.0;
                h_cells->Cells.iskfac[index] = _iskfac;//0.0;
                h_cells->Cells.inafac[index] = _inafac; //2.0;
		h_cells->Cells.skh[index] = _skh;
		h_cells->Cells.nai[index] = _nai;
		h_cells->Cells.nass[index] = _nai;
		h_cells->Cells.ki[index] = _ki;
		h_cells->Cells.kss[index] = _ki;
		h_cells->Cells.ikrfac[index] = _ikrfac;
		h_cells->Cells.iksfac[index] = _iksfac;
		h_cells->Cells.vssfac[index] = _vssfac;
		h_cells->Cells.alphask[index] = _alphask;
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
                h_cells->Cells.VARIABLE1[index] = var1;
                h_cells->Cells.VARIABLE2[index] = var2;
                h_cells->Cells.VARIABLE3[index] = var3;
                h_cells->Cells.VARIABLE4[index] = var4;
#ifdef TT
		printf("%d\t%g\n",index,h_cells->Cells.vcfac[index]);
#endif
		h_cells->pcls[index] = pcl;
#ifdef OHara
		h_cells->Cells.nass[index] = h_cells->Cells.nai[index];
		h_cells->Cells.kss[index] = h_cells->Cells.ki[index];
#endif
        //    }
       // }
   // }
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(APDBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    h_cells->dobif(dt, t);
    
  //  for (int i = 0; i < numpcl; i++) {
  //      for (int j = 0; j < numvar1; j++) {
  //          for (int k = 0; k < numvar2; k++) {
                index = 0;
                fprintf(allbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                for (int l = REMOVEBEATS; l < BEATS; l++) {
                    fprintf(allbifs, "\t%g", h_cells->apds[BEATS*(index) + l]);
                }
                fprintf(allbifs, "\n");
            
            
#ifdef LR1
                fprintf(xrbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(ytosbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(index) + l]);
                fprintf(ytosbifs, "\t%g", h_cells->ytoss[2*BEATS*(index) + l]);
                }
                fprintf(xrbifs, "\n");
                fprintf(ytosbifs, "\n");
#endif
            
#ifdef LR1_taud
                fprintf(xrbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(xrbifs, "\t%g", h_cells->xrs[2*BEATS*(index) + l]);
                }
                fprintf(xrbifs, "\n");
#endif
#ifdef TT
                fprintf(caibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(casrbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(kibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(naibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(cpeakbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
            
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(caibifs, "\t%g", h_cells->cais[2*BEATS*(index) + l]);
                    fprintf(casrbifs, "\t%g", h_cells->casrs[2*BEATS*(index) + l]);
                    fprintf(kibifs, "\t%g", h_cells->kis[2*BEATS*(index) + l]);
                    fprintf(naibifs, "\t%g", h_cells->nais[2*BEATS*(index) + l]);
                }
                for (int l = REMOVEBEATS; l < BEATS; l++) {
                    fprintf(cpeakbifs, "\t%g", h_cells->cpeaks[BEATS*(index) + l]);
                }
                fprintf(caibifs, "\n");
                fprintf(casrbifs, "\n");
                fprintf(kibifs, "\n");
                fprintf(naibifs, "\n");
                fprintf(cpeakbifs, "\n");
#endif
                
#ifdef UCLA
                fprintf(caibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(naibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(cassbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
                    fprintf(caibifs, "\t%g", h_cells->cais[2*BEATS*(index) + l]);
                    fprintf(naibifs, "\t%g", h_cells->nais[2*BEATS*(index) + l]);
                    fprintf(cassbifs, "\t%g", h_cells->cass[2*BEATS*(index) + l]);
                }
                fprintf(caibifs, "\n");
                fprintf(naibifs, "\n");
                fprintf(cassbifs, "\n");
#endif

#ifdef OHara
                fprintf(caibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(naibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(cassbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
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
                fprintf(caibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(kibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(naibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
            
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
                fprintf(caibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(casrbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(kibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
                fprintf(naibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
            
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
 //           }
 //       }
 //   }

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
    fclose(cpeakbifs);
#endif
#ifdef UCLA
    fclose(caibifs);
    fclose(naibifs);
    fclose(cassbifs);
#endif
#ifdef OHara
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


