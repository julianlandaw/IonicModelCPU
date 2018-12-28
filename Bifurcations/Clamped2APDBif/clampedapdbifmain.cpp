//
//  apdbifurcationmain.cu
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
//#define gnaped 0.01
#include "Cells/TTCellIto.cpp"
#elif TTSK
#define TYPECELL TTCellItoSK
#define TYPECELLSTRING "TTSK"
#include "Cells/TTCellItoSK.cpp"
#else
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"
#endif

#ifndef numpcl
#define numpcl 1
#endif

#ifndef numitofac
#define numitofac 16
#endif

#ifndef numvar1
#define numvar1 16
#endif

#ifndef numvar2
#define numvar2 16
#endif

#define NCELLS numpcl*numitofac*numvar1*numvar2

#ifndef BEATS
#define BEATS 600
#endif

#ifndef REMOVEBEATS
#define REMOVEBEATS 400
#endif

#ifndef stimt
#define stimt 100.0L
#endif

#ifndef VARIABLE1
#define VARIABLE1 cai
#endif

#ifndef VARIABLE2
#define VARIABLE2 casr
#endif

#ifndef VARIABLE3
#define VARIABLE3 nai
#endif

#ifndef VARIABLE4
#define VARIABLE4 ki
#endif

#ifndef curfac
#define curfac itofac
#endif

#include "ClampedAPDBif.cpp"

int main(int argc, char *argv[])
{
    const double minpcl = atof(argv[1]);
    const double maxpcl = atof(argv[2]);
    const double minitofac = atof(argv[3]);
    const double maxitofac = atof(argv[4]);
    const double minclamp1 = atof(argv[5]);
    const double maxclamp1 = atof(argv[6]);
    const double minclamp2 = atof(argv[7]);
    const double maxclamp2 = atof(argv[8]);
    const long double dt = atof(argv[9]);

    double bypcl, byitofac, byclamp1, byclamp2;
    bypcl = (numpcl > 1) ? (maxpcl - minpcl)/(numpcl - 1) : 1;
    byitofac = (numitofac > 1) ? (maxitofac - minitofac)/(numitofac - 1) : 1;
    byclamp1 = (numvar1 > 1) ? (maxclamp1 - minclamp1)/(numvar1 - 1) : 1;
    byclamp2 = (numvar2 > 1) ? (maxclamp2 - minclamp2)/(numvar2 - 1) : 1;

    ClampedAPDBif<TYPECELL, NCELLS, BEATS>* h_cells;
    h_cells = new ClampedAPDBif<TYPECELL, NCELLS, BEATS>();
    
    FILE *allbifs;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sbifs_%g_%g_%g_%g.txt", TYPECELLSTRING, minpcl, round(minitofac*10000), minclamp1, minclamp2);
    allbifs = fopen(fileSpec1, "w");

    FILE *var1bif;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%svar1_%g_%g_%g_%g.txt", TYPECELLSTRING, minpcl, round(minitofac*10000),  minclamp1, minclamp2);
    var1bif = fopen(fileSpec2, "w");

    FILE *var2bif;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%svar2_%g_%g_%g_%g.txt", TYPECELLSTRING, minpcl, round(minitofac*10000),  minclamp1, minclamp2);
    var2bif = fopen(fileSpec3, "w");
    
    for (int l = 0; l < numpcl; l++) {
        for (int k = 0; k < numvar1; k++) {
            for (int i = 0; i < numvar2; i++) {
                for (int j = 0; j < numitofac; j++) {
                    int id = numitofac*numvar2*numvar1*l + numitofac*numvar2*k + numitofac*i + j;
#ifdef TTSK
		    h_cells->Cells.iskfac[id] = 5.0;
#endif
                    h_cells->Cells.curfac[id] = byitofac*j + minitofac;
                    h_cells->clampvar3[id] = byclamp1*k + minclamp1;
                    h_cells->clampvar4[id] = byclamp2*i + minclamp2;
                    h_cells->pcls[id] = bypcl*l + minpcl;
#ifdef TT
  //                  h_cells->Cells.ibarcafac[0] = 6e-4/(gcal);
  //                  h_cells->Cells.ikrfac[0] = 0.01/(gkr);
  //                  h_cells->Cells.iksfac[0] = 0.036/(gks);
                    //h_cells->Cells.ibarcafac[id] = 1.05;
		    //h_cells->Cells.nacafac[id] = 1.25;
		    //h_cells->Cells.ikrfac[id] = 0.96875;
#endif
                }
            }
        }
    }
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(ClampedAPDBif<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    // Perform bifurcation on cells in the device
    
    long double t = -stimt;
    int i = 0;
    double t_save = 0;
    while (i < NCELLS) {
        h_cells->iterateall(dt, t);
        t = t + dt;
        if (t > t_save - dt/4.0) {
            while (i < NCELLS && !(h_cells->notdonebifurcation[i])) {
                i = i+1;
            }
            t_save = t_save + 1000;
            printf("%g\t%d\n", (double)t, i);
        }
    }
    
    for (int l = 0; l < numpcl; l++) {
        for (int k = 0; k < numvar1; k++) {
            for (int i = 0; i < numvar2; i++) {
                for (int j = 0; j < numitofac; j++) {
                    fprintf(allbifs, "%g\t%g\t%g\t%g", bypcl*l + minpcl, byclamp1*k + minclamp1, byclamp2*i + minclamp2, byitofac*j + minitofac);
                    fprintf(var1bif, "%g\t%g\t%g\t%g", bypcl*l + minpcl, byclamp1*k + minclamp1, byclamp2*i + minclamp2, byitofac*j + minitofac);
                    fprintf(var2bif, "%g\t%g\t%g\t%g", bypcl*l + minpcl, byclamp1*k + minclamp1, byclamp2*i + minclamp2, byitofac*j + minitofac);
                    int id = numitofac*numvar2*numvar1*l + numitofac*numvar2*k + numitofac*i + j;
                    for (int m = REMOVEBEATS; m < BEATS; m++) {
                        fprintf(allbifs, "\t%g", h_cells->apds[BEATS*id + m]);
                    }
                    for (int m = 2*REMOVEBEATS; m < 2*BEATS; m++) {
                        fprintf(var1bif, "\t%g", h_cells->var1[2*BEATS*id + m]);
                        fprintf(var2bif, "\t%g", h_cells->var2[2*BEATS*id + m]);
                    }
                    fprintf(allbifs, "\n");
                    fprintf(var1bif, "\n");
                    fprintf(var2bif, "\n");
                }
            }
        }
    }

    delete h_cells;
    fclose(allbifs);
    fclose(var1bif);
    fclose(var2bif);
    return 0;
}


