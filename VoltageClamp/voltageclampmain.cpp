//
//  voltageclampmain.cpp
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

#ifdef UCLA
#define TYPECELL UCLACell
#define TYPECELLSTRING "UCLA"
#include "Cells/UCLACell.cpp"

//#define TT
//#ifdef TT
#elif TT
#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
#include "Cells/TTCellIto.cpp"
#else
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"
#endif

#ifndef numitofac
#define numitofac 1
#endif

#ifndef numpcl
#define numpcl 1
#endif

#ifndef numvolt1
#define numvolt1 1
#endif

#ifndef numvolt2
#define numvolt2 1
#endif

#ifndef numapd
#define numapd 1
#endif

#define NCELLS numitofac*numpcl*numvolt1*numvolt2*numapd

#ifndef BEATS
#define BEATS 20
#endif

#ifndef REMOVEBEATS
#define REMOVEBEATS 0
#endif

#ifndef stimt
#define stimt 100.0
#endif

#ifndef VARIABLE
#define VARIABLE cai
#endif

#include "VoltageClamp.cpp"

int main(int argc, char *argv[])
{
    /*
    const double minitofac = atof(argv[1]);
    const double byitofac = atof(argv[2]);
    const double minpcl = atof(argv[3]);
    const double bypcl = atof(argv[4]);
    const double minvolt1 = atof(argv[5]);
    const double byvolt1 = atof(argv[6]);
    const double minvolt2 = atof(argv[7]);
    const double byvolt2 = atof(argv[8]);
    const double minapd = atof(argv[9]);
    const double byapd = atof(argv[10]);
    const double dt = atof(argv[11]);
     */
    const double pcl = atof(argv[1]);
    const double apd1 = atof(argv[2]);
    const double apd2 = atof(argv[3]);
    const long double dt = atof(argv[4]);
    
    //VoltageClamp<TYPECELL, NCELLS, BEATS>* h_cells;
    
    //h_cells = new VoltageClamp<TYPECELL, NCELLS, BEATS>();
    
    VoltageClamp<TYPECELL, 1, BEATS>* h_cells;
    
    h_cells = new VoltageClamp<TYPECELL, 1, BEATS>();
    
    h_cells->commondt = dt;
    /*
    for (int i = 0; i < numpcl; i++) {
        for (int j = 0; j < numitofac; j++) {
            for (int k = 0; k < numvolt1; k++) {
                for (int l = 0; l < numvolt2; l++) {
                    for (int m = 0; m < numapd; m++) {
                        h_cells->pcls[i*numitofac*numvolt1*numvolt2*numapd + j*numvolt1*numvolt2*numapd + k*numvolt2*numapd + l*numapd + m] = bypcl*i + minpcl;
                        h_cells->Cells.itofac[i*numitofac*numvolt1*numvolt2*numapd + j*numvolt1*numvolt2*numapd + k*numvolt2*numapd + l*numapd + m] = byitofac*j + minitofac;
                        h_cells->voltage1[i*numitofac*numvolt1*numvolt2*numapd + j*numvolt1*numvolt2*numapd + k*numvolt2*numapd + l*numapd + m] = byvolt1*k + minvolt1;
                        h_cells->voltage2[i*numitofac*numvolt1*numvolt2*numapd + j*numvolt1*numvolt2*numapd + k*numvolt2*numapd + l*numapd + m] = byvolt2*l + minvolt2;
                        h_cells->apds[i*numitofac*numvolt1*numvolt2*numapd + j*numvolt1*numvolt2*numapd + k*numvolt2*numapd + l*numapd + m] = byapd*m + minapd;
                    }
                }
            }
        }
    }
     */
    h_cells->pcls[0] = pcl;
    h_cells->apds1[0] = apd1;
    h_cells->apds2[0] = apd2;
    
    FILE *allbifs;
    char fileSpec1[100];
    //snprintf(fileSpec1, 100, "%sbifsPCL_%g_ITO_%g_VOLT1_%g_VOLT2_%g_APD_%g.txt", TYPECELLSTRING, minpcl, minitofac, minvolt1, minvolt2, minapd);
    snprintf(fileSpec1, 100, "%sbifsPCL_%g_APD1_%g_APD2_%g.txt", TYPECELLSTRING, pcl, apd1, apd2);
    allbifs = fopen(fileSpec1, "w");

#ifdef TT
    FILE *allcais;
     char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaisPCL_%g_APD1_%g_APD2_%g.txt", TYPECELLSTRING, pcl, apd1, apd2);
    allcais = fopen(fileSpec2, "w");

    FILE *allnais;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%snaisPCL_%g_APD1_%g_APD2_%g.txt", TYPECELLSTRING, pcl, apd1, apd2);
    allnais = fopen(fileSpec3, "w");
#endif
#ifdef UCLA
    FILE *allcais;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaisPCL_%g_APD1_%g_APD2_%g.txt", TYPECELLSTRING, pcl, apd1, apd2);
    allcais = fopen(fileSpec2, "w");

    FILE *allnais;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%snaisPCL_%g_APD1_%g_APD2_%g.txt", TYPECELLSTRING, pcl, apd1, apd2);
    allnais = fopen(fileSpec3, "w");
#endif

    FILE *ap;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%sapPCL_%g_APD1_%g_APD2_%g.txt", TYPECELLSTRING, pcl, apd1, apd2);
    ap = fopen(fileSpec4, "w");
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(VoltageClamp<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    h_cells->dospiketrain(ap);
    
    /*
    for (int i = 0; i < numpcl; i++) {
        for (int j = 0; j < numitofac; j++) {
            for (int k = 0; k < numvolt1; k++) {
                for (int l = 0; l < numvolt2; l++) {
                    for (int m = 0; m < numapd; m++) {
                        fprintf(allbifs, "%g\t%g\t%g\t%g\t%g", bypcl*i + minpcl, byitofac*j + minitofac, byvolt1*k + minvolt1, byvolt2*l + minvolt2, byapd*m + minapd);
                        for (int n = REMOVEBEATS; n < BEATS; n++) {
                            fprintf(allbifs, "\t%g\t%g", h_cells->VARIABLE[2*BEATS*(i*numitofac*numvolt1*numvolt2*numapd + j*numvolt1*numvolt2*numapd + k*numvolt2*numapd + l*numapd + m)+ 2*n], h_cells->VARIABLE[2*BEATS*(i*numitofac*numvolt1*numvolt2*numapd + j*numvolt1*numvolt2*numapd + k*numvolt2*numapd + l*numapd + m)+ 2*n + 1]);
                        }
                        fprintf(allbifs, "\n");
                    }
                }
            }
        }
    }
     */
#ifdef TT
    fprintf(allcais,"%g\t%g\t%g",pcl,apd1,apd2);
    fprintf(allnais,"%g\t%g\t%g",pcl,apd1,apd2);
    for (int n = REMOVEBEATS; n < 2*BEATS; n++) {
        fprintf(allcais, "\t%g\t%g", h_cells->cai[2*n],h_cells->cai[2*n+1]);
	fprintf(allnais, "\t%g\t%g", h_cells->nai[2*n],h_cells->nai[2*n+1]);
    }
    fprintf(allcais,"\n");
    fprintf(allnais,"\n");
#endif

#ifdef UCLA
    fprintf(allcais,"%g\t%g\t%g",pcl,apd1,apd2);
    fprintf(allnais,"%g\t%g\t%g",pcl,apd1,apd2);
    for (int n = REMOVEBEATS; n < 2*BEATS; n++) {
        fprintf(allcais, "\t%g\t%g", h_cells->ci[2*n],h_cells->ci[2*n+1]);
	fprintf(allnais, "\t%g\t%g", h_cells->nai[2*n],h_cells->nai[2*n+1]);
    }
    fprintf(allcais,"\n");
    fprintf(allnais,"\n");
#endif
    

    delete h_cells;
    fclose(allbifs);
#ifdef TT
    fclose(allcais);
    fclose(allnais);
#endif
#ifdef UCLA
    fclose(allcais);
    fclose(allnais);
#endif
    fclose(ap);
    return 0;
}


