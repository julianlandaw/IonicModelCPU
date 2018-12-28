//
//  ttcaibifurcationmain.cpp
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

#define TYPECELL TTCellIto
//#define gnaped 0.01
#define TYPECELLSTRING "TT"
#include "Cells/TTCellIto.cpp"

#ifndef numvar
#define numvar 1
#endif

#ifndef numcai
#define numcai 1
#endif

#define NCELLS numvar*numcai

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

#include "TTCaiBifurcation.cpp"

int main(int argc, char *argv[])
{
    const double minitofac = atof(argv[1]);
    const double maxitofac = atof(argv[2]);
    const double mincai = atof(argv[3]);
    const double maxcai = atof(argv[4]);
    const long double dt = atof(argv[5]);
    
    double byitofac = (numvar > 1) ? (maxitofac - minitofac)/(numvar - 1) : 1;
    double bycai = (numcai > 1) ? (maxcai - mincai)/(numcai - 1) : 1;
    
    TTCaiBifurcation<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new TTCaiBifurcation<TYPECELL, NCELLS, BEATS>();
    
    FILE *allbifs;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sbifs_caithresh_%g_ITO_%g.txt", TYPECELLSTRING, mincai, minitofac);
    allbifs = fopen(fileSpec1, "w");
    
    FILE *allpcls;
    char fileSpec6[100];
    snprintf(fileSpec6, 100, "%spcls_caithresh_%g_ITO_%g.txt", TYPECELLSTRING, mincai, minitofac);
    allpcls = fopen(fileSpec6, "w");
    
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scai_caithresh_%g_ITO_%g.txt", TYPECELLSTRING, mincai, minitofac);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasr_caithresh_%g_ITO_%g.txt", TYPECELLSTRING, mincai, minitofac);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%ski_caithresh_%g_ITO_%g.txt", TYPECELLSTRING, mincai, minitofac);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snai_caithresh_%g_ITO_%g.txt.txt", TYPECELLSTRING, mincai, minitofac);
    naibifs = fopen(fileSpec5, "w");
    
    int index;
    for (int i = 0; i < numcai; i++) {
        for (int j = 0; j < numvar; j++) {
            index = numvar*i + j;
#ifdef clampnai
            h_cells->Cells.nai[index] = 12.0;
#endif
#ifdef clampki
            h_cells->Cells.ki[index] = 138.0;
#endif
            h_cells->Cells.VARIABLE[index] = byitofac*j + minitofac;
            h_cells->caithreshold[index] = bycai*i + mincai;
        }
    }
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(TTCaiBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    h_cells->dobif(dt, t);
    
    for (int i = 0; i < numcai; i++) {
        for (int j = 0; j < numvar; j++) {
            index = numvar*i + j;
            fprintf(allbifs, "%g\t%g", bycai*i + mincai, byitofac*j + minitofac);
            fprintf(allpcls, "%g\t%g", bycai*i + mincai, byitofac*j + minitofac);
            for (int l = REMOVEBEATS; l < BEATS; l++) {
                fprintf(allbifs, "\t%g", h_cells->apds[BEATS*(index) + l]);
                fprintf(allpcls, "\t%g", h_cells->pcls[BEATS*(index) + l]);
            }
            fprintf(allbifs,"\n");
            fprintf(allpcls,"\n");
            fprintf(caibifs, "%g\t%g", bycai*i + mincai, byitofac*j + minitofac);
            fprintf(casrbifs, "%g\t%g", bycai*i + mincai, byitofac*j + minitofac);
            fprintf(kibifs, "%g\t%g", bycai*i + mincai, byitofac*j + minitofac);
            fprintf(naibifs, "%g\t%g", bycai*i + mincai, byitofac*j + minitofac);
            
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
        }
    }
    

    delete h_cells;
    fclose(allbifs);
    fclose(allpcls);
    fclose(caibifs);
    fclose(casrbifs);
    fclose(kibifs);
    fclose(naibifs);
    return 0;
}


