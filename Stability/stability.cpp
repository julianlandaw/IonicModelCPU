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
#elif TTSK
#define TYPECELL TTCellItoSK
#define TYPECELLSTRING "TTSK"
#include "Cells/TTCellItoSK.cpp"
#elif LR1
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"
#elif LR1SK
#define TYPECELL LR1CellItoSK
#define TYPECELLSTRING "LR1SK"
#include "Cells/LR1CellItoSK.cpp"
#elif LR1_taud
#define TYPECELL LR1CellIto_taud
#define TYPECELLSTRING "LR1_taud"
#include "Cells/LR1CellIto_taud.cpp"
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
#define BEATS 100
#endif

#ifndef SAVEBEATS
#define SAVEBEATS 16
#endif

#ifndef stimt
#define stimt 100.0L
#endif

#ifndef curfac
#define curfac itofac
#endif

#include "APDBifurcation.cpp"

int main(int argc, char *argv[])
{
    const double minitofac = atof(argv[1]);
    const double maxitofac = atof(argv[2]);
    const double minpcl = atof(argv[3]);
    const double maxpcl = atof(argv[4]);
    const long double dt = atof(argv[5]);
    
    double byitofac = (numitofac > 1) ? (maxitofac - minitofac)/(numitofac - 1) : 1;
    double bypcl = (numpcl > 1) ? (maxpcl - minpcl)/(numpcl - 1) : 1;
    
    APDBifurcation<TYPECELL, 1, BEATS>* h_cells;
    
    h_cells = new APDBifurcation<TYPECELL, 1, BEATS>();
    
    FILE *stability;
    char fileSpec1[100];
#ifdef TTSK
    snprintf(fileSpec1, 100, "%sstab_PCL_%g_ISK_%g.txt", TYPECELLSTRING, minpcl, minitofac);
#elif LR1SK
    snprintf(fileSpec1, 100, "%sstab_PCL_%g_ISK_%g.txt", TYPECELLSTRING, minpcl, minitofac);
#elif LR1_taud
    snprintf(fileSpec1, 100, "%sstab_PCL_%g_taud_%g.txt", TYPECELLSTRING, minpcl, minitofac);
#else
    snprintf(fileSpec1, 100, "%sstab_PCL_%g_ITO_%g.txt", TYPECELLSTRING, minpcl, minitofac);
#endif
    stability = fopen(fileSpec1, "w");

    //printf("Byte Size of Cells on Device: %lu\n", sizeof(APDBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    double minapd;
    double maxapd;
    double minapd2_1;
    double maxapd2_1;
    double minapd2_2;
    double maxapd2_2;
    double minapd3_1;
    double maxapd3_1;
    double minapd3_2;
    double maxapd3_2;
    double minapd3_3;
    double maxapd3_3;
    double minapd4_1;
    double maxapd4_1;
    double minapd4_2;
    double maxapd4_2;
    double minapd4_3;
    double maxapd4_3;
    double minapd4_4;
    double maxapd4_4;
    
    int j_sim;
    long double t;
    
    for (int i = 0; i < numpcl; i++) {
        for (int j = 0; j < numitofac; j++) {
            h_cells = new APDBifurcation<TYPECELL, 1, BEATS>();
            t = -stimt;
#ifdef LR1_taud
            h_cells->Cells.tauXfac[0] = 5.0;
	    h_cells->Cells.itofac[0] = 1.05;
#endif
            h_cells->Cells.curfac[0] = byitofac*j + minitofac;
            h_cells->pcls[0] = bypcl*i + minpcl;
            h_cells->dobif(dt, t);            
 
            minapd = h_cells->apds[BEATS - SAVEBEATS];
            maxapd = h_cells->apds[BEATS - SAVEBEATS];
            minapd2_1 = h_cells->apds[BEATS - SAVEBEATS];
            maxapd2_1 = h_cells->apds[BEATS - SAVEBEATS];
            minapd2_2 = h_cells->apds[BEATS - SAVEBEATS + 1];
            maxapd2_2 = h_cells->apds[BEATS - SAVEBEATS + 1];
            minapd3_1 = h_cells->apds[BEATS - SAVEBEATS];
            maxapd3_1 = h_cells->apds[BEATS - SAVEBEATS];
            minapd3_2 = h_cells->apds[BEATS - SAVEBEATS + 1];
            maxapd3_2 = h_cells->apds[BEATS - SAVEBEATS + 1];
            minapd3_3 = h_cells->apds[BEATS - SAVEBEATS + 2];
            maxapd3_3 = h_cells->apds[BEATS - SAVEBEATS + 2];
            minapd4_1 = h_cells->apds[BEATS - SAVEBEATS];
            maxapd4_1 = h_cells->apds[BEATS - SAVEBEATS];
            minapd4_2 = h_cells->apds[BEATS - SAVEBEATS + 1];
            maxapd4_2 = h_cells->apds[BEATS - SAVEBEATS + 1];
            minapd4_3 = h_cells->apds[BEATS - SAVEBEATS + 2];
            maxapd4_3 = h_cells->apds[BEATS - SAVEBEATS + 2];
            minapd4_4 = h_cells->apds[BEATS - SAVEBEATS + 3];
            maxapd4_4 = h_cells->apds[BEATS - SAVEBEATS + 3];
            
            for (j_sim = BEATS - SAVEBEATS + 1; j_sim < BEATS; j_sim++) {
                if (h_cells->apds[j_sim] < minapd) {
                    minapd = h_cells->apds[j_sim];
                }
                if (h_cells->apds[j_sim] > maxapd) {
                    maxapd = h_cells->apds[j_sim];
                }
            }
            
            for (j_sim = BEATS - SAVEBEATS + 2; j_sim < BEATS - 1; j_sim = j_sim + 2) {
                if (h_cells->apds[j_sim] < minapd2_1) {
                    minapd2_1 = h_cells->apds[j_sim];
                }
                if (h_cells->apds[j_sim] > maxapd2_1) {
                    maxapd2_1 = h_cells->apds[j_sim];
                }
                if (h_cells->apds[j_sim+1] < minapd2_2) {
                    minapd2_2 = h_cells->apds[j_sim+1];
                }
                if (h_cells->apds[j_sim+1] > maxapd2_2) {
                    maxapd2_2 = h_cells->apds[j_sim+1];
                }
            }
            
            for (j_sim = BEATS - SAVEBEATS + 3; j_sim < BEATS - 2; j_sim = j_sim + 3) {
                if (h_cells->apds[j_sim] < minapd3_1) {
                    minapd3_1 = h_cells->apds[j_sim];
                }
                if (h_cells->apds[j_sim] > maxapd3_1) {
                    maxapd3_1 = h_cells->apds[j_sim];
                }
                if (h_cells->apds[j_sim+1] < minapd3_2) {
                    minapd3_2 = h_cells->apds[j_sim+1];
                }
                if (h_cells->apds[j_sim+1] > maxapd3_2) {
                    maxapd3_2 = h_cells->apds[j_sim+1];
                }
                if (h_cells->apds[j_sim+2] < minapd3_3) {
                    minapd3_3 = h_cells->apds[j_sim+2];
                }
                if (h_cells->apds[j_sim+2] > maxapd3_3) {
                    maxapd3_3 = h_cells->apds[j_sim+2];
                }
            }
            
            for (j_sim = BEATS - SAVEBEATS + 4; j_sim < BEATS - 3; j_sim = j_sim + 4) {
                if (h_cells->apds[j_sim] < minapd4_1) {
                    minapd4_1 = h_cells->apds[j_sim];
                }
                if (h_cells->apds[j_sim] > maxapd4_1) {
                    maxapd4_1 = h_cells->apds[j_sim];
                }
                if (h_cells->apds[j_sim+1] < minapd4_2) {
                    minapd4_2 = h_cells->apds[j_sim+1];
                }
                if (h_cells->apds[j_sim+1] > maxapd4_2) {
                    maxapd4_2 = h_cells->apds[j_sim+1];
                }
                if (h_cells->apds[j_sim+2] < minapd4_3) {
                    minapd4_3 = h_cells->apds[j_sim+2];
                }
                if (h_cells->apds[j_sim+2] > maxapd4_3) {
                    maxapd4_3 = h_cells->apds[j_sim+2];
                }
                if (h_cells->apds[j_sim+3] < minapd4_4) {
                    minapd4_4 = h_cells->apds[j_sim+3];
                }
                if (h_cells->apds[j_sim+3] > maxapd4_4) {
                    maxapd4_4 = h_cells->apds[j_sim+3];
                }
            }
            fprintf(stability,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",h_cells->pcls[0],h_cells->Cells.curfac[0],maxapd,maxapd - minapd, maxapd2_1 - minapd2_1, maxapd2_2 - minapd2_2, maxapd3_1 - minapd3_1,maxapd3_2 - minapd3_2, maxapd3_3 - minapd3_3,maxapd4_1 - minapd4_1,maxapd4_2 - minapd4_2,maxapd4_3-minapd4_3,maxapd4_4 - minapd4_4);
	    printf("%g\t%g\t%g\t%g\n",h_cells->apds[BEATS - SAVEBEATS],h_cells->apds[BEATS - SAVEBEATS+1],h_cells->apds[BEATS - SAVEBEATS+2],h_cells->apds[BEATS - SAVEBEATS+3]);
        }
    }
    
    delete h_cells;
    fclose(stability);
    return 0;
}
