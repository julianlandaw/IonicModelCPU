//
//  voltageclampmain.cpp
//
//  Main file for generating APD bifurcation diagrams.
//  .txt output files are organized as:
//  1st column:         pacing cycle length (PCL)
//  2nd column:         1st APD (for set number of beats BEATS1)
//  3rd column:         2nd APD (for set number of beats BEATS2)
//  The number of remaining columns depends on the macro BEATS
//  where #remaining columns = BEATS + BEATS
//  One later on can modify the code so it's not the same number of beats
//  for each APD
//
//  Created by Julian Landaw on 1/7/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

//#define ibarca (1.2*0.09)

//#define TT
#ifdef TT
#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
#include "Cells/TTCellIto.cpp"
#elif LR2
#define TYPECELL LR2CellIto
#define TYPECELLSTRING "LR2"
#include "Cells/LR2CellIto.cpp"
#else
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"
#endif

#ifndef NCELLS
#define NCELLS 300
#endif

#ifndef stimt
#define stimt 100.0
#endif

#ifndef curfac
#define curfac itofac
#endif

#include "Cable1D.cpp"

int main(int argc, char *argv[])
{
    const double pcl = atof(argv[1]);
    const int numstimulus = atoi(argv[2]);
    const double curfac1 = atof(argv[3]);
    const double curfac2 = atof(argv[4]);
    const double dx1 = atof(argv[5]);
    const double dx2 = atof(argv[6]);
    const double hill1 = atof(argv[7]);
    const double hill2 = atof(argv[8]);
    const double dt = atof(argv[9]);
    
    Cable1D<TYPECELL, NCELLS>* h_cells;
    
    h_cells = new Cable1D<TYPECELL, NCELLS>(pcl, numstimulus);
    
    for (int i = 0; i < NCELLS; i++) {
        //h_cells->Cells.curfac[i] = curfac1 + i*(curfac2 - curfac1)/(NCELLS - 1);
        //h_cells->dx[i] = dx1 + i*(dx2 - dx1)/(NCELLS - 2);
        if (hill1 > 1e-6) {
            h_cells->Cells.curfac[i] = curfac2 + (curfac1 - curfac2)/(1.0 + exp((i + 1.0 - NCELLS/2.0)/hill1));
        }
        else {
            h_cells->Cells.curfac[i] = curfac1 + i*(curfac2 - curfac1)/(NCELLS - 1.0);
        }

       
        if (hill2 > 1e-6) {
            h_cells->dx[i] = dx2 + (dx1 - dx2)/(1.0 + exp((i + 1.0 - NCELLS/2.0)/hill2));
        }
        else {
           h_cells->dx[i] = dx1 + i*(dx2 - dx1)/(NCELLS - 1.0);
        }
        
        if (h_cells->dx[i] < 0) {
            h_cells->dx[i] = 0.0;
        }
#ifdef TT
        h_cells->Cells.ibarcafac[i] = 1.2;
#endif

#ifdef LR1nsca
    h_cells->Cells.tauXfac[i] = 5.0;
    h_cells->Cells.inscafac[i] = 0.0575;
#endif

#ifdef LR2
    h_cells->Cells.icalfac[i] = 0.5;
#endif
    }	

    FILE *cable1d;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sallap_PCL_%g_dx1_%g_dx2_%g_ito1_%g_ito2_%g.txt", TYPECELLSTRING, pcl, dx1, dx2, curfac1, curfac2);
    cable1d = fopen(fileSpec1, "w");
    
#ifdef TT
    FILE *cable1d_cai;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sallca_PCL_%g_dx1_%g_dx2_%g_ito1_%g_ito2_%g.txt", TYPECELLSTRING, pcl, dx1, dx2, curfac1, curfac2);
    cable1d_cai = fopen(fileSpec2, "w");
#endif
#ifdef LR2
    FILE *cable1d_cai;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sallca_PCL_%g_dx1_%g_dx2_%g_ito1_%g_ito2_%g.txt", TYPECELLSTRING, pcl, dx1, dx2, curfac1, curfac2);
    cable1d_cai = fopen(fileSpec2, "w");
#endif

    fprintf(cable1d,"0");
    for (int i = 0; i < NCELLS; i++) {
        fprintf(cable1d,"\t%g",h_cells->Cells.curfac[i]);
#ifdef TT
        fprintf(cable1d_cai,"\t%g",h_cells->Cells.curfac[i]);
#endif
#ifdef LR2
        fprintf(cable1d_cai,"\t%g",h_cells->Cells.curfac[i]);
#endif        
    }
    fprintf(cable1d,"\n");
#ifdef TT
    fprintf(cable1d_cai,"\n");
#endif

    printf("Byte Size of Cells on Device: %lu\n", sizeof(Cable1D<TYPECELL, NCELLS>) );
    //Now run the program

    double t = -stimt;
    double t_save = 0;
    double other_tsave = 0;
    printf("test");
    while (t < pcl*(600) + pcl - 100) {
        h_cells->diffuseall();
        h_cells->iterateall(dt, t);
        t = t + dt;
        if (t > t_save - dt/4.0) {
            fprintf(cable1d, "%g", t);
#ifdef TT
            fprintf(cable1d_cai,"\t%g", t);
#endif
#ifdef LR2
            fprintf(cable1d_cai,"\t%g", t);
#endif
            for (int i = 0; i < NCELLS; i++) {
                fprintf(cable1d, "\t%g", h_cells->Cells.v[i]);
#ifdef TT
                fprintf(cable1d_cai,"\t%g", h_cells->Cells.cai[i]);
#endif
#ifdef LR2
                fprintf(cable1d_cai,"\t%g", h_cells->Cells.cai[i]);
#endif
            }
            fprintf(cable1d, "\n");
#ifdef TT
            fprintf(cable1d_cai, "\n");
#endif
#ifdef LR2
            fprintf(cable1d_cai, "\n");
#endif
            t_save = t_save + 5;
        }
        if (t > other_tsave - dt/4.0) {
            printf("%g\n",t);
            other_tsave = other_tsave + 1;
        }
    }
    
    delete h_cells;
    fclose(cable1d);
#ifdef TT
    fclose(cable1d_cai);
#endif
#ifdef LR2
    fclose(cable1d_cai);
#endif
    return 0;
}


