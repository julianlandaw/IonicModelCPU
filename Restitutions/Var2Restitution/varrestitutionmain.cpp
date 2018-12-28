//
//  restitutionmain.cpp
//
//  Main file for generating S1S2 Restitutition curves.
//  .txt output files are organized as:
//  1st column:         diastolic interval (DI)
//  2nd column:         action potential duration (APD, also called S2)
//  The number of rows depends on the macro NCELLS
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

#ifdef OHara
#define TYPECELL OHaraCell
#define TYPECELLSTRING "OHara"
#include "Cells/OHaraCell.cpp"

#elif UCLA
#define TYPECELL UCLACell
#define TYPECELLSTRING "UCLA"
#include "Cells/UCLACell.cpp"

#elif LR2
#define TYPECELL LR2CellIto
#define TYPECELLSTRING "LR2"
#include "Cells/LR2CellIto.cpp"

#elif TT
#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
#include "Cells/TTCellIto.cpp"

#else
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"
#endif

#ifndef NUMVAR1
#define NUMVAR1 256
#endif

#ifndef NUMVAR2
#define NUMVAR2 256
#endif

#define NCELLS (NUMVAR1*NUMVAR2)

#ifndef BEATS
#define BEATS 1024
#endif

#ifndef stimt
#define stimt 100.0
#endif

#ifndef VARIABLE1
#define VARIABLE1 cai
#endif

#ifndef VARIABLE2
#define VARIABLE2 casr
#endif

#ifndef curfac
#define curfac itofac
#endif

#include "varRestitution.cpp"

int main(int argc, char *argv[])
{
    const double curfac = atof(argv[1]);
    const double pcl = atof(argv[2]);
    const double startvar1 = atof(argv[3]);
    const double endvar1 = atof(argv[4]);
    const double startvar2 = atof(argv[5]);
    const double endvar2 = atof(argv[6]);
    const double dt = atof(argv[7]);
    
    double byvar1 = (endvar1 - startvar1)/(NUMVAR1 - 1);
    double byvar2 = (endvar2 - startvar2)/(NUMVAR2 - 1);
    
    varRestitution<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new varRestitution<TYPECELL, NCELLS, BEATS>();
#ifdef LR1
    h_cells->CommonCell.itofac[0] = 0.0;
    h_cells->CommonCell.yshift[0] = 0.0;
    h_cells->CommonCell.icalfac[0] = 1.0;
    h_cells->CommonCell.ikfac[0] = 1.0;
    h_cells->CommonCell.ikifac[0] = 1.0;
    h_cells->CommonCell.tauXfac[0] = 10.0;
#endif
#ifdef OHara
    h_cells->CommonCell.iskfac[0] = 0.0;
    h_cells->CommonCell.itofac[0] = 0.0;
    h_cells->CommonCell.inafac[0] = 2.0;
#endif
#ifdef LR2
    h_cells->CommonCell.iupfac[0] = 3.0;
    h_cells->CommonCell.itofac[0] = 0.0;
    h_cells->CommonCell.iskfac[0] = 0.0;
    h_cells->CommonCell.skh[0] = 0.0025;
#endif
#ifdef TT
    h_cells->CommonCell.itofac[0] = 0.0;
    h_cells->CommonCell.iskfac[0] = 0.0;
#endif 
    h_cells->CommonCell.curfac[0] = curfac;
    h_cells->basepcl = pcl;
    h_cells->commondt = dt;
    
    FILE *rest;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%srestPCL_%g_ITO_%g.txt", TYPECELLSTRING, pcl, curfac);
    rest = fopen(fileSpec1, "w");
    
    FILE *saveconditions;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sVARS_PCL_%g_ITO_%g.txt", TYPECELLSTRING, pcl, curfac);
    saveconditions = fopen(fileSpec2, "w");
    
    for (int i = 0; i < NUMVAR1; i++) {
        for (int j = 0; j < NUMVAR2; j++) {
            h_cells->VARIABLE1[i*NUMVAR2 + j] = startvar1 + byvar1*i;
            h_cells->VARIABLE2[i*NUMVAR2 + j] = startvar2 + byvar2*j;
        }
    }
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(varRestitution<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    h_cells->dorest();
    
    h_cells->CommonCell.saveconditions(saveconditions,0,true,h_cells->commont);
    
    for (int i = 0; i < NCELLS; i++) {
#ifdef UCLA
        fprintf(rest,"%g\t%g\t%g\t%g\n", h_cells->VARIABLE1[i], h_cells->VARIABLE2[i], h_cells->CommonCell.ci[0], h_cells->S2[i]);
#else
        fprintf(rest,"%g\t%g\t%g\n", h_cells->VARIABLE1[i], h_cells->VARIABLE2[i], h_cells->S2[i]);
#endif
    }

    delete h_cells;
    fclose(rest);
    return 0;
}


