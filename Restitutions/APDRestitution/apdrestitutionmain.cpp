//
//  apdrestitutionmain.cpp
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

//#define TT
#ifdef TT
#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
#define gnaped 0.01
#include "Cells/TTCellIto.cpp"

#elif TP06
#define TYPECELL TP06Cell
#define TYPECELLSTRING "TP06"
#include "Cells/TP06Cell.cpp"

#elif LR1
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"
#elif LR2
#define TYPECELL LR2CellIto
#define TYPECELLSTRING "LR2"
#include "Cells/LR2CellIto.cpp"
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
#elif OHara
#define TYPECELL OHaraCell
#define TYPECELLSTRING "OHara"
#include "Cells/OHaraCell.cpp"
#else
#define TYPECELL TTCellMod
#define TYPECELLSTRING "TTMod"
#include "Cells/TTCellMod.cpp"
#endif

#ifndef NCELLS
#define NCELLS 4000
#endif

#ifndef BEATS
#define BEATS 200
#endif

#ifndef stimt
#define stimt 100.0
#endif

#ifndef curfac
#define curfac itofac
#endif

#include "APDRestitution.cpp"

int main(int argc, char *argv[])
{
    const double itofac = atof(argv[1]);
    const double pcl = atof(argv[2]);
    const double startdi = atof(argv[3]);
    const double maxdi = atof(argv[4]);
    const long double dt = atof(argv[5]);
    
    double bydi = (NCELLS > 1) ? (maxdi - startdi)/(NCELLS - 1) : 1;
    
    APDRestitution<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new APDRestitution<TYPECELL, NCELLS, BEATS>();
#ifndef Fenton
//    h_cells->CommonCell.itofac[0] = itofac;
#endif
    
#ifdef LR2
    h_cells->CommonCell.itofac[0] = 0;
    h_cells->CommonCell.iskfac[0] = 0.0; //20.0;
    h_cells->CommonCell.skh[0] = 0.0025;
    h_cells->CommonCell.iupfac[0] = 1.0; // 3.0;
    h_cells->CommonCell.tauyfac[0] = 1.0;
#endif
    
#ifdef LR1
    h_cells->CommonCell.tauXfac[0] = 1.0;
    h_cells->CommonCell.yshift[0] = -8.0; //-8.0; //0.0;
    h_cells->CommonCell.itofac[0] = 0.0;
    h_cells->CommonCell.icalfac[0] = 1.0;
    h_cells->CommonCell.ikfac[0] = 1.0;
    h_cells->CommonCell.ikifac[0] = 1.0;
    h_cells->CommonCell.itoslowfac[0] = 0.0;
#endif
#ifdef TT
#ifdef clampnai
    h_cells->CommonCell.nai[0] = 12.0;
#endif
#ifdef clampki
    h_cells->CommonCell.ki[0] = 138.0;
#endif
    //h_cells->CommonCell.ibarcafac[0] = 6e-4/(gcal);
    //h_cells->CommonCell.ikrfac[0] = 0.01/(gkr);
    //h_cells->CommonCell.iksfac[0] = 0.036/(gks);
    
#endif
#ifdef OHara
    h_cells->CommonCell.skh[0] = 0.006;
    h_cells->CommonCell.iskfac[0] = 0.0;
    h_cells->CommonCell.inafac[0] = 2.0;
    h_cells->CommonCell.itofac[0] = 0.0;
    h_cells->CommonCell.nai[0] = 9.0;
    h_cells->CommonCell.nass[0] = h_cells->CommonCell.nai[0];
    h_cells->CommonCell.yshift[0] = -8.0;
#endif
        
    h_cells->CommonCell.curfac[0] = itofac;
    printf("%g\n",h_cells->CommonCell.curfac[0]);
    h_cells->basepcl = pcl;
    h_cells->commondt = dt;
    
    FILE *rest;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%srestPCL_%g_ITO_%g_DI1_%g.txt", TYPECELLSTRING, pcl, itofac, startdi);
    rest = fopen(fileSpec1, "w");
    
    FILE *ap;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%saprestPCL_%g_ITO_%g_DI1_%g.txt", TYPECELLSTRING, pcl, itofac, startdi);
    ap = fopen(fileSpec2, "w");
    
    for (int i = 0; i < NCELLS; i++) {
        h_cells->DI[i] = startdi + bydi*i;
    }

    int id1 = 0; 
    int id2 = NCELLS - 1;    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(APDRestitution<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    h_cells->dorest(ap, id1, id2);
    
    for (int i = 0; i < NCELLS; i++) {
        fprintf(rest,"%g\t%g\n", h_cells->DI[i], h_cells->S2[i]);
    }

    delete h_cells;
    fclose(rest);
    fclose(ap);
    return 0;
}


