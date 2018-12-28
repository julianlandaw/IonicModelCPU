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

//#define TT
#ifdef TT
#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
#include "Cells/TTCellIto.cpp"
#else
#define TYPECELL LR1CellIto
#define TYPECELLSTRING "LR1"
#include "Cells/LR1CellIto.cpp"
#endif

#ifndef NUMDIS
#define NUMDIS 64
#endif

#ifndef NUMVARS
#define NUMVARS 64
#endif

#define NCELLS (NUMDIS*NUMVARS)

#ifndef BEATS
#define BEATS 1024
#endif

#ifndef stimt
#define stimt 100.0
#endif

#ifndef VARIABLE
#define VARIABLE xr
#endif

#include "APDvarRestitution.cpp"

int main(int argc, char *argv[])
{
    const double itofac = atof(argv[1]);
    const double pcl = atof(argv[2]);
    const double startdi = atof(argv[3]);
    const double enddi = atof(argv[4]);
    const double startvar = atof(argv[5]);
    const double endvar = atof(argv[6]);
    const double dt = atof(argv[7]);

    double bydi = (NUMDIS > 1) ? (enddi - startdi)/(NUMDIS - 1) : 1.0;
    double byvar = (NUMVARS > 1) ? (endvar - startvar)/(NUMVARS - 1) : 1.0;
    
    APDvarRestitution<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new APDvarRestitution<TYPECELL, NCELLS, BEATS>();
    
    h_cells->CommonCell.itofac[0] = itofac;
    h_cells->basepcl = pcl;
    h_cells->commondt = dt;
    
    FILE *rest;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%srestPCL_%g_ITO_%g_DI_%g.txt", TYPECELLSTRING, pcl, itofac, startdi);
    rest = fopen(fileSpec1, "w");
    
    for (int i = 0; i < NUMDIS; i++) {
        for (int j = 0; j < NUMVARS; j++) {
            h_cells->DI[NUMVARS*i + j] = startdi + i*bydi;
            h_cells->VARIABLE[NUMVARS*i + j] = startvar + j*byvar;
        }
    }
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(APDvarRestitution<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    h_cells->dorest();
    
    for (int i = 0; i < NCELLS; i++) {
        fprintf(rest,"%g\t%g\t%g\n", h_cells->DI[i], h_cells->VARIABLE[i], h_cells->S2[i]);
    }

    delete h_cells;
    fclose(rest);
    return 0;
}


