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

#ifndef NCELLS
#define NCELLS 1001
#endif

#ifndef BEATS
#define BEATS 1024
#endif

#ifndef stimt
#define stimt 100.0
#endif

#ifndef VARIABLE
#define VARIABLE xr
#endif

#include "varRestitution.cpp"

int main(int argc, char *argv[])
{
    const double itofac = atof(argv[1]);
    const double pcl = atof(argv[2]);
    const double startvar = atof(argv[3]);
    const double endvar = atof(argv[4]);
    const long double dt = atof(argv[5]);
    
    double byvar = (NCELLS > 1) ? (endvar - startvar)/(NCELLS - 1) : 1;
    
    varRestitution<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new varRestitution<TYPECELL, NCELLS, BEATS>();
    
    //h_cells->CommonCell.itofac[0] = itofac;
#ifdef LR1
    h_cells->CommonCell.tauXfac[0] = 1.0; //5.0;
    h_cells->CommonCell.yshift[0] = -8.0; //8.0;
    h_cells->CommonCell.icalfac[0] = 1.0; //1.15;
    h_cells->CommonCell.itofac[0] = 0.0; //1.05;
    h_cells->CommonCell.ikfac[0] = 1.0;
    h_cells->CommonCell.ikifac[0] = 1.0; //2.2; 
#endif
#ifdef TT
#ifdef clampki
    h_cells->CommonCell.ki[0] = 138.0;
#endif
#ifdef clampnai
    h_cells->CommonCell.nai[0] = 12.0;
#endif
#endif
    h_cells->CommonCell.itofac[0] = itofac;
    h_cells->basepcl = pcl;
    h_cells->commondt = dt;
    
    FILE *rest;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%srestPCL_%g_ITO_%g.txt", TYPECELLSTRING, pcl, itofac);
    rest = fopen(fileSpec1, "w");
    
    FILE *saveconditions;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%sVARS_PCL_%g_ITO_%g.txt", TYPECELLSTRING, pcl, itofac);
    saveconditions = fopen(fileSpec2, "w");
    
    FILE *ap;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%saprestPCL_%g_ITO_%g.txt", TYPECELLSTRING, pcl, itofac);
    ap = fopen(fileSpec3, "w");
    
    int id1, id2;
    id1 = 0;
    id2 = 0;
    for (int i = 0; i < NCELLS; i++) {
        h_cells->VARIABLE[i] = startvar + byvar*i;
	if (fabs(h_cells->VARIABLE[i] - 12.0) < byvar/4) {
            id2 = i;
	}
#ifdef LR1
        if (fabs(startvar + byvar*i - 0.05) < byvar/4) {
            id2 = i;
        }
#endif
    }
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(varRestitution<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    h_cells->dorest(ap, id1, id2);
    
    h_cells->CommonCell.saveconditions(saveconditions,0,true,h_cells->commont);
    
    for (int i = 0; i < NCELLS; i++) {
        fprintf(rest,"%g\t%g\n", h_cells->VARIABLE[i], h_cells->S2[i]);
    }

    delete h_cells;
    fclose(rest);
    fclose(saveconditions);
    fclose(ap);
    return 0;
}


