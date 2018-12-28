//
//  VoltageClamp.cpp
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "VoltageClamp.h"

#ifndef VoltageClamp_cpp
#define VoltageClamp_cpp

#ifndef threshold
#define threshold -75.0
#endif

#ifndef stimulus
#define stimulus -80.0
#endif

#ifndef stimduration
#define stimduration 0.5
#endif

#ifndef DV_MAX
#define DV_MAX 0.1
#endif

#ifndef ADAPTIVE
#define ADAPTIVE 10
#endif

#ifndef stimt
#define stimt 100.0;
#endif

#ifndef VARIABLE
#define VARIABLE cai
#endif

template <template<int> class typecell, int ncells, int beats>
VoltageClamp<typecell, ncells, beats>::VoltageClamp() {
    for (int i = 0; i < ncells; i++) {
        pcls[i] = 1000.0;
        apds1[i] = 300.0;
        apds2[i] = 300.0;
        voltage1[i] = -85.0;
        voltage2[i] = 0.0;
        curbeat[i] = -1;
        inapd[i] = false;
#ifdef TT
        Cells.itofac[i] = 0.9;
#ifdef clampki
	Cells.ki[i] = 138.0;
#endif
#ifdef clampnai
	Cells.nai[i] = 12.0;
#endif
#endif
        for (int j = 0; j < 2*beats; j++) {
#ifdef TT
            cai[2*beats*i + j] = 0.0;
	    nai[2*beats*i + j] = 0.0;
#endif
#ifdef UCLA
	    ci[2*beats*i + j] = 0.0;
	    nai[2*beats*i + j] = 0.0;
#endif
        }
    }
    donetrain = 0;
    t = -stimt;
    commondt = 0.1;
}

template <template<int> class typecell, int ncells, int beats>
void VoltageClamp<typecell, ncells, beats>::iterate(const int id) {
    if (id < ncells && curbeat[id] < 2*beats) {
        if (t > -commondt/4.0 && curbeat[id] < beats && (fmod(t + commondt/4.0, pcls[id]) < (apds1[id] - commondt/4.0))) {
            Cells.v[id] = voltage2[id];
            if (!inapd[id]) {
                curbeat[id] = curbeat[id] + 1;
                inapd[id] = true;
                //if (curbeat[id] == beats) {
                //    donetrain = donetrain + 1;
                //}
                //else {
#ifdef TT
                    cai[2*(id*2*beats + curbeat[id])] = Cells.cai[id];
		    nai[2*(id*2*beats + curbeat[id])] = Cells.nai[id];
#endif
#ifdef UCLA
		    ci[2*(id*2*beats + curbeat[id])] = Cells.ci[id];
		    nai[2*(id*2*beats + curbeat[id])] = Cells.nai[id];
#endif
                //}
            }
        }
        else if (t > -commondt/4.0 && curbeat[id] >= beats && (fmod(t + commondt/4.0, pcls[id]) < (apds2[id] - commondt/4.0))) {
            Cells.v[id] = voltage2[id];
            if (!inapd[id]) {
                curbeat[id] = curbeat[id] + 1;
                inapd[id] = true;
                if (curbeat[id] == 2*beats) {
                    donetrain = donetrain + 1;
                }
                else {
#ifdef TT
                    cai[2*(id*2*beats + curbeat[id])] = Cells.cai[id];
		    nai[2*(id*2*beats + curbeat[id])] = Cells.nai[id];
#endif
#ifdef UCLA
		    ci[2*(id*2*beats + curbeat[id])] = Cells.ci[id];
		    nai[2*(id*2*beats + curbeat[id])] = Cells.nai[id];
#endif
                }
            }
        }
        else {
            Cells.v[id] = voltage1[id];
            if (inapd[id]) {
                printf("%g\n",(double)t);
#ifdef TT
                cai[2*(id*2*beats + curbeat[id]) + 1] = Cells.cai[id];
		nai[2*(id*2*beats + curbeat[id]) + 1] = Cells.nai[id];
#endif
#ifdef UCLA
		ci[2*(id*2*beats + curbeat[id]) + 1] = Cells.ci[id];
		nai[2*(id*2*beats + curbeat[id]) + 1] = Cells.nai[id];
#endif
                inapd[id] = false;
            }
        }
        //printf("%g\n",t);
//#ifdef TT
//        Cells.ki[id] = 138;
//        Cells.nai[id] = 12;
//#endif
        Cells.stepdt(id, commondt, 0);
    }
}

template <template<int> class typecell, int ncells, int beats>
void VoltageClamp<typecell, ncells, beats>::iterateall(FILE* ap) {
    for (int i = 0; i < ncells; i++) {
        iterate(i);
    }
    t = t + commondt;
}

template <template<int> class typecell, int ncells, int beats>
void VoltageClamp<typecell, ncells, beats>::dospiketrain(FILE* ap) {
    //printf("%g\n",t);
    double t_save = -commondt/2;
    while (donetrain < ncells) {
        iterateall(ap);
        if (t > t_save) {
            //fprintf(ap,"%g\t%g\t%g\n",t,Cells.v[0],Cells.VARIABLE[0]);
            t_save = t_save + 1;
        }
    }
}
#endif //VoltageClamp_cpp
