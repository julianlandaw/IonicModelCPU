//
//  APDRestitution.cpp
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the previous diastolic interval (DI).
//  Key modifiable variable is "basepcl", the prior pacing cycle length (PCL)
//  With default value basepcl = 1000.0;
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "APDRestitution.h"

#ifndef APDRestitution_cpp
#define APDRestitution_cpp

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
#define stimt 100.0
#endif

template <template<int> class typecell, int ncells, int beats>
APDRestitution<typecell, ncells, beats>::APDRestitution() {
    st = 0.0;
    basepcl = 1000.0;
    curbeat = -1;
    commont = -stimt;
    commondt = 0.1;
    for (int i = 0; i < ncells; i++) {
        DI[i] = 10.0;
        S2[i] = -1.0;
    }
    stimtime = 0.0;
    donerestitution = 0;
}

template <template<int> class typecell, int ncells, int beats>
void APDRestitution<typecell, ncells, beats>::dorest(FILE* ap, int id1, int id2)
{
    setup(ap);
    long double t = 0;
    double t_save = commont;
    while (donerestitution < ncells) {
        iterateall(commondt, t);
        t = t + commondt;
        commont = commont + commondt;
        if (commont > t_save) {
            fprintf(ap,"%g\t%g\t%g\n",(double)commont,AllCells.v[id1],AllCells.v[id2]);
            t_save = t_save + 1;
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDRestitution<typecell, ncells, beats>::iterate(const int id, long double dt, long double t) {
    if (S2[id] < 0) {
        double vold = AllCells.v[id];
        double restst = (t > DI[id] -dt/4.0 && t < DI[id] + stimduration - dt/4.0) ? stimulus : 0.0;
        AllCells.stepdt(id, dt, restst);
        double vnew = AllCells.v[id];
        if (vold >= threshold && vnew < threshold) {
            S2[id] = t - DI[id];
            donerestitution = donerestitution + 1;
            printf("%d\n", donerestitution);
        }
    }
    else if (t < S2[id] + DI[id] + 100) {
        AllCells.stepdt(id, dt, 0);
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDRestitution<typecell, ncells, beats>::iterateall(long double dt, long double t) {
    for (int i = 0; i < ncells; i++) {
        iterate(i, dt, t);
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDRestitution<typecell, ncells, beats>::setup(FILE* ap)
{
    double t_save = commont;
    while (curbeat < beats) {
        double vold = CommonCell.v[0];
        st = (commont > -commondt/4.0 && commont > stimtime - commondt/4.0 && commont < stimtime + stimduration - commondt/4.0) ? stimulus : 0.0;
        CommonCell.stepdt(0, commondt, st);
        double vnew = CommonCell.v[0];
        if (vold >= threshold && vnew < threshold) {
            curbeat += 1;
            stimtime += basepcl;
            while (stimtime < commont - commondt/4.0) {
               stimtime += basepcl;
            }
            printf("%d\n", curbeat);
        }
        commont = commont + commondt;
        if (commont > t_save) {
            fprintf(ap,"%g\t%g",(double)commont,CommonCell.v[0]);
#ifdef TT
            fprintf(ap,"%g\n",CommonCell.cai[0]);
#else
            fprintf(ap,"\n");
#endif
            t_save = t_save + 1;
        }
    }
    for (int i = 0; i < ncells; i++) {
        AllCells.setcell(i, &CommonCell);
    }
}
#endif //APDRestitution_cpp
