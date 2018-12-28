//
//  varRestitution.cpp
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the previous diastolic interval (DI).
//  Key modifiable variable is "basepcl", the prior pacing cycle length (PCL)
//  With default value basepcl = 1000.0;
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "varRestitution.h"

#ifndef varRestitution_cpp
#define varRestitution_cpp

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

#ifndef VARIABLE
#define VARIABLE cai
#endif

template <template<int> class typecell, int ncells, int beats>
varRestitution<typecell, ncells, beats>::varRestitution() {
    commonst = 0.0;
    basepcl = 1000.0;
    curbeat = -1;
    commont = -stimt;
    commondt = 0.1;
    for (int i = 0; i < ncells; i++) {
        VARIABLE[i] = 10.0;
        S2[i] = -1.0;
    }
    donerestitution = 0;
}

template <template<int> class typecell, int ncells, int beats>
void varRestitution<typecell, ncells, beats>::dorest(FILE* ap, int id1, int id2)
{
    setup(ap);
    long double startt = commont;
    double t_save = commont;
    while (donerestitution < ncells) {
        iterateall(commondt, commont, startt);
        commont = commont + commondt;
        if (commont > t_save) {
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)commont,AllCells.v[id1],AllCells.VARIABLE[id1],AllCells.v[id2],AllCells.VARIABLE[id2]);
            t_save = t_save + 1;
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void varRestitution<typecell, ncells, beats>::iterate(const int id, long double dt, long double t, long double startt) {
    if (S2[id] < 0) {
        double vold = AllCells.v[id];
        double restst = ((fmod((t + dt/4.0), basepcl) < (stimduration - dt/4.0)) && (t < startt + basepcl/2)) ? stimulus : 0.0;
        AllCells.stepdt(id, dt, restst);
        double vnew = AllCells.v[id];
        if (vold >= threshold && vnew < threshold) {
            S2[id] = t - startt;
            donerestitution = donerestitution + 1;
            printf("%d\n", donerestitution);
        }
    }
    else {
        AllCells.stepdt(id, dt, 0);
    }
}

template <template<int> class typecell, int ncells, int beats>
void varRestitution<typecell, ncells, beats>::iterateall(long double dt, long double t, long double startt) {
    for (int i = 0; i < ncells; i++) {
        iterate(i, dt, t, startt);
    }
}

template <template<int> class typecell, int ncells, int beats>
void varRestitution<typecell, ncells, beats>::setup(FILE* ap)
{
    double t_save = commont;
    while (curbeat < beats) {
        double vold = CommonCell.v[0];
        commonst = (commont > -commondt/4.0 && fmod((commont + commondt/4.0), basepcl) < (stimduration - commondt/4.0) && (CommonCell.v[0] < threshold || commonst < -1.0)) ? stimulus : 0.0;
        CommonCell.stepdt(0, commondt, commonst);
        double vnew = CommonCell.v[0];
        if (vold >= threshold && vnew < threshold) {
            curbeat += 1;
            printf("%d\n", curbeat);
        }
        commont = commont + commondt;
        if (commont > t_save) {
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)commont,CommonCell.v[0],CommonCell.VARIABLE[0],CommonCell.v[0],CommonCell.VARIABLE[0]);
            t_save = t_save + 1;
        }
    } //Goes until the end of the last pre-pacing beat
    while (fmod((commont + commondt/4.0), basepcl) > stimduration - commondt/4.0) {
        CommonCell.stepdt(0, commondt, 0);
        commont = commont + commondt;
        if (commont > t_save) {
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)commont,CommonCell.v[0],CommonCell.VARIABLE[0],CommonCell.v[0],CommonCell.VARIABLE[0]);
            t_save = t_save + 1;
        }
    } //Goes until right when the stimulus is about to occur
    while (CommonCell.v[0] < threshold) {
        CommonCell.stepdt(0, commondt, stimulus);
        commont = commont + commondt;
        if (commont > t_save) {
            fprintf(ap,"%g\t%g\t%g\t%g\t%g\n",(double)commont,CommonCell.v[0],CommonCell.VARIABLE[0],CommonCell.v[0],CommonCell.VARIABLE[0]);
            t_save = t_save + 1;
        }
    } //Goes until slightly above the threshold
    for (int i = 0; i < ncells; i++) {
        AllCells.setcell(i, &CommonCell);
        AllCells.VARIABLE[i] = VARIABLE[i];
    }
    //Sets up the restitution
}
#endif //varRestitution_cpp
