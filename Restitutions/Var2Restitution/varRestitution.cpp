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

#ifndef VARIABLE1
#define VARIABLE1 cai
#endif

#ifndef VARIABLE2
#define VARIABLE2 casr
#endif

template <template<int> class typecell, int ncells, int beats>
varRestitution<typecell, ncells, beats>::varRestitution() {
    commonst = 0.0;
    basepcl = 1000.0;
    curbeat = -1;
    commont = -stimt;
    commondt = 0.1;
    for (int i = 0; i < ncells; i++) {
        VARIABLE1[i] = 10.0;
        VARIABLE2[i] = 10.0;
        S2[i] = -1.0;
    }
    donerestitution = 0;
}

template <template<int> class typecell, int ncells, int beats>
void varRestitution<typecell, ncells, beats>::dorest()
{
    setup();
    double startt = commont;
    while (donerestitution < ncells) {
        iterateall(commondt, commont, startt);
        commont = commont + commondt;
    }
}

template <template<int> class typecell, int ncells, int beats>
void varRestitution<typecell, ncells, beats>::iterate(const int id, double dt, double t, double startt) {
    if (S2[id] < 0) {
        double vold = AllCells.v[id];
        double restst = (fmod((t + dt/4.0), basepcl) < (stimduration - dt/4.0)) ? stimulus : 0.0;
        AllCells.stepdt(id, dt, restst);
        double vnew = AllCells.v[id];
        if (vold >= threshold && vnew < threshold) {
            S2[id] = t - startt;
            donerestitution = donerestitution + 1;
            printf("%d\n", donerestitution);
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void varRestitution<typecell, ncells, beats>::iterateall(double dt, double t, double startt) {
    for (int i = 0; i < ncells; i++) {
        iterate(i, dt, t, startt);
    }
}

template <template<int> class typecell, int ncells, int beats>
void varRestitution<typecell, ncells, beats>::setup()
{
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
    } //Goes until the end of the last pre-pacing beat
    while (fmod((commont + commondt/4.0), basepcl) > stimduration - commondt/4.0) {
        CommonCell.stepdt(0, commondt, 0);
        commont = commont + commondt;
    } //Goes until right when the stimulus is about to occur
    while (CommonCell.v[0] < threshold) {
        CommonCell.stepdt(0, commondt, stimulus);
        commont = commont + commondt;
    } //Goes until slightly above the threshold
    for (int i = 0; i < ncells; i++) {
        AllCells.setcell(i, &CommonCell);
        AllCells.VARIABLE1[i] = VARIABLE1[i];
        AllCells.VARIABLE2[i] = VARIABLE2[i];
    } //Sets up the restitution
}
#endif //varRestitution_cpp
