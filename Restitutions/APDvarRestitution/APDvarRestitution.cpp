//
//  APDvarRestitution.cpp
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the previous diastolic interval (DI).
//  Key modifiable variable is "basepcl", the prior pacing cycle length (PCL)
//  With default value basepcl = 1000.0;
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "APDvarRestitution.h"

#ifndef APDvarRestitution_cpp
#define APDvarRestitution_cpp

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
APDvarRestitution<typecell, ncells, beats>::APDvarRestitution() {
    st = 0.0;
    basepcl = 1000.0;
    curbeat = -1;
    commont = -stimt;
    commondt = 0.1;
    for (int i = 0; i < ncells; i++) {
        notstartrestitution[i] = true;
        VARIABLE[i] = 0.0;
        DI[i] = 10.0;
        S2[i] = -1.0;
    }
    donerestitution = 0;
}

template <template<int> class typecell, int ncells, int beats>
void APDvarRestitution<typecell, ncells, beats>::dorest()
{
    setup();
    double t = 0;
    while (donerestitution < ncells) {
        iterateall(commondt, t);
        t = t + commondt;
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDvarRestitution<typecell, ncells, beats>::iterate(const int id, double dt, double t) {
    if (S2[id] < 0) {
        double vold = AllCells.v[id];
        double restst;
        if (t > DI[id] -dt/4.0 && t < DI[id] + stimduration - dt/4.0) {
            restst = stimulus;
            if (notstartrestitution[id] && vold >= threshold) {
                AllCells.VARIABLE[id] = VARIABLE[id];
                notstartrestitution[id] = false;
            }
        }
        else {restst = 0;}
        AllCells.stepdt(id, dt, restst);
        double vnew = AllCells.v[id];
        if (vold >= threshold && vnew < threshold) {
            S2[id] = t - DI[id];
            donerestitution = donerestitution + 1;
            printf("%d\n", donerestitution);
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDvarRestitution<typecell, ncells, beats>::iterateall(double dt, double t) {
    for (int i = 0; i < ncells; i++) {
        iterate(i, dt, t);
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDvarRestitution<typecell, ncells, beats>::setup()
{
    while (curbeat < beats) {
        double vold = CommonCell.v[0];
        st = (commont > -commondt/4.0 && fmod((commont + commondt/4.0), basepcl) < (stimduration - commondt/4.0) && (CommonCell.v[0] < threshold || st < -1.0)) ? stimulus : 0.0;
        CommonCell.stepdt(0, commondt, st);
        double vnew = CommonCell.v[0];
        if (vold >= threshold && vnew < threshold) {
            curbeat += 1;
            printf("%d\n", curbeat);
        }
        commont = commont + commondt;
    }
    for (int i = 0; i < ncells; i++) {
        AllCells.setcell(i, &CommonCell);
    }
}
#endif //APDvarRestitution_cpp
