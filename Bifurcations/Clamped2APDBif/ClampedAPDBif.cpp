//
//  APDBifurcation.cu
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "ClampedAPDBif.h"

#ifndef ClampedAPDBif_cpp
#define ClampedAPDBif_cpp

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
ClampedAPDBif<typecell, ncells, beats>::ClampedAPDBif() {
    int i;
    for (i = 0; i < ncells*beats; i++) {
        apds[i] = 0.1;
    }
    for (i = 0; i < 2*ncells*beats; i++) {
        var1[i] = 0;
        var2[i] = 0;
    }
    for (i = 0; i < ncells; i++) {
        clampvar3[i] = Cells.VARIABLE3[i];
        clampvar4[i] = Cells.VARIABLE4[i];
        startapds[i] = 0.0;
        st[i] = 0.0;
        curbeat[i] = -1;
        notdonebifurcation[i] = true;
        pcls[i] = 1000.0;
    }
}

template <template<int> class typecell, int ncells, int beats>
void ClampedAPDBif<typecell, ncells, beats>::iterate(const int id, long double dt, long double t) {
    if (notdonebifurcation[id]) {
        double vold = Cells.v[id];
        st[id] = (t > -dt/4.0 && fmod((t + dt/4.0), pcls[id]) < (stimduration - dt/4.0) && (vold < threshold || st[id] < -1.0)) ? stimulus : 0.0;
        Cells.VARIABLE3[id] = clampvar3[id];
        Cells.VARIABLE4[id] = clampvar4[id];
        Cells.stepdt(id, dt, st[id]);
        double vnew = Cells.v[id];
        if (vold < threshold && vnew >= threshold) {
            startapds[id] = t;
            if (curbeat[id] >= -1 && curbeat[id] < beats - 1) {
                var1[2*(beats*id + curbeat[id] + 1)] = Cells.VARIABLE1[id];
                var2[2*(beats*id + curbeat[id] + 1)] = Cells.VARIABLE2[id];
            }
        }
        else if (vold >= threshold && vnew < threshold) {
            curbeat[id] += 1;
            if (curbeat[id] >= 0 && curbeat[id] < beats) {
                var1[2*(beats*id + curbeat[id]) + 1] = Cells.VARIABLE1[id];
                var2[2*(beats*id + curbeat[id]) + 1] = Cells.VARIABLE2[id];
            }
            apds[beats*id + curbeat[id]] = t - startapds[id];
            if (curbeat[id] == beats - 1) {
                notdonebifurcation[id] = false;
            }
            
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void ClampedAPDBif<typecell, ncells, beats>::iterateall(long double dt, long double t) {
    for (int i = 0; i < ncells; i++) {
        iterate(i,dt,t);
    }
}
#endif //ClampedAPDBif
