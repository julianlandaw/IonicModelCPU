//
//  APDBifurcation.cpp
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "APDBifurcation.h"

#ifndef APDBifurcation_cpp
#define APDBifurcation_cpp

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

template <template<int> class typecell, int ncells, int beats>
APDBifurcation<typecell, ncells, beats>::APDBifurcation() {
    int i;
    for (i = 0; i < ncells*beats; i++) {
        apds[i] = 0.1;
    }
#ifdef TT
    for (i = 0; i < 2*ncells*beats; i++) {
        kis[i] = 0;
        nais[i] = 0;
        cais[i] = 0;
        casrs[i] = 0;
    }
#endif
    for (i = 0; i < ncells; i++) {
        pcls[i] = 0.0;
        startapds[i] = 0.0;
        st[i] = 0.0;
        curbeat[i] = -1;
    }
    donebifurcation = 0;
}

template <template<int> class typecell, int ncells, int beats>
void APDBifurcation<typecell, ncells, beats>::iterate(const int id, long double dt, long double t) {
    if (id < ncells && curbeat[id] < beats - 1) {
        double vold = Cells.v[id];
        st[id] = (t > -dt/4.0 && fmod((t + dt/4.0), pcls[id]) < (stimduration - dt/4.0) && (Cells.v[id] < threshold || st[id] < -1.0)) ? stimulus : 0.0;
        Cells.stepdt(id, dt, st[id]);
        double vnew = Cells.v[id];
        if (vold < threshold && vnew >= threshold) {
            startapds[id] = t;
#ifdef TT
            if (curbeat[id] >= -1 && curbeat[id] < beats - 1) {
                kis[2*(beats*id + curbeat[id] + 1)] = Cells.ki[id];
                nais[2*(beats*id + curbeat[id] + 1)] = Cells.nai[id];
                cais[2*(beats*id + curbeat[id] + 1)] = Cells.cai[id];
                casrs[2*(beats*id + curbeat[id] + 1)] = Cells.casr[id];
            }
#endif
        }
        else if (vold >= threshold && vnew < threshold) {
            curbeat[id] += 1;
#ifdef TT
            kis[2*(beats*id + curbeat[id]) + 1] = Cells.ki[id];
            nais[2*(beats*id + curbeat[id]) + 1] = Cells.nai[id];
            cais[2*(beats*id + curbeat[id]) + 1] = Cells.cai[id];
            casrs[2*(beats*id + curbeat[id]) + 1] = Cells.casr[id];
#endif
            apds[beats*id + curbeat[id]] = t - startapds[id];
            if (curbeat[id] == beats - 1) {
                donebifurcation = donebifurcation + 1;
                //printf("%d\n",donebifurcation);
            }
            
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDBifurcation<typecell, ncells, beats>::iterateall(long double dt, long double t) {
    for (int i = 0; i < ncells; i++) {
#ifdef TT
	Cells.ki[i] = 138.0;
	Cells.nai[i] = 12.0;
#endif
	iterate(i, dt, t);
#ifdef TT
	Cells.ki[i] = 138.0;
	Cells.nai[i] = 12.0;
#endif
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDBifurcation<typecell, ncells, beats>::dobif(long double dt, double startt) {
    long double t = startt;
    long double t_save = 0;
    while (donebifurcation < ncells) {
        iterateall(dt, t);
        t = t + dt;
        if (t > t_save - dt/2) {
            //printf("%g\n",(double)t);
            t_save = t_save + 1000;
        }
    }
}
#endif //APDBifurcation_cpp
