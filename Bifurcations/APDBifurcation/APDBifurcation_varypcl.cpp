//
//  APDBifurcation.cpp
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#ifndef PRECTYPE
#define PRECTYPE double
#endif

#include "APDBifurcation_varypcl.h"

#ifndef APDBifurcation_varypcl_cpp
#define APDBifurcation_varypcl_cpp

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
APDBifurcation_varypcl<typecell, ncells, beats>::APDBifurcation_varypcl() {
    int i;
    for (i = 0; i < ncells*beats; i++) {
        apds[i] = 0.1;
    }
#ifdef LR1
    for (i = 0; i < 2*ncells*beats; i++) {
        xrs[i] = 0;
    }
#endif
#ifdef LR1_taud
    for (i = 0; i < 2*ncells*beats; i++) {
        xrs[i] = 0;
    }
#endif
#ifdef TT
    for (i = 0; i < 2*ncells*beats; i++) {
        kis[i] = 0;
        nais[i] = 0;
        cais[i] = 0;
        casrs[i] = 0;
    }
#endif
#ifdef UCLA
    for (i = 0; i < 2*ncells*beats; i++) {
        //kis[i] = 0;
        nais[i] = 0;
        cais[i] = 0;
        cass[i] = 0;
        //casrs[i] = 0;
    }
#endif

#ifdef LR2
    for (i = 0; i < 2*ncells*beats; i++) {
        kis[i] = 0;
        nais[i] = 0;
        cais[i] = 0;
    }
#endif
    for (i = 0; i < ncells; i++) {
        pcls1[i] = 0.0;
        pcls2[i] = 0.0;
        startapds[i] = 0.0;
        st[i] = 0.0;
        curbeat[i] = -1;
        stimtime[i] = 0;
    }
    donebifurcation = 0;
}

template <template<int> class typecell, int ncells, int beats>
void APDBifurcation_varypcl<typecell, ncells, beats>::iterate(const int id, long double dt, long double t) {
    if (id < ncells && curbeat[id] < beats) { //if (id < ncells && curbeat[id] < beats - 1) {
        PRECTYPE vold = Cells.v[id];
        st[id] = (t > -dt/4.0 && t > stimtime[id] - dt/4.0 && t < stimtime[id] + stimduration - dt/4.0) ? stimulus : 0.0;
#ifdef TTMod
        Cells.ki[id] = 140.0;
        Cells.nai[id] = 10.0;
#endif
        Cells.stepdt(id, (PRECTYPE)dt, st[id]);

#ifdef TTMod
        Cells.ki[id] = 140.0;
        Cells.nai[id] = 10.0;
#endif
        PRECTYPE vnew = Cells.v[id];
        if (vold < threshold && vnew >= threshold) {
            startapds[id] = t;
            
#ifdef LR1
            if (curbeat[id] >= -1 && curbeat[id] < beats - 1) {
                xrs[2*(beats*id + curbeat[id] + 1)] = Cells.xr[id];
            }
#endif
#ifdef LR1_taud
            if (curbeat[id] >= -1 && curbeat[id] < beats - 1) {
                xrs[2*(beats*id + curbeat[id] + 1)] = Cells.xr[id];
            }
#endif
#ifdef TT
            if (curbeat[id] >= -1 && curbeat[id] < beats - 1) {
                kis[2*(beats*id + curbeat[id] + 1)] = Cells.ki[id];
                nais[2*(beats*id + curbeat[id] + 1)] = Cells.nai[id];
                cais[2*(beats*id + curbeat[id] + 1)] = Cells.cai[id];
                casrs[2*(beats*id + curbeat[id] + 1)] = Cells.casr[id];
            }
#endif
#ifdef UCLA
            if (curbeat[id] >= -1 && curbeat[id] < beats - 1) {
                //kis[2*(beats*id + curbeat[id] + 1)] = Cells.ki[id];
                nais[2*(beats*id + curbeat[id] + 1)] = Cells.nai[id];
                cais[2*(beats*id + curbeat[id] + 1)] = Cells.ci[id];
                cass[2*(beats*id + curbeat[id] + 1)] = Cells.cs[id];
                //casrs[2*(beats*id + curbeat[id] + 1)] = Cells.casr[id];
            }
#endif
#ifdef LR2
            if (curbeat[id] >= -1 && curbeat[id] < beats - 1) {
                kis[2*(beats*id + curbeat[id] + 1)] = Cells.ki[id];
                nais[2*(beats*id + curbeat[id] + 1)] = Cells.nai[id];
                cais[2*(beats*id + curbeat[id] + 1)] = Cells.cai[id];
            }
#endif
#ifdef TTMod
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
            stimtime[id] += pcls1[id] + curbeat[id]*(pcls2[id] - pcls1[id])/(beats - 1);
            while (stimtime[id] < t - dt/4.0) {
               stimtime[id] += pcls1[id] + curbeat[id]*(pcls2[id] - pcls1[id])/(beats - 1);
            }
#ifdef LR1
            xrs[2*(beats*id + curbeat[id]) + 1] = Cells.xr[id];
#endif
#ifdef LR1_taud
            xrs[2*(beats*id + curbeat[id]) + 1] = Cells.xr[id];
#endif
#ifdef TT
            kis[2*(beats*id + curbeat[id]) + 1] = Cells.ki[id];
            nais[2*(beats*id + curbeat[id]) + 1] = Cells.nai[id];
            cais[2*(beats*id + curbeat[id]) + 1] = Cells.cai[id];
            casrs[2*(beats*id + curbeat[id]) + 1] = Cells.casr[id];
#endif
#ifdef UCLA
            //kis[2*(beats*id + curbeat[id]) + 1] = Cells.ki[id];
            nais[2*(beats*id + curbeat[id]) + 1] = Cells.nai[id];
            cais[2*(beats*id + curbeat[id]) + 1] = Cells.ci[id];
            cass[2*(beats*id + curbeat[id]) + 1] = Cells.cs[id];
            //casrs[2*(beats*id + curbeat[id]) + 1] = Cells.casr[id];
#endif
#ifdef LR2
            kis[2*(beats*id + curbeat[id]) + 1] = Cells.ki[id];
            nais[2*(beats*id + curbeat[id]) + 1] = Cells.nai[id];
            cais[2*(beats*id + curbeat[id]) + 1] = Cells.cai[id];
#endif
#ifdef TTMod
            kis[2*(beats*id + curbeat[id]) + 1] = Cells.ki[id];
            nais[2*(beats*id + curbeat[id]) + 1] = Cells.nai[id];
            cais[2*(beats*id + curbeat[id]) + 1] = Cells.cai[id];
            casrs[2*(beats*id + curbeat[id]) + 1] = Cells.casr[id];
#endif
            apds[beats*id + curbeat[id]] = t - startapds[id];
            if (curbeat[id] == beats - 1) {
                donebifurcation = donebifurcation + 1;
                printf("%d\n",donebifurcation);
            }
            
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDBifurcation_varypcl<typecell, ncells, beats>::iterateall(long double dt, long double t) {
    for (int i = 0; i < ncells; i++) {
        iterate(i, dt, t);
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDBifurcation_varypcl<typecell, ncells, beats>::dobif(long double dt, double startt) {
    long double t = startt;
    while (donebifurcation < ncells) {
        iterateall(dt, t);
        t = t + dt;
        if (t > stimtime[0] - dt/4.0 && t < stimtime[0] + 3*dt/4.0) {
#ifdef LR2
            printf("%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",(double)t,curbeat[0],Cells.v[0],Cells.cai[0],Cells.nai[0],Cells.ki[0],Cells.jsr[0],Cells.nsr[0],Cells.tcicr[0]);
#else
            printf("%g\t%d\t%g\n",(double)t, curbeat[0], Cells.v[0]);
#endif
        }
        //if (curbeat[0] == 225 && t > 226150 && t < 226200) {
        //    printf("%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",(double)t,curbeat[0],Cells.v[0],Cells.cai[0],Cells.nai[0],Cells.ki[0],Cells.jsr[0],Cells.nsr[0],Cells.tcicr[0]);
        //}
    }
}
#endif //APDBifurcation_varypcl_cpp
