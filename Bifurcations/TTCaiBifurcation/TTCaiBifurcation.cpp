//
//  TTCaiBifurcation.cpp
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "TTCaiBifurcation.h"

#ifndef TTCaiBifurcation_cpp
#define TTCaiBifurcation_cpp

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
TTCaiBifurcation<typecell, ncells, beats>::TTCaiBifurcation() {
    int i;
    for (i = 0; i < ncells*beats; i++) {
        apds[i] = 0.1;
        pcls[i] = 0.1;
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
    for (i = 0; i < ncells; i++) {
        //d0[i] = 0.0; // DI = d0 + alpha*APD
        //apdfac[i] = 0.0;
        caithreshold[i] = 4.2e-5;
        startapds[i] = 0.0;
        st[i] = 0.0;
        curbeat[i] = -1;
        stimtime[i] = 0.0;
    }
    donebifurcation = 0;
    //numapdsave = 0;
}

template <template<int> class typecell, int ncells, int beats>
void TTCaiBifurcation<typecell, ncells, beats>::iterate(const int id, long double dt, long double t) {
    if (id < ncells && curbeat[id] < beats) {
        double vold = Cells.v[id];
        double caiold = Cells.cai[id];
        st[id] = (t > -dt/50.0 && t > stimtime[id] - dt/50.0 && t < stimtime[id] + stimduration - dt/50.0) ? stimulus : 0.0;

        Cells.stepdt(id, (double)dt, st[id]);

        double vnew = Cells.v[id];
        double cainew = Cells.cai[id];
        if (caiold >= caithreshold[id] && cainew < caithreshold[id] && t > -dt/50.0 && vnew < threshold && t > stimtime[id] + 5) {
           pcls[beats*id + curbeat[id]] = t + dt - stimtime[id];
           stimtime[id] = t + dt;
        }
        else if (cainew < caithreshold[id] && vnew < threshold && t > stimtime[id] + 5) {
           pcls[beats*id + curbeat[id]] = t + dt - stimtime[id];
           printf("ALREADY BELOW\n");
           stimtime[id] = t + dt;
        }
        if (vold < threshold && vnew >= threshold) {
            startapds[id] = t;
            if (id == 0) {
                printf("%g\t%d\n",(double)t,ncells);
            }

            if (curbeat[id] >= -1 && curbeat[id] < beats - 1) {
                kis[2*(beats*id + curbeat[id] + 1)] = Cells.ki[id];
                nais[2*(beats*id + curbeat[id] + 1)] = Cells.nai[id];
                cais[2*(beats*id + curbeat[id] + 1)] = Cells.cai[id];
                casrs[2*(beats*id + curbeat[id] + 1)] = Cells.casr[id];
            }
        }
        else if (vold >= threshold && vnew < threshold) {
            curbeat[id] += 1;
            double apd = t - startapds[id];

            kis[2*(beats*id + curbeat[id]) + 1] = Cells.ki[id];
            nais[2*(beats*id + curbeat[id]) + 1] = Cells.nai[id];
            cais[2*(beats*id + curbeat[id]) + 1] = Cells.cai[id];
            casrs[2*(beats*id + curbeat[id]) + 1] = Cells.casr[id];

            apds[beats*id + curbeat[id]] = apd;
            if (curbeat[id] == beats - 1) {
                donebifurcation = donebifurcation + 1;
                printf("%d\n",donebifurcation);
            }
            
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void TTCaiBifurcation<typecell, ncells, beats>::iterateall(long double dt, long double t) {
    for (int i = 0; i < ncells; i++) {
        iterate(i, dt, t);
    }
}

template <template<int> class typecell, int ncells, int beats>
void TTCaiBifurcation<typecell, ncells, beats>::dobif(long double dt, long double startt) {
    long double t = startt;
    long double t_save = 0;
    while (donebifurcation < ncells) {
        iterateall(dt, t);
        t = t + dt;
    }
}
#endif //TTCaiBifurcation_cpp
