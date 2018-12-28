//
//  DIBifurcation.cpp
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "VaryDIBifurcation.h"

#ifndef VaryDIBifurcation_cpp
#define VaryDIBifurcation_cpp

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
VaryDIBifurcation<typecell, ncells, beats>::VaryDIBifurcation() {
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
    for (i = 0; i < ncells; i++) {
        d0[i] = 0.0; // DI = d0 + alpha*APD
        apdfac[i] = 0.0;
        startapds[i] = 0.0;
        st[i] = 0.0;
        curbeat[i] = -1;
        stimtime[i] = 0.0;
	apd_star[i] = 0.0;
    }
    donebifurcation = 0;
    numapdsave = 0;
}

template <template<int> class typecell, int ncells, int beats>
void VaryDIBifurcation<typecell, ncells, beats>::iterate(const int id, long double dt, long double t) {
    if (id < ncells && curbeat[id] < beats) {
        double vold = Cells.v[id];
        st[id] = (t > -dt/50.0 && t > stimtime[id] - dt/50.0 && t < stimtime[id] + stimduration - dt/50.0) ? stimulus : 0.0;
#ifdef TTMod
        Cells.ki[id] = 140.0;
        Cells.nai[id] = 10.0;
#endif
        Cells.stepdt(id, (double)dt, st[id]);

#ifdef TTMod
        Cells.ki[id] = 140.0;
        Cells.nai[id] = 10.0;
#endif
        double vnew = Cells.v[id];
        if (vold < threshold && vnew >= threshold) {
            startapds[id] = t;
            if (id == 0) {
                //printf("%g\t%d\t%g\n",(double)t,curbeat[id],(curbeat[id] > -1) ? apds[curbeat[0]] : -1);
            }
            
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
            double apd = t - startapds[id];
	    
            double updatestim;
            if (numapdsave == 0) {
		if (curbeat[id] < (BEATS + REMOVEBEATS)/2) {
                    updatestim = d0[id] + apdfac[id]*apd; //300;
		}	
		else {
                    updatestim = d0[id] + apdfac[id]*apd; //1000; 
		}
/*
#ifdef LR1
		    
		//if (curbeat[id] < 100) {
		if (curbeat[id] < REMOVEBEATS) {
		    apd_star[id] = apd;
		    updatestim = d0[id] + 2.4*apd;
		}			
		else {
                    if (curbeat[id] >= REMOVEBEATS & curbeat[id] < (BEATS + 2*REMOVEBEATS)/3) {
    		        updatestim = d0[id] + 2.4*apd_star[id] + apdfac[id]*(apd - apd_star[id]);
		    }
		    else if (curbeat[id] >= (BEATS + 2*REMOVEBEATS)/3 & curbeat[id] < (2*BEATS + REMOVEBEATS)/3) {
                        updatestim = d0[id] + 2.4*apd_star[id]; //+ 0*(apd - apd_star[id]);
		    }
		    else {
                        updatestim = d0[id] + 2.4*apd_star[id] + apdfac[id]*(apd - apd_star[id]);
		    }
		}
#else
                updatestim = d0[id] + apdfac[id]*apd; // apdfac = 0 implies constant-T pacing, =1 --> constant-DI
#endif
*/
            }
	    else if (curbeat[id] < 1) {
	    //else if (curbeat[id] < REMOVEBEATS) {
#ifdef LR1
		apd_star[id] = apd;	
		//updatestim = d0[id] + 2.4*apd;
		updatestim = d0[id];
#else
		updatestim = d0[id];
#endif
	    }
            else {
                double apdsums = apd;
                for (int i = 0; i < numapdsave; i++) {
                    if (curbeat[id]-i-1 < 1) {
                        apdsums += apd;
                    }
                    else {    
                        apdsums += apds[beats*id + curbeat[id] - i - 1];
                    }
                }
/*		
#ifdef LR1
		if (curbeat[id] >= REMOVEBEATS & curbeat[id] < (BEATS + 2*REMOVEBEATS)/3) {
		    updatestim = d0[id] + 2.4*apd_star[id] + apdfac[id]*(apd - apdsums/(1 + numapdsave));
		}
		else if (curbeat[id] >= (BEATS + 2*REMOVEBEATS)/3 & curbeat[id] < (2*BEATS + REMOVEBEATS)/3) {
                    updatestim = d0[id] + 2.4*apd_star[id];
		}
		else {
                    updatestim = d0[id] + 2.4*apd_star[id] + apdfac[id]*(apd - apdsums/(1 + numapdsave));
		}
#else
                updatestim = d0[id] + apdfac[id]*(apd - apdsums/(1 + numapdsave));
#endif
*/
		updatestim = d0[id] + apdfac[id]*(apd - apdsums/(1 + numapdsave));
            }
            while (updatestim < apd) {
            	updatestim = updatestim + d0[id];
	    } 
            stimtime[id] += updatestim;
	    //stimtime[id] += (updatestim > apd) ? updatestim : apd + 10;
            printf("%g\t%d\t%g\t%g\t%g\n",(double)t,curbeat[id],stimtime[id],d0[id],apdfac[id]);
            
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
            apds[beats*id + curbeat[id]] = apd;
            
            
        }
        if (curbeat[id] == beats - 1 && t > stimtime[id] - stimt - dt - dt/50.0) {
            donebifurcation = donebifurcation + 1;
            printf("%g\t%d\n",(double)t,donebifurcation);
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void VaryDIBifurcation<typecell, ncells, beats>::iterateall(long double dt, long double t) {
    for (int i = 0; i < ncells; i++) {
        iterate(i, dt, t);
    }
}

template <template<int> class typecell, int ncells, int beats>
#ifdef saveap
void VaryDIBifurcation<typecell, ncells, beats>::dobif(long double dt, long double startt, FILE** ap) {
#else
void VaryDIBifurcation<typecell, ncells, beats>::dobif(long double dt, long double startt) {
#endif
    long double t = startt;
    long double t_save = 0;
    while (donebifurcation < ncells) {
        iterateall(dt, t);
        t = t + dt;
#ifdef saveap
	if (t > t_save - dt/4.0) {
            if (curbeat[0] > REMOVEBEATS) {
#ifdef LR1
	    	fprintf(*ap,"%g\t%g\t%g\n",(double)t,Cells.v[0],Cells.xr[0]);
#else
		fprintf(*ap,"%g\t%g\n",(double)t,Cells.v[0]);
#endif
	    }
	    t_save = t_save + 1;
	}
#endif
    }
}
#endif //VaryDIBifurcation_cpp
