//
//  FentonCell.cpp
//
//  Implementation of the LR1 Model, with UCLA Ito
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "FentonCell.h"

#ifndef FentonCell_cpp
#define FentonCell_cpp

#ifndef DV_MAX
#define DV_MAX 0.1
#endif

#ifndef ADAPTIVE
#define ADAPTIVE 10
#endif

#ifndef stimulus
#define stimulus -80.0
#endif

#ifndef stimduration
#define stimduration 0.5
#endif

template <int ncells>
 FentonCell<ncells>::FentonCell()
{
    for (int i = 0; i < ncells; i++) {
        v[i] = -85;
        h[i] = 1.0;
        diffcurrent[i] = 0.0;
    }
}

template <int ncells>
  bool FentonCell<ncells>::iterate(const int id, double dt, double st, double dv_max) {
    double dh;
    double ina, ik;
    
    ina = comp_ina(id, dt, dh);
    ik = comp_ik(id, double dt);
    dv = (diffcurrent[id] - (ik + ina + st))*dt;
    
    //Comment this out to remove adaptive time-stepping
    if (dv_max > 0 && dv*dv > dv_max*dv_max) {return false;}
    
    v[id] += dv;
    h[id] += dh*dt;
    
    return true;
}

template <int ncells>
void FentonCell<ncells>::stepdt (const int id, double dt, double st) {
    if (id > -1 && id < ncells) {
        bool success = iterate(id, dt, st, DV_MAX);
        if (!success) {
            for (int i = 0; i < ADAPTIVE; i++) {
                iterate(id, dt/ADAPTIVE, st, -1);
            }
        }
    }
}

template <int ncells>
  double FentonCell<ncells>::comp_ina (int id, double dt, double& dh) //Fast Sodium Current
{
    double htau, hss, ina;
    double vnew = (v[id] + 85.0)/100.0;
    
    if (vnew >= Vc) {
        hss = 0.0;
        htau = 5.6;
    }
    else {
        hss = 1.0;
        htau = 80.0;
    }
    
    if (vnew >= Vc) {
        ina = -100.0*h*(vnew - Vc)*(1.0 - vnew)/0.4576;
    }
    else {
        ina = 0;
    }

    dh = (hss-(hss-h[id])*exp(-dt/htau) - h[id])/dt;
    
    return ina;
}

template <int ncells>
double FentonCell<ncells>::comp_ik (int id, double dt)
{
    double ik;
    double vnew = (v[id] + 85.0)/100.0;
    
    if (vnew >= Vc) {
        ik = 100.0/360.0;
    }
    else {
        ik = 100.0*vnew/96.4;
    }
    
    return ik;
}


template <int ncells>
void FentonCell<ncells>::setcell (int id, FentonCell<1>* newcell)
{
    v[id] = newcell->v[0];
    h[id] = newcell->h[0];
    diffcurrent[id] = newcell->diffcurrent[0];
}

template <int ncells>
void FentonCell<ncells>::getcell (int id, FentonCell<1>* newcell)
{
    newcell->v[0] = v[id];
    newcell->h[0] = h[id];
    newcell->diffcurrent[0] = diffcurrent[id];
}

template <int ncells>
void FentonCell<ncells>::saveconditions(FILE* file, int id, bool header, double t) {
    if (header) {
        fprintf(file,"t\tv\th\n");
    }
    fprintf(file,"%g\t",t);
    fprintf(file,"%.12f\t",v[id]);
    fprintf(file,"%.12f\t",h[id]);
}
#endif // FentonCell_cpp

