//
//  Cable1D.cpp
//
//  Structure for analyzing the dynamics of variables (e.g. nai, ki)
//  as a function of PCL and APD. This structure allows the user to give
//  two APDs, one APD for set number of beats (BEATS1) and the other APD
//  for another set number of beats (BEATS2)
//
//  Created by Julian Landaw on 1/7/17.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include "Cable1D.h"

#ifndef Cable1D_cpp
#define Cable1D_cpp

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
#define stimt 100.0;
#endif

template <template<int> class typecell, int ncells>
Cable1D<typecell, ncells>::Cable1D() {
    pcl = 1000;
    numstimulus = 10;
    for (int i = 0; i < ncells; i++) {
        dx[i] = 0;
    }
}

template <template<int> class typecell, int ncells>
Cable1D<typecell, ncells>::Cable1D(double _pcl, int _numstimulus) {
    pcl = _pcl;
    numstimulus = _numstimulus;
    for (int i = 0; i < ncells; i++) {
        dx[i] = 0;
    }
}

template <template<int> class typecell, int ncells>
void Cable1D<typecell, ncells>::iterate(const int id, double commondt, double t) {
     if (id < ncells) {
#ifdef TT
         Cells.nai[id] = 12.0;
         Cells.ki[id] = 138.0;
#endif
         if (id < numstimulus && t > -commondt/4.0 && (fmod(t + commondt/4.0, pcl) < stimduration - commondt/4.0)) {
             Cells.stepdt(id, commondt, stimulus);
         }
         else {
             Cells.stepdt(id, commondt, 0);
         }
#ifdef TT
         Cells.nai[id] = 12.0;
         Cells.ki[id] = 138.0;
#endif
     }
}

template <template<int> class typecell, int ncells>
void Cable1D<typecell, ncells>::diffuse(const int id) {
    if (id == 0) {
        Cells.diffcurrent[id] = dx[id]*(Cells.v[id+1] - Cells.v[id]);
    }
    else if (id < ncells - 1) {
        Cells.diffcurrent[id] = dx[id]*(Cells.v[id+1] - Cells.v[id]) + dx[id-1]*(Cells.v[id-1] - Cells.v[id]);
    }
    else if (id == ncells - 1) {
        Cells.diffcurrent[id] = dx[id-1]*(Cells.v[id-1] - Cells.v[id]);
    }
}

template <template<int> class typecell, int ncells>
void Cable1D<typecell, ncells>::iterateall(double commondt, double t) {
    for (int i = 0; i < ncells; i++) {
        iterate(i, commondt, t);
    }
}

template <template<int> class typecell, int ncells>
void Cable1D<typecell, ncells>::diffuseall() {
    for (int i = 0; i < ncells; i++) {
    	diffuse(i);
    }
}

#endif //Cable1D_cpp
