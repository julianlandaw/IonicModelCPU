//
//  VoltageClamp.h
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#ifndef VoltageClamp_h
#define VoltageClamp_h

#ifndef VARIABLE
#def VARIABLE cai
#endif

template <template<int> class typecell, int ncells, int beats>
class VoltageClamp
{
public:
    typecell<ncells> Cells;
    long double t;
    long double commondt;
    int curbeat[ncells];
    double pcls[ncells];
    double apds1[ncells];
    double apds2[ncells];
    double voltage1[ncells];
    double voltage2[ncells];
    bool inapd[ncells];
    int donetrain;
#ifdef TT
    double nai[4*ncells*beats];
    double cai[4*ncells*beats];
#endif
#ifdef UCLA
    double ci[4*ncells*beats];
    double nai[4*ncells*beats];
#endif
    
    VoltageClamp();
    
    void iterate(const int id);
    
    void iterateall(FILE *ap);
    
    void dospiketrain(FILE *ap);
};

#endif //VoltageClamp_h
