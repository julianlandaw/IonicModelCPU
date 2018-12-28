//
//  DIBifurcation.h
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#ifndef VaryDIBifurcation_h
#define VaryDIBifurcation_h

template <template<int> class typecell, int ncells, int beats>
class VaryDIBifurcation
{
public:
    typecell<ncells> Cells;
    double apds[2*ncells*beats];
#ifdef LR1
    double xrs[2*ncells*beats];
#endif
#ifdef LR1_taud
    double xrs[2*ncells*beats];
#endif

#ifdef TT
    double kis[2*ncells*beats];
    double nais[2*ncells*beats];
    double cais[2*ncells*beats];
    double casrs[2*ncells*beats];
#endif
#ifdef UCLA
    //double kis[2*ncells*beats];
    double nais[2*ncells*beats];
    double cais[2*ncells*beats];
    double cass[2*ncells*beats];
    //double casrs[2*ncells*beats];
#endif

#ifdef LR2
    double kis[2*ncells*beats];
    double nais[2*ncells*beats];
    double cais[2*ncells*beats];
#endif
#ifdef TTMod
    double kis[2*ncells*beats];
    double nais[2*ncells*beats];
    double cais[2*ncells*beats];
    double casrs[2*ncells*beats];
#endif

    double d0[ncells]; 
    double apdfac[ncells];
    long double startapds[ncells];
    double st[ncells];
    int curbeat[ncells];
    int donebifurcation;
    int numapdsave;
    double stimtime[ncells];
    double apd_star[ncells];
    
    VaryDIBifurcation();
    
    void iterate(const int id, long double dt, long double t);
    
    void iterateall(long double dt, long double t);
#ifdef saveap 
    void dobif(long double dt, long double startt, FILE** ap);
#else
    void dobif(long double dt, long double startt);
#endif
};

#endif //VaryDIBifurcation_h
