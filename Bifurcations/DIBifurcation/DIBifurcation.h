//
//  DIBifurcation.h
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#ifndef DIBifurcation_h
#define DIBifurcation_h

template <template<int> class typecell, int ncells, int beats>
class DIBifurcation
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
#ifdef TTMod
    double kis[2*ncells*beats];
    double nais[2*ncells*beats];
    double cais[2*ncells*beats];
    double casrs[2*ncells*beats];
#endif
    long double dis[ncells];
    long double startapds[ncells];
    double st[ncells];
    int curbeat[ncells];
    long double t_endap[ncells];
    int donebifurcation;
    
    DIBifurcation();
    
    void iterate(const int id, long double dt, long double t);
    
    void iterateall(long double dt, long double t);
    
    void dobif(long double dt, long double startt);
};

#endif //DIBifurcation_h
