//
//  DIBifurcation.h
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#ifndef TTCaiBifurcation_h
#define TTCaiBifurcation_h

template <template<int> class typecell, int ncells, int beats>
class TTCaiBifurcation
{
public:
    typecell<ncells> Cells;
    double apds[2*ncells*beats];
    double pcls[2*ncells*beats];

    double kis[2*ncells*beats];
    double nais[2*ncells*beats];
    double cais[2*ncells*beats];
    double casrs[2*ncells*beats];

    //long double d0[ncells]; 
    //long double apdfac[ncells];
    long double startapds[ncells];
    double caithreshold[ncells];
    double st[ncells];
    int curbeat[ncells];
    int donebifurcation;
    //int numapdsave;
    double stimtime[ncells];
    
    TTCaiBifurcation();
    
    void iterate(const int id, long double dt, long double t);
    
    void iterateall(long double dt, long double t);
    
    void dobif(long double dt, long double startt);
};

#endif //VaryDIBifurcation_h
