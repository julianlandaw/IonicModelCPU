//
//  APDBifurcation.cuh
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#ifndef ClampedAPDBif_cuh
#define ClampedAPDBif_cuh

#ifndef VARIABLE1
#define VARIABLE1 cai
#endif

#ifndef VARIABLE2
#define VARIABLE2 casr
#endif

#ifndef VARIABLE3
#define VARIABLE3 nai
#endif

#ifndef VARIABLE4
#define VARIABLE4 ki
#endif

template <template<int> class typecell, int ncells, int beats>
class ClampedAPDBif
{
public:
    typecell<ncells> Cells;
    double apds[ncells*beats];
    double var1[2*ncells*beats];
    double var2[2*ncells*beats];
    double pcls[ncells];
    double clampvar3[ncells];
    double clampvar4[ncells];
    double startapds[ncells];
    double st[ncells];
    int curbeat[ncells];
    bool notdonebifurcation[ncells];
    
    ClampedAPDBif();
    
    void iterate(const int id, long double dt, long double t);
    
    void iterateall(long double dt, long double t);
};

#endif //ClampedAPDBif_cuh
