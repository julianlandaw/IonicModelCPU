//
//  varRestitution.h
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the previous diastolic interval (DI).
//  Key modifiable variable is "basepcl", the prior pacing cycle length (PCL)
//  With default value basepcl = 1000.0;
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#ifndef varRestitution_h
#define varRestitution_h

#ifndef VARIABLE1
#define VARIABLE1 cai
#endif

#ifndef VARIABLE2
#define VARIABLE2 casr
#endif

template <template<int> class typecell, int ncells, int beats>
class varRestitution
{
public:
    typecell<ncells> AllCells;
    typecell<1> CommonCell;
    double commonst;
    double basepcl;
    int curbeat;
    double commont;
    double commondt;
    double VARIABLE1[ncells];
    double VARIABLE2[ncells];
    double S2[ncells];
    int donerestitution;
    
    varRestitution();
    
    void dorest();
    
    void iterate(const int id, double dt, double t, double startt);
    
    void iterateall(double dt, double t, double startt);
    
    void setup();
};

#endif //varRestitution_h
