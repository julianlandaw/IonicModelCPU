//
//  APDRestitution.h
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the previous diastolic interval (DI).
//  Key modifiable variable is "basepcl", the prior pacing cycle length (PCL)
//  With default value basepcl = 1000.0;
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#ifndef APDRestitution_h
#define APDRestitution_h

template <template<int> class typecell, int ncells, int beats>
class APDRestitution
{
public:
    typecell<ncells> AllCells;
    typecell<1> CommonCell;
    double st;
    double basepcl;
    int curbeat;
    long double commont;
    long double commondt;
    double DI[ncells];
    double S2[ncells];
    int donerestitution;
    double stimtime;
    
    APDRestitution();
    
    void dorest(FILE* ap, int id1, int id2);
    
    void iterate(const int id, long double dt, long double t);
    
    void iterateall(long double dt, long double t);
    
    void setup(FILE* ap);
};

#endif //APDRestitution_h
