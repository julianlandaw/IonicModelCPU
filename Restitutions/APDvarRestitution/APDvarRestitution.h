//
//  APDvarRestitution.h
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the previous diastolic interval (DI).
//  Key modifiable variable is "basepcl", the prior pacing cycle length (PCL)
//  With default value basepcl = 1000.0;
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#ifndef APDvarRestitution_h
#define APDvarRestitution_h

#ifndef VARIABLE
#define VARIABLE cai
#endif

template <template<int> class typecell, int ncells, int beats>
class APDvarRestitution
{
public:
    typecell<ncells> AllCells;
    typecell<1> CommonCell;
    double st;
    double basepcl;
    int curbeat;
    double commont;
    double commondt;
    double VARIABLE[ncells];
    double DI[ncells];
    double S2[ncells];
    bool notstartrestitution[ncells];
    int donerestitution;
    
    APDvarRestitution();
    
    void dorest();
    
    void iterate(const int id, double dt, double t);
    
    void iterateall(double dt, double t);
    
    void setup();
};

#endif //APDvarRestitution_h
