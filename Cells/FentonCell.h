//
//  FentonCell.h
//
//  Implementation of the LR1 model, with UCLA Ito
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#ifndef FentonCell_h
#define FentonCell_h

template <int ncells>
class FentonCell
{
public:
    double v[ncells]; //Membrane voltage
    double diffcurrent[ncells];
    double h[ncells]; // Na inactivation
    
    FentonCell();
    
    bool iterate(const int id, double dt, double st, double dv_max);
    
    void stepdt(const int id, double dt, double st);
    
    double comp_ina (int id, double dt, double& dh); //Fast Sodium Current
    
    double comp_ik (int id, double dt);
    
    void setcell (int id, FentonCell<1>* newcell);
    
    void getcell (int id, FentonCell<1>* newcell);
    
    void saveconditions(FILE* file, int id, bool header, double t=0.0);
};

#endif // FentonCell_h

