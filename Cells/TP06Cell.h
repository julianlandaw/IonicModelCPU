//
//  TP06Cell.h
//
//  Implementation of the ten Tusscher model, with UCLA Ito
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#ifndef TP06Cell_h
#define TP06Cell_h

#include <fstream>

template <int ncells>
class TP06Cell
{
public:
    double v[ncells]; //Membrane voltage
    double diffcurrent[ncells];
    
    double nai[ncells];
    double ki[ncells];
    double cai[ncells];
    double casr[ncells];
    double cass[ncells];
    
    /* Fast Sodium Current (time dependant) */
    double m[ncells]; // Na activation
    double h[ncells]; // Na inactivation
    double j[ncells]; // Na inactivation
    
    /* L-type Calcium Window Current ([Ca] dependent) */
    double d[ncells];
    double f[ncells];
    double f2[ncells];
    double fca[ncells];
    
    /* Ito Transient Outward Current (Dumaine et al. Circ Res 1999;85:803-809) */
    double zdv[ncells]; // Ito activation
    double ydv[ncells]; // Ito inactivation
    double itofac[ncells];
    double nacafac[ncells];
    double nakfac[ncells];
    double ibarcafac[ncells];
    double ikrfac[ncells];
    double iksfac[ncells];
    
    double xs[ncells];
    
    double xr1[ncells];
    double xr2[ncells];
    
    double rbar[ncells]; //19 arrays
    
    // SK Variables
    double iskfac[ncells];
    double skh[ncells];
    double skn[ncells];
    
    #ifndef UCLAito
    double rinfshift[ncells];
    double sinfshift[ncells];
    double taurshift[ncells];
    double tausshift[ncells];
    #endif
    
    // Na and K clamp variables
    bool naiclamped[ncells];
    bool kiclamped[ncells];
    
    TP06Cell();
    
    bool iterate(const int id, double dt, double st, double dv_max);
    
    void stepdt(const int id, double dt, double st);
    
    double comp_ina (int id, double dt, double& dm, double& dh, double& dj); //Fast Sodium Current
    
    double comp_ical (int id, double dt, double& dd, double& df, double& df2, double& dfca); // L-type Calcium Current
    
    double comp_ito (int id, double dt, double& dzdv, double& dydv);
    
    double comp_isk (int id);  // SK current
    
    double comp_iks (int id, double dt, double& dxs);
    
    double comp_ikr (int id, double dt, double& dxr1, double& dxr2);
    
    double comp_ik1 (int id);
    
    double comp_inaca (int id);
    
    double comp_inak (int id);
    
    double comp_ipca (int id);
    
    double comp_ipk (int id);
    
    double comp_ibna (int id);
    
    double comp_ibca (int id);
    
    void comp_calcdyn (int id, double dt, double ical, double ibca, double ipca, double inaca, double& drbar, double& dcasr, double& dcai, double& dcass); // MAKE SURE THIS IS COMPUTED AFTER ALL OTHER CURRENTS
    
    void setcell (int id, TP06Cell<1>* newcell);
    
    void getcell (int id, TP06Cell<1>* newcell);
    
    void saveconditions(FILE* file, int id, bool header, double t=0.0);
    
    void write(std::fstream& file);
    
    void read(std::fstream& file);
    
};

#endif // TP06Cell_h

