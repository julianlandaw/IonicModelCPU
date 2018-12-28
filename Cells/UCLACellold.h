//
//  UCLACell.h
//
//  Implementation of the UCLA model
//  Created by Julian Landaw on 5/24/18.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#ifndef UCLACell_h
#define UCLACell_h

template <int ncells>
class UCLACell
{
public:
    double xm[ncells];// sodium m-gate
    double xh[ncells];// sodium h-gate
    double xj[ncells];// soiumj-gate

    double xr[ncells];// ikr gate variable 
    double xs1[ncells];// iks gate variable
    double xs2[ncells];// iks gate varaible 

    double xtos[ncells];// ito slow activation
    double ytos[ncells];// ito slow inactivation
    
    double v[ncells]; // voltage

    double cp[ncells];// averaged dyadic space con.
    double cs[ncells];// averaged submembrane conc.
        
    double ci[ncells];// myoplasm conc.

    double cj[ncells];// NSR load
    double cjp[ncells];// average JSR load

    double xir[ncells];// SR current flux

        // Markov gate variables 

    double c1[ncells];// C1
    double c2[ncells];// C2
    double xi1ca[ncells];// I1_Ca
    double xi1ba[ncells];// I1_Ba
    double xi2ca[ncells];// I2_Ca
    double xi2ba[ncells];// I2_Ba

    double nai[ncells];// internal Na conc.

    double xtof[ncells];// ito fast activation
    double ytof[ncells];// ito slow inactivation

    double tropi[ncells];// time dependent buffers in myplasm (troponin)
    double trops[ncells];// time dependent buffers in submembrane (troponin)

    double vold[ncells];
    double jparam[ncells];
    
    double diffcurrent[ncells];
    
    double itofac[ncells];
    
    UCLACell();
    
    bool iterate(const int id, double dt, double st, double dv_max);
    void stepdt(const int id, double dt, double st);
    double comp_ina(const int id, double dt, double& dxm, double& dxh, double& dxj);
    double comp_ikr(const int id, double dt, double& dxr);
    double comp_iks(const int id, double dt, double& dxs1, double& dxs2);
    double comp_ik1(const int id);
    double comp_ito(const int id, double dt, double& dxtos, double& dytos, double& dxtof, double& dytof);
    double comp_inak(const int id);
    double comp_inaca(const int id, double csm);
    double comp_icalpo(const int id, double& dc1, double& dc2, double& dxi1ca, double& dxi2ca, double& dxi1ba, double& dxi2ba);
    double comp_iup(const int id);
    double comp_ileak(const int id);
    double comp_inst_buffer(double c);
    
    double comp_rxa(const int id, double csm);
    double comp_Q(const int id);
    double comp_dir(const int id, double po, double Qr, double rxa, double dcj);
    double comp_dcp(const int id, double po, double Qr, double rxa);
    
    void setcell (int id, UCLACell<1>* newcell);
    
    void getcell (int id, UCLACell<1>* newcell);
    
    void saveconditions(FILE* file, int id, bool header, double t=0.0);
};

#endif // UCLACell_h