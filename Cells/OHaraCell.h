//
//  UCLACell.h
//
//  Implementation of the UCLA model
//  Created by Julian Landaw on 5/24/18.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#ifndef OHaraCell_h
#define OHaraCell_h

#define nao 140.0
#define cao 1.8
#define ko 5.4
#define BSRmax 0.047
#define KmBSR 0.00087
#define BSLmax 1.124
#define KmBSL 0.0087
#define cmdnmax 0.05
#define kmcmdn 0.00238
#define trpnmax 0.07
#define kmtrpn 0.0005
#define csqnmax 10.0
#define kmcsqn 0.8

#define aCaMK 0.05
#define bCaMK 0.00068
#define CaMKo 0.05
#define KmCaM 0.0015
#define KmCaMK 0.15

#define R 8314.0
#define T 310.0
#define F 96485.0
#define frt (F/(R*T))
//added frt myself
#define L 0.01
#define rad 0.0011
#define pi 3.14159265359
#define vcell (1000*pi*rad*rad*L)
#define Ageo (2*pi*rad*rad+2*pi*rad*L)
#define Acap (2*Ageo)
#define vmyo (0.68*vcell)
#define vmito (0.26*vcell)
#define vsr (0.06*vcell)
#define vnsr (0.0552*vcell)
#define vjsr (0.0048*vcell)
#define vss (0.02*vcell)


template <int ncells>
class OHaraCell
{
public:
    /*
    //constants
    double const nao=140.0;//extracellular sodium in mM
    double const cao=1.8;//extracellular calcium in mM
    double const ko=5.4;//extracellular potassium in mM

    //buffer paramaters
    double const BSRmax=0.047;
    double const KmBSR=0.00087;
    double const BSLmax=1.124;
    double const KmBSL=0.0087;
    double const cmdnmax=0.05;
    double const kmcmdn=0.00238;
    double const trpnmax=0.07;
    double const kmtrpn=0.0005;
    double const csqnmax=10.0;
    double const kmcsqn=0.8;

    //CaMK paramaters
    double const aCaMK=0.05;
    double const bCaMK=0.00068;
    double const CaMKo=0.05;
    double const KmCaM=0.0015;
    double const KmCaMK=0.15;

    //physical constants
    double const R=8314.0;
    double const T=310.0;
    double const F=96485.0;
    double const frt = F/(R*T); // added this myself

    //cell geometry
    double const L=0.01;
    double const rad=0.0011;
    double const vcell=1000*3.14*rad*rad*L;
    double const Ageo=2*3.14*rad*rad+2*3.14*rad*L;
    double const Acap=2*Ageo;
    double const vmyo=0.68*vcell;
    double const vmito=0.26*vcell;
    double const vsr=0.06*vcell;
    double const vnsr=0.0552*vcell;
    double const vjsr=0.0048*vcell;
    double const vss=0.02*vcell;
    */
    
    // state variables
    double v[ncells];
    double nai[ncells];
    double nass[ncells];
    double ki[ncells];
    double kss[ncells];
    double cai[ncells];
    double cass[ncells];
    double cansr[ncells];
    double cajsr[ncells];
    double m[ncells];
    double hf[ncells];
    double hs[ncells];
    double j[ncells];
    double hsp[ncells];
    double jp[ncells];
    double mL[ncells];
    double hL[ncells];
    double hLp[ncells];
    double a[ncells];
    double iF[ncells];
    double iS[ncells];
    double ap[ncells];
    double iFp[ncells];
    double iSp[ncells];
    double d[ncells];
    double ff[ncells];
    double fs[ncells];
    double fcaf[ncells];
    double fcas[ncells];
    double jca[ncells];
    double nca[ncells];
    double ffp[ncells];
    double fcafp[ncells];
    double xrf[ncells];
    double xrs[ncells];
    double xs1[ncells];
    double xs2[ncells];
    double xk1[ncells];
    double Jrelnp[ncells];
    double Jrelp[ncells];
    double CaMKt[ncells];
    
    int celltype[ncells];
    
    double diffcurrent[ncells];
    
    double itofac[ncells];
    double inafac[ncells];
    double ikrfac[ncells];
    double iksfac[ncells];
    
    double iskfac[ncells];
    double skh[ncells];
    double skn[ncells];
    
    OHaraCell();
    
    bool iterate(const int id, double dt, double st, double dv_max);
    void stepdt(const int id, double dt, double st);
    
    void comp_ina(const int id, double dt, double& INa, double& INaL, double& dm, double& dhf, double& dhs, double& dj, double& dhsp, double& djp, double& dmL, double& dhL, double& dhLp);
    double comp_ito(const int id, double dt, double& da, double& diF, double& diS, double& dap, double& diFp, double& diSp);
    void comp_ical(const int id, double dt, double& ICaL, double& ICaNa, double& ICaK, double& dd, double& dff, double& dfs, double& dfcaf, double& dfcas, double& djca, double& dffp, double& dfcafp, double& dnca);
    double comp_ikr(const int id, double dt, double& dxrf, double& dxrs);
    double comp_iks(const int id, double dt, double& dxs1, double& dxs2);
    double comp_ik1(const int id, double dt, double& dxk1);
    void comp_inaca(const int id, double& INaCa_i, double& INaCa_ss);
    double comp_inak(const int id);
    double comp_ikb(const int id);
    double comp_inab(const int id);
    double comp_icab(const int id);
    double comp_ipca(const int id);
    
    double comp_rxa(const int id, double csm);
    double comp_Q(const int id);
    double comp_dir(const int id, double po, double Qr, double rxa, double dcj);
    double comp_dcp(const int id, double po, double Qr, double rxa);
    
    // SK
    double comp_isk(const int id);
    
    void setcell (int id, OHaraCell<1>* newcell);
    
    void getcell (int id, OHaraCell<1>* newcell);
    
    void saveconditions(FILE* file, int id, bool header, double t);
    
};

#endif // OHaraCell_h