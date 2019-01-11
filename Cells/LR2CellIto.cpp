//
//  LR2CellIto.cpp
//
//  Implementation of the LRd Model
//  Created by Julian Landaw on 2/26/18.
//  Copyright © 2018 Julian Landaw. All rights reserved.
//
#include <math.h>
#include <fstream>
#include "LR2CellIto.h"

#ifndef LR2CellIto_cpp
#define LR2CellIto_cpp

#ifndef ADAPTIVE
#define ADAPTIVE 10
#endif

#ifndef DV_MAX
#define DV_MAX 0.1
#endif

#ifndef stimulus
#define stimulus -80.0
#endif

#ifndef stimduration
#define stimduration 0.5
#endif

template <int ncells>
LR2CellIto<ncells>::LR2CellIto()
{
    #define l 0.01  
    #define a 0.0011 
    #define pi 3.141592
    #define vcell (1000.0*pi*a*a*l) 
    #define ageo (2.0*pi*a*a+2.0*pi*a*l)
    #define acap (ageo*2.0)
    #define vmyo (vcell*0.68) 
    #define vmito (vcell*0.26) 
    #define vsr (vcell*0.06) 
    #define vnsr (vcell*0.0552)
    #define vjsr (vcell*0.0048)
    #define vcleft (vcell*0.12/0.88)
    
    /* Terms for Solution of Conductance and Reversal Potential */ 
    #define R 8314.0
    #define frdy 96485.0 
    #define temp 310.0
    
    /* Ion Valences */ 
    #define zna 1.0
    #define zk 1.0 
    #define zca 2.0
    
    /* Ion Concentrations */ 
    #define taudiff 1000.0 
    
    /* NSR Ca Ion Concentration Changes */ 
    #define kmup 0.00092
    #define iupbar 0.00875 
    #define nsrbar 15.0
    
    /* JSR Ca Ion Concentration Changes */ 
    #define tauon 0.5
    #define tauoff 0.5
    #define csqnth 8.75
    #define gmaxrel 150.0 
    #define csqnbar 10.0
    #define kmcsqn 0.8
    
    
    /* Translocation of Ca Ions from NSR to JSR */ 
    #define tautr 180.0

    /* Myoplasmic Ca Ion Concentration Changes */ 
    #define cmdnbar 0.050 
    #define trpnbar 0.070
    #define kmcmdn 0.00238
    #define kmtrpn 0.0005
    
    /* Fast Sodium Current (time dependant) */ 
    #define gna 16.0
    
    /* Current through L-type Ca Channel */ 
    #define kmca 0.0006 
    #define pca 0.00054
    #define gacai 1.0
    #define gacao 0.341
    #define pna 0.000000675 
    #define ganai 0.75 
    #define ganao 0.75 
    #define pk 0.000000193
    #define gaki 0.75
    #define gako 0.75
    
    /* Current through T-type Ca Channel */ 
    #define gcat 0.05
    
    /* Rapidly Activating Potassium Current */ 
    // Constant gkr embedded in the formula as 0.02614*sqrt(ko[id]/5.4)
    
    /* Slowly Activating Potassium Current */ 
    // Constant gks embedded in the formula as 0.433*(1+0.6/(1+pow((0.000038/cai[id]),1.4))); 
    #define prnak 0.01833
    
    /* Potassium Current (time-independent) */ 
    // Constant gki embedded in the formula as gki = 0.75*(sqrt(ko[id]/5.4)); 
    
    /* Plateau Potassium Current */ 
    #define gkp 0.00552
    
    /* Na-Activated K Channel */ 
    #define gkna 0.12848
    #define nkna 2.8
    #define kdkna 66.0
    
    /* ATP-Sensitive K Channel */ 
    #define natp 0.24
    #define nicholsarea 0.00005
    #define atpi 3.0
    #define hatp 2.0
    #define katp 0.250
    
    /* Ito Transient Outward Current (Dumaine et al. Circ Res 1999;85:803-809) */  
    #define gitodv 0.5
    
    /* Isk calcium-dependent potassium current */
    #define gsk 0.005
    //#define skh 0.0006
      //'normal' is 0.0006
    //#define skn 4 
      

    /* Sodium-Calcium Exchanger V-S */ 
    #define c1 0.00025 
    #define c2 0.0001
    #define gammas 0.15 
    
    /* Sodium-Potassium Pump */ 
    #define ibarnak 2.25 
    #define kmnai 10.0
    #define kmko 1.5
    
    /* Nonspecific Ca-activated Current */ 
    #define pnsca 0.000000175 
    #define kmnsca 0.0012
    
    /* Sarcolemmal Ca Pump */ 
    #define ibarpca 1.15 
    #define kmpca 0.0005
    
    /* Ca Background Current */ 
    #define gcab 0.003016
    
    /* Na Background Current */ 
    #define gnab 0.004
    
    /* Cleft Space */
    #define nabm 140.0
    #define kbm 4.5
    #define cabm 1.8
    
    
    /*
    nai = 12.236437; // Initial Intracellular Na (mM) 
    nao = 140; // Initial Extracellular Na (mM) 
    nabm = 140; // Initial Bulk Medium Na (mM) 
    ki = 136.89149; // Initial Intracellular K (mM) 
    ko = 4.5; // Initial Extracellular K (mM) 
    kbm = 4.5; // Initial Bulk Medium K (mM) 
    cai = 0.000079; // Initial Intracellular Ca (mM) 
    cao = 1.8; // Initial Extracellular Ca (mM) 
    cabm = 1.8; // Initial Bulk Medium Ca (mM) 

    // Initial Gate Conditions //
    m = 0.000838; 
    h = 0.993336; 
    j = 0.995484; 
    d = 0.000003; 
    f = 0.999745; 
    xs1 = 0.004503; 
    xs2 = 0.004503; 
    xr = 0.000129; 
    b = 0.000994; 
    g = 0.994041; 
    zdv = 0.0120892; 
    ydv = 0.999978;

    // Initial Conditions //
    grelbarjsrol = 0; 
    tjsrol = 1000; 
    tcicr = 1000; 
    jsr = 1.179991; 
    nsr = 1.179991; 
    trpn = 0.0143923; 
    cmdn = 0.00257849; 
    csqn = 6.97978; 
    flag = 0; 
    dt = udt; 
    utsc = 50; 
    dcaiont = 0; 
    i=-1;
    */
    for (int i = 0; i < ncells; i++) {
        v[i] = -88.654973; 
        nai[i] = 12.236437; // Initial Intracellular Na (mM) 
        nao[i] = 140; // Initial Extracellular Na (mM) 
        //nabm[i] = 140; // Initial Bulk Medium Na (mM) 
        ki[i] = 136.89149; // Initial Intracellular K (mM) 
        ko[i] = 4.5; // Initial Extracellular K (mM) 
        //kbm[i] = 4.5; // Initial Bulk Medium K (mM) 
        cai[i] = 0.000079; // Initial Intracellular Ca (mM) 
        cao[i] = 1.8; // Initial Extracellular Ca (mM) 
        //cabm[i] = 1.8; // Initial Bulk Medium Ca (mM) 
        
        m[i] = 0.000838; 
        h[i] = 0.993336; 
        j[i] = 0.995484; 
        d[i] = 0.000003; 
        f[i] = 0.999745; 
        xs1[i] = 0.004503; 
        xs2[i] = 0.004503; 
        xr[i] = 0.000129; 
        b[i] = 0.000994; 
        g[i] = 0.994041; 
        zdv[i] = 0.0120892; 
        ydv[i] = 0.999978;
        
        //grelbarjsrol[i] = 0; 
        //tjsrol[i] = 1000; 
        tcicr[i] = 1000.0; 
        jsr[i] = 1.179991; 
        nsr[i] = 1.179991; 
        flag[i] = 0; 
        dcaiontold[i] = 0.0;
        caiontold[i] = 0.0;
        
        icalfac[i] = 1.0;
        itofac[i] = 1.0;
        iskfac[i] = 0.0;
        skh[i] = 0.0006;
        skn[i] = 4;
        tauXfac[i] = 1.0;
        iupfac[i] = 1.0;
        ikrfac[i] = 1.0;
        iksfac[i] = 1.0;

		nacafac[i] = 1.0;
    }
    
}

template <int ncells>
bool LR2CellIto<ncells>::iterate(const int id, double dt, double st, double dv_max)
{    
    double ina, inab, ilcana, inak, inaca; //sodium currents
    double ikr, iks, iki, ikp, ilcak, ito; //potassium currents
    double ilca, icab, ipca, icat; //calcium currents
    double naiont, kiont, caiont, it, dv; //total currents
    double dm, dh, dj; //sodium currents gate changes
    double dd, df; //ical gate changes
    double db, dg; //icat gate changes
    double dxr; //ikr gate changes
    double dxs1, dxs2; //iks gate changes
    double dzdv, dydv; //ito gate changes
    double dnai, dki, dcai; //ion concentration changes
    double dcaiont, dtcicr, djsr, dnsr; //calcium dynamics changes
    //double dnao, dko, dcao; //cleft concentration changes
    double isk;
    
    comp_ina(id, dt, dm, dh, dj, ina);
    comp_ical(id, dt, dd, df, ilca, ilcana, ilcak);
    comp_icat(id, dt, db, dg, icat);
    comp_ikr(id, dt, dxr, ikr);
    comp_iks(id, dt, dxs1, dxs2, iks);
    comp_iki(id, iki);
    comp_ikp(id, ikp);
    //comp_ikna(id, ikna); 
    //comp_ikatp(id, ikatp);
    comp_ito(id, dt, dzdv, dydv, ito);
    comp_isk(id, isk);
    comp_inaca(id, inaca);
    comp_inak(id, inak);
    //comp_insca(id, insna, insk); 
    comp_ipca(id, ipca);
    comp_icab(id, icab);
    comp_inab(id, inab);
    
    naiont = ina+inab+ilcana+3*inak+3*inaca; 
    kiont = ikr+iks+iki+ikp+ilcak-2*inak+ito + isk;
    caiont = ilca+icab+ipca-2*inaca+icat; 
    
    it = st+naiont+kiont+caiont - diffcurrent[id];
    dv = -it*dt; //since current is oriented outward from the cell, by convention
    
    if (dv_max > 0 && dv*dv > dv_max*dv_max) {return false;}
    
    conc_nai (id, naiont, dnai);
    conc_ki (id, kiont, st - diffcurrent[id], dki);
    calc_dyn(id, dt, caiont, dcaiont, dtcicr, djsr, dnsr, dcai);
    //conc_cleft(id, dt, st - diffcurrent[id], naiont, kiont, caiont, dnao, dko, dcao);
    
    if(v[id] > -35.0 && dcaiont>dcaiontold[id] && flag[id] == 0) 
    {
        flag[id] = 1; 
        tcicr[id] = 0;
    }
    
    else if (it < -10.0 && tcicr[id] > 10.0 && flag[id] == 1) 
    {
        flag[id] = 0;
        tcicr[id] += dtcicr*dt;
    }
    else {
        tcicr[id] += dtcicr*dt;
    }
    dcaiontold[id] = dcaiont;
    caiontold[id] = caiont;
    
    v[id] += dv;
    m[id] += dm*dt;
    h[id] += dh*dt;
    j[id] += dj*dt;
    d[id] += dd*dt;
    f[id] += df*dt;
    b[id] += db*dt;
    g[id] += dg*dt;
    xr[id] += dxr*dt;
    xs1[id] += dxs1*dt;
    xs2[id] += dxs2*dt;
    zdv[id] += dzdv*dt;
    ydv[id] += dydv*dt;
    nai[id] += dnai*dt;
    ki[id] += dki*dt;
    cai[id] += dcai*dt;
    jsr[id] += djsr*dt;
    nsr[id] += dnsr*dt;
    //nao[id] += dnao*dt;
    //ko[id] += dko*dt;
    //cao[id] += dcao*dt;
    
    return true;
}

template <int ncells>
void LR2CellIto<ncells>::stepdt (int id, double dt, double st) {
    if (id > -1 && id < ncells) {
        bool success = iterate(id, dt, st, DV_MAX);
        if (!success) {
            for (int i = 0; i < ADAPTIVE; i++) {
                iterate(id, dt/ADAPTIVE, st, -1);
            }
        }
    }
}

template <int ncells>
void LR2CellIto<ncells>::comp_ina (int id, double dt, double& dm, double& dh, double& dj, double &ina) 
{ 
    double ena, am, bm, ah, bh, aj, bj, mtau, htau, jtau, mss, hss, jss; 
    ena = ((R*temp)/frdy)*log(nao[id]/nai[id]); 

    if (fabs(v[id] + 47.13) < 1e-2) {
        am = 0.32/0.1;    
    }
    else {
        am = 0.32*(v[id]+47.13)/(1.0-exp(-0.1*(v[id]+47.13))); 
    }
    bm = 0.08*exp(-v[id]/11.0); 
    if (v[id] < -40.0) 
    {
        ah = 0.135*exp((80.0+v[id])/-6.8); 
        bh = 3.56*exp(0.079*v[id])+310000.0*exp(0.35*v[id]); 
        aj = (-127140.0*exp(0.2444*v[id])-0.00003474*exp(-0.04391*v[id]))*((v[id]+37.78)/(1.0+exp(0.311*(v[id]+79.23)))); 
        bj = (0.1212*exp(-0.01052*v[id]))/(1.0+exp(-0.1378*(v[id]+40.14)));
    }
    else 
    {
        ah = 0.0; 
        bh = 1.0/(0.13*(1.0+exp((v[id]+10.66)/-11.1))); 
        aj = 0.0; 
        bj = (0.3*exp(-0.0000002535*v[id]))/(1.0+exp(-0.1*(v[id]+32.0)));
    }
    mtau = 1.0/(am+bm); 
    htau = 1.0/(ah+bh); 
    jtau = 1.0/(aj+bj); 

    mss = am*mtau; 
    hss = ah*htau; 
    jss = aj*jtau; 

    dm = (mss-(mss-m[id])*exp(-dt/mtau) - m[id])/dt; 
    dh = (hss-(hss-h[id])*exp(-dt/htau) - h[id])/dt; 
    dj = (jss-(jss-j[id])*exp(-dt/jtau) - j[id])/dt; 

    ina = gna*m[id]*m[id]*m[id]*h[id]*j[id]*(v[id]-ena); 
}

template <int ncells>
void LR2CellIto<ncells>::comp_ical (int id, double dt, double& dd, double& df, double& ilca, double& ilcana, double& ilcak) 
{ 
    double dss, taud, fss, tauf, fca, ibarca, ibarna, ibark, vfrt;
#ifndef Miyoshi
    dss = 1.0/(1.0+exp(-(v[id]+10.0)/6.24)); 
    if (fabs(v[id] + 10.0) < 1e-2) {
        taud = dss*1.0/(6.24*0.035);    
    }
    else {
        taud = dss*(1.0-exp(-(v[id]+10.0)/6.24))/(0.035*(v[id]+10.0)); 
    }
    //taud = dss*(1.0-exp(-(v[id]+10.0)/6.24))/(0.035*(v[id]+10.0)); 

    fss = (1.0/(1.0+exp((v[id]+32.0)/8.0)))+(0.6/(1.0+exp((50.0-v[id])/20.0))); 
    tauf = 1.0/(0.0197*exp(-pow(0.0337*(v[id]+10.0),2.0))+0.02); 
#else
    dss = 1.0/(1.0 + exp(-(v[id]+7.0)/6.0));
    taud = 0.29 + 1.1*exp(0.052*(v[id]+13.0))/(1.0 + exp(0.132*(v[id]+13.0)));
    
    fss = 0.99/(1.0 + exp((v[id]+28.9)/4.9)) + 0.23/(1.0 + exp(-(v[id]-40.0)/32.0));
    tauf = 22.0 + 280.0*exp(0.062*(v[id]+28.3))/(1.0 + exp(0.25*(v[id]+28.3)));
#endif

    dd = (dss-(dss-d[id])*exp(-dt/taud) - d[id])/dt; 
    df = (fss-(fss-f[id])*exp(-dt/tauf) - f[id])/dt; 
    
    vfrt = v[id]*frdy/(R*temp);
    
    if (fabs(vfrt) < 1e-3) {
        ibarca = pca*zca*frdy*(gacai*cai[id]*exp(zca*vfrt) - gacao*cao[id])/(1.0 + (zca*vfrt)/2.0 + (zca*vfrt)*(zca*vfrt)/6.0 + (zca*vfrt)*(zca*vfrt)*(zca*vfrt)/24.0);
        ibarna = pna*zna*frdy*(ganai*nai[id]*exp(zna*vfrt) - ganao*nao[id])/(1.0 + (zna*vfrt)/2.0 + (zna*vfrt)*(zna*vfrt)/6.0 + (zna*vfrt)*(zna*vfrt)*(zna*vfrt)/24.0);
        ibark = pk*zk*frdy*(gaki*ki[id]*exp(zk*vfrt) - gako*ko[id])/(1.0 + (zk*vfrt)/2.0 + (zk*vfrt)*(zk*vfrt)/6.0 + (zk*vfrt)*(zk*vfrt)*(zk*vfrt)/24.0);
    }
    else {
        ibarca = pca*zca*(zca*vfrt*frdy)
            *((gacai*cai[id]*exp(zca*vfrt)-gacao*cao[id])/(exp(zca*vfrt)-1.0)); 
        ibarna = pna*zna*(zna*vfrt*frdy)
            *((ganai*nai[id]*exp(zna*vfrt)-ganao*nao[id])/(exp(zna*vfrt)-1.0)); 
        ibark = pk*zk*(zk*vfrt*frdy)
            *((gaki*ki[id]*exp(zk*vfrt)-gako*ko[id])/(exp(zk*vfrt)-1.0));
    }
    /*
    ibarca = pca*zca*zca*((v[id]*frdy*frdy)/(R*temp))
        *((gacai*cai[id]*exp((zca*v[id]*frdy)/(R*temp))-gacao*cao[id])/(exp((zca*v[id]*frdy)/(R*temp))-1.0)); 
    ibarna = pna*zna*zna*((v[id]*frdy*frdy)/(R*temp))
        *((ganai*nai[id]*exp((zna*v[id]*frdy)/(R*temp))-ganao*nao[id])/(exp((zna*v[id]*frdy)/(R*temp))-1.0)); 
    ibark = pk*zk*zk*((v[id]*frdy*frdy)/(R*temp))
        *((gaki*ki[id]*exp((zk*v[id]*frdy)/(R*temp))-gako*ko[id])/(exp((zk*v[id]*frdy)/(R*temp))-1.0));
    */
    fca = 1.0/(1.0+cai[id]/kmca); 
    
    ilca = icalfac[id]*d[id]*f[id]*fca*ibarca; 
    ilcana = icalfac[id]*d[id]*f[id]*fca*ibarna; 
    ilcak = icalfac[id]*d[id]*f[id]*fca*ibark;
}

template <int ncells>
void LR2CellIto<ncells>::comp_icat (int id, double dt, double& db, double& dg, double& icat) 
{ 
    double bss, taub, gss, taug, eca;
    eca = (R*temp/(zca*frdy))*log(cao[id]/cai[id]);
    bss = 1.0/(1.0+exp(-(v[id]+14.0)/10.8)); 
    taub = 3.7+6.1/(1.0+exp((v[id]+25.0)/4.5)); 

    gss = 1.0/(1.0+exp((v[id]+60.0)/5.6)); 
    if (v[id] <= 0.0) {
        taug = -0.875*v[id]+12.0; 
    }
    
    else {
        taug = 12.0; 
    }

    db = (bss-(bss-b[id])*exp(-dt/taub) - b[id])/dt; 
    dg = (gss-(gss-g[id])*exp(-dt/taug) - g[id])/dt; 
    
    icat = gcat*b[id]*b[id]*g[id]*(v[id]-eca); 
}

template <int ncells>
void LR2CellIto<ncells>::comp_ikr (int id, double dt, double& dxr, double& ikr) 
{ 
    double gkr, ekr, xrss, tauxr, r;
    gkr = 0.02614*sqrt(ko[id]/5.4); 
    ekr = ((R*temp)/frdy)*log(ko[id]/ki[id]); 

    xrss = 1.0/(1.0+exp(-(v[id]+21.5)/7.5)); 
    tauxr = 0.0;
    if (fabs(v[id] + 14.2) < 1e-2) {
        tauxr = tauxr + 1.0/(0.00138/0.123);    
    }
    else {
        tauxr = tauxr + 1.0/(0.00138*(v[id]+14.2)/(1.0-exp(-0.123*(v[id]+14.2))));
    }
    if (fabs(v[id] + 38.9) < 1e-2) {
        tauxr = tauxr + 0.00061/0.145;    
    }
    else {
        tauxr = tauxr + 0.00061*(v[id]+38.9)/(exp(0.145*(v[id]+38.9))-1.0);    
    }
    
    //tauxr = 1.0/(0.00138*(v[id]+14.2)/(1.0-exp(-0.123*(v[id]+14.2)))+0.00061*(v[id]+38.9)/(exp(0.145*(v[id]+38.9))-1.0)); 
    
    
    tauxr = tauXfac[id]*tauxr;
    dxr = (xrss-(xrss-xr[id])*exp(-dt/tauxr) - xr[id])/dt; 
    r = 1.0/(1.0+exp((v[id]+9.0)/22.4)); 
    ikr = ikrfac[id]*gkr*xr[id]*r*(v[id]-ekr); 
}

template <int ncells>
void LR2CellIto<ncells>::comp_iks (int id, double dt, double& dxs1, double& dxs2, double& iks) 
{ 
    double gks, eks, xs1ss, xs2ss, tauxs1, tauxs2;
    
    gks = 0.433*(1.0+0.6/(1.0+pow((0.000038/cai[id]),1.4))); 
    eks = ((R*temp)/frdy)*log((ko[id]+prnak*nao[id])/(ki[id]+prnak*nai[id])); 

    xs1ss = 1.0/(1.0+exp(-(v[id]-1.5)/16.7)); 
    xs2ss = xs1ss; 
    if (fabs(v[id] + 30.0) < 1e-2) {
        tauxs1 = 1.0/(0.0000719/0.148 + 0.000131/0.0687);    
    }
    else {
        tauxs1 = 1.0/(0.0000719*(v[id]+30.0)/(1.0-exp(-0.148*(v[id]+30.0)))+0.000131*(v[id]+30.0)/(exp(0.0687*(v[id]+30.0))-1.0)); 
    }
    //tauxs1 = 1.0/(0.0000719*(v[id]+30.0)/(1.0-exp(-0.148*(v[id]+30.0)))+0.000131*(v[id]+30.0)/(exp(0.0687*(v[id]+30.0))-1.0)); 
    tauxs2 = 4.0*tauxs1; 
    dxs1 = (xs1ss-(xs1ss-xs1[id])*exp(-dt/tauxs1) - xs1[id])/dt; 
    dxs2 = (xs2ss-(xs2ss-xs2[id])*exp(-dt/tauxs2) - xs2[id])/dt; 
    iks = iksfac[id]*gks*xs1[id]*xs2[id]*(v[id]-eks); 
} 

template <int ncells>
void LR2CellIto<ncells>::comp_iki (int id, double& iki) 
{ 
    double gki, eki, aki, bki, kin;
    gki = 0.75*(sqrt(ko[id]/5.4)); 
    eki = ((R*temp)/frdy)*log(ko[id]/ki[id]); 

    aki = 1.02/(1.0+exp(0.2385*(v[id]-eki-59.215))); 
    bki = (0.49124*exp(0.08032*(v[id]-eki+5.476))+exp(0.06175*(v[id]-eki-594.31)))/(1.0+exp(-0.5143*(v[id]-eki+4.753))); 
    kin = aki/(aki+bki); 
    
    iki = gki*kin*(v[id]-eki); 
} 

template <int ncells>
void LR2CellIto<ncells>::comp_ikp (int id, double& ikp) 
{ 
    double ekp, kp;
    ekp = ((R*temp)/frdy)*log(ko[id]/ki[id]);

    kp = 1.0/(1.0+exp((7.488-v[id])/5.98)); 

    ikp = gkp*kp*(v[id]-ekp); 
}

template <int ncells>
void LR2CellIto<ncells>::comp_ikna (int id, double& ikna) 
{ 
    double ekna, pona, pov;
    ekna = ((R*temp)/frdy)*log(ko[id]/ki[id]); 
    pona = 0.85/(1.0+pow((kdkna/nai[id]),2.8)); 
    pov = 0.8-0.65/(1+exp((v[id]+125.0)/15.0)); 
    
    ikna = gkna*pona*pov*(v[id]-ekna); 
}

template <int ncells>
void LR2CellIto<ncells>::comp_ikatp (int id, double& ikatp) 
{ 
/* Note: If you wish to use this current in your simulations, there are additional */ 
/* changes which must be made to the code as detailed in Cardiovasc Res 1997;35:256-272 */ 
    double ekatp, gkatp, patp, gkbaratp;
    
    ekatp = ((R*temp)/frdy)*log(ko[id]/ki[id]); 
    gkatp = 0.000195/nicholsarea; 
    patp = 1.0/(1.0+(pow((atpi/katp),hatp))); 
    gkbaratp = gkatp*patp*(pow((ko[id]/4.0),natp)); 

    ikatp = gkbaratp*(v[id]-ekatp); 
} 

#ifdef UCLAito

#define gtof 0.11

template <int ncells>
void LR2CellIto<ncells>::comp_ito (int id, double dt, double& dzdv, double& dydv, double& ito)
{
    double ek, rt1, rt2;
    double xtof_inf, ytof_inf, rt4, rt5, txf, tyf;
    
    //ek = (1.0/frt)*log(xko/xki);
    ek = ((R*temp)/frdy)*log((ko[id])/(ki[id]));
    rt1 = -(v[id] + 3.0)/15.0;
    rt2 = (v[id] + 33.5)/10.0;
    //rt3 = (v[id] + 60.0)/10.0;
    //xtos_inf = 1.0/(1.0 + exp(rt1));
    //ytos_inf = 1.0/(1.0 + exp(rt2));
    //rs_inf = 1.0/(1.0 + exp(rt2));
    //txs = 9.0/(1.0 + exp(-rt1)) + 0.5;
    //tys = 3000.0/(1.0 + exp(rt3)) + 30.0;
    //xitos = gtos*xtos[id]*(ytos[id] + 0.5*rs_inf)*(v[id] - ek);
    //dxtos = xtos_inf - (xtos_inf - xtos[id])*exp(-dt/txs);
    //dytos = ytos_inf - (ytos_inf - ytos[id])*exp(-dt/tys);
    
    xtof_inf = 1.0/(1.0 + exp(rt1));
    ytof_inf = 1.0/(1.0 + exp(rt2));
    rt4 = -(v[id]/30.0)*(v[id]/30.0);
    rt5 = (v[id] + 33.5)/10.0;
    txf = 3.5*exp(rt4) + 1.5;
    tyf = 20.0/(1.0 + exp(rt5)) + 20.0;
    //xitof = itofac[id]*gtof*xtof[id]*ytof[id]*(v[id] - ek);
    ito = itofac[id]*gtof*zdv[id]*ydv[id]*(v[id] - ek);
    dzdv = (xtof_inf - (xtof_inf - zdv[id])*exp(-dt/txf) - zdv[id])/dt; //z = x
    dydv = (ytof_inf - (ytof_inf - ydv[id])*exp(-dt/tyf) - ydv[id])/dt;
    
    //return xitof; 
}

#else

template <int ncells>
void LR2CellIto<ncells>::comp_ito (int id, double dt, double& dzdv, double& dydv, double& ito) 
{ 
    double ekdv, rvdv, azdv, bzdv, tauzdv, zssdv, aydv, bydv, tauydv, yssdv;
    
    ekdv = ((R*temp)/frdy)*log((ko[id])/(ki[id])); 
    rvdv = exp(v[id]/100.0); 
    azdv = (10.0*exp((v[id]-40.0)/25.0))/(1.0+exp((v[id]-40.0)/25.0)); 
    bzdv = (10.0*exp(-(v[id]+90.0)/25.0))/(1.0+exp(-(v[id]+90.0)/25.0)); 
    tauzdv = 1.0/(azdv+bzdv); 
    zssdv = azdv/(azdv+bzdv); 
    dzdv = (zssdv-(zssdv-zdv[id])*exp(-dt/tauzdv) - zdv[id])/dt; 

    aydv = 0.015/(1.0+exp((v[id]+60.0)/5.0)); 
    bydv = (0.1*exp((v[id]+25.0)/5.0))/(1.0+exp((v[id]+25.0)/5.0)); 
    tauydv = 1.0/(aydv+bydv); 
    yssdv = aydv/(aydv+bydv); 
    dydv = (yssdv-(yssdv-ydv[id])*exp(-dt/tauydv) - ydv[id])/dt; 
    ito = itofac[id]*gitodv*zdv[id]*zdv[id]*zdv[id]*ydv[id]*rvdv*(v[id]-ekdv); 
}

#endif

template <int ncells>
void LR2CellIto<ncells>::comp_isk (int id, double& isk)
{
    double ek, z;
    double rat = pow(skh[id]/cai[id],skn[id]);
    z = 1.0/(1.0 + rat);
    ek = ((R*temp)/frdy)*log((ko[id])/(ki[id]));
    
    isk = iskfac[id]*gsk*z*(v[id] - ek);
}

template <int ncells>
void LR2CellIto<ncells>::comp_inaca (int id, double& inaca) 
{ 
    inaca = nacafac[id]*c1*exp((gammas-1.0)*v[id]*frdy/(R*temp))
*((exp(v[id]*frdy/(R*temp))*nai[id]*nai[id]*nai[id]*cao[id]-nao[id]*nao[id]*nao[id]*cai[id])
/(1.0+c2*exp((gammas-1)*v[id]*frdy/(R*temp))*(exp(v[id]*frdy/(R*temp))*nai[id]*nai[id]*nai[id]*cao[id]+nao[id]*nao[id]*nao[id]*cai[id])));  
}

template <int ncells>
void LR2CellIto<ncells>::comp_inak (int id, double& inak) 
{ 
    double sigma, fnak;
    sigma = (exp(nao[id]/67.3)-1.0)/7; 

    fnak = 1.0/(1.0 + 0.1245*exp((-0.1*v[id]*frdy)/(R*temp))+0.0365*sigma*exp((-v[id]*frdy)/(R*temp))); 

    inak = ibarnak*fnak*(1.0/(1.0+pow(kmnai/nai[id],2.0)))*(ko[id]/(ko[id]+kmko)); 
}

template <int ncells>
void LR2CellIto<ncells>::comp_insca (int id, double& insna, double& insk) 
{ 
    double ibarnsna, ibarnsk;
    ibarnsna = pnsca*zna*zna*((v[id]*frdy*frdy)/(R*temp))
*((ganai*nai[id]*exp((zna*v[id]*frdy)/(R*temp))-ganao*nao[id])/(exp((zna*v[id]*frdy)/(R*temp))-1.0)); 
    ibarnsk = pnsca*zk*zk*((v[id]*frdy*frdy)/(R*temp))
*((gaki*ki[id]*exp((zk*v[id]*frdy)/(R*temp))-gako*ko[id])/(exp((zk*v[id]*frdy)/(R*temp))-1.0)); 

    insna = ibarnsna/(1.0+pow(kmnsca/cai[id],3.0)); 
    insk = ibarnsk/(1.0+pow(kmnsca/cai[id],3.0)); 
}

template <int ncells>
void LR2CellIto<ncells>::comp_ipca (int id, double& ipca) 
{ 
    ipca = (ibarpca*cai[id])/(kmpca+cai[id]); 
} 

template <int ncells>
void LR2CellIto<ncells>::comp_icab (int id, double& icab) 
{ 
    double ecan; 
    ecan = ((R*temp)/(2.0*frdy))*log(cao[id]/cai[id]); 
    icab = gcab*(v[id]-ecan); 
}

template <int ncells>
void LR2CellIto<ncells>::comp_inab (int id, double& inab) 
{ 
    double enan; 
    enan = ((R*temp)/frdy)*log(nao[id]/nai[id]); 
    inab = gnab*(v[id]-enan); 
} 

template <int ncells>
void LR2CellIto<ncells>::conc_nai (int id, double naiont, double& dnai) 
{ 
// The units of dnai is in mM. Note that naiont should be multiplied by the 
// cell capacitance to get the correct units. Since cell capacitance = 1 uF/cm^2, 
// it doesn't explicitly appear in the equation below. 
// This holds true for the calculation of dki and dcai. */ 

    dnai = -(naiont*acap)/(vmyo*zna*frdy); 
    //nai = dnai + nai; 
} 

template <int ncells>
void LR2CellIto<ncells>::conc_ki (int id, double kiont, double st, double& dki) 
{ 
    dki = -((kiont+st)*acap)/(vmyo*zk*frdy); 
} 

template <int ncells>
void LR2CellIto<ncells>::calc_dyn (int id, double dt, double caiont, double& dcaiont, double& dtcicr, double& djsr, double& dnsr, double& dcai) 
{ 
    double itr, kleak, ileak, iup, on, off, magrel, irelcicr;
    double csqn, bjsr, cjsr, trpn, cmdn, catotal, bmyo, cmyo, dmyo, gpig;
    itr = (nsr[id]-jsr[id])/tautr; 
    
    kleak = iupbar/nsrbar; 
    ileak = kleak*nsr[id]; 

    iup = iupfac[id]*iupbar*cai[id]/(cai[id]+kmup); //THIS IS SERCA PUMP, INCREASE THIS

    dcaiont = (caiont-caiontold[id])/dt; 

    /*
    if(v[id] >-35 && dcaiont>dcaiontold[id] && flag[id]==0) 
    {
        flag[id] = 1; 
        tcicr[id] = 0;
    } 
    */
    on = 1.0/(1.0+exp((-tcicr[id]+4)/tauon)); 
    off = (1.0-1.0/(1.0+exp((-tcicr[id]+4)/tauoff))); 
    magrel = 1.0/(1.0+exp((caiont+5)/0.9)); 
    irelcicr = gmaxrel*on*off*magrel*(jsr[id]-cai[id]); 

    dtcicr = 1.0; 

    /*
    greljsrol = grelbarjsrol*(1.0-exp(-tjsrol[id]/tauon))*exp(-tjsrol[id]/tauoff); 
    ireljsrol = greljsrol*(jsr[id]-cai[id]); 
    
    dtjsrol = 1; 
    */

    csqn = csqnbar*(jsr[id]/(jsr[id]+kmcsqn)); 
    djsr = dt*(itr-irelcicr); //djsr = dt*(itr-irelcicr-ireljsrol); 
    bjsr = csqnbar-csqn-djsr-jsr[id]+kmcsqn; 
    cjsr = kmcsqn*(csqn+djsr+jsr[id]); 

    djsr = ((sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2.0 - jsr[id])/dt;
    
    dnsr = (iup-ileak-itr*vjsr/vnsr); 
    //nsr = nsr+dnsr;
    
    dcai = -dt*(((caiont*acap)/(vmyo*zca*frdy))+((iup-ileak)*vnsr/vmyo)
    -(irelcicr*vjsr/vmyo)); //dcai = -dt*(((caiont*acap)/(vmyo*zca*frdy))+((iup-ileak)*vnsr/vmyo)-(irelcicr*vjsr/vmyo)-(ireljsrol*vjsr/vmyo)); 
    trpn = trpnbar*(cai[id]/(cai[id]+kmtrpn)); 
    cmdn = cmdnbar*(cai[id]/(cai[id]+kmcmdn)); 

    catotal = trpn+cmdn+dcai+cai[id]; 
    bmyo = cmdnbar+trpnbar-catotal+kmtrpn+kmcmdn; 
    cmyo = (kmcmdn*kmtrpn)-(catotal*(kmtrpn+kmcmdn))+(trpnbar*kmcmdn)+(cmdnbar*kmtrpn); 
    dmyo = -kmtrpn*kmcmdn*catotal; 
    gpig = sqrt(bmyo*bmyo-3*cmyo); 

    dcai = ((2.0*gpig/3.0)*cos(acos((9.0*bmyo*cmyo-2.0*bmyo*bmyo*bmyo-27.0*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-(bmyo/3.0) - cai[id])/dt;
} 

template <int ncells>
void LR2CellIto<ncells>::conc_cleft(int id, double dt, double st, double naiont, double kiont, double caiont, double& dnao, double& dko, double& dcao) 
{ 
    dnao = ((nabm-nao[id])/taudiff+naiont*acap/(vcleft*frdy)); 
    //nao = dnao+nao;

    dko = ((kbm-ko[id])/taudiff+(kiont+st)*acap/(vcleft*frdy));
    //ko = dko+ko;
    
    dcao = ((cabm-cao[id])/taudiff+caiont*acap/(vcleft*frdy*2)); 
    //cao = dcao+cao; 
} 

template <int ncells>
void LR2CellIto<ncells>::setcell (int id, LR2CellIto<1>* newcell)
{
    v[id] = newcell->v[0];
    m[id] = newcell->m[0];
    h[id] = newcell->h[0];
    j[id] = newcell->j[0];
    d[id] = newcell->d[0];
    f[id] = newcell->f[0];
    b[id] = newcell->b[0];
    g[id] = newcell->g[0];
    xr[id] = newcell->xr[0];
    xs1[id] = newcell->xs1[0];
    xs2[id] = newcell->xs2[0];
    zdv[id] = newcell->zdv[0];
    ydv[id] = newcell->ydv[0];
    nai[id] = newcell->nai[0];
    ki[id] = newcell->ki[0];
    cai[id] = newcell->cai[0];
    jsr[id] = newcell->jsr[0];
    nsr[id] = newcell->nsr[0];
    nao[id] = newcell->nao[0];
    ko[id] = newcell->ko[0];
    cao[id] = newcell->cao[0];
    flag[id] = newcell->flag[0];
    
    tcicr[id] = newcell->tcicr[0];
    dcaiontold[id] = newcell->dcaiontold[0];
    caiontold[id] = newcell->caiontold[0];
    
    icalfac[id] = newcell->icalfac[0];
    itofac[id] = newcell->itofac[0];
    iskfac[id] = newcell->iskfac[0];
    skh[id] = newcell->skh[0];
    tauXfac[id] = newcell->tauXfac[0];
    iupfac[id] = newcell->iupfac[0];
    ikrfac[id] = newcell->ikrfac[0];
    iksfac[id] = newcell->iksfac[0];
    nacafac[id] = newcell->nacafac[0];
}

template <int ncells>
void LR2CellIto<ncells>::getcell (int id, LR2CellIto<1>* newcell)
{
    newcell->v[0] = v[id];
    newcell->m[0] = m[id];
    newcell->h[0] = h[id];
    newcell->j[0] = j[id];
    newcell->d[0] = d[id];
    newcell->f[0] = f[id];
    newcell->b[0] = b[id];
    newcell->g[0] = g[id];
    newcell->xr[0] = xr[id];
    newcell->xs1[0] = xs1[id];
    newcell->xs2[0] = xs2[id];
    newcell->zdv[0] = zdv[id];
    newcell->ydv[0] = ydv[id];
    newcell->nai[0] = nai[id];
    newcell->ki[0] = ki[id];
    newcell->cai[0] = cai[id];
    newcell->jsr[0] = jsr[id];
    newcell->nsr[0] = nsr[id];
    newcell->nao[0] = nao[id];
    newcell->ko[0] = ko[id];
    newcell->cao[0] = cao[id];
    newcell->flag[0] = flag[id];
    
    newcell->tcicr[0] = tcicr[id];
    newcell->dcaiontold[0] = dcaiontold[id];
    newcell->caiontold[0] = caiontold[id];
    
    newcell->icalfac[0] = icalfac[id];
    newcell->itofac[0] = itofac[id];
    newcell->iskfac[0] = iskfac[id];
    newcell->skh[0] = skh[id];
    newcell->tauXfac[0] = tauXfac[id];
    newcell->iupfac[0] = iupfac[id];
    newcell->ikrfac[0] = ikrfac[id];
    newcell->iksfac[0] = iksfac[id];
    newcell->nacafac[0] = nacafac[id];
}

template <int ncells>
void LR2CellIto<ncells>::saveconditions(FILE* file, int id, bool header, double t) {
    if (header) {
        fprintf(file,"t\tv\tm\th\tj\td\tf\tb\tg\txr\txs1\txs2\tzdv\tydv\tnai\tki\tcai\tjsr\tnsr\tnao\tko\tcao\tflag\ttcicr\n");
    }
    fprintf(file,"%g\t",t);
    fprintf(file,"%.12f\t",v[id]);
    fprintf(file,"%.12f\t",m[id]);
    fprintf(file,"%.12f\t",h[id]);
    fprintf(file,"%.12f\t",j[id]);
    fprintf(file,"%.12f\t",d[id]);
    fprintf(file,"%.12f\t",f[id]);
    fprintf(file,"%.12f\t",b[id]);
    fprintf(file,"%.12f\t",g[id]);
    fprintf(file,"%.12f\t",xr[id]);
    fprintf(file,"%.12f\t",xs1[id]);
    fprintf(file,"%.12f\t",xs2[id]);
    fprintf(file,"%.12f\t",zdv[id]);
    fprintf(file,"%.12f\t",ydv[id]);
    fprintf(file,"%.12f\t",nai[id]);
    fprintf(file,"%.12f\t",ki[id]);
    fprintf(file,"%.12f\t",cai[id]);
    fprintf(file,"%.12f\n",jsr[id]);
    fprintf(file,"%.12f\t",nsr[id]);
    fprintf(file,"%.12f\t",nao[id]);
    fprintf(file,"%.12f\t",ko[id]);
    fprintf(file,"%.12f\t",cao[id]);
    fprintf(file,"%d\t",flag[id]);
    fprintf(file,"%.12f\n",tcicr[id]);
}

template <int ncells>
void LR2CellIto<ncells>::write(std::fstream& file) {
    file.write((char*)&this, sizeof(this) );
}

template <int ncells>
void LR2CellIto<ncells>::read(std::fstream& file) {
    file.read((char*)&this, sizeof(this) );
}

    #undef l   
    #undef a 
    #undef pi 
    #undef vcell  
    #undef ageo 
    #undef acap 
    #undef vmyo  
    #undef vmito  
    #undef vsr  
    #undef vnsr 
    #undef vjsr 
    #undef vcleft 
    
    /* Terms for Solution of Conductance and Reversal Potential */ 
    #undef R 
    #undef frdy  
    #undef temp 
    
    /* Ion Valences */ 
    #undef zna 
    #undef zk  
    #undef zca 
    
    /* Ion Concentrations */ 
    #undef taudiff  
    
    /* NSR Ca Ion Concentration Changes */ 
    #undef kmup 
    #undef iupbar  
    #undef nsrbar 
    
    /* JSR Ca Ion Concentration Changes */ 
    #undef tauon 
    #undef tauoff 
    #undef csqnth
    #undef gmaxrel  
    #undef csqnbar 
    #undef kmcsqn 
    
    
    /* Translocation of Ca Ions from NSR to JSR */ 
    #undef tautr 

    /* Myoplasmic Ca Ion Concentration Changes */ 
    #undef cmdnbar  
    #undef trpnbar 
    #undef kmcmdn 
    #undef kmtrpn 
    
    /* Fast Sodium Current (time dependant) */ 
    #undef gna 
    
    /* Current through L-type Ca Channel */ 
    #undef kmca  
    #undef pca 
    #undef gacai 
    #undef gacao 
    #undef pna  
    #undef ganai  
    #undef ganao  
    #undef pk 
    #undef gaki 
    #undef gako 
    
    /* Current through T-type Ca Channel */ 
    #undef gcat 
    
    /* Rapidly Activating Potassium Current */ 
    // Constant gkr embedded in the formula as 0.02614*sqrt(ko[id]/5.4)
    
    /* Slowly Activating Potassium Current */ 
    // Constant gks embedded in the formula as 0.433*(1+0.6/(1+pow((0.000038/cai[id]),1.4))); 
    #undef prnak 
    
    /* Potassium Current (time-independent) */ 
    // Constant gki embedded in the formula as gki = 0.75*(sqrt(ko[id]/5.4)); 
    
    /* Plateau Potassium Current */ 
    #undef gkp 
    
    /* Na-Activated K Channel */ 
    #undef gkna 
    #undef nkna 
    #undef kdkna 
    
    /* ATP-Sensitive K Channel */ 
    #undef natp 
    #undef nicholsarea 
    #undef atpi 
    #undef hatp 
    #undef katp 
    
    /* Ito Transient Outward Current (Dumaine et al. Circ Res 1999;85:803-809) */  
    #undef gitodv 

    /* Isk calcium-dependent potassium current */
    #undef gsk 
    #undef skh 
    #undef skn

    /* Sodium-Calcium Exchanger V-S */ 
    #undef c1  
    #undef c2 
    #undef gammas  
    
    /* Sodium-Potassium Pump */ 
    #undef ibarnak  
    #undef kmnai 
    #undef kmko
    
    /* Nonspecific Ca-activated Current */ 
    #undef pnsca  
    #undef kmnsca 
    
    /* Sarcolemmal Ca Pump */ 
    #undef ibarpca  
    #undef kmpca 
    
    /* Ca Background Current */ 
    #undef gcab 
    
    /* Na Background Current */ 
    #undef gnab
    
    /* Cleft Space */
    #undef nabm 
    #undef kbm 
    #undef cabm 



#endif
