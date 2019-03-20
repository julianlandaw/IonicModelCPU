//
//  UCLACell.cpp
//
//  Implementation of the UCLA Model
//  Created by Julian Landaw on 5/30/18.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#include <math.h>
#include <fstream>
#include "UCLACell.h"

#ifndef UCLACell_cpp
#define UCLACell_cpp

#define N 16
#define vc -80.0
#define temp 308.0
#define xxr 8.314
#define xf 96.485
#define frt (xf/(xxr*temp))
#define xnao 136.0
#define xki 140.0
#define xko 5.40
#define cao 1.8
//#define ek ((1.0/frt)*log(xko/xki))

#ifdef xiaodong
#define gca 40.0
#define gnaca 1.0
#define gkr 0.0
#define gks 0.50
#else
#define gca 182.0
#define gnaca 0.84
#define gkr 0.0125
#define gks 0.1386
#endif

#define gtos 0.04
#define gtof 0.11
#define vup 0.4
#define gna 12.0
#define gkix 0.3
#define gnak 1.5

#define taur 30.0
#define taud 4.0
#define taua 100.0
#define av 11.3
#define cstar 90.0
#define gleak 0.00002069
#define kj 50.0
#define cup 0.5

#ifndef DV_MAX
#define DV_MAX 0.1
#endif

#ifndef ADAPTIVE
#define ADAPTIVE 10
#endif

#ifndef stimulus
#define stimulus -80.0
#endif

#ifndef stimduration
#define stimduration 0.5
#endif

template <int ncells>
UCLACell<ncells>::UCLACell()
{
    for (int i = 0; i < ncells; i++) {
        xm[i] = 0.001145222753;// sodium m-gate
        xh[i] = 0.9898351676;// sodium h-gate
        xj[i] = 0.9930817518;// soiumj-gate

        xr[i] = 1.0; //0.008709989976;// ikr gate variable 
        xs1[i] = 1.0; //0.08433669901;// iks gate variable
        xs2[i] = 1.0; //0.1412866149;// iks gate varaible 

        xtos[i] = 0.003757746357;// ito slow activation
        ytos[i] = 0.1553336368;// ito slow inactivation
    
        v[i] = -86.79545769; // voltage

        cp[i] = 1.682601371;// averaged dyadic space con.
        cs[i] = 0.3205609256;// averaged submembrane conc.
        ci[i] = 0.3863687451;// myoplasm conc.

        cj[i] = 107.0388739;// NSR load
        cjp[i] = 95.76256179;// average JSR load

        xir[i] = 0.006462569526;// SR current flux

        // Markov gate variables 

        #ifdef xiaodong
        d[i] = 0.0;
        f[i] = 1.0;
        #else
        c1[i] = 1.925580885e-05;// C1
        c2[i] = 0.9535940241;// C2
        xi1ca[i] = 0.007052299702;// I1_Ca
        xi1ba[i] = 3.629261123e-05;// I1_Ba
        xi2ca[i] = 0.02316349806;// I2_Ca
        xi2ba[i] = 0.01613268649;// I2_Ba
        #endif

        nai[i] = 14.01807252;// internal Na conc.

        xtof[i] = 0.003737842131;// ito fast activation
        ytof[i] = 0.9823715315;// ito slow inactivation

        tropi[i] = 29.64807803;// time dependent buffers in myplasm (troponin)
        trops[i] = 26.37726416;// time dependent buffers in submembrane (troponin)

        vold[i] = v[i];
        jparam[i]=1.0;
        
        diffcurrent[i] = 0.0;
        
        inafac[i] = 1.0;
        itofac[i] = 1.0;
        itoslowfac[i] = 1.0;
        ikrfac[i] = 1.0;
        iksfac[i] = 1.0;
        nacafac[i] = 1.0;
        icalfac[i] = 1.0;
        // SK variables. default is SK current = 0
        iskfac[i] = 0.0;
        skh[i] = 0.6;
        skn[i] = 4;

#ifdef clampnai
        naiclamped[i] = true;
#else
        naiclamped[i] = false;
#endif
    }
}

template <int ncells>
bool UCLACell<ncells>::iterate(const int id, double dt, double st, double dv_max) {
    double xik1, xito, xinak, csm, jnaca, jd, xkon, xkoff, btrop, xbi, xbs, jup, jleak;
    double po, rxa, jca, dcs, dci, dcj, dcjp, Qr, dir, dcp, xina, xikr, xiks, wca;
    double xinaca, xica, xrr, dv;
    double dxtos, dytos, dxtof, dytof;
#ifdef xiaodong
    double dd, df;
#else
    double dc1, dc2, dxi1ca, dxi2ca, dxi1ba, dxi2ba;
#endif
    double dxm, dxh, dxj;
    double dxr, dxs1, dxs2;
    double xisk;
    
    xisk = comp_isk(id);
    xik1 = comp_ik1(id);
    xito = comp_ito(id, dt, dxtos, dytos, dxtof, dytof);
    xinak = comp_inak(id);
    csm = cs[id]/1000.0;
    jnaca = comp_inaca(id, csm);
    jd = (cs[id] - ci[id])/taud;
    
    xkon = 0.0327;
    xkoff = 0.0196;
    btrop = 70.0;
    xbi = xkon*ci[id]*(btrop - tropi[id]) - xkoff*tropi[id];
    xbs = xkon*cs[id]*(btrop - trops[id]) - xkoff*trops[id];
    
    jup = comp_iup(id);
    jleak = comp_ileak(id);
    
#ifdef xiaodong
    po = comp_icalpo(id, dt, dd, df);
#else
    po = comp_icalpo(id, dc1, dc2, dxi1ca, dxi2ca, dxi1ba, dxi2ba);
#endif
    rxa = comp_rxa(id, csm);
    jca = icalfac[id]*gca*po*rxa;
    
    dcs = comp_inst_buffer(cs[id])*(50.0*(xir[id] - jd - jca + jnaca)-xbs);
    dci = comp_inst_buffer(ci[id])*(jd - jup + jleak - xbi);
    
    dcj = -xir[id] + jup - jleak;
    dcjp = (cj[id] - cjp[id])/taua;
    Qr = comp_Q(id);
    dir = comp_dir(id, po, Qr, rxa, dcj);
    dcp = comp_dcp(id, po, Qr, rxa);
    
    xina = comp_ina(id, dt, dxm, dxh, dxj);
    xikr = comp_ikr(id, dt, dxr);
    xiks = comp_iks(id, dt, dxs1, dxs2);
    
    wca = 8.0;
    xinaca = wca*jnaca;
    xica = 2.0*wca*jca;
    xrr = (1.0/wca)/1000.0;
    //nai += (-xrr*(xina+3.0*xinak+3.0*xinaca))*dt;
    dv = (diffcurrent[id] - (xina + xik1 + xikr + xiks + xisk + xito + xinaca + xica + xinak + st))*dt;
    
    //Comment this out to remove adaptive time-stepping
    if (dv_max > 0 && dv*dv > dv_max*dv_max) {return false;}
//#ifndef clampnai 
    if (!naiclamped[id]) {
        nai[id] += (-xrr*(xina + 3.0*xinak + 3.0*xinaca))*dt;
    }
//#endif
    vold[id] = v[id]; // Save the old voltage
    v[id] += dv;
    xm[id] = dxm;
    xh[id] = dxh;
    xj[id] = dxj;
    xr[id] = dxr;
    xs1[id] = dxs1;
    xs2[id] = dxs2;
    xtos[id] = dxtos;
    ytos[id] = dytos;
    xtof[id] = dxtof;
    ytof[id] = dytof;
    cp[id] += dcp*dt;
    cs[id] += dcs*dt;
    ci[id] += dci*dt;
    cj[id] += dcj*dt;
    xir[id] += dir*dt;
    cjp[id] += dcjp*dt;

    tropi[id] += xbi*dt;
    trops[id] += xbs*dt;
    
#ifdef xiaodong
    d[id] = dd;
    f[id] = df;
#else
    c1[id] += dc1*dt;
    c2[id] += dc2*dt;
    xi1ca[id] += dxi1ca*dt;
    xi2ca[id] += dxi2ca*dt;
    xi1ba[id] += dxi1ba*dt;
    xi2ba[id] += dxi2ba*dt;
#endif
    
    return true;
}

template <int ncells>
void UCLACell<ncells>::stepdt (const int id, double dt, double st) {
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
double UCLACell<ncells>::comp_ina (int id, double dt, double& dxm, double& dxh, double& dxj) //Fast Sodium Current
{
    double ena, taum, tauh, tauj, mss, hss, jss, xina;
    double am, bm, ah, bh, aj, bj;
    
    ena = (1.0/frt)*log(xnao/nai[id]);
    am = 3.2;
    if (fabs(v[id] + 47.13) > 0.001) {
        am = 0.32*(v[id]+47.13)/(1.0-exp(-0.1*(v[id]+47.13)));
    }
    
    bm = 0.08*exp(-v[id]/11.0);
    
    if (v[id] < -40.0) {
        ah = 0.135*exp((80+v[id])/-6.8);
        bh = 3.56*exp(0.079*v[id])+310000.0*exp(0.35*v[id]);
        aj = (-127140.0*exp(0.2444*v[id])-0.00003474*exp(-0.04391*v[id]))*((v[id]+37.78)/(1.0+exp(0.311*(v[id]+79.23))));
        bj = (0.1212*exp(-0.01052*v[id]))/(1.0+exp(-0.1378*(v[id]+40.14)));
    }
    else {
        ah = 0.0;
        bh = 1.0/(0.13*(1.0+exp((v[id]+10.66)/-11.1)));
        aj = 0.0;
        bj = (0.3*exp(-0.0000002535*v[id]))/(1.0+exp(-0.1*(v[id]+32.0)));
    }
    
    taum = 1.0/(am+bm);
    tauh = 1.0/(ah+bh);
    tauj = 1.0/(aj+bj)*jparam[id];
    
    mss = am*taum;
    hss = ah*tauh;
    jss = aj*tauj;
    
    xina = inafac[id]*gna*xm[id]*xm[id]*xm[id]*xh[id]*xj[id]*(v[id]-ena);
    
    dxm = mss-(mss-xm[id])*exp(-dt/taum);
    dxh = hss-(hss-xh[id])*exp(-dt/tauh);
    dxj = jss-(jss-xj[id])*exp(-dt/tauj);
    
    return xina;
}

template <int ncells>
double UCLACell<ncells>::comp_ito (int id, double dt, double& dxtos, double& dytos, double& dxtof, double& dytof)
{
    double ek, rt1, rt2, rt3, xtos_inf, ytos_inf, rs_inf, txs, tys, xitos;
    double xtof_inf, ytof_inf, rt4, rt5, txf, tyf, xitof;
    
    ek = (1.0/frt)*log(xko/xki);
    rt1 = -(v[id] + 3.0)/15.0;
    rt2 = (v[id] + 33.5)/10.0;
    rt3 = (v[id] + 60.0)/10.0;
    xtos_inf = 1.0/(1.0 + exp(rt1));
    ytos_inf = 1.0/(1.0 + exp(rt2));
    rs_inf = 1.0/(1.0 + exp(rt2));
    txs = 9.0/(1.0 + exp(-rt1)) + 0.5;
    tys = 3000.0/(1.0 + exp(rt3)) + 30.0;
    xitos = itoslowfac[id]*gtos*xtos[id]*(ytos[id] + 0.5*rs_inf)*(v[id] - ek);
    dxtos = xtos_inf - (xtos_inf - xtos[id])*exp(-dt/txs);
    dytos = ytos_inf - (ytos_inf - ytos[id])*exp(-dt/tys);
    
    xtof_inf = xtos_inf;
    ytof_inf = ytos_inf;
    rt4 = -(v[id]/30.0)*(v[id]/30.0);
    rt5 = (v[id] + 33.5)/10.0;
    txf = 3.5*exp(rt4) + 1.5;
    tyf = 20.0/(1.0 + exp(rt5)) + 20.0;
    xitof = itofac[id]*gtof*xtof[id]*ytof[id]*(v[id] - ek);
    dxtof = xtof_inf - (xtof_inf - xtof[id])*exp(-dt/txf);
    dytof = ytof_inf - (ytof_inf - ytof[id])*exp(-dt/tyf);
    
    #ifdef noslowito
    return xitof;
    #else
    return xitos + xitof;
    #endif
    
}

template <int ncells>
double UCLACell<ncells>::comp_ikr (int id, double dt, double& dxr)
{
    double ek, gss, xkrv1, xkrv2, taukr, xkrinf, rg, xikr;
    
    ek = (1.0/frt)*log(xko/xki);
    
    gss = sqrt(xko/5.4);
    xkrv1 = 0.00138/0.123;
    if (fabs(v[id] + 7.0) > 0.001) {
        xkrv1 = 0.00138*(v[id] + 7.0)/(1.0 - exp(-0.123*(v[id] + 7.0)));
    }
    
    xkrv2 = 0.00061/0.145;
    if (fabs(v[id] + 10.0) > 0.001) {
        xkrv2 = 0.00061*(v[id] + 10.0)/(exp(0.145*(v[id] + 10.0)) - 1.0);
    }
    
    taukr = 1.0/(xkrv1 + xkrv2);
    xkrinf = 1.0/(1.0 + exp(-(v[id] + 50.0)/7.5));
    rg = 1.0/(1.0 + exp((v[id] + 33.0)/22.4));
    
    xikr = ikrfac[id]*gkr*gss*xr[id]*rg*(v[id] - ek);
    
    dxr = xkrinf - (xkrinf - xr[id])*exp(-dt/taukr);
    
    return xikr;
}

#define prnak 0.018330

template <int ncells>
double UCLACell<ncells>::comp_iks(int id, double dt, double& dxs1, double& dxs2)
{
    double eks, xs1ss, xs2ss, tauxs1, tauxs2, gksx, xiks;
    //prnak = 0.018330;
    eks = (1.0/frt)*log((xko + prnak*xnao)/(xki + prnak*nai[id]));
    xs1ss = 1.0/(1.0 + exp(-(v[id] - 1.50)/16.70));
    xs2ss = xs1ss;
    tauxs1 = 1.0/(0.0000719/0.148 + 0.000131/0.0687);
    if (fabs(v[id] + 30.0) > 0.001) {
        tauxs1 = 1.0/(0.0000719*(v[id] + 30.0)/(1.0 - exp(-0.148*(v[id] + 30.0))) + 0.000131*(v[id] + 30.0)/(exp(0.0687*(v[id] + 30.0)) - 1.0));
    }
    tauxs2 = 4*tauxs1;
    gksx = (1.0 + 0.8/(1.0 + (0.5/ci[id])*(0.5/ci[id])*(0.5/ci[id])));
    xiks = iksfac[id]*gks*gksx*xs1[id]*xs2[id]*(v[id] - eks);
    
    dxs1 = xs1ss - (xs1ss - xs1[id])*exp(-dt/tauxs1);
    dxs2 = xs2ss - (xs2ss - xs2[id])*exp(-dt/tauxs2);
    
    return xiks;
}

#define gki (sqrt(xko/5.4))

template <int ncells>
double UCLACell<ncells>::comp_ik1 (int id)
{
    double ek, aki, bki, xkin, xik1;
    
    ek = (1.0/frt)*log(xko/xki);
    //gki = sqrt(xko/5.4);
    aki = 1.02/(1.0 + exp(0.2385*(v[id] - ek - 59.215)));
    bki = (0.49124*exp(0.08032*(v[id] - ek + 5.465)) + exp(0.061750*(v[id] - ek - 594.31)))/(1.0 + exp(-0.5143*(v[id] - ek + 4.753)));
    xkin = aki/(aki + bki);
    xik1 = gkix*gki*xkin*(v[id] - ek);
    return xik1;
}

#define xkmko 1.5
#define xkmnai 12.0
#define sigma ((exp(xnao/67.3)-1.0)/7.0)

template <int ncells>
double UCLACell<ncells>::comp_inak(int id) {
    double fnak, xinak;
    //xkmko = 1.5;
    //xkmnai = 12.0;
    //sigma = (exp(xnao/67.3)-1.0)/7.0;
    fnak = 1.0/(1.0 + 0.1245*exp(-0.1*v[id]*frt) + 0.0365*sigma*exp(-v[id]*frt));
    xinak = gnak*fnak*(1.0/(1.0 + (xkmnai/nai[id])))*xko/(xko+xkmko);
    return xinak;
}

#define xkdna 0.3
#define xmcao 1.3
#define xmnao 87.5
#define xmnai 12.3
#define xmcai 0.0036

template <int ncells>
double UCLACell<ncells>::comp_inaca(int id, double csm) {
    double zw3, zw4, aloss, yz1, yz2, yz3, yz4, zw8, xinacaq;
    zw3 = (nai[id]*nai[id]*nai[id])*cao*exp(v[id]*0.35*frt) - (xnao*xnao*xnao)*csm*exp(v[id]*(0.35 - 1.0)*frt);
    zw4 = 1.0 + 0.2*exp(v[id]*(0.35 - 1.0)*frt);
    //xkdna = 0.3;
    aloss = 1.0/(1.0 + (xkdna/cs[id])*(xkdna/cs[id])*(xkdna/cs[id]));
    //xmcao = 1.3;
    //xmnao = 87.5;
    //xmnai = 12.3;
    //xmcai = 0.0036;
    yz1 = xmcao*(nai[id]*nai[id]*nai[id]) + (xmnao*xmnao*xmnao)*csm;
    yz2 = (xmnai*xmnai*xmnai)*cao*(1.0+csm/xmcai);
    yz3 = xmcai*(xnao*xnao*xnao)*(1.0 + (nai[id]/xmnai)*(nai[id]/xmnai)*(nai[id]/xmnai));
    yz4 = (nai[id]*nai[id]*nai[id])*cao + (xnao*xnao*xnao)*csm;
    zw8 = yz1 + yz2 + yz3 + yz4;
    xinacaq = nacafac[id]*gnaca*aloss*zw3/(zw4*zw8);
    return xinacaq;
}

#define pca 0.00054

template <int ncells>
double UCLACell<ncells>::comp_rxa(int id, double csm) 
{
    double za, factor1, factor, rxa;
    //pca = 0.00054;
    za = v[id]*2.0*frt;
    factor1 = 4.0*pca*xf*frt;
    factor = v[id]*factor1;
    if (fabs(za) < 0.001) {
        rxa = factor1*(csm*exp(za)-0.341*(cao))/(2.0*frt);
    }
    else {
        rxa = factor*(csm*exp(za)-0.341*(cao))/(exp(za)-1.0);
    }
    return rxa;
}

#ifdef xiaodong
template <int ncells>
double UCLACell<ncells>::comp_icalpo(int id, double dt, double& dd, double& df) {
    double dss, taudx, fss, tauf, fca, po;
    dss = 1.0/(1.0+exp(-(v[id]+10.0)/6.24)); 
    taudx = dss*1.0/(6.24*0.035);
    if (fabs(v[id] + 10.0) > 0.001) {
        taudx = dss*(1.0-exp(-(v[id]+10.0)/6.24))/(0.035*(v[id]+10.0)); 
    }
    //taudx = dss*(1.0-exp(-(v[id]+10.0)/6.24))/(0.035*(v[id]+10.0)); 

    fss = (1.0/(1.0+exp((v[id]+32.0)/8.0)))+(0.6/(1.0+exp((50.0-v[id])/20.0))); 
    tauf = 1.0/(0.0197*exp(-(0.0337*(v[id]+10.0)*0.0337*(v[id]+10.0))) + 0.02);
    //tauf = 1.0/(0.0197*exp(-pow(0.0337*(v[id]+10.0),2.0))  +0.02);
    
    dd = dss-(dss-d[id])*exp(-dt/taudx); 
    df = fss-(fss-f[id])*exp(-dt/tauf);
    
    fca = 1.0;
    if (cs[id] > 0) {
        fca = 1.0/(1.0 + cs[id]/0.6);
    }
    
    po = d[id]*f[id]*fca;
    
    return po;
}
#else
template <int ncells>
double UCLACell<ncells>::comp_icalpo(int id, double& dc1, double& dc2, double& dxi1ca, double& dxi2ca, double& dxi1ba, double& dxi2ba) 
{
    double vth, s6, taupo, poinf, alpha, beta, r1, r2, cat, fca, s1, s1t, k1, k2, k1t, k2t, s2, s2t;
    double vx, sx, poi, tau3, k3, k3t, vy, sy, Pr, recov, tca, cpt, tau_ca, tauca, tauba, vyr, syr, Ps;
    double k6, k5, k6t, k5t, k4, k4t, po; 

    vth=0.0;
    s6 = 8.0;
    taupo=1.0;
    poinf=1.0/(1.0+exp(-(v[id] - vth)/s6));
    
    alpha = poinf/taupo;
    beta = (1.0-poinf)/taupo;
    
    r1 = 0.30;
    r2 = 3.0;
    cat = 3.0;
    fca = 1.0/(1.0 + (cat/cp[id])*(cat/cp[id])*(cat/cp[id]));
    s1 = 0.0182688*fca;
    s1t = 0.00195;
    k1 = 0.024168*fca;
    k2 = 1.03615e-4;
    k1t = 0.00413;
    k2t = 0.00224;
    s2 = s1*(r1/r2)*(k2/k1);
    s2t = s1t*(r1/r2)*(k2t/k1t);
    
    vx = -40.0;
    sx = 3.0;
    poi = 1.0/(1.0 + exp(-(v[id] - vx)/sx));
    tau3 = 3.0;
    
    k3 = (1.0-poi)/tau3;
    k3t = k3;
    
    vy = -40.0;
    sy = 4.0;
    Pr = 1.0-1.0/(1.0 + exp(-(v[id] - vy)/sy));
    
    recov = 10.0+4954.0*exp(v[id]/15.6);
    
    tca = 78.0329;
    cpt = 6.09365;
    tau_ca = tca/(1.0 + (cp[id]/cpt)*(cp[id]/cpt)*(cp[id]/cpt)*(cp[id]/cpt)) + 0.1;
    
    tauca = (recov-tau_ca)*Pr + tau_ca;
    tauba = (recov-450.0)*Pr + 450.0;
    
    vyr = -40.0;
    syr = 11.32;
    Ps = 1.0/(1.0 + exp(-(v[id]-vyr)/syr));
    
    k6 = fca*Ps/tauca;
    k5 = (1.0-Ps)/tauca;
    
    k6t = Ps/tauba;
    k5t = (1.0-Ps)/tauba;
    
    k4=k3*(alpha/beta)*(k1/k2)*(k5/k6);
    k4t=k3t*(alpha/beta)*(k1t/k2t)*(k5t/k6t);
    
    po=1.0-xi1ca[id]-xi2ca[id]-xi1ba[id]-xi2ba[id]-c1[id]-c2[id];
    
    dc2 = beta*c1[id] + k5*xi2ca[id] + k5t*xi2ba[id] - (k6 + k6t + alpha)*c2[id];
    dc1 = alpha*c2[id] + k2*xi1ca[id] + k2t*xi1ba[id] + r2*po - (beta + r1 + k1t + k1)*c1[id];
    dxi1ca = k1*c1[id] + k4*xi2ca[id] + s1*po - (k3 + k2 + s2)*xi1ca[id];
    dxi2ca = k3*xi1ca[id] + k6*c2[id] - (k5 + k4)*xi2ca[id];
    dxi1ba = k1t*c1[id] + k4t*xi2ba[id] + s1t*po - (k3t + k2t + s2t)*xi1ba[id];
    dxi2ba = k3t*xi1ba[id] + k6t*c2[id] - (k5t + k4t)*xi2ba[id];
    
    return po;
}
#endif

template <int ncells>
double UCLACell<ncells>::comp_iup(int id) 
{
    return vup*ci[id]*ci[id]/(ci[id]*ci[id] + cup*cup);
}

template <int ncells>
double UCLACell<ncells>::comp_ileak(int id) 
{
    return gleak*(cj[id]*cj[id]/(cj[id]*cj[id] + kj*kj))*(cj[id]*16.667-ci[id]); //vsr/vcell = 0.06
}

#define bcal 24.0
#define xkcal 7.0
#define srmax 47.0
#define srkd 0.6
#define bmem 15.0
#define kmem 0.3
#define bsar 42.0
#define ksar 13.0

template <int ncells>
double UCLACell<ncells>::comp_inst_buffer(double c) 
{
    double bpx, spx, mempx, sarpx;
    //bcal = 24.0;
    //xkcal = 7.0;
    //srmax = 47.0;
    //srkd = 0.6;
    //bmem = 15.0;
    //kmem = 0.3;
    //bsar = 42.0;
    //ksar = 13.0;
    bpx=bcal*xkcal/((xkcal+c)*(xkcal+c));
    spx=srmax*srkd/((srkd+c)*(srkd+c));
    mempx=bmem*kmem/((kmem+c)*(kmem+c));
    sarpx=bsar*ksar/((ksar+c)*(ksar+c));
    
    return 1.0/(1.0+bpx+spx+mempx+sarpx);
}

template <int ncells>
double UCLACell<ncells>::comp_Q(int id) 
{
    double bv, Qr;
    bv = (cstar - 50.0) - av*cstar;
    if (cjp[id] > 50.0 && cjp[id] < cstar) {
        Qr = cjp[id] - 50.0;
    }
    else {
        Qr = av*cjp[id] + bv;
    }
    return cj[id]*Qr/cstar;
}

#define ay 0.05
#define gryr 2.58079

template <int ncells>
double UCLACell<ncells>::comp_dir(int id, double po, double Qr, double rxa, double dcj) 
{
    double sparkV=exp(-ay*(v[id]+30))/(1.0+exp(-ay*(v[id]+30)));
    double spark_rate=gryr*po*fabs(rxa)*sparkV;
    
    return spark_rate*Qr - xir[id]*(1.0 - taur*dcj/cj[id])/taur;
}

#define gbarsr 26841.8
#define ax 0.3576
#define gdyad 9000.0
#define taups 0.5

template <int ncells>
double UCLACell<ncells>::comp_dcp(int id, double po, double Qr, double rxa) 
{
    //const double gbarsr=26841.8;// m mol/(cm C) 
    //const double ax=0.3576;
    //const double gdyad=9000.0;// m mol/(cm C) 
    double gsr = gbarsr*exp(-ax*(v[id] + 30.0))/(1.0 + exp(-ax*(v[id] + 30.0)));
    double xirp = po*Qr*fabs(rxa)*gsr;

    double xicap = po*gdyad*fabs(rxa);
    //const double taups = 0.5;
    return xirp + xicap - (cp[id] - cs[id])/taups;
}

#define gsk 0.005

template <int ncells>
double UCLACell<ncells>::comp_isk (int id)
{
    double ek, z, isk;
    double rat = pow(skh[id]/ci[id],skn[id]);
    //double rat = pow(skh[id]/cs[id],skn[id]);
    z = 1.0/(1.0 + rat);
    ek = (1.0/frt)*log(xko/xki);
    
    isk = iskfac[id]*gsk*z*(v[id] - ek);
    
    return isk;
}

template <int ncells>
void UCLACell<ncells>::saveconditions(FILE* file, int id, bool header, double t) {
    if (header) {
        fprintf(file,"t\tv\txm\txh\txj\txr\txs1\txs2\txtos\tytos\txir\tc1\tc2\txi1ca\txi1ba\txi2ca\txi2ba\txtof\tytof\ttropi\ttrops\tcj\tcjp\tci\tcs\tcp\tnai\n");
    }
    fprintf(file,"%g\t",t);
    fprintf(file,"%.12f\t",v[id]); //2
    fprintf(file,"%.12f\t",xm[id]); //3
    fprintf(file,"%.12f\t",xh[id]); //4
    fprintf(file,"%.12f\t",xj[id]); //5
    fprintf(file,"%.12f\t",xr[id]); //6
    fprintf(file,"%.12f\t",xs1[id]); //7
    fprintf(file,"%.12f\t",xs2[id]); //8
    fprintf(file,"%.12f\t",xtos[id]); //9
    fprintf(file,"%.12f\t",ytos[id]); //10
    fprintf(file,"%.12f\t",xir[id]); //11
    #ifdef xiaodong
    fprintf(file,"%.12f\t",d[id]); //12
    fprintf(file,"%.12f\t",f[id]); //13
    #else
    fprintf(file,"%.12f\t",c1[id]); //12
    fprintf(file,"%.12f\t",c2[id]); //13 
    fprintf(file,"%.12f\t",xi1ca[id]); //14
    fprintf(file,"%.12f\t",xi1ba[id]); //15
    fprintf(file,"%.12f\t",xi2ca[id]); //16
    fprintf(file,"%.12f\t",xi2ba[id]); //17
    #endif
    fprintf(file,"%.12f\t",xtof[id]); //18
    fprintf(file,"%.12f\t",ytof[id]); //19
    fprintf(file,"%.12f\t",tropi[id]); //20
    fprintf(file,"%.12f\t",trops[id]); //21
    fprintf(file,"%.12f\t",cj[id]); //22
    fprintf(file,"%.12f\t",cjp[id]); //23
    fprintf(file,"%.12f\t",ci[id]); //24
    fprintf(file,"%.12f\t",cs[id]); //25
    fprintf(file,"%.12f\t",cp[id]); //26
    fprintf(file,"%.12f\n",nai[id]); //27
}

template <int ncells>
void UCLACell<ncells>::setcell (int id, UCLACell<1>* newcell)
{
    xm[id] = newcell->xm[0];
    xh[id] = newcell->xh[0];
    xj[id] = newcell->xj[0];
    xr[id] = newcell->xr[0];
    xs1[id] = newcell->xs1[0];
    xs2[id] = newcell->xs2[0];
    xtos[id] = newcell->xtos[0];
    ytos[id] = newcell->ytos[0];
    v[id] = newcell->v[0];
    cp[id] = newcell->cp[0];
    cs[id] = newcell->cs[0];
    ci[id] = newcell->ci[0];
    cj[id] = newcell->cj[0];
    cjp[id] = newcell->cjp[0];
    xir[id] = newcell->xir[0];
    #ifdef xiaodong
    d[id] = newcell->d[0];
    f[id] = newcell->f[0];
    #else
    c1[id] = newcell->c1[0];
    c2[id] = newcell->c2[0];
    xi1ca[id] = newcell->xi1ca[0];
    xi1ba[id] = newcell->xi1ba[0];
    xi2ca[id] = newcell->xi2ca[0];
    xi2ba[id] = newcell->xi2ba[0];
    #endif
    nai[id] = newcell->nai[0];
    xtof[id] = newcell->xtof[0];
    ytof[id] = newcell->ytof[0];
    tropi[id] = newcell->tropi[0];
    trops[id] = newcell->trops[0];
    vold[id] = newcell->vold[0];
    jparam[id] = newcell->jparam[0];
    diffcurrent[id] = newcell->diffcurrent[0];
    
    inafac[id] = newcell->inafac[0]; //
    itoslowfac[id] = newcell->itoslowfac[0]; //
    itofac[id] = newcell->itofac[0];
    ikrfac[id] = newcell->ikrfac[0];
    iksfac[id] = newcell->iksfac[0];
    nacafac[id] = newcell->nacafac[0];
    icalfac[id] = newcell->icalfac[0]; //

    iskfac[id] = newcell->iskfac[0];
    skh[id] = newcell->skh[0];
    skn[id] = newcell->skn[0];
    
    naiclamped[id] = newcell->naiclamped[0]; //need 42
}

template <int ncells>
void UCLACell<ncells>::getcell (int id, UCLACell<1>* newcell)
{
    newcell->xm[0] = xm[id];
    newcell->xh[0] = xh[id];
    newcell->xj[0] = xj[id];
    newcell->xr[0] = xr[id];
    newcell->xs1[0] = xs1[id];
    newcell->xs2[0] = xs2[id];
    newcell->xtos[0] = xtos[id];
    newcell->ytos[0] = ytos[id];
    newcell->v[0] = v[id];
    newcell->cp[0] = cp[id];
    newcell->cs[0] = cs[id];
    newcell->ci[0] = ci[id];
    newcell->cj[0] = cj[id];
    newcell->cjp[0] = cjp[id];
    newcell->xir[0] = xir[id];
    #ifdef xiaodong
    newcell->d[0] = d[id];
    newcell->f[0] = f[id];
    #else
    newcell->c1[0] = c1[id];
    newcell->c2[0] = c2[id];
    newcell->xi1ca[0] = xi1ca[id];
    newcell->xi1ba[0] = xi1ba[id];
    newcell->xi2ca[0] = xi2ca[id];
    newcell->xi2ba[0] = xi2ba[id];
    #endif
    newcell->nai[0] = nai[id];
    newcell->xtof[0] = xtof[id];
    newcell->ytof[0] = ytof[id];
    newcell->tropi[0] = tropi[id];
    newcell->trops[0] = trops[id];
    newcell->vold[0] = vold[id];
    newcell->jparam[0] = jparam[id];
    newcell->diffcurrent[0] = diffcurrent[id];
    
    newcell->inafac[0] = inafac[id];
    newcell->itofac[0] = itofac[id];
    newcell->itoslowfac[0] = itoslowfac[id];
    newcell->ikrfac[0] = ikrfac[id];
    newcell->iksfac[0] = iksfac[id];
    newcell->nacafac[0] = nacafac[id];
    newcell->icalfac[0] = icalfac[id];

    newcell->iskfac[0] = iskfac[id];
    newcell->skh[0] = skh[id];
    newcell->skn[0] = skn[id];
    
    newcell->naiclamped[0] = naiclamped[id];
}



#endif // UCLACell_cpp
