//
//  HundCell.cpp
//
//  Implementation of the Hund-Rudy Model
//  Created by Julian Landaw on 2/26/18.
//  Copyright © 2018 Julian Landaw. All rights reserved.
//
#include <math.h>
#include <fstream>
#include "HundCell.h"

#ifndef HundCell_cpp
#define HundCell_cpp

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
HundCell<ncells>::HundCell()
{
    #define     zna 		1.0
    #define     zk 		1.0
    #define     zcl  		-1.0
    #define     zca      	2.0
    #define     ganai     	0.75
    #define     ganao           0.75
    #define     gaki        	0.75
    #define     gako          	0.75
    #define     gacao         	0.341
    #define     gacai        	1.0
    #define     csqnbar        	10.0
    #define     kmcsqn     	0.8
    #define     kmtrpn    	0.5e-3
    #define     kmcmdn      	2.38e-3
    #define     trpnbar 	70e-3
    #define     cmdnbar 	50e-3
    #define     alphacamk  	0.05
    #define     betacamk   	0.00068
    #define     bsrbar    	0.047
    #define     kmbsr   	0.00087
    #define     bslbar  	1.124
    #define     kmbsl     	0.0087
    #define     sstau    	0.2
    #define     temp     	310.0
    #define     frdy      	96485.0
    #define     R         	8314.0
    #define     pi         	3.14
    #define     radius       	0.0011
    #define     length     	0.01
    #define     rcg         	2.0
    #define     vcell     	(1000.0*pi*radius*radius*length)
    #define     ageo      	(2.0*pi*radius*radius + 2.0*pi*radius*length)
    #define     acap       	(rcg*ageo)  
    #define     vmyo         	(vcell * 0.68)
    #define     vnsr       	(vcell * 0.0552)
    #define     vjsr        	(vcell * 0.0048)
    #define     vss         	(vcell * 0.02)
    #define     camk0      	0.05
    #define     kmcam       	0.0015
    #define     kmcamk      	0.15
    #define     nao       	140
    #define     ko             	5.4
    #define     cao           	1.8
    #define     clo            	100.0

    for (int k = 0; k < ncells; k++) {
        v[k]                  =       -86.9351482;
        m[k]               =       0.00111859;
        h[k]               =       0.990199813;
        j[k]               =       0.993630289;
        d[k]               =       1.38e-06;
        f[k]               =       0.999967371;
        f2[k]              =       0.999964108;
        fca[k]             =       0.985118263;
        fca2[k]            =       1.0;
        xs1[k]             =       0.018988645    ;
        xs2[k]             =       0.018988645    ;
        xr[k]              =       1.40e-08;
        cli[k]             =       18.90029866;
        camktrap[k]        =       0.000632046;
        nai[k]             =       9.875593321;
        ki[k]              =       142.0142622;
        cai[k]             =       8.34e-05;
        cansr[k]           =       1.271695903;
        cajsr[k]           =       1.27169598;
        car[k]             =       8.34e-05;
        a[k]               =       0.01303108;
        i[k]               =       0.999972058;
        i2[k]              =       0.999813749;
        aa[k]              =       0.000555155;
        ml[k]              =       0.00111859;
        hl[k]              =       0.339310414;
        ri[k]              =       0.999950186;
        ro[k]              =       0.0;
        dpow[k]          =       8.987394846;
        camkactive[k]      =       0.003264868;
        icasave[k]          =       0;

	diffcurrent[k] = 0.0;
	itofac[k] = 1.0;
	iksfac[k] = 1.0;
	ikrfac[k] = 1.0;
	nacafac[k] = 1.0;
	nakfac[k] = 1.0;
	inafac[k] = 1.0;

	naiclamp[k] = false;
	kiclamp[k] = false;
    }
    
}

template <int ncells>
bool HundCell<ncells>::iterate(const int id, double dt, double st, double dv_max)
{    
	double dm, dh, dj, ina;
	double dml, dhl, inal;    
	double dd, df, df2, dfca, dfca2, ddpow, ica;
	double ik1;
	double dxs1, dxs2, iks;
	double ikp;
	double icab;
	double iclb;
	double dxr, ikr;
	double inaca;
	double inak;
	double ipca;
	double da, di, di2, ito1;
	double daa, ito2;
	double ctkcl;
	double ctnacl;
	double icatot, iktot, inatot, icltot, itot, dv;
	double qleak, qup, qrel, dro, dri, qtr;

	comp_ina (id, dt, dm, dh, dj, ina);
	comp_inal (id, dt, dml, dhl, inal);
	comp_ical (id, dt, dd, df, df2, dfca, dfca2, ddpow, ica);
	comp_ik1(id, ik1);
	comp_iks(id, dt, dxs1, dxs2, iks);
	comp_ikp(id, ikp);
	comp_icab(id, icab);
	comp_iclb(id, iclb);
	comp_ikr(id, dt, dxr, ikr);
    	comp_inaca(id, inaca);
	comp_inak(id, inak);
	comp_ipca(id, ipca);
	comp_ito1(id, dt, da, di, di2, ito1);
	comp_ito2(id, dt, daa, ito2);
	comp_kcl(id, ctkcl);
	comp_nacl(id, ctnacl);

 	icatot          = ica + icab + ipca - 2.0*inaca;
    	iktot           = ikr + iks + ikp + ik1 - 2.0*inak + ito1 + 0.5*(st - diffcurrent[id]);
    	inatot          = 3.0*inak + ina + 3.0*inaca + inal;
    	icltot          = ito2 + iclb + 0.5*(st - diffcurrent[id]);
    	itot            = icatot  + iktot + inatot + icltot;
    
    	dv = -itot*dt; //since current is oriented outward from the cell, by convention
    
    	if (dv_max > 0 && dv*dv > dv_max*dv_max) {return false;}

	comp_qleak(id, qleak);
	comp_qup(id, qup);
	comp_qrel(id, dt, dro, dri, qrel);
	comp_qtr(id, qtr);

// Calculate concentration changes here
	
	double qdiff           =       (car[id]-cai[id])/sstau;
        double dcar            =       dt*(-(ica)*acap/(vss*2.0*frdy)+qrel*vjsr/vss-qdiff);
 
    	double bsr      =       bsrbar*(car[id]/(car[id]+kmbsr));
    	double bsl      =       bslbar*(car[id]/(car[id]+kmbsl));
 
    	double cartot   =       car[id]       +       bsr     +       bsl     +       dcar;
 
    	double b1      =       bsrbar   +       bslbar   -       cartot  +       kmbsr    +       kmbsl;
    	double c1      =       (kmbsr*kmbsl)-(cartot*(kmbsr+kmbsl))
                                        +bsrbar*kmbsl+bslbar*kmbsr;
    	double d1      =       -1.0*kmbsr*kmbsl*cartot;
 
    	car[id]       =       (2.0/3.0)*sqrt(b1*b1-3.0*c1)*
                                cos(acos((9.0*b1*c1-2*b1*b1*b1-27*d1)/(2.0*pow((b1*b1-3.0*c1),1.5)))/3.0)-b1/3.0;
 
 
    //    JSR
    	double dcajsr          =       dt*(qtr-qrel);
    	double cajsrtot        =       dcajsr / (1.0+csqnbar*kmcsqn/pow(kmcsqn+cajsr[id],2));
    	cajsr[id]     =       cajsr[id] + cajsrtot;
    //    NSR
    	double dcansr            =       dt * (qup - qleak - qtr * vjsr/vnsr);
    	cansr[id]        =       cansr[id]     + dcansr;
    //    NAI,KI,CLI
    	double dnai              =       dt*((-inatot*acap)/(vmyo*zna*frdy) + ctnacl);
    	double dki              =       dt*((-iktot *acap)/(vmyo*zk*frdy)  + ctkcl);
    	double dcli              =       dt*((-icltot*acap)/(vmyo*zcl*frdy) + ctnacl + ctkcl);     
    	if (kiclamp[id] == false) {
		ki[id]          =       dki     + ki[id];
	}
	if (naiclamp[id] == false) {
    		nai[id]          =       dnai    + nai[id];
	}
    	cli[id]          =       dcli    + cli[id];
   
    //          NOTE:  USING ANALTYICAL SOLUTION FOR MYOPLASMIC CALCIUM BUFFERING AS OF 9/29/2005     
    	double dcai        =           dt*-((icab + ipca -2.0*inaca )*acap/
                        (zca*frdy*vmyo)+(qup-qleak)*vnsr/vmyo-qdiff*vss/vmyo);
 
    	double trpn            =       trpnbar*(cai[id]/(cai[id] + kmtrpn));
    	double cmdn            =       cmdnbar*(cai[id]/(cai[id] + kmcmdn));
    	double catotal         =       trpn + cmdn + dcai + cai[id];
    	double bmyo            =       cmdnbar + trpnbar - catotal + kmtrpn + kmcmdn;
    	double cmyo            =       kmcmdn*kmtrpn - catotal*(kmtrpn+kmcmdn)
                                        + (trpnbar*kmcmdn) + cmdnbar*kmtrpn;
    	double dmyo            =       -kmtrpn*kmcmdn*catotal;
 
 
    	cai[id]
                        =       (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)
                                *cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;    //    CAMKINASE
 
    	double camkbound               =       camk0*(1.0-camktrap[id])*1.0/(1.0 + (kmcam/car[id]));
    	camktrap[id]          =       dt*(alphacamk*camkbound*(camkbound+camktrap[id])-betacamk*camktrap[id]) + camktrap[id];
 
    	camkactive[id]        =       camkbound       +       camktrap[id];
// Concentration updates complete

	icasave[id] = ica;

	v[id] += dv;
	m[id] += dm*dt;
	h[id] += dh*dt;
	j[id] += dj*dt;
	ml[id] += dml*dt;
	hl[id] += dhl*dt;
	d[id] += dd*dt;
	f[id] += df*dt;
	f2[id] += df2*dt;
	fca[id] += dfca*dt;
	fca2[id] += dfca2*dt;
	dpow[id] += ddpow*dt;
	xs1[id] += dxs1*dt;
	xs2[id] += dxs2*dt;
	xr[id] += dxr*dt;
	a[id] += da*dt;
	i[id] += di*dt;
	i2[id] += di2*dt;
	aa[id] += daa*dt;
    
    return true;
}

template <int ncells>
void HundCell<ncells>::stepdt (int id, double dt, double st) {
    if (id > -1 && id < ncells) {
        bool success = iterate(id, dt, st, DV_MAX);
        if (!success) {
            for (int i = 0; i < ADAPTIVE; i++) {
                iterate(id, dt/ADAPTIVE, st, -1);
            }
        }
    }
}

#define MaxGNa 8.25

template <int ncells>
void HundCell<ncells>::comp_ina (int id, double dt, double& dm, double& dh, double& dj, double &ina) 
{ 
    	double bj1a,bj2a,bj2b,bj,aj1a,aj1b,aj1c,aj, ah, bh, am, bm;
    	double ms,tm,hs,th,js,tj;
    	//double camfact;
	//double KMCAM=0.06;  
    	//double vshift; 
	//double deltaalpha = -.18;
    	//double MaxGNa, ENa;
	double ENa;

	//MaxGNa = 8.25;

    	//camfact= 1.0/(1.0 + (KMCAM/caM[id])*(KMCAM/caM[id])*(KMCAM/caM[id])*(KMCAM/caM[id]));
    	//vshift = -3.25*camfact;

    	ENa=(R*temp/frdy)*log(nao/nai[id]);
	
    	am = 0.32*(v[id]+47.13)/(1.0-exp(-0.1*(v[id]+47.13)));  
    	bm = 0.08*exp(-v[id]/11.0);				 	

	if((v[id])>=-40.0)
	{
	        ah=0.0;	
		bh=1.0/(0.13*(1.0 + exp((v[id] + 10.66)/(-11.1))));  	
	}
	else
	{
	        ah = 0.135*exp((80.0+v[id])/(-6.8));
		bh = 3.56*exp(.079*(v[id]))+(3.1e5)*exp(0.35*(v[id]));
        }

	if((v[id])>=-40.0)
	{
		aj=0.0;	

		bj1a = exp((-2.535e-7)*(v[id]));	
		bj = (0.3)*bj1a/(exp(-0.1*(v[id]+32.0))+1.0);  
	}
	
	else
	{				 		
		aj1a=(-1.2714e5)*exp(0.2444*(v[id]));
		aj1b=(3.474e-5)*exp(-0.04391*(v[id]));
		aj1c=(v[id]+37.78)/(1+exp(.311*(v[id]+79.23)));
		
		//aj=(1+camfact*deltaalpha)*(aj1a-aj1b)*aj1c;
		aj = (aj1a - aj1b)*aj1c;    

		bj2a = 0.1212*exp(-0.01052*(v[id]));
		bj2b = 1.0+exp(0-.1378*(v[id] + 40.14));
		bj = bj2a/bj2b;
	}
	
	ms=am/(am+bm);	
	tm=1.0/(am+bm);
	dm = (ms-(ms-m[id])*exp(-dt/tm) - m[id])/dt;
	
	hs=ah/(ah+bh);
	th=1.0/(ah+bh);
	dh = (hs-(hs-h[id])*exp(-dt/th) - h[id])/dt;
	
	js=aj/(aj+bj);
	tj=1.0/(aj+bj);
	dj = (js-(js-j[id])*exp(-dt/tj) - j[id])/dt;

	ina = inafac[id]*(MaxGNa*m[id]*m[id]*m[id]*h[id]*j[id])*(v[id]-ENa); 
}

template <int ncells>
void HundCell<ncells>::comp_inal (int id, double dt, double& dml, double& dhl, double &inal) {

	double ENa=(R*temp/frdy)*log(nao/nai[id]);
	
	double aml = 0.32*(v[id]+47.13)/(1.0-exp(-0.1*(v[id]+47.13)));  
	double bml = 0.08*exp(-v[id]/11.0);				 
	
	double hlinf = 1.0/(1.0 + exp((v[id]+91)/6.1));  

	double ms = aml/(aml+bml);	
	double tml = 1.0/(aml+bml);
	dml = (ms-(ms - ml[id])*exp(-dt/tml) - ml[id])/dt;
	
	double thl = 600.0;
	dhl = (hlinf - (hlinf-hl[id])*exp(-dt/thl) - hl[id])/dt;
		
	inal = inafac[id]*(0.0065)*ml[id]*ml[id]*ml[id]*hl[id]*(v[id]-ENa); 
}

#define pca 2.43e-4
#define fca_dtaucamkbar 10.0

template <int ncells>
void HundCell<ncells>::comp_ical (int id, double dt, double& dd, double& df, double& df2, double& dfca, double& dfca2, double& ddpow, double& ica) 
{ 
        double ibarca                  =       pca*zca*zca*(((v[id]-15)*frdy*frdy)/(R*temp))
                                *((gacai*car[id]*exp((zca*(v[id]-15.0)*frdy)/
                                (R*temp))-gacao*cao)/
                                (exp((zca*(v[id]-15.0)*frdy)/(R*temp))-1.0));
        double dss             = (1.0/(1.0+exp(-(v[id]-4.0)/6.74)));                         //  activation
        double dtau            = (0.59+.8*exp(0.052*(v[id]+13))/(1.0+exp(0.132*(v[id]+13))));
        double fss             = 0.30 + 0.7/(1.0 + exp((v[id]+17.12)/7.0));        //  fast voltage dependent inactivation
        double ftau            = 1.0/(0.2411*exp( -(0.045*(v[id]-9.6914))*(0.045*(v[id]-9.6914 )))+0.0529);
        double f2ss            = 0.23 + 0.77/(1.0 + exp((v[id]+17.12)/7.0));
        double f2tau           = 1.0/(0.0423*exp( -(0.059*(v[id]-18.5726))*(0.059*(v[id]-18.5726)))+0.0054);
        double fcass           = 0.3/(1.0 - icasave[id]/0.05) + 0.55/(1.0 + car[id]/0.003) + 0.15;
        double fca2ss          = 1.0/(1.0 - icasave[id]/0.01);
        double fca_dtaucamk    = fca_dtaucamkbar*camkactive[id]/(kmcamk + camkactive[id]);
        double fcatau          = fca_dtaucamk  + 0.5 + 1.0/(1.0 + car[id]/0.003);
        double fca2tau         = 300.0/( 1.0 + exp((-icasave[id] - 0.175)/0.04)) + 125.0;
        double powss           = 9.0 - 8.0/(1.0+exp(-(v[id]+65)/3.4));
        double powtau          = 10;
        dd                = (dss-(dss-d[id])*exp(-dt/dtau) - d[id])/dt;
        df                = (fss-(fss-f[id])*exp(-dt/ftau) - f[id])/dt;
        df2               = (f2ss-(f2ss-f2[id])*exp(-dt/f2tau) - f2[id])/dt;
        dfca              = (fcass-(fcass-fca[id])*exp(-dt/fcatau) - fca[id])/dt;
        dfca2             = (fca2ss-(fca2ss-fca2[id])*exp(-dt/fca2tau) - fca2[id])/dt;
        ddpow           = (powss-(powss-dpow[id])*exp(-dt/powtau) - dpow[id])/dt;
        ica   =       icalfac[id]*pow(d[id],dpow[id])*f[id]*f2[id]*fca[id]*fca2[id]*ibarca;
}

#define pcab (0.1225*0.003016*0.00054)

template <int ncells>
void HundCell<ncells>::comp_icab(int id, double& icab) {
	icab    = pcab*zca*zca*((v[id]*frdy*frdy)/(R*temp))
                *((gacai*cai[id]*exp((zca*v[id]*frdy)/
                  (R*temp))-gacao*cao)
                  /(exp((zca*v[id]*frdy)/(R*temp))-1.0));
}

#define gbarclb 0.000225

template <int ncells>
void HundCell<ncells>::comp_iclb(int id, double& iclb) {
	double ecl     =     -(R*temp/frdy)*log(clo/cli[id]);
	iclb    =       gbarclb*(v[id] - ecl);
}

#define gbark1 0.5

template <int ncells>
void HundCell<ncells>::comp_ik1(int id,double& ik1) {
	double ek1     = ((R*temp)/frdy)*log(ko/ki[id]);
        double k1a     = 1.02/(1.0+exp(0.2385*(v[id] - ek1-59.215)));
        double k1b     = (0.49124*exp(0.08032*(v[id] - ek1+5.476))+exp(0.06175*(v[id] - ek1-594.31)))/(1.0+exp(-0.5143*(v[id] - ek1 + 4.753)));
        double k1ss    = k1a/(k1a+k1b);
      double gk1     = gbark1*sqrt(ko/5.4);
      ik1     = gk1*k1ss*(v[id] - ek1);
}

#define gkp 0.00276

template <int ncells>
void HundCell<ncells>::comp_ikp(int id, double& ikp) {
	double ekp     =       ((R*temp)/frdy)*log(ko/ki[id]);
        double kp      =       1.0/(1.0+exp((7.488-v[id])/5.98));
        ikp     =       gkp*kp*(v[id] - ekp);
}

#define gbarkr 0.0138542

template <int ncells>
void HundCell<ncells>::comp_ikr(int id, double dt, double& dxr, double& ikr) {
	double gkr     =       gbarkr*sqrt(ko/5.4);
        double ekr     =       ((R*temp)/frdy)*log(ko/ki[id]);
        double xrss    =       1.0/(1.0+exp(-(v[id] + 10.085)/4.25));
        double xrtau   =       1.0/(0.0006*(v[id] - 1.7384)/(1.0-exp(-0.136*(v[id] - 1.7384)))+ 0.0003*(v[id] + 38.3608)/(exp(0.1522*(v[id] + 38.3608))-1.0));
        double rkr       =       1.0/(1.0+exp((v[id] + 10)/15.4));
        dxr           =       (xrss-(xrss-xr[id])*exp(-dt/xrtau) - xr[id])/dt;
        ikr             =       ikrfac[id]*gkr*xr[id]*rkr*(v[id]-ekr);
}

#define prnak 0.01833

template <int ncells>
void HundCell<ncells>::comp_iks(int id, double dt, double& dxs1, double& dxs2, double& iks) {
	double eks     =       ((R*temp)/frdy)*log((ko+prnak*nao)/(ki[id]+prnak*nai[id]));
        double gksbar  = 0.0575*(0.433*(1.0+0.6/(1.0+pow((0.000038/cai[id]),1.4))));
        double xsss    = 1.0/(1.0+exp(-(v[id]-10.5)/24.7));
        double xs1tau  = 1.0/(0.0000761*(v[id]+44.6)/(1.0-exp(-9.97*(v[id]+44.6))) + 0.00036 *(v[id] - 0.55)/(exp(0.128*(v[id] - 0.55))-1.0));
        double xs2tau  = 2*xs1tau;
        dxs1 = (xsss-(xsss-xs1[id])*exp(-dt/xs1tau) - xs1[id])/dt;
        dxs2 = (xsss-(xsss-xs2[id])*exp(-dt/xs2tau) - xs2[id])/dt;
        iks     = iksfac[id]*gksbar * xs1[id] * xs2[id] *(v[id] - eks);
}

#define vmax 4.5
#define kmcaact 0.000125
#define kmnai_inaca 12.3
#define kmnao 87.5 
#define kmcai 0.0036
#define kmcao 1.3
#define nu 0.35
#define ksat 0.27

template <int ncells>
void HundCell<ncells>::comp_inaca(int id, double& inaca) {
  double allo            = 1.0 / (1.0 + (kmcaact/(1.5*cai[id]))*(kmcaact/ (1.5*cai[id])));
  double num             = vmax*(nai[id]*nai[id]*nai[id]*cao*exp(nu*v[id]*frdy/
                            (R*temp)) - nao*nao*nao*1.5*cai[id]*exp((nu-1)*v[id]*frdy/(R*temp)));
        double denommult       = 1.0 + ksat*exp((nu-1)*v[id]*frdy/(R*temp));
        double denomterm1      = kmcao*nai[id]*nai[id]*nai[id] + kmnao*kmnao*kmnao*1.5*cai[id] + kmnai_inaca*kmnai_inaca*kmnai_inaca*cao*(1.0+1.5*cai[id]/kmcai);
        double denomterm2      = kmcai*nao*nao*nao*(1.0 + (nai[id]/kmnai_inaca)*(nai[id]/kmnai_inaca)*(nai[id]/kmnai_inaca))+nai[id]*nai[id]*nai[id]*cao + nao*nao*nao*1.5*cai[id];
        double deltaE          = num/(denommult*(denomterm1+denomterm2));
        inaca           = nacafac[id]*allo*deltaE;
}

#define gnakbar (0.275*2.25)
#define kmnai_inak 10.0
#define kmko 1.5

template <int ncells>
void HundCell<ncells>::comp_inak(int id, double& inak) {
	double sigma   = (1.0/7.0)*(exp(nao/67.3)-1.0);
        double fnak    = 1.0/(1.0+0.1245*exp(-0.1*v[id]*frdy/(R*temp))+0.0365*sigma*exp(-v[id]*frdy/(R*temp)));
        inak    = nakfac[id]*gnakbar*fnak*  1.0/(1.0 + (kmnai_inak/nai[id])*(kmnai_inak/nai[id]))  *  (ko/(ko+kmko));
}

#define ipcabar 0.0575
#define kmpca 0.0005

template <int ncells>
void HundCell<ncells>::comp_ipca(int id, double& ipca) {
	ipca    =       ipcabar*cai[id]/(kmpca+cai[id]);
}

#define gbarto1 0.19

template <int ncells>
void HundCell<ncells>::comp_ito1(int id, double dt, double& da, double& di, double& di2, double& ito1) {
	double alphaa    =     25.0*exp((v[id]-40.0)/25.0)/(1.0 + exp((v[id]-40.0)/25.0));
        double betaa     =     25.0*exp(-(v[id]+90.0)/25.0)/(1.0 + exp(-(v[id]+90.0)/25.0));
        double alphai          =       0.03/(1.0 + exp((v[id]+60.0)/5.0));
        double betai           =       0.2*exp((v[id]+25)/5.0)/(1.0 + exp((v[id]+25)/5));
        double alphai2         =       0.00225/(1.0 + exp((v[id]+60.0)/5.0));
        double betai2          =       0.1*exp((v[id]+25.0)/5.0)/(1.0 + exp((v[id]+25.0)/5.0));
        double atau      =     1.0/(alphaa   +   betaa);
        double itau      =     1.0/(alphai   +   betai);
        double i2tau     =     1.0/(alphai2  +   betai2);
        double ass       =     alphaa    /(alphaa  + betaa);
        double iss       =     alphai    /(alphai  + betai);
        double i2ss      =     alphai2    /(alphai2  + betai2);
        da         = (ass   -(ass-a[id])    *exp(-dt/atau) - a[id])/dt;
        di         = (iss   -(iss-i[id])    *exp(-dt/itau) - i[id])/dt;
        di2        = (i2ss  -(i2ss-i2[id])  *exp(-dt/i2tau) - i2[id])/dt;
        double rto1    =       exp(v[id]/300);
        double eto1    =       ((R*temp)/frdy)*log(ko/ki[id]);
        ito1    =       itofac[id]*gbarto1 *a[id]*a[id]*a[id]*i[id]*i2[id]*rto1*(v[id]-eto1);
}

#define pcl__ 0.0000004
#define kmto2 0.1502

template <int ncells>
void HundCell<ncells>::comp_ito2(int id, double dt, double& daa, double& ito2) {
	double aass        =     1.0/(1.0 + kmto2/car[id]);
        double aatau       =     1.0;
        daa          =     (aass-(aass-aa[id])*exp(-dt/aatau) - aa[id])/dt;
        double ibarto2     =       pcl__*zcl*zcl*((v[id]*frdy*frdy)/(R*temp))*
                            (cli[id] - clo*exp(-zcl*v[id]*frdy/(R*temp)))/
                            (1.0 - exp(-zcl*v[id]*frdy/(R*temp)));
        ito2        =       itofac[id]*ibarto2*aa[id];
}

#define ctkclbar 0.0000070756

template <int ncells>
void HundCell<ncells>::comp_kcl(int id, double& ctkcl) {
	double ecl     =       -(R*temp/frdy)*log(clo/cli[id]);
        double ek      =       (R*temp/frdy)*log(ko/ki[id]);
        ctkcl   =       ctkclbar*(ek - ecl)/((ek - ecl) + 87.8251);
}

#define ctnaclbar 0.0000098443

template <int ncells>
void HundCell<ncells>::comp_nacl(int id, double& ctnacl) {
	double ena     =       (R*temp/frdy)*log(nao/nai[id]);
        double ecl     =       -(R*temp/frdy)*log(clo/cli[id]);
        ctnacl  =       ctnaclbar*(ena-ecl)*(ena-ecl)*(ena-ecl)*(ena-ecl)/((ena-ecl)*(ena-ecl)*(ena-ecl)*(ena-ecl) + 87.8251*87.8251*87.8251*87.8251);
}

#define qreldtaucamkbar 10.0

template <int ncells>
void HundCell<ncells>::comp_qrel(int id, double dt, double& dro, double& dri, double& qrel) {
	double ibarca                  =       pca*zca*zca*(((v[id]-15)*frdy*frdy)/(R*temp))
                                *((gacai*car[id]*exp((zca*(v[id]-15.0)*frdy)/
                                (R*temp))-gacao*cao)/
                                (exp((zca*(v[id]-15.0)*frdy)/(R*temp))-1.0));	

	double qreldtaucamk            =       qreldtaucamkbar*(camkactive[id])/(kmcamk + camkactive[id]);
        double rossjsr                 =       pow(cajsr[id],1.9)/
                                        ( pow(cajsr[id],1.9) + pow((49.28*car[id])/(car[id] + 0.0028),1.9) );
        double ross                    =       rossjsr*icasave[id]*icasave[id]/(icasave[id]*icasave[id] + 1.0);
        double rotau                   =       3.0;
        dro                =       (ross - (ross - ro[id])*exp(-dt/rotau) - ro[id])/dt;
        double cafac                   =       1.0/(1.0 + exp((icasave[id] + 0.05)/0.015));
        double riss                    =       1.0/(1.0 + exp((car[id] - 0.0004 + 0.002*cafac)/0.000025));
        double ritau                   =       (350.0-qreldtaucamk)/(1.0 + exp((car[id] - 0.003 + 0.003*cafac)/0.0002)) + 3.0 + qreldtaucamk;
        dri                =       (riss - (riss - ri[id])*exp(-dt/ritau) - ri[id])/dt;
        double vg                      =       1.0/(1.0 + exp((ibarca + 13.0)/5.0));
        double grelbar                 =       3000.0*vg;
        qrel                    =       grelbar*ro[id]*ri[id]*(cajsr[id] - car[id]);
}

#define qleakbar 0.004375
#define nsrbar 15.0

template <int ncells>
void HundCell<ncells>::comp_qleak(int id, double& qleak) {
	qleak   =       qleakbar*cansr[id]/nsrbar;
}

#define qupbar 0.004375
#define kmup 0.00092
#define dqupcamkbar 0.75
#define dkmplbbar 0.00017

template <int ncells>
void HundCell<ncells>::comp_qup(int id, double& qup) {
	double dkmplb  =       dkmplbbar*camkactive[id]/(kmcamk + camkactive[id]);
        double dqupcamk=       dqupcamkbar*camkactive[id]/(kmcamk + camkactive[id]);
        qup     =       (dqupcamk + 1.0)*qupbar*cai[id]/(cai[id] + kmup - dkmplb);
}

#define tautr 120.0

template <int ncells>
void HundCell<ncells>::comp_qtr(int id, double& qtr) {
	qtr             =       (cansr[id] - cajsr[id])/tautr;
}

template <int ncells>
void HundCell<ncells>::setcell (int id, HundCell<1>* newcell)
{
    	v[id] = newcell->v[0];
	m[id] = newcell->m[0];
	h[id] = newcell->h[0];
	j[id] = newcell->j[0];
	d[id] = newcell->d[0];
	f[id] = newcell->f[0];
	f2[id] = newcell->f2[0];
	fca[id] = newcell->fca[0];
	fca2[id] = newcell->fca2[0];
	xs1[id] = newcell->xs1[0];
	xs2[id] = newcell->xs2[0];
	xr[id] = newcell->xr[0];
	cli[id] = newcell->cli[0];
	camktrap[id] = newcell->camktrap[0];
	nai[id] = newcell->nai[0];
	ki[id] = newcell->ki[0];
	cai[id] = newcell->cai[0];
	cansr[id] = newcell->cansr[0];
	cajsr[id] = newcell->cajsr[0];
	car[id] = newcell->car[0];
	a[id] = newcell->a[0];
	i[id] = newcell->i[0];
	i2[id] = newcell->i2[0];
	aa[id] = newcell->aa[0];
	ml[id] = newcell->ml[0];
	hl[id] = newcell->hl[0];
	ri[id] = newcell->ri[0];
	ro[id] = newcell->ro[0];
	dpow[id] = newcell->dpow[0];
	camkactive[id] = newcell->camkactive[0];
	icasave[id] = newcell->icasave[0];
	diffcurrent[id] = newcell->diffcurrent[0];
	
	naiclamp[id] = newcell->naiclamp[0];
	kiclamp[id] = newcell->kiclamp[0];
	inafac[id] = newcell->inafac[0];
	icalfac[id] = newcell->icalfac[0];
	ikrfac[id] = newcell->ikrfac[0];
	iksfac[id] = newcell->iksfac[0];
	itofac[id] = newcell->itofac[0];
	nacafac[id] = newcell->nacafac[0];
	nakfac[id] = newcell->nakfac[0];
}

template <int ncells>
void HundCell<ncells>::getcell (int id, HundCell<1>* newcell) {
	newcell->v[0] = v[id];
 	newcell->m[0] = m[id];
	newcell->h[0] = h[id];
	newcell->j[0] = j[id];
	newcell->d[0] = d[id];
	newcell->f[0] = f[id];
	newcell->f2[0] = f2[id];
	newcell->fca[0] = fca[id];
	newcell->fca2[0] = fca2[id];
	newcell->xs1[0] = xs1[id];
	newcell->xs2[0] = xs2[id];
	newcell->xr[0] = xr[id];
	newcell->cli[0] = cli[id];
	newcell->camktrap[0] = camktrap[id];
	newcell->nai[0] = nai[id];
	newcell->ki[0] = ki[id];
	newcell->cai[0] = cai[id];
	newcell->cansr[0] = cansr[id];
	newcell->cajsr[0] = cajsr[id];
	newcell->car[0] = car[id];
	newcell->a[0] = a[id];
	newcell->i[0] = i[id];
	newcell->i2[0] = i2[id];
	newcell->aa[0] = aa[id];
	newcell->ml[0] = ml[id];
	newcell->hl[0] = hl[id];
	newcell->ri[0] = ri[id];
	newcell->ro[0] = ro[id];
	newcell->dpow[0] = dpow[id];
	newcell->camkactive[0] = camkactive[id];
	newcell->icasave[0] = icasave[id];
	newcell->diffcurrent[0] = diffcurrent[id];
	newcell->naiclamp[0] = naiclamp[id];
	newcell->kiclamp[0] = kiclamp[id];
	newcell->inafac[0] = inafac[id];
	newcell->ikrfac[0] = ikrfac[id];
	newcell->iksfac[0] = iksfac[id];
	newcell->itofac[0] = itofac[id];
	newcell->nacafac[0] = nacafac[id];
	newcell->nakfac[0] = nakfac[id];
	newcell->icalfac[0] = icalfac[id];
}


#endif //HundCell_cpp
