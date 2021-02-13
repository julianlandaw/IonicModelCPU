#ifndef LR2CellIto_h
#define LR2CellIto_h

template <int ncells>
class LR2CellIto
{
public:
    
    double v[ncells];
    double nai[ncells]; // Initial Intracellular Na (mM) 
    double nao[ncells]; // Initial Extracellular Na (mM) 
    //double nabm[ncells]; // Initial Bulk Medium Na (mM) 
    double ki[ncells]; // Initial Intracellular K (mM) 
    double ko[ncells]; // Initial Extracellular K (mM) 
    //double kbm[ncells]; // Initial Bulk Medium K (mM) 
    double cai[ncells]; // Initial Intracellular Ca (mM) 
    double cao[ncells]; // Initial Extracellular Ca (mM) 
    //double cabm[ncells]; // Initial Bulk Medium Ca (mM) 
        
    double m[ncells]; 
    double h[ncells]; 
    double j[ncells]; 
    double d[ncells]; 
    double f[ncells]; 
    double xs1[ncells]; 
    double xs2[ncells]; 
    double xr[ncells]; 
    double b[ncells]; 
    double g[ncells]; 
    double zdv[ncells]; 
    double ydv[ncells];
        
        //grelbarjsrol[i] = 0; 
        //tjsrol[i] = 1000; 
    double tcicr[ncells]; 
    double jsr[ncells]; 
    double nsr[ncells]; 
    int flag[ncells]; 
    
    double dcaiontold[ncells];
    double caiontold[ncells];
        //dt = udt; 
        //utsc = 50; 
        //i=-1;
        
    double icalfac[ncells];
    double itofac[ncells];
    double tauXfac[ncells];
    double iupfac[ncells]; //SERCA
    double ikrfac[ncells];
    double iksfac[ncells];
    double ik1fac[ncells];
    double nacafac[ncells];
    double ikatpfac[ncells];

    double iskfac[ncells];
    double skh[ncells];
    double skn[ncells];
    double xsk[ncells];
    
    double tauyfac[ncells];
    
    double diffcurrent[ncells];
    
    LR2CellIto();
    
    bool iterate(const int id, double dt, double st, double dv_max);
    
    void stepdt(const int id, double dt, double st);
    
    void comp_ina (int id, double dt, double& dm, double& dh, double& dj, double &ina);
    void comp_ical (int id, double dt, double& dd, double& df, double& ilca, double& ilcana, double& ilcak);
    void comp_icat (int id, double dt, double& db, double& dg, double& icat);
    void comp_ikr (int id, double dt, double& dxr, double& ikr);
    void comp_iks (int id, double dt, double& dxs1, double& dxs2, double& iks);
    void comp_iki (int id, double& iki);
    void comp_ikp (int id, double& ikp);
    void comp_ikna (int id, double& ikna);
    void comp_ikatp (int id, double& ikatp);
    void comp_ito (int id, double dt, double& dzdv, double& dydv, double& ito);
    void comp_isk(int id, double dt, double& dxsk, double& isk);
    void comp_inaca (int id, double& inaca);
    void comp_inak (int id, double& inak);
    void comp_insca (int id, double& insna, double& insk);
    void comp_ipca (int id, double& ipca);
    void comp_icab (int id, double& icab);
    void comp_inab (int id, double& inab);
    void conc_nai (int id, double naiont, double& dnai);
    void conc_ki (int id, double kiont, double st, double& dki); 
    void calc_dyn (int id, double dt, double caiont, double& dcaiont, double& dtcicr, double& djsr, double& dnsr, double& dcai);
    void conc_cleft(int id, double dt, double st, double naiont, double kiont, double caiont, double& dnao, double& dko, double& dcao); 
    
    void setcell (int id, LR2CellIto<1>* newcell);
    
    void getcell (int id, LR2CellIto<1>* newcell);
    
    void saveconditions(FILE* saveconditions, int id, bool header, double t);
    
    void write(std::fstream &file);
    void read(std::fstream &file);
};

#endif // LR2CellIto_h
