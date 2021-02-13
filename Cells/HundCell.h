#ifndef HundCell_h
#define HundCell_h

template <int ncells>
class HundCell
{
public:
    
    double v[ncells];
    double m[ncells];
    double h[ncells];
    double j[ncells];
    double d[ncells];
    double f[ncells];
    double f2[ncells];
    double fca[ncells];
    double fca2[ncells];
    double xs1[ncells];
    double xs2[ncells];
    double xr[ncells];
    double cli[ncells];
    double camktrap[ncells];
    double nai[ncells];
    double ki[ncells];
    double cai[ncells];
    double cansr[ncells];
    double cajsr[ncells];
    double car[ncells];
    double a[ncells];
    double i[ncells];
    double i2[ncells];
    double aa[ncells];
    double ml[ncells];
    double hl[ncells];
    double ri[ncells];
    double ro[ncells];
    double dpow[ncells];
    double camkactive[ncells];
    double icasave[ncells];
    
    double diffcurrent[ncells];

    bool naiclamp[ncells];
    bool kiclamp[ncells];
  
  double inafac[ncells];
  double ikrfac[ncells];
  double iksfac[ncells];
  double itofac[ncells];
  double icalfac[ncells];
  double nacafac[ncells];
  double nakfac[ncells];
    
    HundCell();
    
    bool iterate(const int id, double dt, double st, double dv_max);
    
    void stepdt(const int id, double dt, double st);
    
    void comp_ina (int id, double dt, double& dm, double& dh, double& dj, double &ina);
    void comp_inal (int id, double dt, double& dml, double& dhl, double &inal);
    void comp_ical (int id, double dt, double& dd, double& df, double& df2, double& dfca, double& dfca2, double& ddpow, double& ica);
    void comp_icab (int id, double& icab);
    void comp_iclb (int id, double& iclb);
    void comp_ik1 (int id, double& ik1);
    void comp_ikp (int id, double& ikp);
    void comp_ikr (int id, double dt, double& dxr, double& ikr);
    void comp_iks (int id, double dt, double& dxs1, double& dxs2, double& iks);
    void comp_inaca (int id, double& inaca);
    void comp_inak (int id, double& inak);
    void comp_ipca (int id, double& ipca);
    void comp_ito1 (int id, double dt, double& da, double& di, double& di2, double& ito1);
    void comp_ito2 (int id, double dt, double& daa, double& ito2);
    void comp_kcl (int id, double& ctkcl);
    void comp_nacl (int id, double& ctnacl);
    void comp_qrel (int id, double dt, double& dro, double& dri, double& qrel);
    void comp_qleak (int id, double& qleak);
    void comp_qup (int id, double& qup);
    void comp_qtr (int id, double& qtr);
    
    void setcell (int id, HundCell<1>* newcell);
    
    void getcell (int id, HundCell<1>* newcell);
    
    //void saveconditions(FILE* saveconditions, int id, bool header, double t);
    
    //void write(std::fstream &file);
    //void read(std::fstream &file);
};

#endif // HundCell_h
