#ifndef CELLS
#define CELLS 300
#endif
//#define slowito
#include "LR1CellIto.cpp"
//#undef slowito

extern "C"
{
    LR1CellIto<CELLS>* LR1CellIto_new() {
        return new LR1CellIto<CELLS>();
    }
    LR1CellIto<1>* LR1CellIto_newsing() {
        return new LR1CellIto<1>();
    }
    void LR1_stepdt(LR1CellIto<CELLS>* LR1Cell,int id, double dt, double st) {
        LR1Cell->stepdt(id,dt,st);
    }
    void LR1_stepdtsing(LR1CellIto<1>* LR1Cell, double dt, double st) {
        LR1Cell->stepdt(0,dt,st);
    }
    double LR1_getv(LR1CellIto<CELLS>* LR1Cell, int id) {
        if (id > -1 && id < CELLS) {  
             return LR1Cell->v[id];
        }
        else {return LR1Cell->v[0];}
    }
    double LR1_getcai(LR1CellIto<CELLS>* LR1Cell, int id) {
        if (id > -1 && id < CELLS) {  
             return LR1Cell->cai[id];
        }
        else {return LR1Cell->cai[0];}
    }
    double LR1_getvsing(LR1CellIto<1>* LR1Cell) {
        return LR1Cell->v[0];
    }
    double LR1_getcaising(LR1CellIto<1>* LR1Cell) {
        return LR1Cell->cai[0];
    }
    void LR1_setv(LR1CellIto<CELLS>* LR1Cell, int id, double v) {
        if (id > -1 && id < CELLS) {  
            LR1Cell->v[id] = v;
        }
    }
    void LR1_setvsing(LR1CellIto<1>* LR1Cell, double v) {
        LR1Cell->v[0] = v;
    }
    void LR1_setxr(LR1CellIto<CELLS>* LR1Cell, int id, double xr) {
        if (id > -1 && id < CELLS) {  
            LR1Cell->xr[id] = xr;
        }
    }
    void LR1_setxrsing(LR1CellIto<1>* LR1Cell, double xr) {
        LR1Cell->xr[0] = xr;
    }
    double LR1_getxr(LR1CellIto<CELLS>* LR1Cell, int id) {
        if (id > -1 && id < CELLS) {  
             return LR1Cell->xr[id];
        }
        else {return LR1Cell->xr[0];}
    }
    double LR1_getxrsing(LR1CellIto<1>* LR1Cell) {
        return LR1Cell->xr[0];
    }
    void LR1_setito(LR1CellIto<CELLS>* LR1Cell, int id, double itofac) {
        if (id > -1 && id < CELLS) {  
            LR1Cell->itofac[id] = itofac;
        }
    }
    void LR1_setitosing(LR1CellIto<1>* LR1Cell, double itofac) {
        LR1Cell->itofac[0] = itofac;
    }
    void LR1_settauxfac(LR1CellIto<CELLS>* LR1Cell, int id, double tauxf) {
         if (id > -1 && id < CELLS) {
             LR1Cell->tauXfac[id] = tauxf;
         }
    }
    void LR1_settauxfacsing(LR1CellIto<1>* LR1Cell, double tauxf) {
         LR1Cell->tauXfac[0] = tauxf;
    }
    void LR1_seticalfac(LR1CellIto<CELLS>* LR1Cell, int id, double ibarca) {
         if (id > -1 && id < CELLS) {
             LR1Cell->icalfac[id] = ibarca;
         }
    }
    void LR1_seticalfacsing(LR1CellIto<1>* LR1Cell, double ibarca) {
         LR1Cell->icalfac[0] = ibarca;
    }
    
    
    void LR1_copycell(LR1CellIto<CELLS>* LR1Cell, int id, LR1CellIto<1>* LR1Cell2) {
        LR1Cell->setcell(id,LR1Cell2);
    }
    void LR1_copycellsing(LR1CellIto<1>* LR1Cell, LR1CellIto<1>* LR1Cell2) {
        LR1Cell->setcell(0,LR1Cell2);
    }
    void LR1_diffuse1D(LR1CellIto<CELLS>* LR1Cell, double D, int id) {
        if (CELLS > 1) {
            if (id == 0){
                LR1Cell->diffcurrent[id] = D*(LR1Cell->v[id+1] - LR1Cell->v[id]);
            }
            else if (id < CELLS - 1) {
                LR1Cell->diffcurrent[id] = D*(LR1Cell->v[id+1] - LR1Cell->v[id]) + D*(LR1Cell->v[id-1] - LR1Cell->v[id]);
            }
            else if (id == CELLS - 1) {
                LR1Cell->diffcurrent[id] = D*(LR1Cell->v[id-1] - LR1Cell->v[id]);
            }
        }
    }
    void LR1_diffuse1Dall(LR1CellIto<CELLS>* LR1Cell, double D) {
        for (int i = 0; i < CELLS; i++) {
            LR1_diffuse1D(LR1Cell, D, i);
        }
    }
}

#define SK
#include "LR2CellIto.cpp"

extern "C"
{
    LR2CellIto<CELLS>* LR2CellIto_new() {
        return new LR2CellIto<CELLS>();
    }
    LR2CellIto<1>* LR2CellIto_newsing() {
        return new LR2CellIto<1>();
    }
    void LR2_stepdt(LR2CellIto<CELLS>* LR2Cell,int id, double dt, double st) {
        LR2Cell->stepdt(id,dt,st);
    }
    void LR2_stepdtsing(LR2CellIto<1>* LR2Cell, double dt, double st) {
        LR2Cell->stepdt(0,dt,st);
    }
    double LR2_getv(LR2CellIto<CELLS>* LR2Cell, int id) {
        if (id > -1 && id < CELLS) {  
             return LR2Cell->v[id];
        }
        else {return LR2Cell->v[0];}
    }
    double LR2_getcai(LR2CellIto<CELLS>* LR2Cell, int id) {
        if (id > -1 && id < CELLS) {  
             return LR2Cell->cai[id];
        }
        else {return LR2Cell->cai[0];}
    }
    double LR2_getvsing(LR2CellIto<1>* LR2Cell) {
        return LR2Cell->v[0];
    }
    double LR2_getcaising(LR2CellIto<1>* LR2Cell) {
        return LR2Cell->cai[0];
    }
    double LR2_getnai(LR2CellIto<CELLS>* LR2Cell, int id) {
        if (id > -1 && id < CELLS) {
            return LR2Cell->nai[id];
        }
        else {return LR2Cell->nai[0];}
    }
    double LR2_getki(LR2CellIto<CELLS>* LR2Cell, int id) {
        if (id > -1 && id < CELLS) {
            return LR2Cell->ki[id];
        }
        else {return LR2Cell->ki[0];}
    }
    double LR2_getnaising(LR2CellIto<1>* LR2Cell) {
        return LR2Cell->nai[0];
    }
    double LR2_getkising(LR2CellIto<1>* LR2Cell) {
        return LR2Cell->ki[0];
    }    
    void LR2_setv(LR2CellIto<CELLS>* LR2Cell, int id, double v) {
        if (id > -1 && id < CELLS) {  
            LR2Cell->v[id] = v;
        }
    }
    void LR2_setvsing(LR2CellIto<1>* LR2Cell, double v) {
        LR2Cell->v[0] = v;
    }
    void LR2_setxr(LR2CellIto<CELLS>* LR2Cell, int id, double xr) {
        if (id > -1 && id < CELLS) {  
            LR2Cell->xr[id] = xr;
        }
    }
    void LR2_setxrsing(LR2CellIto<1>* LR2Cell, double xr) {
        LR2Cell->xr[0] = xr;
    }
    double LR2_getxr(LR2CellIto<CELLS>* LR2Cell, int id) {
        if (id > -1 && id < CELLS) {  
             return LR2Cell->xr[id];
        }
        else {return LR2Cell->xr[0];}
    }
    double LR2_getxrsing(LR2CellIto<1>* LR2Cell) {
        return LR2Cell->xr[0];
    }
    void LR2_setito(LR2CellIto<CELLS>* LR2Cell, int id, double itofac) {
        if (id > -1 && id < CELLS) {  
            LR2Cell->itofac[id] = itofac;
        }
    }
    void LR2_setitosing(LR2CellIto<1>* LR2Cell, double itofac) {
        LR2Cell->itofac[0] = itofac;
    }
    
    #ifdef SK
    void LR2_setisk(LR2CellIto<CELLS>* LR2Cell, int id, double iskfac) {
        if (id > -1 && id < CELLS) {  
            LR2Cell->iskfac[id] = iskfac;
        }
    }
    void LR2_setisksing(LR2CellIto<1>* LR2Cell, double iskfac) {
        LR2Cell->iskfac[0] = iskfac;
    }

    void LR2_setskh(LR2CellIto<CELLS>* LR2Cell, int id, double skh)
    {
        if (id > -1 && id < CELLS) {
            LR2Cell->skh[id] = skh;
        }
    }

    void LR2_setskhsing(LR2CellIto<1>* LR2Cell, double skh)
    {
        LR2Cell->skh[0] = skh;
    }
    #endif
    
    void LR2_settauxfac(LR2CellIto<CELLS>* LR2Cell, int id, double tauxf) {
         if (id > -1 && id < CELLS) {
             LR2Cell->tauXfac[id] = tauxf;
         }
    }
    void LR2_settauxfacsing(LR2CellIto<1>* LR2Cell, double tauxf) {
         LR2Cell->tauXfac[0] = tauxf;
    }
    void LR2_seticalfac(LR2CellIto<CELLS>* LR2Cell, int id, double ibarca) {
         if (id > -1 && id < CELLS) {
             LR2Cell->icalfac[id] = ibarca;
         }
    }
    void LR2_seticalfacsing(LR2CellIto<1>* LR2Cell, double ibarca) {
         LR2Cell->icalfac[0] = ibarca;
    }
    
    void LR2_setiupfac(LR2CellIto<CELLS>* LR2Cell, int id, double iup) {
        if (id > -1 && id < CELLS) { 
            LR2Cell->iupfac[id] = iup;
        }
    }
    
    void LR2_setiupfacsing(LR2CellIto<1>* LR2Cell, double iup) {
         LR2Cell->iupfac[0] = iup;
    }
    
    void LR2_copycell(LR2CellIto<CELLS>* LR2Cell, int id, LR2CellIto<1>* LR2Cell2) {
        LR2Cell->setcell(id,LR2Cell2);
    }
    void LR2_copycellsing(LR2CellIto<1>* LR2Cell, LR2CellIto<1>* LR2Cell2) {
        LR2Cell->setcell(0,LR2Cell2);
    }
    void LR2_diffuse1D(LR2CellIto<CELLS>* LR2Cell, double D, int id) {
        if (CELLS > 1) {
            if (id == 0){
                LR2Cell->diffcurrent[id] = D*(LR2Cell->v[id+1] - LR2Cell->v[id]);
            }
            else if (id < CELLS - 1) {
                LR2Cell->diffcurrent[id] = D*(LR2Cell->v[id+1] - LR2Cell->v[id]) + D*(LR2Cell->v[id-1] - LR2Cell->v[id]);
            }
            else if (id == CELLS - 1) {
                LR2Cell->diffcurrent[id] = D*(LR2Cell->v[id-1] - LR2Cell->v[id]);
            }
        }
    }
    void LR2_diffuse1Dall(LR2CellIto<CELLS>* LR2Cell, double D) {
        for (int i = 0; i < CELLS; i++) {
            LR2_diffuse1D(LR2Cell, D, i);
        }
    }
}

#include "LR1CellIto_nsca.cpp"

extern "C"
{
    LR1CellIto_nsca<CELLS>* LR1CellIto_nsca_new() {
        return new LR1CellIto_nsca<CELLS>();
    }
    LR1CellIto_nsca<1>* LR1CellIto_nsca_newsing() {
        return new LR1CellIto_nsca<1>();
    }
    void LR1_nsca_stepdt(LR1CellIto_nsca<CELLS>* LR1Cell,int id, double dt, double st) {
        LR1Cell->stepdt(id,dt,st);
    }
    void LR1_nsca_stepdtsing(LR1CellIto_nsca<1>* LR1Cell, double dt, double st) {
        LR1Cell->stepdt(0,dt,st);
    }
    double LR1_nsca_getv(LR1CellIto_nsca<CELLS>* LR1Cell, int id) {
        if (id > -1 && id < CELLS) {  
             return LR1Cell->v[id];
        }
        else {return LR1Cell->v[0];}
    }
    double LR1_nsca_getcai(LR1CellIto_nsca<CELLS>* LR1Cell, int id) {
        if (id > -1 && id < CELLS) {  
             return LR1Cell->cai[id];
        }
        else {return LR1Cell->cai[0];}
    }
    double LR1_nsca_getvsing(LR1CellIto_nsca<1>* LR1Cell) {
        return LR1Cell->v[0];
    }
    double LR1_nsca_getcaising(LR1CellIto_nsca<1>* LR1Cell) {
        return LR1Cell->cai[0];
    }
    void LR1_nsca_setv(LR1CellIto_nsca<CELLS>* LR1Cell, int id, double v) {
        if (id > -1 && id < CELLS) {  
            LR1Cell->v[id] = v;
        }
    }
    void LR1_nsca_setvsing(LR1CellIto_nsca<1>* LR1Cell, double v) {
        LR1Cell->v[0] = v;
    }
    void LR1_nsca_setito(LR1CellIto_nsca<CELLS>* LR1Cell, int id, double itofac) {
        if (id > -1 && id < CELLS) {  
            LR1Cell->itofac[id] = itofac;
        }
    }
    void LR1_nsca_setitosing(LR1CellIto_nsca<1>* LR1Cell, double itofac) {
        LR1Cell->itofac[0] = itofac;
    }
    void LR1_nsca_setinsca(LR1CellIto_nsca<CELLS>* LR1Cell, int id, double inscafac) {
        if (id > -1 && id < CELLS) {
            LR1Cell->inscafac[id] = inscafac;
        }
    }
    void LR1_nsca_setinscasing(LR1CellIto_nsca<1>* LR1Cell, double inscafac) {
        LR1Cell->inscafac[0] = inscafac;
    }
    void LR1_nsca_setxr(LR1CellIto_nsca<CELLS>* LR1Cell, int id, double xr) {
        if (id > -1 && id < CELLS) {  
            LR1Cell->xr[id] = xr;
        }
    }
    void LR1_nsca_setxrsing(LR1CellIto_nsca<1>* LR1Cell, double xr) {
        LR1Cell->xr[0] = xr;
    }
    double LR1_nsca_getxr(LR1CellIto_nsca<CELLS>* LR1Cell, int id) {
        if (id > -1 && id < CELLS) {  
             return LR1Cell->xr[id];
        }
        else {return LR1Cell->xr[0];}
    }
    double LR1_nsca_getxrsing(LR1CellIto_nsca<1>* LR1Cell) {
        return LR1Cell->xr[0];
    }
    void LR1_nsca_settauxfac(LR1CellIto_nsca<CELLS>* LR1Cell, int id, double tauxf) {
         if (id > -1 && id < CELLS) {
             LR1Cell->tauXfac[id] = tauxf;
         }
    }
    void LR1_nsca_settauxfacsing(LR1CellIto_nsca<1>* LR1Cell, double tauxf) {
         LR1Cell->tauXfac[0] = tauxf;
    }
    void LR1_nsca_seticalfac(LR1CellIto_nsca<CELLS>* LR1Cell, int id, double ibarca) {
         if (id > -1 && id < CELLS) {
             LR1Cell->icalfac[id] = ibarca;
         }
    }
    void LR1_nsca_seticalfacsing(LR1CellIto_nsca<1>* LR1Cell, double ibarca) {
         LR1Cell->icalfac[0] = ibarca;
    }
    
    void LR1_nsca_copycell(LR1CellIto_nsca<CELLS>* LR1Cell, int id, LR1CellIto_nsca<1>* LR1Cell2) {
        LR1Cell->setcell(id,LR1Cell2);
    }
    void LR1_nsca_copycellsing(LR1CellIto_nsca<1>* LR1Cell, LR1CellIto_nsca<1>* LR1Cell2) {
        LR1Cell->setcell(0,LR1Cell2);
    }
    void LR1_nsca_diffuse1D(LR1CellIto_nsca<CELLS>* LR1Cell, double D, int id) {
        if (CELLS > 1) {
            if (id == 0){
                LR1Cell->diffcurrent[id] = D*(LR1Cell->v[id+1] - LR1Cell->v[id]);
            }
            else if (id < CELLS - 1) {
                LR1Cell->diffcurrent[id] = D*(LR1Cell->v[id+1] - LR1Cell->v[id]) + D*(LR1Cell->v[id-1] - LR1Cell->v[id]);
            }
            else if (id == CELLS - 1) {
                LR1Cell->diffcurrent[id] = D*(LR1Cell->v[id-1] - LR1Cell->v[id]);
            }
        }
    }
    void LR1_nsca_diffuse1Dall(LR1CellIto_nsca<CELLS>* LR1Cell, double D) {
        for (int i = 0; i < CELLS; i++) {
            LR1_nsca_diffuse1D(LR1Cell, D, i);
        }
    }
}

#define EPI
#include "TTCellIto.cpp"
#undef EPI

extern "C"
{
    TTCellIto<CELLS>* TTCellIto_new() {
        return new TTCellIto<CELLS>();
    }
    TTCellIto<1>* TTCellIto_newsing() {
        return new TTCellIto<1>();
    }
    void TT_stepdt(TTCellIto<CELLS>* TTCell,int id, double dt, double st) {
        TTCell->stepdt(id,dt,st);
    }
    void TT_stepdtsing(TTCellIto<1>* TTCell, double dt, double st) {
        TTCell->stepdt(0,dt,st);
    }
    double TT_getv(TTCellIto<CELLS>* TTCell, int id) {
        if (id > -1 && id < CELLS) {  
             return TTCell->v[id];
        }
        else {return TTCell->v[0];}
    }
    double TT_getcai(TTCellIto<CELLS>* TTCell, int id) {
        if (id > -1 && id < CELLS) {  
             return TTCell->cai[id];
        }
        else {return TTCell->cai[0];}
    }
    double TT_getvsing(TTCellIto<1>* TTCell) {
        return TTCell->v[0];
    }
    double TT_getcaising(TTCellIto<1>* TTCell) {
        return TTCell->cai[0];
    }
    void TT_setv(TTCellIto<CELLS>* TTCell, int id, double v) {
        if (id > -1 && id < CELLS) {  
            TTCell->v[id] = v;
        }
    }
    void TT_setvsing(TTCellIto<1>* TTCell, double v) {
        TTCell->v[0] = v;
    }
    void TT_setito(TTCellIto<CELLS>* TTCell, int id, double itofac) {
        if (id > -1 && id < CELLS) {  
            TTCell->itofac[id] = itofac;
        }
    }
    void TT_setitosing(TTCellIto<1>* TTCell, double itofac) {
        TTCell->itofac[0] = itofac;
    }
    void TT_copycell(TTCellIto<CELLS>* TTCell, int id, TTCellIto<1>* TTCell2) {
        TTCell->setcell(id,TTCell2);
    }
    void TT_copycellsing(TTCellIto<1>* TTCell, TTCellIto<1>* TTCell2) {
        TTCell->setcell(0,TTCell2);
    }
    void TT_diffuse1D(TTCellIto<CELLS>* TTCell, double D, int id) {
        if (CELLS > 1) {
            if (id == 0){
                TTCell->diffcurrent[id] = D*(TTCell->v[id+1] - TTCell->v[id]);
            }
            else if (id < CELLS - 1) {
                TTCell->diffcurrent[id] = D*(TTCell->v[id+1] - TTCell->v[id]) + D*(TTCell->v[id-1] - TTCell->v[id]);
            }
            else if (id == CELLS - 1) {
                TTCell->diffcurrent[id] = D*(TTCell->v[id-1] - TTCell->v[id]);
            }
        }
    }
    void TT_diffuse1Dall(TTCellIto<CELLS>* TTCell, double D) {
        for (int i = 0; i < CELLS; i++) {
            TT_diffuse1D(TTCell, D, i);
        }
    }
}
