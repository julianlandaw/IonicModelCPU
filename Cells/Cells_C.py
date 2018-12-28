import ctypes
import numpy as np
import matplotlib.pyplot as plt

lib = ctypes.cdll.LoadLibrary('./libCells_C.so')

class LR1CellIto(object):
    stimduration = 0.5
    stim = -80.0
    onecell = True
    typecell = 1
    def __init__(self, onecell = True):
        self.onecell = onecell
        
        lib.LR1CellIto_new.argtypes = None
        lib.LR1CellIto_new.restype = ctypes.c_void_p

        lib.LR1_stepdt.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_double]
        lib.LR1_stepdt.restype = ctypes.c_void_p
        
        lib.LR1_getv.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR1_getv.restype = ctypes.c_double
        
        lib.LR1_getcai.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR1_getcai.restype = ctypes.c_double
        
        lib.LR1_setv.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_setv.restype = ctypes.c_void_p
        
        lib.LR1_setito.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_setito.restype = ctypes.c_void_p
        
        lib.LR1_setxr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_setxr.restype = ctypes.c_void_p
        
        lib.LR1_getxr.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR1_getxr.restype = ctypes.c_double
        
        lib.LR1_settauxfac.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_settauxfac.restype = ctypes.c_void_p
        
        lib.LR1_setibarcafac.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_setibarcafac.restype = ctypes.c_void_p
        
        lib.LR1CellIto_newsing.argtypes = None
        lib.LR1CellIto_newsing.restype = ctypes.c_void_p

        lib.LR1_stepdtsing.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double]
        lib.LR1_stepdtsing.restype = ctypes.c_void_p
        
        lib.LR1_getvsing.argtypes = [ctypes.c_void_p]
        lib.LR1_getvsing.restype = ctypes.c_double
        
        lib.LR1_getcaising.argtypes = [ctypes.c_void_p]
        lib.LR1_getcaising.restype = ctypes.c_double
        
        lib.LR1_setvsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_setvsing.restype = ctypes.c_void_p
        
        lib.LR1_setitosing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_setitosing.restype = ctypes.c_void_p
        
        lib.LR1_setxrsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_setxrsing.restype = ctypes.c_void_p
        
        lib.LR1_getxrsing.argtypes = [ctypes.c_void_p]
        lib.LR1_getxrsing.restype = ctypes.c_double
        
        lib.LR1_settauxfacsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_settauxfacsing.restype = ctypes.c_void_p
        
        lib.LR1_setibarcafacsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_setibarcafacsing.restype = ctypes.c_void_p
        
        lib.LR1_copycell.argtypes = [ctypes.c_void_p,ctypes.c_int, ctypes.c_void_p]
        lib.LR1_copycell.restype = ctypes.c_void_p
        
        lib.LR1_copycellsing.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        lib.LR1_copycellsing.restype = ctypes.c_void_p
        
        lib.LR1_diffuse1D.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_int]
        lib.LR1_diffuse1D.restype = ctypes.c_void_p
        
        lib.LR1_diffuse1Dall.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_diffuse1Dall.restype = ctypes.c_void_p
        
        if (onecell):
            self.obj = lib.LR1CellIto_newsing()
        else:
            self.obj = lib.LR1CellIto_new()
       
    def stepdt(self,dt,st,ind=0):
        if (self.onecell is False):
            lib.LR1_stepdt(self.obj,ctypes.c_int(ind),ctypes.c_double(dt),ctypes.c_double(st))
        else:
            lib.LR1_stepdtsing(self.obj,ctypes.c_double(dt),ctypes.c_double(st))
        
    def getv(self,ind=0):
        if (self.onecell is False):
            return lib.LR1_getv(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR1_getvsing(self.obj)
        
    def getcai(self,ind=0):
        if (self.onecell is False):
            return lib.LR1_getcai(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR1_getcaising(self.obj)    
    
    def setv(self,v,ind=0):
        if (self.onecell is False):
            lib.LR1_setv(self.obj,ctypes.c_int(ind),ctypes.c_double(v))
        else:
            lib.LR1_setvsing(self.obj,ctypes.c_double(v))
        
    def setito(self,itofac,ind=0):
        if (self.onecell is False):
            lib.LR1_setito(self.obj,ctypes.c_int(ind),ctypes.c_double(itofac))
        else:
            lib.LR1_setitosing(self.obj,ctypes.c_double(itofac))
    def setxr(self,xr,ind=0):
        if (self.onecell is False):
            lib.LR1_setxr(self.obj,ctypes.c_int(ind),ctypes.c_double(xr))
        else:
            lib.LR1_setxrsing(self.obj,ctypes.c_double(xr))
    def getxr(self,ind=0):
        if (self.onecell is False):
            return lib.LR1_getxr(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR1_getxrsing(self.obj)
    
    def settauxfac(self,tauxf,ind=0):
        if (self.onecell is False):
            lib.LR1_settauxfac(self.obj,ctypes.c_int(ind),ctypes.c_double(tauxf))
        else:
            lib.LR1_settauxfacsing(self.obj,ctypes.c_double(tauxf))
            
    def setibarcafac(self,ibarca,ind=0):
        if (self.onecell is False):
            lib.LR1_setibarcafac(self.obj,ctypes.c_int(ind),ctypes.c_double(ibarca))
        else:
            lib.LR1_setibarcafacsing(self.obj,ctypes.c_double(ibarca))
    
    def copycell(self,newcell,ind=0):
        if (self.onecell is False):
            lib.LR1_copycell(self.obj,ctypes.c_int(ind),newcell.obj)
        else:
            lib.LR1_copycellsing(self.obj,newcell.obj)
    
    def diffuse(self, D, ind):
        if (self.onecell is False):
            lib.LR1_diffuse1D(self.obj, ctypes.c_double(D), ctypes.c_int(ind))
    
    def diffuseall(self, D):
        if (self.onecell is False):
            lib.LR1_diffuse1Dall(self.obj, ctypes.c_double(D))
            
class LR2CellIto(object):
    stimduration = 0.5
    stim = -80.0
    onecell = True
    typecell = 2
    def __init__(self, onecell = True):
        self.onecell = onecell
        
        lib.LR2CellIto_new.argtypes = None
        lib.LR2CellIto_new.restype = ctypes.c_void_p

        lib.LR2_stepdt.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_double]
        lib.LR2_stepdt.restype = ctypes.c_void_p
        
        lib.LR2_getv.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR2_getv.restype = ctypes.c_double
        
        lib.LR2_getcai.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR2_getcai.restype = ctypes.c_double
        
        lib.LR2_getnai.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR2_getnai.restype = ctypes.c_double
        
        lib.LR2_getki.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR2_getki.restype = ctypes.c_double
        
        lib.LR2_setv.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR2_setv.restype = ctypes.c_void_p
        
        lib.LR2_setito.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR2_setito.restype = ctypes.c_void_p
        
        lib.LR2_setisk.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR2_setisk.restype = ctypes.c_void_p

        lib.LR2_setskh.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR2_setskh.restype = ctypes.c_void_p
        
        lib.LR2_setxr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR2_setxr.restype = ctypes.c_void_p
        
        lib.LR2_getxr.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR2_getxr.restype = ctypes.c_double
        
        lib.LR2_settauxfac.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR2_settauxfac.restype = ctypes.c_void_p
        
        lib.LR2_setibarcafac.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR2_setibarcafac.restype = ctypes.c_void_p
        
        lib.LR2_setiupfac.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR2_setiupfac.restype = ctypes.c_void_p
        
        lib.LR2CellIto_newsing.argtypes = None
        lib.LR2CellIto_newsing.restype = ctypes.c_void_p

        lib.LR2_stepdtsing.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double]
        lib.LR2_stepdtsing.restype = ctypes.c_void_p
        
        lib.LR2_getvsing.argtypes = [ctypes.c_void_p]
        lib.LR2_getvsing.restype = ctypes.c_double
        
        lib.LR2_getcaising.argtypes = [ctypes.c_void_p]
        lib.LR2_getcaising.restype = ctypes.c_double
        
        lib.LR2_getnaising.argtypes = [ctypes.c_void_p]
        lib.LR2_getnaising.restype = ctypes.c_double
        
        lib.LR2_getkising.argtypes = [ctypes.c_void_p]
        lib.LR2_getkising.restype = ctypes.c_double
        
        lib.LR2_setvsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR2_setvsing.restype = ctypes.c_void_p
        
        lib.LR2_setitosing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR2_setitosing.restype = ctypes.c_void_p
        
        lib.LR2_setisksing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR2_setisksing.restype = ctypes.c_void_p
        
        lib.LR2_setskhsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR2_setskhsing.restype = ctypes.c_void_p

        lib.LR2_setxrsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR2_setxrsing.restype = ctypes.c_void_p
        
        lib.LR2_getxrsing.argtypes = [ctypes.c_void_p]
        lib.LR2_getxrsing.restype = ctypes.c_double
        
        lib.LR2_settauxfacsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR2_settauxfacsing.restype = ctypes.c_void_p
        
        lib.LR2_setiupfacsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR2_setiupfacsing.restype = ctypes.c_void_p
        
        lib.LR2_setibarcafacsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR2_setibarcafacsing.restype = ctypes.c_void_p
        
        lib.LR2_copycell.argtypes = [ctypes.c_void_p,ctypes.c_int, ctypes.c_void_p]
        lib.LR2_copycell.restype = ctypes.c_void_p
        
        lib.LR2_copycellsing.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        lib.LR2_copycellsing.restype = ctypes.c_void_p
        
        lib.LR2_diffuse1D.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_int]
        lib.LR2_diffuse1D.restype = ctypes.c_void_p
        
        lib.LR2_diffuse1Dall.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR2_diffuse1Dall.restype = ctypes.c_void_p
        
        if (onecell):
            self.obj = lib.LR2CellIto_newsing()
        else:
            self.obj = lib.LR2CellIto_new()
       
    def stepdt(self,dt,st,ind=0):
        if (self.onecell is False):
            lib.LR2_stepdt(self.obj,ctypes.c_int(ind),ctypes.c_double(dt),ctypes.c_double(st))
        else:
            lib.LR2_stepdtsing(self.obj,ctypes.c_double(dt),ctypes.c_double(st))
        
    def getv(self,ind=0):
        if (self.onecell is False):
            return lib.LR2_getv(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR2_getvsing(self.obj)
        
    def getcai(self,ind=0):
        if (self.onecell is False):
            return lib.LR2_getcai(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR2_getcaising(self.obj) 
        
    def getnai(self,ind=0):
        if (self.onecell is False):
            return lib.LR2_getnai(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR2_getnaising(self.obj)
        
    def getki(self,ind=0):
        if (self.onecell is False):
            return lib.LR2_getki(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR2_getkising(self.obj)    
    
    def setv(self,v,ind=0):
        if (self.onecell is False):
            lib.LR2_setv(self.obj,ctypes.c_int(ind),ctypes.c_double(v))
        else:
            lib.LR2_setvsing(self.obj,ctypes.c_double(v))
        
    def setito(self,itofac,ind=0):
        if (self.onecell is False):
            lib.LR2_setito(self.obj,ctypes.c_int(ind),ctypes.c_double(itofac))
        else:
            lib.LR2_setitosing(self.obj,ctypes.c_double(itofac))
            
    def setisk(self,iskfac,ind=0):
        if (self.onecell is False):
            lib.LR2_setisk(self.obj,ctypes.c_int(ind),ctypes.c_double(iskfac))
        else:
            lib.LR2_setisksing(self.obj,ctypes.c_double(iskfac))

    def setskh(self,skh,ind=0):
        if (self.onecell is False):
            lib.LR2_setskh(self.obj,ctypes.c_int(ind),ctypes.c_double(skh))
        else:
            lib.LR2_setskhsing(self.obj,ctypes.c_double(skh))

    def setxr(self,xr,ind=0):
        if (self.onecell is False):
            lib.LR2_setxr(self.obj,ctypes.c_int(ind),ctypes.c_double(xr))
        else:
            lib.LR2_setxrsing(self.obj,ctypes.c_double(xr))
    def getxr(self,ind=0):
        if (self.onecell is False):
            return lib.LR2_getxr(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR2_getxrsing(self.obj)
    
    def settauxfac(self,tauxf,ind=0):
        if (self.onecell is False):
            lib.LR2_settauxfac(self.obj,ctypes.c_int(ind),ctypes.c_double(tauxf))
        else:
            lib.LR2_settauxfacsing(self.obj,ctypes.c_double(tauxf))
            
    def setibarcafac(self,ibarca,ind=0):
        if (self.onecell is False):
            lib.LR2_setibarcafac(self.obj,ctypes.c_int(ind),ctypes.c_double(ibarca))
        else:
            lib.LR2_setibarcafacsing(self.obj,ctypes.c_double(ibarca))
            
    def setiup(self,iup,ind=0):
        if (self.onecell is False):
            lib.LR2_setiupfac(self.obj,ctypes.c_int(ind),ctypes.c_double(iup))
        else:
            lib.LR2_setiupfacsing(self.obj,ctypes.c_double(iup))        
    
    def copycell(self,newcell,ind=0):
        if (self.onecell is False):
            lib.LR2_copycell(self.obj,ctypes.c_int(ind),newcell.obj)
        else:
            lib.LR2_copycellsing(self.obj,newcell.obj)
    
    def diffuse(self, D, ind):
        if (self.onecell is False):
            lib.LR2_diffuse1D(self.obj, ctypes.c_double(D), ctypes.c_int(ind))
    
    def diffuseall(self, D):
        if (self.onecell is False):
            lib.LR2_diffuse1Dall(self.obj, ctypes.c_double(D))            
            

class LR1CellIto_nsca(object):
    stimduration = 0.5
    stim = -80.0
    onecell = True
    typecell = 3
    def __init__(self, onecell = True):
        self.onecell = onecell
        
        lib.LR1CellIto_nsca_new.argtypes = None
        lib.LR1CellIto_nsca_new.restype = ctypes.c_void_p

        lib.LR1_nsca_stepdt.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_double]
        lib.LR1_nsca_stepdt.restype = ctypes.c_void_p
        
        lib.LR1_nsca_getv.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR1_nsca_getv.restype = ctypes.c_double
        
        lib.LR1_nsca_getcai.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR1_nsca_getcai.restype = ctypes.c_double
        
        lib.LR1_nsca_setv.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_nsca_setv.restype = ctypes.c_void_p
        
        lib.LR1_nsca_setito.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_nsca_setito.restype = ctypes.c_void_p
        
        lib.LR1_nsca_setinsca.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_nsca_setinsca.restype = ctypes.c_void_p
        
        lib.LR1_nsca_setxr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_nsca_setxr.restype = ctypes.c_void_p   
        
        lib.LR1_nsca_getxr.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.LR1_nsca_getxr.restype = ctypes.c_double
        
        lib.LR1_nsca_settauxfac.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_nsca_settauxfac.restype = ctypes.c_void_p
        
        lib.LR1_nsca_setibarcafac.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.LR1_nsca_setibarcafac.restype = ctypes.c_void_p
        
        lib.LR1CellIto_nsca_newsing.argtypes = None
        lib.LR1CellIto_nsca_newsing.restype = ctypes.c_void_p

        lib.LR1_nsca_stepdtsing.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double]
        lib.LR1_nsca_stepdtsing.restype = ctypes.c_void_p
        
        lib.LR1_nsca_getvsing.argtypes = [ctypes.c_void_p]
        lib.LR1_nsca_getvsing.restype = ctypes.c_double
        
        lib.LR1_nsca_getcaising.argtypes = [ctypes.c_void_p]
        lib.LR1_nsca_getcaising.restype = ctypes.c_double
        
        lib.LR1_nsca_setvsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_nsca_setvsing.restype = ctypes.c_void_p
        
        lib.LR1_nsca_setitosing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_nsca_setitosing.restype = ctypes.c_void_p
        
        lib.LR1_nsca_setinscasing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_nsca_setinscasing.restype = ctypes.c_void_p
        
        lib.LR1_nsca_setxrsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_nsca_setxrsing.restype = ctypes.c_void_p   
        
        lib.LR1_nsca_getxrsing.argtypes = [ctypes.c_void_p]
        lib.LR1_nsca_getxrsing.restype = ctypes.c_double
        
        lib.LR1_nsca_settauxfacsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_nsca_settauxfacsing.restype = ctypes.c_void_p
        
        lib.LR1_nsca_setibarcafacsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_nsca_setibarcafacsing.restype = ctypes.c_void_p
        
        lib.LR1_nsca_copycell.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
        lib.LR1_nsca_copycell.restype = ctypes.c_void_p 
        
        lib.LR1_nsca_copycellsing.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        lib.LR1_nsca_copycellsing.restype = ctypes.c_void_p
        
        lib.LR1_nsca_diffuse1D.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_int]
        lib.LR1_nsca_diffuse1D.restype = ctypes.c_void_p
        
        lib.LR1_nsca_diffuse1Dall.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.LR1_nsca_diffuse1Dall.restype = ctypes.c_void_p
        
        if (onecell):
            self.obj = lib.LR1CellIto_nsca_newsing()
        else:
            self.obj = lib.LR1CellIto_nsca_new()
            
    def stepdt(self,dt,st,ind=0):
        if (self.onecell is False):
            lib.LR1_nsca_stepdt(self.obj,ctypes.c_int(ind),ctypes.c_double(dt),ctypes.c_double(st))
        else:
            lib.LR1_nsca_stepdtsing(self.obj,ctypes.c_double(dt),ctypes.c_double(st))
        
    def getv(self,ind=0):
        if (self.onecell is False):
            return lib.LR1_nsca_getv(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR1_nsca_getvsing(self.obj)
    
    def getcai(self,ind=0):
        if (self.onecell is False):
            return lib.LR1_nsca_getcai(self.obj,ctypes.c_int(ind))
        else:
            return lib.LR1_nsca_getcaising(self.obj)
        
    def setv(self,v,ind=0):
        if (self.onecell is False):
            lib.LR1_nsca_setv(self.obj,ctypes.c_int(ind),ctypes.c_double(v))
        else:
            lib.LR1_nsca_setvsing(self.obj,ctypes.c_double(v))
            
    def setito(self,itofac,ind=0):
        if (self.onecell is False):
            lib.LR1_nsca_setito(self.obj,ctypes.c_int(ind),ctypes.c_double(itofac))        
        else:
            lib.LR1_nsca_setitosing(self.obj,ctypes.c_double(itofac)) 
            
    def setinsca(self,inscafac,ind=0):
        if (self.onecell is False):
            lib.LR1_nsca_setinsca(self.obj,ctypes.c_int(ind),ctypes.c_double(inscafac))
        else:
            lib.LR1_nsca_setinscasing(self.obj,ctypes.c_double(inscafac))
    
    def setxr(self,xr,ind=0):
        if (self.onecell is False):
            lib.LR1_nsca_setxr(self.obj,ctypes.c_int(ind),ctypes.c_double(xr))
        else:
            lib.LR1_nsca_setxrsing(self.obj,ctypes.c_double(xr))
        
    def getxr(self,ind=0):
        if (self.onecell is False):
            return lib.LR1_nsca_getxr(self.obj,ctypes.c_int(ind)) 
        else:
            return lib.LR1_nsca_getxrsing(self.obj)
    
    def settauxfac(self,tauxf,ind=0):
        if (self.onecell is False):
            lib.LR1_nsca_settauxfac(self.obj,ctypes.c_int(ind),ctypes.c_double(tauxf))
        else:
            lib.LR1_nsca_settauxfacsing(self.obj,ctypes.c_double(tauxf))
            
    def setibarcafac(self,ibarca,ind=0):
        if (self.onecell is False):
            lib.LR1_nsca_setibarcafac(self.obj,ctypes.c_int(ind),ctypes.c_double(ibarca))
        else:
            lib.LR1_nsca_setibarcafacsing(self.obj,ctypes.c_double(ibarca))        
    
    def copycell(self,newcell,ind=0):
        if (self.onecell is False):
            lib.LR1_nsca_copycell(self.obj,ctypes.c_int(ind),newcell.obj)    
        else:
            lib.LR1_nsca_copycellsing(self.obj,newcell.obj)
           
    def diffuse(self, D, ind):
        if (self.onecell is False):
            lib.LR1_nsca_diffuse1D(self.obj, ctypes.c_double(D), ctypes.c_int(ind))
    
    def diffuseall(self, D):
        if (self.onecell is False):
            lib.LR1_nsca-diffuse1Dall(self.obj, ctypes.c_double(D))    
        
class TTCellIto(object):
    stimduration = 1.0
    stim = -52.0
    onecell = True
    typecell = 4
    def __init__(self, onecell = True):
        self.onecell = onecell
        
        lib.TTCellIto_new.argtypes = None
        lib.TTCellIto_new.restype = ctypes.c_void_p

        lib.TT_stepdt.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_double]
        lib.TT_stepdt.restype = ctypes.c_void_p
        
        lib.TT_getv.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.TT_getv.restype = ctypes.c_double
        
        lib.TT_getcai.argtypes = [ctypes.c_void_p,ctypes.c_int]
        lib.TT_getcai.restype = ctypes.c_double
        
        lib.TT_setv.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.TT_setv.restype = ctypes.c_void_p
        
        lib.TT_setito.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.TT_setito.restype = ctypes.c_void_p
        
        lib.TTCellIto_newsing.argtypes = None
        lib.TTCellIto_newsing.restype = ctypes.c_void_p

        lib.TT_stepdtsing.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double]
        lib.TT_stepdtsing.restype = ctypes.c_void_p
        
        lib.TT_getvsing.argtypes = [ctypes.c_void_p]
        lib.TT_getvsing.restype = ctypes.c_double
        
        lib.TT_getcaising.argtypes = [ctypes.c_void_p]
        lib.TT_getcaising.restype = ctypes.c_double
        
        lib.TT_setvsing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.TT_setvsing.restype = ctypes.c_void_p
        
        lib.TT_setitosing.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.TT_setitosing.restype = ctypes.c_void_p
        
        lib.TT_copycell.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
        lib.TT_copycell.restype = ctypes.c_void_p
        
        lib.TT_copycellsing.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        lib.TT_copycellsing.restype = ctypes.c_void_p
        
        lib.TT_diffuse1D.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_int]
        lib.TT_diffuse1D.restype = ctypes.c_void_p
        
        lib.TT_diffuse1Dall.argtypes = [ctypes.c_void_p, ctypes.c_double]
        lib.TT_diffuse1Dall.restype = ctypes.c_void_p
        
        if (onecell):
            self.obj = lib.TTCellIto_newsing()
        else:
            self.obj = lib.TTCellIto_new()
            
    def stepdt(self,dt,st,ind=0):
        if (self.onecell is False):
            lib.TT_stepdt(self.obj,ctypes.c_int(ind),ctypes.c_double(dt),ctypes.c_double(st))
        else:
            lib.TT_stepdtsing(self.obj,ctypes.c_double(dt),ctypes.c_double(st))
        
    def getv(self,ind=0):
        if (self.onecell is False):
            return lib.TT_getv(self.obj,ctypes.c_int(ind))
        else:
            return lib.TT_getvsing(self.obj)
        
    def getcai(self,ind=0):
        if (self.onecell is False):
            return lib.TT_getcai(self.obj,ctypes.c_int(ind))
        else:
            return lib.TT_getcaising(self.obj)    
    
    def setv(self,v,ind=0):
        if (self.onecell is False):
            lib.TT_setv(self.obj,ctypes.c_int(ind),ctypes.c_double(v))
        else:
            lib.TT_setvsing(self.obj,ctypes.c_double(v))
        
    def setito(self,itofac,ind=0):
        if (self.onecell is False):
            lib.TT_setito(self.obj,ctypes.c_int(ind),ctypes.c_double(itofac))
        else:
            lib.TT_setitosing(self.obj,ctypes.c_double(itofac))
        
    def copycell(self,newcell,ind=0):
        if (self.onecell is False):
            lib.TT_copycell(self.obj,ctypes.c_int(ind),newcell.obj)    
        else:
            lib.TT_copycellsing(self.obj,newcell.obj)
            
    def diffuse(self, D, ind):
        if (self.onecell is False):
            lib.TT_diffuse1D(self.obj, ctypes.c_double(D), ctypes.c_int(ind))
    
    def diffuseall(self, D):
        if (self.onecell is False):
            lib.TT_diffuse1Dall(self.obj, ctypes.c_double(D))        
            
def pacecell(x,pcl=1000,beats=50,beatsave=10,dt=0.05):
    globalt = -1000
    tsave = -200
    savets = []
    savevs = []
    savecais = []
    savenais = []
    savekis = []
        
    while (globalt < -dt/4):
        x.stepdt(dt,0)
        globalt = globalt + dt
        if (globalt > tsave - dt/4):
            savets.append(globalt)
            savevs.append(x.getv())
            savecais.append(x.getcai())
            if (x.typecell == 2 or x.typecell == 4):
                savenais.append(x.getnai())
                savekis.append(x.getki())
            tsave = tsave + 1

    tstim = 0.0
    inapd = False
    apds = []
    startapd = 0
    while (globalt < pcl*beats + 500):
        if (tstim > -dt/4 and tstim < x.stimduration - dt/4):
            if (inapd is False and x.getv() > -75):
                x.stepdt(dt,0)
                tstim = tstim + dt - pcl
                globalt = globalt + dt
                if (globalt > tsave - dt/4):
                    savets.append(globalt)
                    savevs.append(x.getv())
                    savecais.append(x.getcai())
                    if (x.typecell == 2 or x.typecell == 4):
                        savenais.append(x.getnai())
                        savekis.append(x.getki())
                    tsave = tsave + 1
            else:    
                if (inapd is False):
                    startapd = globalt
                
                inapd = True
                x.stepdt(dt,x.stim)
                globalt = globalt + dt
                tstim = tstim + dt
                if (globalt > tsave - dt/4):
                    savets.append(globalt)
                    savevs.append(x.getv())
                    savecais.append(x.getcai())
                    if (x.typecell == 2 or x.typecell == 4):
                        savenais.append(x.getnai())
                        savekis.append(x.getki())
                    tsave = tsave + 1
        else:
            inapd = False
            vold = x.getv()
            x.stepdt(dt,0)
            vnew = x.getv()
            globalt = globalt + dt
            if (vnew <= -75 and vold > -75):
                apds.append(globalt - startapd)
            tstim = tstim + dt
            while (tstim > x.stimduration):
                tstim = tstim - pcl
            if (globalt > tsave - dt/4):
                savets.append(globalt)
                savevs.append(x.getv())
                savecais.append(x.getcai())
                if (x.typecell == 2 or x.typecell == 4):
                    savenais.append(x.getnai())
                    savekis.append(x.getki())
                tsave = tsave + 1
    apds = np.asarray(apds)
    savets = np.asarray(savets)
    savevs = np.asarray(savevs)
    savecais = np.asarray(savecais)
    if (x.typecell == 2 or x.typecell == 4):
        savenais = np.asarray(savenais)
        savekis = np.asarray(savekis)
    tstart = 0 #pcl*(beats - beatsave + 1)
    apdsaved = len(apds)
    return (apds[np.arange(max(0,apdsaved -beatsave),apdsaved)],savets[tstart:],savevs[tstart:],savecais[tstart:],savenais[tstart:],savekis[tstart:])

def restitution(x,dirange,pcl,prebeats,dt):
    globalt = -1000
    while (globalt < -dt/4):
        x.stepdt(dt,0)
        globalt = globalt + dt
    tstim = 0.0
    inapd = False
    numprebeats = 0
    while (numprebeats < prebeats and globalt < pcl*(4*prebeats)):
        if (tstim > -dt/4 and tstim < x.stimduration - dt/4):
            if (inapd is False and x.getv() > -75):
                x.stepdt(dt,0)
                tstim = tstim + dt - pcl
                globalt = globalt + dt
            else:    
                if (inapd is False):
                    startapd = globalt
                
                inapd = True
                x.stepdt(dt,x.stim)
                globalt = globalt + dt
                tstim = tstim + dt
        else:
            inapd = False
            vold = x.getv()
            x.stepdt(dt,0)
            vnew = x.getv()
            globalt = globalt + dt
            if (vnew <= -75 and vold > -75):
                numprebeats = numprebeats + 1
            tstim = tstim + dt
            while (tstim > x.stimduration):
                tstim = tstim - pcl
    apds = []
    xnew = type(x)(onecell = True)
    for di in dirange:
        tnew = 0            
        xnew.copycell(x)
        while (tnew < di - dt/4):
            xnew.stepdt(dt,0)
            tnew = tnew + dt
        tstim = 0
        while (tstim < xnew.stimduration - dt/4):
            xnew.stepdt(dt,xnew.stim)
            tstim = tstim + dt
        while (xnew.getv() > -75):
            xnew.stepdt(dt,0)
            tstim = tstim + dt
        apds.append(tstim)
        #print(di)
    return np.asarray(apds)

def Xrestitution(x,xrange_,pcl,prebeats,dt):
    globalt = -1000
    while (globalt < -dt/4):
        x.stepdt(dt,0)
        globalt = globalt + dt
    tstim = 0.0
    inapd = False
    numprebeats = 0
    while (numprebeats < prebeats and globalt < pcl*(4*prebeats)):
        if (tstim > -dt/4 and tstim < x.stimduration - dt/4):
            if (inapd is False and x.getv() > -75):
                x.stepdt(dt,0)
                tstim = tstim + dt - pcl
                globalt = globalt + dt
            else:    
                if (inapd is False):
                    startapd = globalt
                
                inapd = True
                x.stepdt(dt,x.stim)
                globalt = globalt + dt
                tstim = tstim + dt
        else:
            inapd = False
            vold = x.getv()
            x.stepdt(dt,0)
            vnew = x.getv()
            globalt = globalt + dt
            if (vnew <= -75 and vold > -75):
                numprebeats = numprebeats + 1
            tstim = tstim + dt
            while (tstim > x.stimduration):
                tstim = tstim - pcl
    while (tstim < -dt/4):
        x.stepdt(dt,0)
        globalt = globalt + dt
        tstim = tstim + dt
    print(tstim,globalt,x.getv(),x.getxr())
    apds = []
    xnew = type(x)(onecell = True)
    for xr in xrange_:      
        tstim = 0
        xnew.copycell(x)
        xnew.setxr(xr)
        while (tstim < xnew.stimduration - dt/4):
            xnew.stepdt(dt,xnew.stim)
            tstim = tstim + dt
        while (xnew.getv() > -75):
            xnew.stepdt(dt,0)
            tstim = tstim + dt
        apds.append(tstim)
    return np.asarray(apds)

def cable1D(x,D,pcl=1000,beats=50,beatsave=10,dt=0.05):
    globalt = -1000
    tsave = -200
    savets = []
    savevsall = []
    savecaisall = []
    while (globalt < -dt/4):
        for j in np.arange(0,300):
            x.stepdt(dt,0,j)      
        
        x.diffuseall(D)
        globalt = globalt + dt
        if (globalt > tsave - dt/4):
            savets.append(globalt)
            vs = []
            cais = []
            for j in np.arange(0,300):
                vs.append(x.getv(j))
                cais.append(x.getcai(j))
            savevsall.append(vs)
            savecaisall.append(cais)
            tsave = tsave + 1

    tstim = 0.0
    inapd = False
    apds = []
    startapd = 0
    while (globalt < pcl*beats + 500):
        if (tstim > -dt/4 and tstim < x.stimduration - dt/4):
            if (inapd is False and x.getv(9) > -75):
                for j in np.arange(0,300):
                    x.stepdt(dt,0,j)
                x.diffuseall(D)    
                tstim = tstim + dt - pcl
                globalt = globalt + dt
                if (globalt > tsave - dt/4):
                    vs = []
                    cais = []
                    savets.append(globalt)
                    for j in np.arange(0,300):
                        vs.append(x.getv(j))
                        cais.append(x.getcai(j))
                    savevsall.append(vs)
                    savecaisall.append(cais)
                    tsave = tsave + 1
            else:    
                if (inapd is False):
                    startapd = globalt
                
                inapd = True
                for j in np.arange(0,10):
                    x.stepdt(dt,x.stim,j)
                    
                for j in np.arange(10,300):
                    x.stepdt(dt,0,j)
                    
                x.diffuseall(D)    
                    
                globalt = globalt + dt
                tstim = tstim + dt
                if (globalt > tsave - dt/4):
                    vs = []
                    cais = []
                    savets.append(globalt)
                    for j in np.arange(0,300):
                        vs.append(x.getv(j))
                        cais.append(x.getcai(j))
                    savevsall.append(vs)
                    savecaisall.append(cais)
                    tsave = tsave + 1
        else:
            inapd = False
            vold = x.getv(9)
            for j in np.arange(0,300):
                x.stepdt(dt,0,j)
                
            x.diffuseall(D)    
            vnew = x.getv(9)
            globalt = globalt + dt
            if (vnew <= -75 and vold > -75):
                apds.append(globalt - startapd)
            tstim = tstim + dt
            while (tstim > x.stimduration):
                tstim = tstim - pcl
            if (globalt > tsave - dt/4):
                vs = []
                cais = []
                savets.append(globalt)
                for j in np.arange(0,300):
                    vs.append(x.getv(j))
                    cais.append(x.getcai(j))
                    
                savevsall.append(vs)
                savecaisall.append(cais)
                tsave = tsave + 1
                
    apds = np.asarray(apds)
    savets = np.asarray(savets)
    savevsall = np.asmatrix(savevsall)
    tstart = pcl*(beats - beatsave + 1)
    apdsaved = len(apds)
    return (apds[np.arange(max(0,apdsaved - beatsave),apdsaved)],savets[tstart:],savevsall[tstart:,:],savecaisall[tstart:,:])

        
