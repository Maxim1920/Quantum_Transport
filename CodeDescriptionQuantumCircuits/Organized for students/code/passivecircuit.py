# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 18:03:39 2026
implementation of passive linear circuit
contains conductors, capacitors, inductors and transformers
Version 0.1 28-01-2026
@author: dhyul
"""
#reserved variable. Conventionally, j= imaginary unit. Our definition of FT: j= - imaginary unit
J=complex(0,1)
import numpy as np
#for Hall elements
Tri=np.array([[-1.0,0,1.0],[1.0,-1.0,0],[0,1.0,-1.0]]) #current goes from/to 0->1, 1->2, 2->0
Four=np.array([[-1.0,0,0,1],[1.0,-1.0,0,0],[0,1.0,-1.0,0],[0,0,1.0,-1.0]]) #current goes from/to 0->1, 1->2, 2->3, 3->0



import numpy as np

from lingencthy import LinGenCircuit

class PassiveCircuit(LinGenCircuit):
    def __init__(self,omega=0.001):
        LinGenCircuit.__init__(self,dt=np.complex128) #complex numbers 
        self.omega=omega #has frequncy
    #####
    #adding functions
    #####    
    #define terminals and nodes: no properties, just dimension 1        
    def addNode(self,name):
        self._addNode(name,1,None)
    def addTerminal(self,name):
        self._addTerminal(name, 1, [1], None)
    #add conductor
    def addConductor(self,name,code,name1,name2,conductance=1.0):
        assert(conductance > 0),"Attempt to set a negative conductance"
        self._addConnector(name,code,name1,name2,conductance)
        #add type
        self.dics.connectors[name].type='R'
    def addInductor(self,name, code, name1,name2,inductance=1.0):
        assert(inductance > 0),"Attempt to set a negative inductance"
        self._addConnector(name,code,name1,name2,inductance)
        self.dics.connectors[name].type='L'
    def addCapacitor(self,name, code, name1,name2,capacitance=1.0):
        assert(capacitance > 0),"Attempt to set a negative capacitance"
        self._addConnector(name,code,name1,name2,capacitance)
        self.dics.connectors[name].type='C'
        #explanation of the transformer parameters
        #two branches: a - nodes 0, 1 b - nodes 2, 3
        #inductances of the branches: La and Lb
        #very small resistance of the branches: Ra, Rb
        #transinductance = alpha*polarity*sqrt(La*Lb)
        #polarity = \pm 1, alpha from 0 to 1, default = 1.0 -> ideal transformer            
    def addTransformer(self,name,code,name0,name1,name2,name3,La=1.0,Lb=1.0,Ra=0,Rb=0,alpha=0.999,polarity=1.0):
        #checkups
        for item in [La,Lb,Ra,Rb]:
            assert item >=0, "Attempt to set a negative parameter for "+name+" value ="+str(item)
        assert La !=-0, "zero inductance does not work for transformer "+name
        assert Lb !=-0, "zero inductance does not work for transformer"
        assert np.abs(alpha) <= 1.0, "alpha should be in (-1,1) range for the transformer "+name
        assert np.abs(polarity) == 1.0, "polarity should be plus-minus 1 for the transformer "+name
        assert (not ((Ra*Rb ==0) and  np.abs(alpha) ==1.0)), "Unfortunately, you cannot set zero resistances for an ideal transformer"
        self._addMultipole(name,code,[name0,name1,name2,name3],[La,Lb,Ra,Rb,alpha,polarity])
        self.dics.connectors[name].type='T'
    def addHallTripole(self,name,code,name0,name1,name2,conductance):
        assert conductance >0, "Conductance must be positive in "+name
        self._addMultipole(name,code,[name0,name1,name2],conductance)
        self.dics.connectors[name].type='H3'
    def addHallFourpole(self,name,code,name0,name1,name2,name3,conductance):
        assert conductance >0, "Conductance must be positive in "+name
        self._addMultipole(name,code,[name0,name1,name2,name3],conductance)
        self.dics.connectors[name].type='H4'
    #####
    #setting functions part
    #####
    def setOmega(self,omega):
        assert omega !=0, "Can not set zero omega."
        self.omega=omega
        for c in self.dics.connectors.values():
            c.storeConductances(self.computeConductances(c))
            self._putConductances(c)
        self.compute()
    def setConductor(self,name,conductance):
        assert(conductance > 0),"Attempt to set a negative conductance in "+name
        assert (self.connector(name).type == 'R'), "Attempt to set conductance to non-conductor"
        self._setConnector(name,conductance)
    def setInductor(self,name,inductance):
        assert(inductance > 0),"Attempt to set a negative inducance"
        assert (self.connector(name).type == 'L'), "Attempt to set inductance to non-inductor"
        self._setConnector(name,inductance)
    def setCapacitor(self,name,capacitance):
        assert(capacitance > 0),"Attempt to set a negative capacitance"
        assert (self.connector(name).type == 'C'), "Attempt to set capacitance to non-capacitor"
        self._setConnector(name,capacitance)
    def setTransformer(self,name,La=None,Lb=None,Ra=None,Rb=None,alpha=None,polarity=None):
        #understand if the parameters are set
        tr=self.connector(name)
        li=[La,Lb,Ra,Rb,alpha,polarity]
        for i in range(6):
            if li[i] == None:
                li[i]=tr.pars[i]    
        print(li)        
        for item in li[0:4]:
            assert item >=0, "Attempt to set a negative parameter for "+name+" value ="+str(item)
        for item in li [0:2]:
            assert item>0, "Transformer "+name+":No zero inductance allowed, sorry."
        assert np.abs(li[4]) <= 1.0, "alpha should be in (-1,1) range for the transformer "+name
        assert np.abs(li[5]) == 1.0, "polarity should be plus-minus 1 for the transformer "+name
        self._setConnector(name,li)
    def setHallTripole(self,name,conductance):
        assert(conductance > 0),"Attempt to set a negative conductance in "+name
        assert(self.connector(name).type =='H3'), "Not a Hall Tripole"
        self._setConnector(name,conductance)
    def setHallFourpole(self,name,conductance):
        assert(conductance > 0),"Attempt to set a negative conductance in "+name
        assert(self.connector(name).type =='H4'), "Not a Hall Fourpole"
        self._setConnector(name,conductance)
        ####
        #compute conductances part
        ####
    def computeConductances(self,c):
        co=c.pars
        ty=c.type
        if(ty=='R'):
            return -co,co,co,-co
        elif (ty=='C'):
            v=J*co*self.omega
            return -v,v,v,-v
        elif(ty=='L'):
            v=1.0/(J*co*self.omega)
            return -v,v,v,-v
        elif(ty=='T'): #must be 4 x 4 matrix
            La,Lb,Ra,Rb,alpha,p=co
            D=-La*Lb*(1-alpha**2)*(self.omega)**2+J*self.omega*(La*Rb+Ra*Lb)+Ra*Rb
            xa,xb,xc=Ra+J*self.omega*La, Rb+J*self.omega*Lb, -J*p*alpha*np.sqrt(La*Lb)*self.omega
            return((1.0/D)*np.array([[-xb,xb,-xc,xc],[xb,-xb,xc,-xc],[-xc,xc,-xa,xa],[xc,-xc,xa,-xa]]))
        elif(ty=='H3'): #3 x 3 matrix
            return(co*Tri)
        elif(ty=='H4'): #3 x 3 matrix
            return(co*Four)      
        else:
            assert False, "Unknown conductor type "+c.type+" for conductor "+str(c)
            
    
        