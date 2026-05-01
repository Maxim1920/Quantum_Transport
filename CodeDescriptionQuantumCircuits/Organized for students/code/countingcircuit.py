# -*- coding: utf-8 -*-
"""
Started on Mon Feb  2 12:31:17 2026

Counting circuits. Statistics in normal semiclassical nanostructures
made upon the publication
-----
Nazarov, YV., & Bagrets, DA. (2002). 
Circuit theory for full counting statistics in multiterminal circuits. 
Physical Review Letters, 88(19), 196801-1-196801-4.
----
Version 0.2 17-02-2026
"""

import numpy as np
from gencithy import Circuit
from scipy.optimize import minimize
from scipy.differentiate import jacobian
#fermi distibution

def Fermi (e, v, t):
    if (t == 0):
        return(0.5*(np.sign(e-v)+1) )
    else:
        return(1.0/(1.0+np.exp((e-v)/t)))

# tau3 just in case    
tau3=np.array([[1.0,0],[0,-1.0]])    
#simple derivative
def derivative(f,x,h=1e-2):
    return((f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h))/(12*h))
#Green function in a terminal

def G(f,cf):
    return(np.array([[1-2*f,-2*f*np.exp((1j)*cf)],[-2*(1-f)*np.exp(-(1j)*cf),2*f-1]])) 
def GT(energy,temperature,voltage,cf):
    f=Fermi(energy,voltage,temperature)
    return(G(f,cf))

#projection to find the counting current from the matrix current
def pro(a):
    return(-0.125*(a[0,0]-a[1,1]))
    
#main class
class CountingCircuit(Circuit):
    #initiation
    def __init__(self, energy=0, temperature=0.01):
        self.Ci = Circuit.__init__(self)
        self.energy = energy
        self.globaltemperature = temperature
        self.mask = False
        self.integrationparameters=[1.0,0.0,200] #default integration parameters
######
#Section: helper functions
######        
    def setIntegrationParameters(self, width=1.0,center=0.0, num = 200)    :
        self.integrationparameters=[width,center,num]
######
#Section: adding functions
#####
    #terminal
    def addTerminal(self, name, voltage=0.0, temperature=None,cf=0.0):
        if (temperature == None):
            te = self.globaltemperature
        else:
            te = temperature
        Circuit.addTerminal(self,name, GT(self.energy, te,voltage,cf))
        # attach properties to the terminal
        self.terminal(name).temperature = te
        self.terminal(name).voltage = voltage
        self.terminal(name).cf=cf
     #node
    def addNode(self, name, G=None):
         if (G == None):  # adding normal G
             Circuit.addNode(self, name, GT(
                 self.energy, self.globaltemperature, 0, 0))
         else:
             Circuit.addNode(self, name, G)
    # the connectors are added as usual
#######
# Section: Setting the elements functions
#######
    # this changes the energy for ALL terminals. 
    #Makes sure that the energy is always the same for all terminals
    def setEnergy(self, en):
        assert (not self.mask),"Cannot set energy when the mask is applied! Do nothing. Please remove the mask"          
        self.energy = en
        for t in self.Terminals.values():
            t.setG(GT(self.energy, t.temperature,t.voltage,t.cf))
     #does not recompute              
    def setGlobalTemperature(self, te):
        for t in self.Terminals.values():
            if (t.temperature == self.globaltemperature):
                t.setG(GT(self.energy, t.voltage, te, t.cf))
                t.temperature=te
        self.globaltemperature=te
     #set the terminal
    def setTerminal(self, name, temperature=None, voltage=None, cf=None):
         di=locals()
         t = self.terminal(name)
         #print(di)
         for s in di.keys():
             if (not s=='name' and not s=='self'): 
                 if di[s] == None:
                     pass
                 else:
                     if (self.mask and (s !='cf')):
                         print("The mask is set. Cannot change "+s+" of the terminal "+name)
                     else:
                         setattr(t, s, di[s])
         #setting the terminal
         if (self.mask):
             if (name in self.masklist) :
                 t.setG(G(1.0,t.cf))
             else:
                 t.setG(G(0.0,t.cf))
         else:         
             t.setG(GT(self.energy, t.voltage, t.temperature, t.cf))
    #sets the mask: group of terminals where f = 1, for all other terminals f=0
    def putMask(self, listofterminals):
        self.mask=True
        self.masklist=listofterminals
        for tname,t in self.Terminals.items():            
            if (tname in self.masklist):
                t.setG(G(1.0,t.cf))
            else:
                t.setG(G(0.0,t.cf))
    #unsets the mask
    def removeMask(self):
        if (self.mask):
            self.mask = False
            self.masklist=[]
            for tname in self.Terminals:
                self.setTerminal(tname)
################
# Bad events section
################
    # is called when verbose is on and iterations do not converge. can be expanded
    def badmoment(self):
        print("Energy: ",self.energy,"Global T", self.globaltemperature)
    # use this to adjust the iteration params for all nodes
    def changeAdjust(self,newadjust):
        for n in self.Nodes.values():
            n.adjust=newadjust            
###############
# Hidden function section
###############
    #returns sorted dic of the terminals
    def _sortedterminals(self):
        d={}
        for tn,tp in self.Terminals.items():
            d[tn]=tp.voltage
        d=list(sorted(d.items(),key=lambda item: item[1],reverse=True))
        d1={}
        for i in range(len(d)-1):
           thi,thiv=d[i]
           ne,nev=d[i+1]
           d1[thi]=thiv-nev
        return(d1) 
    #checks terminal list for concistencies
    def _checkTerminalList(self,li):
        assert (len(li)==len(set(li))), "Duplicates in the list of terminals"
        for t in li:
            assert( t in self.Terminals), t+" is not a terminal name"
    #this is the engine for zero-T calculations
    def _Z(self,fu,arg=None):
        a,li=0.0,[]        
        for ter,width in self._sortedterminals().items():
            li.append(ter)
            self.putMask(li)
            if arg == None: 
                s=fu()
            else:
                s=fu(arg)
            a+=s*(width)
        return(a)
    #this is the same engive for 2 returning values
    def _Z2(self,fu,arg=None):
        a,b,li=0.0,0.0,[]        
        for ter,width in self._sortedterminals().items():
            li.append(ter)
            self.putMask(li)
            if arg == None: 
                s,s1=fu()
            else:
                s,s1=fu(arg)
            a+=s*(width)
            b+=s1*(width)
        return(a,b)
    #now we can define various integration methods
    #we need integration parameters incorporated in the class...
    #this is the engine to integrate over energies
    def _INT(self,fu,check=False,arg=None): #if check = True, integrates a given function of energy...
        oldmask=self._savemask()
        self.removeMask()
        oldenergy=self.energy
        w,c,num=self.integrationparameters #get integration parameters
        # t = tanh((e-c)/2w)
        seq=[[2*w*np.arctanh(t)+c,4*w/(num*(1-t**2))] for t in np.linspace(-1+1.0/num,1-1.0/num,num)]
        a=0.0
        for en, weight in seq: #do to: store existing energy
            self.setEnergy(en)
            if (check==True):
                a+=weight*fu(en)
            else:
                if (arg==None):
                    a+=weight*fu()
                else:
                    a+=weight*fu(arg)
        self.setEnergy(oldenergy)
        if (oldmask !=None):
            self.putMask(oldmask)
        return(a)    
     #same engine for two returning values
    def _INT2(self,fu,check=False,arg=None): #if check = True, integrates a given function of energy...
         oldenergy=self.energy
         oldmask=self._savemask()
         self.removeMask()
         w,c,num=self.integrationparameters #get integration parameters
         # t = tanh((e-c)/2w)
         seq=[[2*w*np.arctanh(t)+c,4*w/(num*(1-t**2))] for t in np.linspace(-1+1.0/num,1-1.0/num,num)]
         a,b=0.0,0.0
         for en, weight in seq: #do to: store existing energy
             self.setEnergy(en)
             if (check==True):
                 s,s1=fu(en)
             else:
                 if (arg==None):
                     s,s1=fu()
                 else:
                     s,s1=fu(arg)
             a+=weight*s
             b+=weight*s1
         self.setEnergy(oldenergy)
         if (oldmask !=None):
             self.putMask(oldmask)
         return(a,b)  
     #save/restore functions
    def _savecf(self):
         self.oldcf = [ x.cf for x in self.Terminals.values()] 
    def _restorecf(self):
        for t,cf in zip(self.Terminals.keys(),self.oldcf):
            self.setTerminal(t,cf=cf)
    def _savemask(self):
        if (self.mask):
            oldmasklist=self.masklist.copy()
        else:
            oldmasklist=None
        return(oldmasklist)
    
#########
#Section: getting the results
#########
    #returns the action per interval at given settings
    def actionPI(self):
        self.solve()
        return(self.getTotalAction())
    #returns the action assuming zero temperature
    def actionZ(self):
        return(self._Z(self.actionPI))
    #now we can define various integration methods
    #we need integration parameters incorporated in the class...
    def actionINT(self):  
        return(self._INT(self.actionPI))
    #computes currents per interval
    def _currents(self,li):
        self.solve()
        return(np.real(np.array([pro(self.getTotalCurrentTerminal(d)) for d in li])))
    #takes care of the cf configuration and checks the list
    def currents(self,li,method='PI'): 
        self._savecf()
        self._checkTerminalList(li)
        for  t in self.Terminals:
            self.setTerminal(t,cf=0.0) #list of terminals
        if (method=='PI'):
            res=self._currents(li)
        elif(method=='Z'):
            res=self._Z(self._currents,arg=li)
        elif(method=='INT'):
            res=self._INT(self._currents,arg=li)
        else:
            assert False, "wrong method!"
        self._restorecf()
        return(res)
    #computes noises
    def _noises(self,li): #list of terminals
        def curs(cf): #function for chi-dependent currents. Here need its imaginary part
            for c, t in zip(cf,li):
                self.setTerminal(t,cf=c)
            self.solve()
            a=-np.imag(np.array([ pro(self.getTotalCurrentTerminal(d)) for d in li]))
            return(a)
        def f(x):
            return np.apply_along_axis(curs, axis=0, arr=x)
        return(jacobian(f,np.zeros(len(li))).df)
    #elaborations with methods
    def noises(self,li,method='PI'):
        self._savecf()
        self._checkTerminalList(li)
        for t in self.Terminals:
            if not (t in li):
                self.setTerminal(t,cf=0.0)
        if (method=='PI'):
            res=self._noises(li)
        elif(method=='Z'):
            res=self._Z(self._noises,arg=li)
        elif(method=='INT'):
            res=self._INT(self._noises,arg=li)
        else:
            assert False, "wrong method!"
        self._restorecf()
        return(res)
    #third cumulants
    def _thirdcumulants(self,li):
        def curs(cf):
            for c, t in zip(cf,li):
                self.setTerminal(t,cf=c)
            self.solve()
            return(-np.real(np.array([ pro(self.getTotalCurrentTerminal(d)) for d in li])))
        def f(x):
            return np.apply_along_axis(curs, axis=0, arr=x)
        def noises(cf):
            return(jacobian(f,cf,tolerances={'atol':1e-12,'rtol':1e-4}).df)
        def f1(x):
            return np.apply_along_axis(noises, axis=0, arr=x)
        return(jacobian(f1,np.zeros(len(li)),tolerances={'atol':1e-12,'rtol':1e-4}).df)
    def thirdcumulants(self,li,method='PI'):
        self._savecf()
        self._checkTerminalList(li)
        for t in self.Terminals:
            if not (t in li):
                self.setTerminal(t,cf=0.0)
        if (method=='PI'):
            res=self._thirdcumulants(li)
        elif(method=='Z'):
            res=self._Z(self._thirdcumulants,arg=li)
        elif(method=='INT'):
            res=self._INT(self._thirdcumulants,arg=li)
        else:
            assert False, "wrong method!"
        self._restorecf()
        return(res)
    #fourth cumulants
    def _fourthcumulants(self,li):
        
        def curs(cf):
            for c, t in zip(cf,li):
                self.setTerminal(t,cf=c)
            self.solve()
            re=[ pro(self.getTotalCurrentTerminal(d)) for d in li]
            return(np.imag(np.array(re)))
        def f(x):
            return np.apply_along_axis(curs, axis=0, arr=x)
        def noises(cf):
            return(jacobian(f,cf,tolerances={'atol':1e-12,'rtol':1e-4}).df)
        def f1(x):
            return np.apply_along_axis(noises, axis=0, arr=x)  
        def thirds(cf):
            return(jacobian(f1,cf,tolerances={'atol':1e-12,'rtol':1e-4}).df)
        def f2(x):
            return np.apply_along_axis(thirds, axis=0, arr=x) 
        return(jacobian(f2,np.zeros(len(li)),tolerances={'atol':1e-12,'rtol':1e-4}).df)
    def fourthcumulants(self,li,method='PI'):
        self._savecf()
        self._checkTerminalList(li)
        for t in self.Terminals:
            if not (t in li):
                self.setTerminal(t,cf=0.0)                 
        if (method=='PI'):
            res=self._fourthcumulants(li)
        elif(method=='Z'):
            res=self._Z(self._fourthcumulants,arg=li)
        elif(method=='INT'):
            res=self._INT(self._fourthcumulants,arg=li)
        else:
            assert False, "wrong method!"
        self._restorecf()
        return(res)
    
    #computes transmission distribution density for the list of transmissions
    def fullTransmissionDistribution(self,list1,listofTransmissions,integratedd=False):
        list2=list(filter(lambda x: (not(x in list1)),self.Terminals.keys()))
        return(self.partialTransmissionDistribution(list1,list2,listofTransmissions,integratedd=integratedd))
    #more general transmission distribution tool
    def partialTransmissionDistribution(self,list1,list2,listofTransmissions,integratedd=False):
        #check the lists
        for t in list1:
            assert( not(t in list2)), "partialTransmissionDistibution: inconsistent lists!"
        self._savecf()
        #stores old mask if present
        oldmasklist=self._savemask()
        self.putMask(list1)
        for t in self.Terminals:
            self.setTerminal(t,cf=0.0)
        def act(xi): #care about sign of xi: opposite now?
            for t in list2:
                self.setTerminal(t,cf=-xi)
            return(self.actionPI())
        eps=1e-3 #to change in future versions
        def integrated(T):
                return((1.0/np.pi)*np.imag(act(np.pi-eps-1j*np.log(1.0/T-1))))
        def density(T):
            return(-derivative(integrated,T))    
        res=[]
        #main cycle
        for t in listofTransmissions:
            print("Transmission: ",t)
            if(integratedd==True):
                res.append(integrated(t))
            elif(integratedd==False):
                res.append(density(t))
            elif(integratedd=='both'):
                res.append([density(t),integrated(t)])
        #restore stored values
        self.removeMask()
        if (oldmasklist != None):
           self.putMask(oldmasklist)
        self._restorecf()
        #returning
        return(res)
    #computes action and currents at given setting of real chi
    def _actioncurrents(self):
        a=np.real(self.actionPI())
        re=[ pro(self.getTotalCurrentTerminal(d)) for d in self.Terminals]
        return(a,np.array(re))  
    #computes action and currents at given setting of real chi
    def _gfchi(self,listofrealchi,method='PI'):
        for t,s in zip(self.Terminals,listofrealchi):
            self.setTerminal(t,cf=-1j*s)
        if (method=='PI'):
            res=self._actioncurrents()
        elif(method=='Z'):
            res=self._Z2(self._actioncurrents)
        elif(method=='INT'):
            res=self._INT2(self._actioncurrents)
        self._restorecf()
        return(res)
    # gigantic fluctuations along a ray.
    # arguments: list of chi, ray, method
    def gfRay(self,chis,ray,method='Z'): 
        self._savecf()
        res1,res2=[],[]
        for chi in chis:
            s,currs=self._gfchi(chi*ray,method=method)
            res1.append(s+np.dot(chi*ray,currs))
            res2.append(currs)
        self._restorecf()
        return(res1,res2) 
    #returns list of actions, list of currents.
    #gigantic fluctuations at given currents
    def gfIs(self,dic,method='PI',bound=5.0): #dictionary of terminals and currents
        self._checkTerminalList(dic.keys())
        for t in self.Terminals:
            if not (t in dic.keys()):
                self.setTerminal(t,cf=0.0)
        self._savecf()
        if method=='PI':
            fu=self.actionPI
        elif method=='Z':
            fu=self.actionZ
        elif(method=='INT'):
            fu=self.actionINT
        def tomini(chis):
            extra=0
            for t,chi in zip(dic.keys(),chis):
                extra+=chi*dic[t]
                self.setTerminal(t,cf=-1j*chi)
            return(np.real(fu())+extra)    
        res=minimize(tomini,np.zeros(len(dic)),bounds=[(-bound,bound) for x in range(len(dic))]) 
        s,curr=self._gfchi(np.array( [-np.imag(t.cf) for t in self.Terminals.values()]),method=method)
        self._restorecf()
        return(res.fun,res.x,np.real(curr))
    #node inspection ADDED 18-02-2026
    def inspectNode(self, name):
        G=self.node(name).getG()
        f=0.5*(1.0-G[0,0])
        if (np.abs(f) > 1e-2):
            c=-1j*np.log(G[0,1]/(-2*f))
        else:
            c=1j*np.log(G[1,0]/(2*(f-1)))
        return(f,c)