# -*- coding: utf-8 -*-
"""
defines Keldysh Superconducting Matrices
matrices and functions specific for superconducting circuits
SuperCircuit class
version 1.0 1-1-2026
@author: Yuli
"""
import numpy as np

from gencithy import normalize, Circuit
from scipy.integrate import quad_vec #needed for adaptive integration only

#######
# Keldysh equilibrium functions and useful matrices
#######
#advanced G  Dynes parameter
def G_A(energy, phase=0, Delta=1.0, Gamma=1E-4):
    # this is two-by-two complex matrix
    en = (energy-1j*Gamma)
    de = Delta*np.exp(1j*phase)
    return (normalize(-1j*np.array([[en, np.conj(de)], [-de, -en]])))

#retarded
def G_R(energy, phase=0, Delta=1.0, Gamma=1E-4):
    # this is two-by-two complex matrix
    en = energy+1j*Gamma
    de = Delta*np.exp(1j*phase)
    return (normalize(-1j*np.array([[en, np.conj(de)], [-de, -en]])))

#Keldysh component
def G_K(energy, temperature=0, phase=0, Delta=1.0, Gamma=1E-4):
    RminA = G_R(energy, phase=phase, Delta=Delta, Gamma=Gamma) - \
        G_A(energy, phase=phase, Delta=Delta, Gamma=Gamma)
    if (temperature == 0):
        return (RminA*np.sign(energy))
    else:
        return (RminA*np.tanh(energy/(2*temperature)))
# 4x4 superconducting
def G_S(energy, temperature=0, phase=0, Delta=1.0, Gamma=1E-4):
    R = G_R(energy, phase=phase, Delta=Delta, Gamma=Gamma)
    A = G_A(energy, phase=phase, Delta=Delta, Gamma=Gamma)
    K = G_K(energy, temperature=temperature,
            phase=phase, Delta=Delta, Gamma=Gamma)
    Z = np.zeros((2, 2), dtype=complex)
    return (np.block([[R, K], [Z, A]]))
# 2x2 tau3
def tau3():
    return (np.array([[1.0+0j, 0], [0, -1.0+0j]]))
# 4x4 tau3
def tau3four():
    Z = np.zeros((2, 2), dtype=complex)
    return (np.block([[tau3(), Z], [Z, tau3()]]))
# 4x4 normal
def G_N(energy, voltage=0, temperature=0):
    R = tau3()
    A = -tau3()
    Z = np.zeros((2, 2), dtype=complex)
    if (temperature == 0):
        h = np.array([[np.sign(energy+voltage)+0j, 0],
                     [0, np.sign(energy-voltage)]])
    else:
        h = np.array([[np.tanh((energy+voltage)/(2*temperature))+0j, 0],
                     [0, np.tanh((energy-voltage)/(2*temperature))]])
    return (np.block([[R, np.dot(R, h)-np.dot(h, A)], [Z, A]]))
#el. current from matrix current
def elcurrent(G):
    return (float(0.125*np.real(G[0, 2]-G[1, 3])))
# Matrices used for differential conductance
def G_d(h1, h2):
    R = tau3()
    A = -tau3()
    Z = np.zeros((2, 2), dtype=complex)
    h = np.array([[h1+0j, 0], [0, h2]])
    return (np.block([[R, np.dot(R, h)-np.dot(h, A)], [Z, A]]))
#this function inspects the Green function. Returns density of states, pseudophase, h1, h2, ...
def inspectG(G):
    dos=float(np.real(G[0,0])) #from BZ
    pseudophase =0.5*complex(1j*np.log(-G[0,1]/G[1,0])) #chi from BZ
    h1=float(np.real(G[0,2]/(G[0,0]-G[2,2]))) #h1= 1-2fel
    h2=float(np.real(G[1,3]/(G[1,1]-G[3,3]))) #h2= 2fh -1
    return (dos,pseudophase,0.5*(1-h1),0.5*(h2+1))
########
# main class
########
class SuperCircuit(Circuit):
    def __init__(self, energy=0, temperature=0.01):
        self.Ci = Circuit.__init__(self)
        self.energy = energy
        self.globaltemperature = temperature
        #adds the leakage terminal
        self.addTerminal("Leakage",tau3four())
        #adds a dictionary for leakage Connectors
        self.LeakageConnectors={}
######
# Section: building functions
##### 
    # adding the terminals
    # superconducting
    def addSuperTerminal(self, name, Delta=1.0, phase=0.0, temperature=None, Gamma=1e-4):
        if (temperature == None):
            te = self.globaltemperature
        else:
            te = temperature
        self.addTerminal(name, G_S(self.energy, temperature=te,
                         phase=phase, Delta=Delta, Gamma=Gamma))
        # attach properties to the terminal
        self.terminal(name).temperature = te
        self.terminal(name).Delta = Delta
        self.terminal(name).phase = phase
        self.terminal(name).Gamma = Gamma
    # normal
    def addNormalTerminal(self, name, voltage=0.0, temperature=None):
        if (temperature == None):
            te = self.globaltemperature
        else:
            te = temperature
        self.addTerminal(
            name, G_N(self.energy, temperature=te, voltage=voltage))
        # attach parameters to the terminal
        self.terminal(name).temperature = te
        self.terminal(name).voltage = voltage
    # adding a Node: rewrite the function from Circuit
    def addNode(self, name, G=None):
        if (G == None):  # adding normal G
            Circuit.addNode(self, name, G_N(
                self.energy, voltage=0, temperature=self.globaltemperature))
            # Circuit.addNode(self, name, G_S(self.energy,voltage=0,temperature=self.globaltemperature))
        else:
            Circuit.addNode(self, name, G)
    # adding common connectors goes through the Circuit interface
    #
    # adding leakage connector
    def addLeakageConnector(self, name, node, eth): #need eth
        self.addConnector(name,'NT','Leakage',node, 'Leakage',conductance=1j)
        c=self.connector(name)
        c.eth=eth
        self.setConductance(name, -1j*self.energy/eth)
        self.LeakageConnectors[name]=c
    # adding decoherence connector: not implemented in the version 1.0
    def addDecoherenceConnector(self, node, DefRate):
        pass
#######
# Section: Setting the elements functions
#######
    # this changes energy for ALL terminals. 
    #Makes sure that the energy is always the same for all terminals
    def setEnergy(self, en):
        self.energy = en
        for k,t in self.Terminals.items():
            if (k!="Leakage"): #skips leakage terminal
                if hasattr(t, 'phase'):  # superconducting terminal
                    t.setG(G_S(self.energy, temperature=t.temperature,
                           phase=t.phase, Delta=t.Delta, Gamma=t.Gamma))
                else:  # normal terminal
                    t.setG(G_N(self.energy, temperature=t.temperature, voltage=t.voltage))
        for c,v in self.LeakageConnectors.items():
            self.setConductance(c, -1j*self.energy/v.eth)
    # this changes temperature in all terminals PROVIDED it is set to the global temperature
    # if this is not the case, set the terminal temperature individualyy
    def setGlobalTemperature(self, te):
        for n,t in self.Terminals.items():
            if (not n=="Leakage"):
                if (t.temperature == self.globaltemperature):
                    if hasattr(t, 'phase'):  # superconducting terminal
                        t.setG(G_S(self.energy, temperature=te,
                               phase=t.phase, Delta=t.Delta, Gamma=t.Gamma))
                    else:  # normal terminal
                        t.setG(G_N(self.energy, temperature=te, voltage=t.voltage))
                    t.temperature=te
        self.globaltemperature=te
    #Changes the parameters of supterminal
    def setSuperTerminal(self, name, Delta=None, Gamma=None, temperature=None, phase=None):
        t = self.terminal(name)
        if (Delta == None):
            pass
        else:
            t.Delta = Delta
        if (Gamma == None):
            pass
        else:
            t.Gamma = Gamma
        if (temperature == None):
            pass  # its own temperature, NOT the global temperature
        else:
            t.temperature = temperature
        if (phase == None):
            pass  # its own temperature, NOT the global temperature
        else:
            t.phase = phase
        t.setG(G_S(self.energy, temperature=t.temperature,
               phase=t.phase, Delta=t.Delta, Gamma=t.Gamma))
    #Changes the parameters of supterminal
    def setNormalTerminal(self, name, temperature=None, voltage=None):
        t = self.terminal(name)
        if (voltage == None):
            pass
        else:
            t.voltage = voltage
        if (temperature == None):
            pass  # its own temperature, NOT the global temperature
        else:
            t.tempetarure = temperature
        t.setG(G_N(self.energy, temperature=t.temperature, voltage=t.voltage))
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
################
# Getting results section
################
    #inspects a node, returns dos, pseudophase, h1,h2 in a tuple
    def inspectNode(self,name):
        G=self.node(name).getG()
        return(inspectG(G))
    #returns current densities at certain energy as np array
    def currentsdensities(self, energy):
        self.setEnergy(energy)
        self.solve()
        di = []
        for t in self.Terminals.keys():
            if (t!="Leakage"):
                di.append( elcurrent(self.getTotalCurrentTerminal(t)))
        return (np.array(di))
    #this is a function used in fast dc
    def cdHfixed(self, name, sg, h1, h2):
        t = self.terminal(name)
        assert (hasattr(t, 'voltage')), "This is not a normal terminal"
        self.setEnergy(sg*np.abs(t.voltage))
        t.setG(G_d(h1, h2))
        self.solve()
        li = []
        for te in self.Terminals.keys():
            if (te!="Leakage"):
                li.append(elcurrent(self.getTotalCurrentTerminal(te)))
        return np.array(li)
# this returns changes of currents upon changing the voltage at terminal 'name'
# assuming this terminal has zero temperature
    def fastDifferentialConductances(self, name):
        if (self.terminal(name).voltage > 0):
            dcs = self.cdHfixed(name,1,1,-1)+self.cdHfixed(name,-1,1,-1)-self.cdHfixed(name,1,1,1)-self.cdHfixed(name,-1,-1,-1)
        else:
            dcs =-self.cdHfixed(name,1,-1,1)-self.cdHfixed(name,-1,-1,1)+self.cdHfixed(name,1,1,1)+self.cdHfixed(name,-1,-1,-1)
        return(list(0.5*dcs))
# here are integration methods
    #all together
    def FullCurrents(self,method='trans',scale=1.0,numpoints=200,minen=0,maxen=5.0,acc=1e-2):
        if(method=='trans'):
            return(self.FullCurrents1(minen,numpoints,scale))
        elif (method=='straight'):
            return(self.FullCurrents0(minen,numpoints,maxen))
        elif (method=='adaptive'):
            return(self.FullCurrents2(scale,acc))
        else:
            print("Unknown method in FullCurrents!!! Use trans, straight, adaptive")
    #simplest addition
    def FullCurrents0(self,minen, numpoints,maxen):
        energies=np.linspace(minen,maxen,numpoints)
        de=maxen/(numpoints-1)
        ans=np.zeros(len(self.Terminals)-1)
        for en in energies:
            ans+=self.currentsdensities(en)
        return(ans*de)
    # integrates with transformation of energy axis
    # t(e) = (e/scale)/1+(e/scale)
    # switches direction of integration at any call - to accelerate convergence
    def FullCurrents1(self,minen,numpoints,scale):
        if (not hasattr(self.FullCurrents1.__func__,'dir')):
            self.FullCurrents1.__func__.dir=1
        if(self.FullCurrents1.__func__.dir ==1):  #direction of integration  
            ts=np.linspace(minen,1.0,numpoints,endpoint=False)
        else:
            ts=np.linspace(1.0-1/numpoints,minen,numpoints)
        self.FullCurrents1.__func__.dir = - self.FullCurrents1.__func__.dir
        #print("dir ",self.FullCurrents1.__func__.dir)
        dt=1.0/(numpoints-1) 
        ans=np.zeros(len(self.Terminals)-1)
        for t in ts:
            ans+=(scale/(1-t)**2)*self.currentsdensities(scale*t/(1-t))
        return(ans*dt)   
    #adaptive method using scipy. Usually much slower... Prone to divergencies
    def FullCurrents2(self, minen,scale,acc):
        def f(t):
            if (t==1.0):
                return(np.zeros(len(self.Terminals)-1))
            else: 
                return(scale/(1-t)**2*np.array(list(self.currentsdensities(scale*t/(1-t)).values())))
        return(quad_vec(f,minen,1.0,quadrature='trapezoid',epsrel=acc)[0])   
######
#end of the SuperCircuit class
######