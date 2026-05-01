# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 12:11:51 2026
Linear general circuit theory
Zero version 
Pre-code notes:
    ----------------
1. Nodes and terminals: characterized by voltage vectors of various dimensions
2. Connectors: different currents on both ends, proportional to voltages on both ends
the proportionality matrices are (re)computed using  the "pointer" to the connector
3. When conductances are (re)computed, they are rewritten to a big matrix
4. This matrix is used to compute the matrix of general conductances
5. And that of electric conductances
6. Has to be a function to return currents at given terminal voltages
7. Next, the parameters of connectors and terminals are changed, and the cycle 
repeats
8-01-2026 look at array linking. No definite solution yet. Better if
only one connector per pair of nodes/terminals. But this seems restrictive. 16:09
added draft element classes
9-01-2026 did some work - skeleton is kinda ready 15:47 
10-01-2026 21:40 all right, I need to decide something...
computeConductances( connector) - returns the conductance matrices
it is natural to make it a part of connector class?
storeConductances: backs up the old values puts new values 
putConductances: puts the conductances into the input matrices 
11-01-2026 new code 'N' for connectors
Need to add connector lists to terminals: only for update! 21:59
now need to write set functions: they take care of updates
22:26 ok, I guess I'm done with the draft class
next step: subclass trivial circuit
12-01-2026 afrer thinking a while, I think I was wrong and need to remade the 
class. 1. function computeConductance must be in the main class: otherwise subclassing does not work
2. put _ for the functions not called by the user: easy to subclass
12-01-2026 ok, I think I'm done with this 12:48 subclassing...
14:16 many errors. messed up with indexing, do not understand... Will do something else.
15:44 finally, got it
13-01-2026 adding addressing functions
22-01-2026 adding calculations of currents, gencurrents, genvoltages
and helper functions
23-01-2026 wiil try complex or real types, 
and multipole
Is this like version 2.0? 
Need to check the multipole 
I checheck the mutlipole, it works: dissapointing result, should think of alternative Hall circuit.
Meanwhile, I need to correct implementation of multipole
28-2-2026 fixed a bug in one-node connectors
2-3-2026 added a connectivity check
@author: dhyul
"""
import numpy as np
from connectivitycheck import connectivityReport
# a class for convenient syntax
class AttrDict(dict):
    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        self[key] = value
######
#
#Section: functions to work with the dictionaries
#
######
def addtodic(dic,name,thing):
    assert (name not in dic.keys()), "Name "+str(name)+" is already taken"
    dic[name]=thing
def checkdic(dic,name,dicname):
    assert (name in dic.keys()), "No "+str(name)+" in dictionary "+str(dicname)
def getfromdic(dic,name,dicname):
    checkdic(dic,name,dicname)  
    return(dic[name]) 
#########
#
#Section: micsellaneous functions
########
#for resistance calculation
def minor(arr,i,j):
    # returns an array with ith row, jth column removed needed for resistance
    return arr[np.array(list(range(i))+list(range(i+1,arr.shape[0])))[:,np.newaxis],
               np.array(list(range(j))+list(range(j+1,arr.shape[1])))]
########
#
#Section: element classes
#
########
class Node:
    def __init__(self,dim):
        self.dim=dim
        self.offset=0 #to put data to
        self.connectors=[]
class Terminal:
    def __init__(self,dim,mask):
        self.dim=dim
        self.offset=0
        self.mask=mask
        self.connectors=[]
class HalfConnector:
    def __init__(self,Head,Tail,dt):
        self.head=Head # this is surely present
        #place for conductance matrices
        self.GHH=np.zeros((self.head.dim, self.head.dim),dtype=dt)
        self.GHHold=np.zeros((self.head.dim, self.head.dim),dtype=dt) 
        self.tail=Tail
        if not self.tail == None: #if there is a tail
            self.GHT=np.zeros((self.head.dim, self.tail.dim),dtype=dt)
            self.GHTold=np.zeros((self.head.dim, self.tail.dim),dtype=dt)
    def backup(self):
        self.GHHold[0:,0:]=self.GHH
        if not self.tail == None:
            self.GHTold[0:,0:]=self.GHT
#Connector        
class Connector:
    def __init__(self,code,A,B,pars,dt):
        self.code=code
        if code =='N':
            self.AA=HalfConnector(A,None,dt)
            A.connectors.append(self) #perhaps do not need this if relaxation does not depend on the node state
        else:
            self.AB=HalfConnector(A,B,dt)
            self.BA=HalfConnector(B,A,dt)
            A.connectors.append(self)
            B.connectors.append(self)
        self.pars=pars
    #def computeConductances(self):
    #    pass
    #returns four matrices
    def storeConductances(self,tu):
        #compute conductances
        #AA,AB,BA,BB=self.computeConductances()
        #backup
        AA,AB,BA,BB=tu
        if self.code =='N':
            self.AA.backup()
            self.AA.GHH[0:,0:]=AA
        else:
            self.AB.backup()
            self.BA.backup()
            self.AB.GHH[0:,0:]=AA
            self.AB.GHT[0:,0:]=AB
            self.BA.GHH[0:,0:]=BB
            self.BA.GHT[0:,0:]=BA
            
class Multipole():
    def __init__(self,code,poles,pars,dim,dt): #dim is computed when adding the multipole
        self.code=code
        self.poles=poles
        self.pars=pars
        #reserve the space for matrices
        self.matrix=np.zeros((dim,dim),dtype=dt)
        self.oldmatrix=np.zeros((dim,dim),dtype=dt)
 ###################       
    def storeConductances(self,newmatrix):
        self.oldmatrix[0:,0:]=self.matrix
        self.matrix[0:,0:]=newmatrix
        
############        
# main class
############
class LinGenCircuit():
    #initialization
    def __init__(self,dt=np.float64):
        self.numbers,self.dics=AttrDict(),AttrDict()
        sorts=['nodes','terminals','connectors']
        for s in sorts:
            setattr(self.numbers, s, 0)
            setattr(self.dics,s,{})
        self.ndim,self.tdim=0,0
        self.finalized=False
        self.dt=dt
    #adding functions
    def _addTerminal(self,name,dim,mask,pars):
        assert(len(mask)==dim),"Terminal "+name+": wrong mask size!"
        t=Terminal(dim,mask)
        t.pars=pars
        t.offset=self.tdim
        self.tdim+=dim
        t.num=self.numbers.terminals
        self.numbers.terminals+=1
        #add to dic
        addtodic(self.dics.terminals,name,t)
    def _addNode(self,name,dim,pars):
        n=Node(dim)
        setattr(n,'pars',pars)
        n.offset=self.ndim
        self.ndim+=dim
        n.num=self.numbers.nodes
        self.numbers.nodes+=1
        #add to dic
        addtodic(self.dics.nodes,name,n)        
    def _addConnector(self,name,code,nameA,nameB,pars):
        #check if nodes/terminals are there
        if code =='NN':
            a,b=self.dics.nodes,self.dics.nodes
            an,bn= 'Nodes',"Nodes"
        elif code == 'NT':
            a,b=self.dics.nodes,self.dics.terminals
            an,bn= 'Nodes',"Terminals"
        elif code == 'TT':
            a,b=self.dics.terminals,self.dics.terminals
            an,bn= 'Terminals',"Terminals"
        elif code == 'N':
            a,an=self.dics.nodes,"Nodes"
        else:
            assert False, "No such code"
        A=getfromdic(a, nameA, an)
        if not code=='N':
            B=getfromdic(b,nameB, bn)
            c=Connector(code, A, B, pars,self.dt)
        else:
            c=Connector(code, A, None, pars,self.dt)
        c.num=self.numbers.connectors
        self.numbers.connectors+=1
        #add to dic
        addtodic(self.dics.connectors,name,c)
    #adds
    def _addMultipole(self,name,code,elements,pars):
        assert(len(code)==len(code)),"Cannot add multpole "+name+" mismatch of the number of poles"
        poles=[]
        offs=[]
        dim=0
        for i in range(len(code)):
            if code[i]=='N':
                a,an=self.dics.nodes,"Nodes"
            elif code[i]=='T':
                a,an=self.dics.terminals,"Terminals"
            else:
                assert False, "Cannot add multpole "+name+". Wrong code!"
            E=getfromdic(a,elements[i], an)
            offs.append(dim)
            poles.append(E)
            dim+=E.dim
        #create
        m=Multipole(code,poles,pars,dim,self.dt) #not sure if code is needed
        m.offs=offs
        #multipole is counted as a connector        
        m.num=self.numbers.connectors
        self.numbers.connectors+=1
        #add to dic
        addtodic(self.dics.connectors,name,m)
    #######   
    #operation functions
    #######
    #finalize
    def finalize(self):
        if self.finalized:
            print("the circuit is already finalized!")
            return(0)
        #reserve space for matrices - now we know the dimensions
        #input matrices
        self.GNN=np.zeros((self.ndim,self.ndim),dtype=self.dt)
        self.GNT=np.zeros((self.ndim,self.tdim),dtype=self.dt)
        self.GTN=np.zeros((self.tdim,self.ndim),dtype=self.dt)
        self.GTT=np.zeros((self.tdim,self.tdim),dtype=self.dt)
        self.mask=np.zeros((self.numbers.terminals,self.tdim),dtype=self.dt)
        #now make the mask
        for t in self.dics.terminals.values():
            self.mask[t.num,t.offset:(t.offset+t.dim)]=t.mask
        #output matrices
        self.GeneralConductanceMatrix=np.zeros((self.tdim,self.tdim),dtype=self.dt)
        self.ConductanceMatrix=np.zeros((self.numbers.terminals,self.numbers.terminals),dtype=self.dt)
        #reserves space for voltages, currents, nodevoltages
        self.Voltages=np.zeros(self.numbers.terminals,dtype=self.dt)
        self.Voltages[0]=1.0
        self.Currents=np.zeros(self.numbers.terminals,dtype=self.dt)
        self.GenCurrents=np.zeros(self.tdim,dtype=self.dt)
        self.NodeVoltages=np.zeros(self.ndim,dtype=self.dt)       
        #fills the input matrices for the first time
        #print(self.dics.connectors)
        for n,c in self.dics.connectors.items():
            c.storeConductances(self.computeConductances(c))
            self._putConductances(c)
        #marks finalized
        self.finalized=True
        #computes the output matrices 
        self.compute()
    #compute
    def compute(self):
        assert (self.finalized),"You have to finalize the circuit first!"
        if (self.ndim >0):
            if np.linalg.cond(self.GNN) < 1e10:
                inv = np.linalg.inv(self.GNN)
            else:
                inv= np.linalg.inv(self.GNN+1e-12*np.diag(np.ones(self.ndim)))
        #handle it
            self.GeneralConductanceMatrix[0:,0:]= self.GTT-self.GTN@inv@self.GNT
        else:
            self.GeneralConductanceMatrix[0:,0:]= self.GTT
        self.ConductanceMatrix[0:,0:]=self.mask@self.GeneralConductanceMatrix@np.transpose(self.mask)        
    #computes currents given the voltages
    def computeCurrents(self):
        self.Currents[0:]=np.dot(self.ConductanceMatrix,self.Voltages)
        self.GenCurrents[0:]=np.dot(self.GeneralConductanceMatrix @ np.transpose(self.mask) ,self.Voltages)
    def computeNodeVoltages(self):
        self.NodeVoltages[0:]=-np.dot(np.linalg.inv(self.GNN)@ self.GNT @ np.transpose(self.mask),self.Voltages)
    #fills the input matrices for the first time
    def _putConductances(self,c):
        if type(c)==Connector:
            if c.code=='N':
                a,ad=c.AA.head.offset,c.AA.head.dim
            else:
                a,ad=c.AB.head.offset,c.AB.head.dim
                b,bd=c.BA.head.offset,c.BA.head.dim     
            if c.code =='TT':
                #GTT is the target matrix
                self.GTT[a:(a+ad),a:(a+ad)]+=c.AB.GHH-c.AB.GHHold
                self.GTT[b:(b+bd),b:(b+bd)]+=c.BA.GHH-c.BA.GHHold
                self.GTT[a:(a+ad),b:(b+bd)]+=c.AB.GHT-c.AB.GHTold
                self.GTT[b:(b+bd),a:(a+ad)]+=c.BA.GHT-c.BA.GHTold            
            elif c.code == 'NT':
                # a is a node, b is a terminal
               self.GNN[a:(a+ad),a:(a+ad)]+=c.AB.GHH-c.AB.GHHold
               self.GTT[b:(b+bd),b:(b+bd)]+=c.BA.GHH-c.BA.GHHold
               self.GNT[a:(a+ad),b:(b+bd)]+=c.AB.GHT-c.AB.GHTold
               self.GTN[b:(b+bd),a:(a+ad)]+=c.BA.GHT-c.BA.GHTold      
            elif c.code == 'NN':
                #GNN is the target matrix
                self.GNN[a:(a+ad),a:(a+ad)]+=c.AB.GHH-c.AB.GHHold
                self.GNN[b:(b+bd),b:(b+bd)]+=c.BA.GHH-c.BA.GHHold
                self.GNN[a:(a+ad),b:(b+bd)]+=c.AB.GHT-c.AB.GHTold
                self.GNN[b:(b+bd),a:(a+ad)]+=c.BA.GHT-c.BA.GHTold  
            elif c.code=='N': #only GNN
                self.GNN[a:(a+ad),a:(a+ad)]+=c.AA.GHH-c.AA.GHHold
            else:
                assert False, "No such code" #have to add the code
        elif type(c)==Multipole:
            for i in range(len(c.code)):
                for j in range(len(c.code)):
                    localcode=c.code[i]+c.code[j] 
                    if(localcode=='NN'):
                        toadd=self.GNN
                    elif (localcode=='NT'):
                        toadd=self.GNT
                    elif (localcode=='TN'):
                        toadd=self.GTN
                    elif (localcode=='TT'):
                        toadd=self.GTT
                    else:
                        assert False, "Problem with Multipole code"
                    a,ad,ase=c.poles[i].offset,c.poles[i].dim,c.offs[i]
                    b,bd,bse=c.poles[j].offset,c.poles[j].dim,c.offs[j]
                    toadd[a:(a+ad),b:(b+bd)]+=c.matrix[ase:(ase+ad),bse:(bse+bd)]-c.oldmatrix[ase:(ase+ad),bse:(bse+bd)]
        #print("GNN",self.GNN)
           
    def computeConductances(self,connector):
        pass
    ####
    #setting functions
    ####
    def _setConnector(self,name,pars):
        c=getfromdic(self.dics.connectors,name,"Connectors")
        c.pars=pars
        c.storeConductances(self.computeConductances(c)) #update
        self._putConductances(c)
    def _setNode(self,name,pars):
        n=getfromdic(self.dics.nodes,name,"Nodes")
        n.pars=pars
        for c in n.connectors:
           c.storeConductances(self.computeConductances(c)) #update
           self._putConductances(c)
    def _setTerminal(self,name,pars):
        t=getfromdic(self.dics.terminals,name,"Terminals")
        t.pars=pars
        for c in t.connectors:
           c.storeConductances(self.computeConductances(c)) #update
           self._putConductances(c)
    ##this is somewhat special: sets voltages. accepts sequencies
    def setVoltage(self,names,values):
        if(type(names)==str):
            self.Voltages[self.terminal(names).num]=values
        else: #assuming names = sequence
            assert(len(names)==len(values)), "in setVoltage: wrong number of values"
            for i in range(len(names)):
                self.Voltages[self.terminal(names[i]).num]=values[i]
                
    ####
    #addressing functions:return pointer by name
    ####        
    def terminal(self,name):
         return(getfromdic(self.dics.terminals,name,"Terminals"))
    def node(self,name):
         return(getfromdic(self.dics.nodes,name,"Nodes"))
    def connector(self,name):
         return(getfromdic(self.dics.connectors,name,"Connectors"))       
    ###
    #helper functions
    ###        
       #element of conductance matrix by names
    def CMel(self,name1,name2):
        t1,t2=self.terminal(name1),self.terminal(name2)
        return(self.ConductanceMatrix[t1.num,t2.num])
        #submatrix of general conductance matrix by names
    def GCMsub(self,name1,name2):
        t1,t2=self.terminal(name1),self.terminal(name2)
        return(self.GeneralConductanceMatrix[t1.offset:(t1.offset+t1.dim),t2.offset:(t2.offset+t2.dim)])           
    def current(self,name):
        return(self.Currents[self.terminal(name).num])
    def gencurrent(self,name):
        t=self.terminal(name)
        return(self.GenCurrents[t.offset:(t.offset+t.dim)])
    def nodevoltage(self,name):
        n=self.node(name)
        return(self.NodeVoltages[n.offset:(n.offset+n.dim)])
    def resistanceMatrix(self,ground='Ground'):
        assert (ground in self.dics.terminals), "There is no terminal "+ground+"Are you trying to compute resistance?" 
        i=self.terminal(ground).num
        return(-np.linalg.inv(minor(self.ConductanceMatrix, i, i)))
    #this arranges the dics for the connectivitycheck
    def _arrangethedics(self):
        reversedConnectors= {v: k for k, v in self.dics.connectors.items()}
        reversedNodes= {v: k for k, v in self.dics.nodes.items()}
        reversedTerminals={v: k for k, v in self.dics.terminals.items()}
        def lists(c): #returns names of nodes and terminals connected to a connector
            if (type(c)==Connector):
                if c.code =='N':
                    return ([],reversedNodes[c.AA.head])
                elif c.code=='NN':
                    return ([],[reversedNodes[c.AB.head],reversedNodes[c.BA.head]])
                elif c.code=='TT':
                    return ([reversedTerminals[c.AB.head],reversedTerminals[c.BA.head]],[])
                elif c.code=='NT':
                    return ([reversedTerminals[c.BA.head]],[reversedNodes[c.AB.head]])
                else:
                    assert False, 'Wrong code connector '+reversedConnectors[c]
            elif (type(c)==Multipole):
                t,n=[],[]
                for ch, p in zip(c.code,c.poles):
                    if ch == 'T':
                        t.append(p)
                    elif ch == 'N':
                        n.append(p)
                return(t,n)
        
            
        T = {reversedTerminals[t]: [reversedConnectors[c] for c in t.connectors] for t in reversedTerminals.keys()}
        N = {reversedNodes[n]: [reversedConnectors[c] for c in n.connectors] for n in reversedNodes.keys()}
        C = {reversedConnectors[c]: [*(lists(c))] for c in reversedConnectors.keys()}
        return(T,N,C)
    def reportConnectivity(self):
        print(connectivityReport(*(self._arrangethedics())))

    
