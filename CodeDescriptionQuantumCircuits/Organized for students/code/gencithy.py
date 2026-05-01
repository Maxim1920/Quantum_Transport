# -*- coding: utf-8 -*-
"""
@author: Yuli
version 1.0 as of 1-1-2026
2-02-2026 added action for the connector.
25-02-2026 changed definition of conductance of an arbitrary connector
"""
import time #needed for time management
import numpy as np
from  scipy.linalg import eig
from connectivitycheck import connectivityReport

#a small class for solve parameters
class SolvePars:
    def __init__(self):
        self.verbose=False
        self.maxiter=1000
        self.goal=1E-8
        self.miniter=3
#alternative function of matrix to suppress myfunm warnings
def myfunm(A,f):
    ei,left,right=eig(A,left=True,right=True)
    newei=f(ei)
    return(right @ np.diag(newei) @ np.linalg.inv(right))
#########
#
#The section of usefull functions
#
########
def tripleMmul(A,B,C): #product of three matrices
    return (np.dot(A,np.dot(B,C)))
def fornormalize(x): #needed for normalization
    return(np.sign(x.real))
def normalize(M): #normalization of G's
    return(myfunm(M,fornormalize))
def commutator(A,B): #commutator of two matrices
    return(np.dot(A,B)-np.dot(B,A))
def anticommutator(A,B):#anticommutator of two matrices
    return(np.dot(A,B)+np.dot(B,A))
######
#
#Section functions for connector types: arbirtrary, diffusive, tunnel, 
#ballistic  
#
#######
def tcurrent(t,x): #needed for an arbitrary connector 
    return(t/(1+t*(x-2)/4.0)) 
def taction(t,x): #needed for an arbitrary connector
    return 0.5*np.log(1+0.25*t*(x-2)) 
def Arbitraryfunctie(x,p): #current for an arbitrary connector
    #return(np.sum(tcurrent(np.array(p.transmissions),x)))
    #does not work for some reasons making it plain 9-7-2024
    s=0
    for t in p.transmissions:
        s+=tcurrent(t,x)
    return(s)
def Arbitraryactie(x,p): #action for an arbitrary connector
        #return(0.5*np.sum(vtaction(p.transmissions,x)))
    s=0
    for t in p.transmissions:
        s+=taction(t,x)
    return(s)
def Ballisticfunctie(x,p): #current for a ballistic connector
    return(4.0/(x+2.0)) 
def Ballisticactie(x,p): #action for a ballistic connector
    return(0.5*np.log(1+0.25*(x-2.0)))

def dddf(x):
    if (x==2.0):
        return(1.0)
    if(x==-2.0):
        return(0.0)
    return(np.arccos(x/2)/np.sqrt(1-x*x/4.0))

ddf=np.vectorize(dddf)

def Diffusivefunctie(x,p): #current for diffusive connector: changed 25-12-2025
    return(ddf(x))

def Diffusiveactie(x,p): #action for diffusive connector added minus 03-02-2026
    return(-0.125*(np.arccos(x/2))**2)
def Tunnelactie(x,p): #action for tunnel connector
        return(0.125*(x-2.0))
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

#######
#
#Section: classes of circuit elements: HalfConnector, Connector, Terminal, Node and their
#derivatives.
#
#######

#this is the most general class incorporating all functionality of the connectors
class HalfConnector:
    def __init__(self,Hnode,Tnode,conductance=1.0,tunnel=True,U=None,V=None,functie=None,actie=None,transmissions=None):
        self.Hnode=Hnode
        self.Tnode=Tnode
        self.conductance=conductance
        self.transmissions=transmissions
        self.tunnel=tunnel #need to treate a tunnel connector separately
        self.U=U
        self.V=V
        self.functie=functie
        self.actie=actie
        #now connect to the Node
        Hnode.HCList.append(self)
        #####
    def tailG(self): #this returns tail Green function, if phase factor - incoporates it
        if (self.U is not None):
            return(tripleMmul(self.U,self.Tnode.G,self.V))
        else:
            return(self.Tnode.G)
    def AC(self): #makes relevant anticommutator
        return(self.conductance*(myfunm(anticommutator(self.Hnode.G,self.tailG()),lambda x: self.functie(x,self))))
    def give(self): #main function, most important for iterations
        if self.tunnel:
            return(self.conductance*self.tailG())
        else:
            return(np.dot(self.tailG(),self.AC()))
    def matrixcurrent(self): #returns the matrix current
        if self.tunnel:
            return(self.conductance*commutator(self.Hnode.G,self.tailG()))
        else:
            return(np.dot(commutator(self.Hnode.G,self.tailG()),self.AC()))
    def action(self): #computes the action. it is eventually The FULL action of the connector...
        x=np.linalg.eigvals(anticommutator(self.Hnode.G,self.tailG()))
        return(self.conductance*np.sum(self.actie(x,self)))
    
#Connector: mostly to keep the track of half-connectors, and define connector types    
class Connector: 
    def __init__(self,contype,Hnode,Tnode,conductance=1.0,U=None,V=None,transmissions=None):
        if (contype=='Dephasing'): #here we make only one half-connector. Tnode ignored
            assert(U is not None),"U should be defined"
            assert(V is not None),"V should be defined"
            self.H=HalfConnector(Hnode,Hnode,conductance=conductance,U=U,V=V,actie=Tunnelactie)
            self.T=None
        elif (contype=='Leakage'): #should be connected to a leakage terminal
            self.H=HalfConnector(Hnode,Tnode,conductance=conductance,actie=Tunnelactie)
            self.T=None
        elif (contype=='Tunnel'): #tunnel contact. can contain the phase
            self.H=HalfConnector(Hnode,Tnode,conductance=conductance,U=U,V=V,actie=Tunnelactie)
            self.T=HalfConnector(Tnode,Hnode,conductance=conductance,U=V,V=U,actie=Tunnelactie)
        elif (contype=='Ballistic'):
            self.H=HalfConnector(Hnode,Tnode,conductance=conductance,tunnel=False,U=U,V=V,functie=Ballisticfunctie,actie=Ballisticactie)
            self.T=HalfConnector(Tnode,Hnode,conductance=conductance,tunnel=False,U=V,V=U,functie=Ballisticfunctie,actie=Ballisticactie)           
        elif (contype=='Diffusive'):
            self.H=HalfConnector(Hnode,Tnode,conductance=conductance,tunnel=False,U=U,V=V,functie=Diffusivefunctie,actie=Diffusiveactie)
            self.T=HalfConnector(Tnode,Hnode,conductance=conductance,tunnel=False,U=V,V=U,functie=Diffusivefunctie,actie=Diffusiveactie)           
        elif (contype=='Arbitrary'):
            assert (np.any(transmissions)), "A set of transmissions should be supplied for an arbitrary contact"
            #insert 25-02-2026
            c=conductance/(np.sum(np.array(transmissions)))
            self.H=HalfConnector(Hnode,Tnode,conductance=c,tunnel=False,U=U,V=V,functie=Arbitraryfunctie,actie=Arbitraryactie,transmissions=transmissions)
            self.T=HalfConnector(Tnode,Hnode,conductance=c,tunnel=False,U=V,V=U,functie=Arbitraryfunctie,actie=Arbitraryactie,transmissions=transmissions)              
        else: 
            assert False, str(contype)+": Unknown connector type!"       
    def setconductance(self,cond):
        #insert 25-02-2026
        if (np.any(self.H.transmissions)):
            c=cond/np.sum(np.array(self.H.transmissions))
        else:
            c=cond
        self.H.conductance=c
        if(self.T):
            self.T.conductance=c
    def matrixcurrent(self):
            if (self.T):
                return( self.H.matrixcurrent(),self.T.matrixcurrent())
            else:
                return(self.H.matrixcurrent())
    def action(self):
        return(self.H.action())  #02-02-2026 not quite sure about it in general...
        
#class Terminal: storing G not subject to iterations
class Terminal:
    def __init__(self,G):
        self.setG(G)
        self.HCList=[]
    def setG(self,G):
        self.G=G.copy()
    def totalcurrent(self):
        return(sum([x.matrixcurrent() for x in self.HCList]))
    def getG(self):
        return(self.G)
#class Node: storing G, making the iterations 
#adjustment coefficient is individual for each node. 
# change if there are problems with convergence   
class Node(Terminal):
    def __init__(self,G):
        Terminal.__init__(self, G)
        self.M=G.copy()
        self.oldaccu=1.0
        self.accu=1.0
        self.adjust=0.5
    def makeM(self):
        self.M=normalize(sum([x.give() for x in self.HCList]))
        self.oldaccu=self.accu 
        self.accu = np.linalg.norm(commutator(self.M,self.G))
        return(self.accu)
    def iterate(self):
        self.G=normalize((1-self.adjust)*self.G+self.adjust*self.M)
#######
#        
# The main class: Circuit        
# 
#######       
class Circuit:
    def __init__(self):
        self.Nodes={}
        self.Terminals={}
        self.Connectors={}
        #parameters for solving
        self.solvepars=SolvePars()
        
#building functions
    def addNode(self,name,G):
       addtodic(self.Nodes,name,Node(G))
    def addTerminal(self,name,G):
        addtodic(self.Terminals,name,Terminal(G))
    def addConnector(self,name,code,contype,Hname,Tname,conductance=1.0,U=None,V=None,transmissions=None):
        if code =='NN':
            a,b=self.Nodes,self.Nodes
            an,bn= 'Nodes',"Nodes"
        elif code == 'NT':
            a,b=self.Nodes,self.Terminals
            an,bn= 'Nodes',"Terminals"
        elif code == 'TT':
            a,b=self.Terminals,self.Terminals
            an,bn= 'Terminals',"Terminals"
        else:
            assert False, "No such code"
        Hnode=getfromdic(a, Hname, an)
        Tnode=getfromdic(b,Tname, bn)
        addtodic(self.Connectors,name,Connector(contype,Hnode,Tnode,conductance=conductance,U=U,V=V,transmissions=transmissions))
#checking function    
    def check(self): 
        #need to implement a better checking procedure...
        #Find disconnected nodes
        for key in self.Nodes:
            if len(self.Nodes[key].HCList)==0 :
                print("Node "+str(key)+" is disconnected!")
        #Find disconnected terminals
        for key in self.Terminals:
            if len(self.Terminals[key].HCList)==0 :
                print("Terminal "+str(key)+" is disconnected!")        
     #formats the dictionaries into strings
    def _formatthedics(self):
         reversedTerminals = {v:k for k, v in self.Terminals.items()}
         reversedNodes= {v:k for k, v in self.Nodes.items()}
         reversedConnectors= {v:k for k, v in self.Connectors.items()}
         T={reversedTerminals[x]:[] for x in reversedTerminals.keys()}
         N={reversedNodes[x]:[] for x in reversedNodes.keys()}
         C={reversedConnectors[x]:[] for x in reversedConnectors.keys()}
         for c in reversedConnectors.keys():
             head = c.H.Hnode
             tail = c.H.Tnode
             te,ne =[],[]
             if (head in reversedTerminals):
                 T[reversedTerminals[head]].append(reversedConnectors[c])
                 te.append(reversedTerminals[head])
             elif (head in reversedNodes):
                 N[reversedNodes[head]].append(reversedConnectors[c])
                 ne.append(reversedNodes[head])
             else:
                 assert False, "Mistake in the dictionaries, connector "+str(reversedConnectors[c])
             if (head !=tail):
                 if (tail in reversedTerminals):
                     T[reversedTerminals[tail]].append(reversedConnectors[c])
                     te.append(reversedTerminals[tail])
                 elif (tail in reversedNodes):
                     N[reversedNodes[tail]].append(reversedConnectors[c])
                     ne.append(reversedNodes[tail])
                 else:
                     assert False, "Mistake in the dictionaries, connector "+str(reversedConnectors[c])
             C[reversedConnectors[c]]=[te,ne].copy()
         return(T,N,C)
    def reportConnectivity(self):
        print(connectivityReport(*(self._formatthedics())))
             
            
#the functions for solving the circuit 
    def iteration(self):
        for key in self.Nodes:
            self.Nodes[key].iterate()
    def check_consistence(self):
        a=0.0
        for key in self.Nodes:
           a+= self.Nodes[key].makeM()
        return(a)
    def solve(self): #change 24-12: does not rum if no nodes 
        #verbose: levels
        # 0: no warning 1: warning if not converged 2: reports results of each run
        #3 or more: reports results of each iteration
        if (len(self.Nodes)==0):
            return(None)
         #time management, number of iterations, acc achieved
        if(self.solvepars.verbose>0):
            start_time=time.time()
            if(self.solvepars.verbose>2):
                print("Start iterations")
        for i in range(self.solvepars.maxiter):
            a = self.check_consistence()
            if(self.solvepars.verbose>2):
                print("Iteration: "+str(i)+" Achieved: "+str(a)+" Elapsed: "+str(time.time()-start_time) )                   
            if ((a>self.solvepars.goal)or (i<self.solvepars.miniter)): #goal not reached
             #here it would be time to adjust the adjustable parameters
                self.iteration()
            else: 
                    break
                
        if(self.solvepars.verbose>0):
            if(i==(self.solvepars.maxiter-1)):
                print("Solve report: bad. Goal not achieved. Accuracy :"+str(a))
                self.badmoment() #to be implemented in subclasses
            else:
                if(self.solvepars.verbose>1):
                    print("Solve report: good.Goal is achieved.")
            if(self.solvepars.verbose>1):
                print("Iterations: "+str(i)+" Time:"+ str(time.time()-start_time))
# the function to handle non-convergence. Implemented in subclasses
    def badmomemt(self):
         pass            
#addressing functions: return a pointer to an element
    def terminal(self,name):
        return(getfromdic(self.Terminals,name,"Terminals"))
    def node(self,name):
        return(getfromdic(self.Nodes,name,"Nodes"))
    def connector(self,name):
        return(getfromdic(self.Connectors,name,"Connectors"))
#setting functions
    def setGTerminal(self,name,G):
        self.terminal(name).setG(G)
    def setGNode(self,name,G):
        self.node(name).setG(G)
    def setConductance(self,name,c):
        self.connector(name).setconductance(c)
#getting functions 
    def getTotalCurrentTerminal(self,name):
        return(self.terminal(name).totalcurrent())
    def getGTerminal(self,name):
        return(self.terminal(name).getG())
    def getGNode(self,name):
        return(self.node(name).getG())
    def getMatrixCurrent(self,name):
        return(self.connector(name).matrixcurrent())
    def getTotalAction(self):
        a=0
        for cname in self.Connectors:
            a+=self.Connectors[cname].action()
        return(a)
        
#####
#end of the circuit class
####    
    