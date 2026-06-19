# -*- coding: utf-8 -*-
"""
07-03-2026 
Spintronics version 2.0
"""
import numpy as np
from lingencthy import LinGenCircuit

def normagnetization(m):
    length=np.linalg.norm(np.array(m))
    assert length>0, "zero vector lenght!"
    return(np.array(m)/length)
def le(a):
    if (type(a)==float):
        return (1)
    else:
        return(len(a))

#print(normagnetization([1.0,1.0,1.0]))

#dictionary for number of parameters in a conductor:
howManyPars={'NN':1, 'NF':4, 'FN':4, 'HN':3, 'NH':3, 'FF':6, 'HH': 1, 'FH':3, 'HF':3, 'RF':1, 'RN':1, 'P':4 }

##################
#Section: functions for conductance
##################

#contuctance matrices of FF connector
def coFF(m1,m2,pars):
    #m1,m2 = 3d vectors of magnetization
    cosinus=np.dot(m1,m2)
    # pars:
    # gp: conductance in parallel configuration
    # gap: conductance in antiparallel configuration
    # pp: polarization parallell
    # pa: polarization unparallel
    # grl: relaxation conductance left
    # grr: relaxation conductance right
    #extracting:
    gp,gap,pp,pa,grl,grr=pars
    #GB,GS=(1+cosinus)*gp*(1+pp)/4,(1+cosinus)*gp*(1-pp)/4
    #GSB,GBS=(1-cosinus)*gap*(1+pa)/4,(1-cosinus)*gap*(1-pa)/4
    #GRL,GRR=grl*(1-cosinus**2)/2,grl*(1-cosinus**2)/2
    GP,GAP=gp*(1.0+cosinus)/2.0,gap*(1-cosinus)/2.0
    GRL,GRR=grl*(1-cosinus**2)/2,grl*(1-cosinus**2)/2
    G=GP+GAP
    #left = A, right =B
    AA=-np.array([[G,GP*pp+GAP*pa],[GP*pp+GAP*pa,G+GRL]])
    AB=np.array([[G,GP*pp-GAP*pa],[GP*pp+GAP*pa,GP-GAP]])
    BA=np.array([[G,GP*pp+GAP*pa],[GP*pp-GAP*pa,GP-GAP]])
    BB=-np.array([[G,GP*pp-GAP*pa],[GP*pp-GAP*pa,G+GRR]])  
    return(AA,AB,BA,BB)
#conductance matrices of HH connector
def coHH(m1,m2,pars):
    #pars: just one parameter, conductance in parallel configuration
    cosinus=np.dot(m1,m2)
    cond=pars*(1+cosinus)/2
    return(-cond,cond,cond,-cond)
#conductance matrix of NN connector
def coNN(pars):
    #pars: just one parameter, conductance
    a=pars*np.diag(np.ones(4))
    return(-a,a,a,-a)
# need to add at least three more. Even afwachten... 
def coFH(m1,m2,pars):
    #pars: 
    # gp - conductance in parallel configuration
    # gap - conductance in anti-parallel configuration
    # gr - relaxation at ferro side
    cosinus= np.dot(m1,m2)
    gp,gap,gr=pars
    GP,GAP=gp*(1+cosinus)/2,gap*(1-cosinus)/2
    GR=gr*(1-cosinus**2)/2
    G=GP+GAP
    AA=-np.array([[G,GP-GAP],[GP-GAP,G+GR]])
    AB=np.array([[G],[GP-GAP]])
    BA=np.array([[G,GP-GAP]])
    BB=-np.array([G])
    return(AA,AB,BA,BB)
def coHF(m1,m2,pars): #exchamge of m1 and m2 does not matter
    aa,ab,ba,bb=coFH(m1,m2,pars)
    return(bb,ba,ab,aa)
#normal-ferro
def coFN(m,pars):
    G,P,reMG,imMG=pars
    #Gu,Gd=G*(1+P)/2,G*(1-P)/2
    crossM=np.array([[0,0,0,0],[0,0,-m[2],m[1]],[0,m[2],0,-m[0]],[0,-m[1],m[0],0]])
    quaM=np.array([[0,0,0,0],[0,m[0]**2,m[0]*m[1],m[0]*m[2]],[0,m[0]*m[1],m[1]**2,m[1]*m[2]],[0,m[0]*m[2],m[1]*m[2],m[2]**2]])
    eens=np.diag([1.0,0,0,0])
    eenv=np.diag([0,1.0,1.0,1.0])
    veM=np.array([[0,m[0],m[1],m[2]],[m[0],0,0,0],[m[1],0,0,0],[m[2],0,0,0]])
    AA=-G*np.array([[1.0,P],[P,1.0]])
    AB=G*np.array([[1.0,P*m[0],P*m[1],P*m[2]],[P,m[0],m[1],m[2]]])
    BA=G*np.array([[1.0,P],[P*m[0],m[0]],[P*m[1],m[1]],[P*m[2],m[2]]])
    BB=-G*(eens+quaM) - reMG*(eenv-quaM) -P*G*veM +imMG*crossM
    return(AA,AB,BA,BB)
def coNF(m,pars):
    aa,ab,ba,bb=coFN(m,pars)
    return(bb,ba,ab,aa)
#normal- halfmetal
def coHN(m,pars):
    #same as above, polarization =1
    G,reMG,imMG=pars
    crossM=np.array([[0,0,0,0],[0,0,-m[2],m[1]],[0,m[2],0,-m[0]],[0,-m[1],m[0],0]])
    quaM=np.array([[0,0,0,0],[0,m[0]**2,m[0]*m[1],m[0]*m[2]],[0,m[0]*m[1],m[1]**2,m[1]*m[2]],[0,m[0]*m[2],m[1]*m[2],m[2]**2]])
    eens=np.diag([1.0,0,0,0])
    eenv=np.diag([0,1.0,1.0,1.0])
    veM=np.array([[0,m[0],m[1],m[2]],[m[0],0,0,0],[m[1],0,0,0],[m[2],0,0,0]])
    AA=-G*np.array([[1]])
    AB=G*np.array([[1,m[0],m[1],m[2]]])
    BA=G*np.array([[1.0],[m[0]],[m[1]],[m[2]]])
    BB=-G*(eens+quaM) - reMG*(eenv-quaM) -G*veM +imMG*crossM   
    return(AA,AB,BA,BB)

def coNH(m,pars):
    aa,ab,ba,bb=coHN(m,pars)
    return(bb,ba,ab,aa)
#precession
def coP(pars):
    co,m=pars[0],pars[1:]
    crossM=np.array([[0,0,0,0],[0,0,-m[2],m[1]],[0,m[2],0,-m[0]],[0,-m[1],m[0],0]])
    AA=co*crossM 
    return(AA,None,None,None)
#relaxation
def coRN(pars):
    co=pars
    AA=-co*np.diag([0,1.0,1.0,1.0])
    return(AA,None,None,None)
def coRF(pars):
    co=pars
    AA=-co*np.array([[0,0],[0,1.0]])
    return(AA,None,None,None)
###########################
    
    


class SpintronicCircuit(LinGenCircuit):
    #no special init
    #####
    #adding functions
    ####
    def addFerroNode(self, name,magnetization=[0.0,0.0,1.0]):
        self._addNode(name,2,normagnetization(magnetization))
        self.node(name).type='F'
    def addHMNode(self, name,magnetization=[0.0,0.0,1.0]):
        self._addNode(name,1,normagnetization(magnetization))
        self.node(name).type='H'
    def addNormalNode(self,name):
        self._addNode(name,4,None)
        self.node(name).type='N'
    def addFerroTerminal(self, name,magnetization=[0.0,0.0,1.0]):
        self._addTerminal(name, 2, [1.0,0], normagnetization(magnetization)) #changed 6-3-2026
        self.terminal(name).type='F'
    def addHMTerminal(self, name,magnetization=[0.0,0.0,1.0]):
        self._addTerminal(name, 1, [1], normagnetization(magnetization))
        self.terminal(name).type='H'
    def addNormalTerminal(self, name):
        self._addTerminal(name, 4, [1.0,0.0,0.0,0.0], None)
        self.terminal(name).type='N'
    def addConductor(self, name,code,name1,name2,pars):
        #we first add the connector 
        self._addConnector(name,code,name1,name2,pars)
        #and them the check the number of parameters
        added=self.connector(name)
        added.type=added.AB.head.type+added.BA.head.type #define connector type
        #look into howManyPars dictionary
        assert (howManyPars[added.type]==le(added.pars)),"Wrong number of parameters for connector "+name+". For the type "+added.type+" "+str(howManyPars[added.type])+" is required, "+str(le(added.pars))+" is provided."      

    def addRelaxation(self,name,nodename,pars): #only to F and N nodes
        nodetype=self.node(nodename).type
        assert (not nodetype=='H'), "Cannot add relaxation to a half-metal node!"
        self._addConnector(name, 'N', nodename, None, pars)
        self.connector(name).type="R"+nodetype
        
    def addPrecession(self,name,nodename,conductance,pars): #only to normal nodes
        assert (le(pars) == 4), "Wrong number of parameters for precession"
        nodetype=self.node(nodename).type
        assert (nodetype=='N'), "Precession can be added to normal nodes only."
        m=normagnetization(pars[1:])
        self._addConnector(name,'N',nodename,None,[pars[0],m[0],m[1],m[2]])
    ####
    #changing functions
    ####
    #this sets the magtetization of a given node
    def setMagnetization(self,name,code,magnetization): #code: N or T
        m=normagnetization(magnetization)
        if (code=='N'): #node
            assert (not self.node(name).type == 'N'),'Cannot set magnetization of a normal node!'
            self._setNode(name, m)
        elif (code=='T'): #terminal
            assert (not self.terminal(name).type == 'N'),'Cannot set magnetization of a normal terminal!'
            self._setTerminal(name, m)
        else: 
            assert False, 'Wrong code: neither node nor terminal'
    def setConductance(self, name, pars):
        c=self.connector(name)
        #check the number of pars
        assert (howManyPars[c.type]==le(pars)),"Wrong number of parameters for connector "+name+". For the type "+c.type+" "+str(howManyPars[c.type])+" is required, "+str(le(pars))+" is provided."      
        self._setConnector(name, pars)
        
    ####
    #compute conductance
    ####
    def computeConductances(self,connector):
        #this is just the dispatcher for conductance functions
        typ=connector.type
        code=connector.code
        #a,b - pointers to nodes/terminals
        if(code=='N'): #no tail
            a,b=connector.AA.head, None
        else:
            a,b=connector.AB.head, connector.BA.head
        if(typ=='FF'):
            return(coFF(a.pars,b.pars,connector.pars))
        elif(typ=='HH'):
            return(coHH(a.pars,b.pars,connector.pars))
        elif(typ=='NN'):
            return(coNN(connector.pars))
        elif(typ=='FN'):
            return(coFN(a.pars,connector.pars))
        elif(typ=='NF'):
            return(coNF(b.pars,connector.pars))
        elif(typ=='HN'):
            return(coHN(a.pars,connector.pars))
        elif(typ=='NH'):
            return(coNH(b.pars,connector.pars))
        elif(typ=='FH'):
            return(coFH(a.pars, b.pars, connector.pars))
        elif(typ=='HF'):
            return(coHF(a.pars, b.pars, connector.pars))
        elif(typ=='P'):
            return(coP(a.pars, connector.pars))
        elif(typ=='RF'):
            return(coRF(connector.pars))
        elif(typ=='RN'):
            return(coRN(connector.pars))
        else: 
            print("Unknown connector type!")
        
        
        