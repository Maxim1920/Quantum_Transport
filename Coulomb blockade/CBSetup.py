# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:03:35 2023
@author: Yuli
Simulator for Coulomb blockade circuits
---------------
Version 0.0.1 Jan 31 2024
24-2-2024 correcting
implementing debug printing
25-2-2024 sorphisticating debug.
25-2-2024 implementing an alternative stopping function
So it becomes 
Version 0.0.2 Feb 25 2024
some bugs fixed. sorted_cache function added. movie_theater in the simulator
So it becomes
Version 0.0.3 March 16 2024 
Fixed voltage sign bug.
Version 0.0.4 Jan 16 2025
"""

import simulator as sim
import numpy as np
import time

def addtodic(dic,key,value):
    if(key in dic.keys()):
        dic[key].append(value)
    else:
        dic[key]=[value]

#class ChargingState(sim.State):
#    def transfertofrom(self,to,fro):
#        self.id[to]+=1
#        self.id[fro]-=1
#    def addto(self,to,value):
#        self.id+=value
        
class JunctionBTI(sim.Actor):
    def __init__(self,state,conductance,capacitance,firstisland, secondisland):
        self.conductance=conductance
        self.capacitance=capacitance
        self.first=firstisland
        self.second=secondisland
        self.cb=0 #must be set upon finalising
        sim.Actor.__init__(self,2,state)
    def rates(self):
        dv=self.state.voltage[self.first]-self.state.voltage[self.second]
        return(self.conductance*np.array([self.state.prerate(dv+self.cb),self.state.prerate(-dv+self.cb)]))
    def work(self,event):
        if (event==0):
            self.state.transfertofrom(self.first,self.second)
        else:
            self.state.transfertofrom(self.second,self.first)

class JunctionBIL(sim.Actor):
    def __init__(self,state,conductance,capacitance,island, lead):
        self.conductance=conductance
        self.capacitance=capacitance
        self.island=island
        self.lead=lead
        self.cb=0 #must be set upon finalising
        sim.Actor.__init__(self,2,state)
    def rates(self):
        dv=self.state.voltage[self.island]-self.state.biasvoltage[self.lead]
        return(self.conductance*np.array([self.state.prerate(dv+self.cb),self.state.prerate(-dv+self.cb)]))
    def work(self,event):
        self.state.addto(self.island,1-2*event) #0 adds an electron to the island

    
class CBSetup(sim.State):
    #######
    #initialization
    #######
    def __init__(self): #creates nessesary elements and lists
        self.temperature=0.1 #change if required
        self.finished=False
        self.islanddic={} #all dics: name->ordernumber
        self.leaddic={}
        self.gatedic={}
        self.junctionBTIdic={}
        self.junctionBILdic={}
        self.capacitanceBTIdic={}
        self.capacitanceBILdic={}
        #junctions
        self.JunctionsBTI=[]
        self.JunctionsBIL=[]
        #connectivity information
        self.islandstoleads={}
        self.islandstogates={}
        self.islandstoislands={}
        self.islandstoislandscap={}
        self.islandstoleadscap={}
        #reference to the charging state. Cannot set it at this point yet
        #self.S=ChargingState() setup is now the state
    ###########    
    #adding element functions
    ##########
    def add_island(self,name):
        assert(not self.finished), "Error. The setup is already finished"
        assert(not (name in self.islanddic.keys())), "Error. The island "+name+" is already added."
        l=len(self.islanddic)
        self.islanddic[name]=l
        pass
    def add_lead(self,name):
        assert(not self.finished), "Error. The setup is already finished"
        assert(not (name in self.leaddic.keys())), "Error. The lead "+name+" is already added."
        l=len(self.leaddic)
        self.leaddic[name]=l
        pass
    def add_gate(self,name,island,capacitance):
        assert(not self.finished), "Error. The setup is already finished"
        assert(not (name in self.gatedic.keys())), "Error. The gate "+name+" is already added."
        assert((island in self.islanddic.keys())), "Error. The island "+island+" is not here."
        l=len(self.gatedic)
        isl=self.islanddic[island]
        self.gatedic[name]=(l,isl, capacitance)
        addtodic(self.islandstogates,isl,[l,capacitance])
        pass
    def add_junctionBTI(self,name,firstisland,secondisland,conductance,capacitance):
        assert(not self.finished), "Error. The setup is already finished"
        assert(not (name in self.junctionBTIdic.keys())), "Error. The junction "+name+" is already added."
        assert( (firstisland in self.islanddic.keys()) and (secondisland in self.islanddic.keys()) ), "Error. These islands are not present!"
        assert( not (firstisland == secondisland)), "Error. It's a perversion to connect an island to itself."
        l=len(self.junctionBTIdic)
        self.junctionBTIdic[name]=l
        #connectivity
        addtodic(self.islandstoislands,self.islanddic[firstisland],[self.islanddic[secondisland],capacitance])
        addtodic(self.islandstoislands,self.islanddic[secondisland],[self.islanddic[firstisland],capacitance])
        #making an actor
        self.JunctionsBTI.append(JunctionBTI(self,conductance,capacitance,self.islanddic[firstisland],self.islanddic[secondisland]))
        pass
    def add_junctionBIL(self,name,island,lead,conductance,capacitance):
        assert(not self.finished), "Error. The setup is already finished"
        assert(not (name in self.junctionBILdic)), "Error. The junction "+name+" is already added."
        assert( (island in self.islanddic.keys()) and (lead in self.leaddic.keys()) ), "Error. Either island or gate is not present!"
        l=len(self.junctionBILdic)
        self.junctionBILdic[name]=l
        #connectivity
        addtodic(self.islandstoleads,self.islanddic[island],[self.leaddic[lead],capacitance])
        #making an actor
        self.JunctionsBIL.append(JunctionBIL(self,conductance,capacitance,self.islanddic[island],self.leaddic[lead]))
        pass    
    def add_capacitanceBTI(self,name,firstisland,secondisland,capacitance):
        assert(not self.finished), "Error. The setup is already finished"
        assert(not (name in self.capacitanceBTIdic.keys())), "Error. The capacitance "+name+" is already added."
        assert( (firstisland in self.islanddic.keys()) and (secondisland in self.islanddic.keys()) ), "Error. These islands are not present!"
        assert( not (firstisland == secondisland)), "Error. It's a perversion to connect an island to itself."
        l=len(self.capacitanceBTIdic)
        self.capacitanceBTIdic[name]=l
        #connectivity
        addtodic(self.islandstoislandscap,self.islanddic[firstisland],[self.islanddic[secondisland],capacitance])
        addtodic(self.islandstoislandscap,self.islanddic[secondisland],[self.islanddic[firstisland],capacitance])
    def add_capacitanceBIL(self,name,island,lead,capacitance):
        assert(not self.finished), "Error. The setup is already finished"
        assert(not (name in self.capacitanceBILdic.keys())), "Error. The capacitance "+name+" is already added."
        assert( (island in self.islanddic.keys()) and (lead in self.gatedic.keys()) ), "Error. Either island or gate is not present!"
        l=len(self.capacitanceBILdic)
        self.capacitanceBIL[name]=l
        #connectivity
        addtodic(self.islandstoleadscap,self.islanddic[island],[self.leaddic[lead],capacitance])
    #####################
    #finalising the setup
    ####################
    def finish(self):
        #set up dimensions
        self.dim=len(self.islanddic)
        self.numgates=len(self.gatedic)
        self.numleads=len(self.leaddic)
        ###################
        #here we keep the charging state
        ###################
        self.id=np.zeros(self.dim,dtype=np.int8)
        #build capacitance matrix:loook up connectivity lists
        self.capmatrix=np.zeros((self.dim,self.dim))
        #first island to island junctions and island to island capacitances
        for di in (self.islandstoislands, self.islandstoislandscap):        
            for i in range(self.dim):
                if i in di.keys():
                    li=di[i]
                    for c in li:
                        self.capmatrix[i,i]+=c[1] 
                        self.capmatrix[i,c[0]]-=c[1]
        #second capacitancies to the leads and to the gates
        for di in (self.islandstogates, self.islandstoleads):
             for i in range(self.dim):
                 if i in di.keys():
                     li=di[i]
                     for c in li:
                         self.capmatrix[i,i]+=c[1]             
        self.invcapmatrix=np.linalg.inv(self.capmatrix) #now we have inverse capmatrix
        #set cbs in the junctions
        for ju in self.JunctionsBIL:
            ju.cb=0.5*self.invcapmatrix[ju.island,ju.island]
        for ju in self.JunctionsBTI:
            ju.cb=0.5*(self.invcapmatrix[ju.first,ju.first]+self.invcapmatrix[ju.second,ju.second])-self.invcapmatrix[ju.first,ju.second]
        #reserve space for voltages and quties        
        self.biasvoltage=np.zeros(self.numleads) #lead voltages
        self.gatevoltage=np.zeros(self.numgates) #gate voltages
        self.q=np.zeros(self.dim) #this will be for induced charges in each island
        self.voltage=np.zeros(self.dim) #for voltages
        self.qvoltage=np.zeros(self.dim)
        #create the simulator
        self.simulator=sim.Simulator(self.JunctionsBIL+self.JunctionsBTI, self)
        #it would be good to check the connectivities here... Would do so later
        self.finished=True #done
        #Make desciption of the setup:
        self.info="Number of junctions: "+str(len(self.JunctionsBIL)+len(self.JunctionsBTI))+"\n"
        self.info+="Number of leads: "+str(self.numleads)+"\n"
        self.info+="Number of gates: "+str(self.numgates)+"\n"
        self.info+="Number of islands: "+str(self.dim)+"\n"
        self.changed=False #this shows if the qs need to be recomputed
    #################    
    #print setup info
    ################
    def print_info(self):
        print(self.info)
    ##################    
    # for setting parameters
    ##################
    def set_lead_voltage(self,name,value):
        self.biasvoltage[self.leaddic[name]]=value
        self.changed=True
    def set_gate_voltage(self,name,value):
        self.gatevoltage[self.gatedic[name][0]]=value
        self.changed=True
    def set_lead_voltages(self,varray): #alternative way. requires ids of the leads
        self.biasvoltage[:]=varray
        self.changed=True
    def set_gate_voltages(self,varray): #alternative way. requires ids of the gates
        self.gatevoltage[:]=varray
        self.changed=True
    ##################
    #internal functions required for computing the rates
    ##################
    def prerate (self,E):
        # Gamma = G*E/(exp(E/T)-1)  if E<0  if E>0 Gamma=G*E*exp(-E/T) if E=0 G*T
        # so we need temperature and energy differences...
        # we get them from the state: it is simpler than to reference to the setup...
        if (abs(E)<1e-6):
            return(self.temperature)
        if (E<0):
            return(E/(np.exp(E/self.temperature)-1))
        if (E>0):
            return(E*np.exp(-E/self.temperature)/(1-np.exp(-E/self.temperature)))
        #-------------------------------
    def _compute_q_and_vq(self):
        self.q[:]=0 #zeroit
        #first gates
        for isl in self.islandstogates.keys():
            for item in self.islandstogates[isl]:
                self.q[isl]+=item[1]*self.gatevoltage[item[0]]
        #second junctions to leads,capacitanes to leads
        for di in (self.islandstoleads,self.islandstoleadscap):
            for isl in di.keys():
                for item in di[isl]:
                    self.q[isl]+=item[1]*self.biasvoltage[item[0]]
        #qs are ready, we compute the corresponding voltages
        self.qvoltage=self.invcapmatrix @ self.q #14-01-2025 wrong sign? Changing sign
        
    def _compute_voltages(self):
        if(self.changed):
            self._compute_q_and_vq()
        self.voltage=(self.invcapmatrix @ (self.id))+self.qvoltage
    def guess_state(self):
        self.id=np.floor(self.q +0.5)
    def prepare_to_compute(self): #this is called by the simulator    
        self._compute_voltages()
    ##################
    #state changin functions
    ##################
    def transfertofrom(self,to,fro):
        self.id[to]+=1
        self.id[fro]-=1
    def addto(self,to,value):
        self.id[to]+=value
    #chache investigation functions
    def sorted_cache(self,num=5): #returns num most probable configurations and their probabilities
        cache_dic=self.simulator.statetimes
        sorted_chache = sorted(cache_dic.items(), key=lambda x:x[1], reverse=True)
        totaltime=sum(x[1] for x in sorted_chache)
        result=[(x[0],x[1]/totaltime) for x in sorted_chache]
        if (num >0):
            return(result[:num])
        else:
            return(result)
        
    #simulation 
    def give_a_run(self,stop='time', timeinterval=1000.0,maxsteps=1000, charge = ("0",100),debug=False): 
        if(debug):
            starttime=time.time()
        if (stop=='charge'): #figure out what junction. Its name charge[0]
            num=self.junctionBILdic[charge[0]]
            difference=(2*num,2*num+1,charge[1])
            stop='difference'
        else: 
            difference =(0,0,100) #not importan anyway
        #returns list of currents from each lead
        expired,steps=self.simulator.simcycle(stop=stop,timeinterval=timeinterval,maxsteps=maxsteps,difference=difference)
        #extractcurrent
        currents=[]
        for i in range(len(self.JunctionsBIL)):
            currents.append((self.simulator.eventsoccured[2*i]-self.simulator.eventsoccured[2*i+1])/expired)
        if(debug):
            print("CPU time: "+"{:.1f}".format(time.time()-starttime)+" Simulated time: "+"{:.1f}".format(expired)+" Events: "+str(steps)+" Cache: "+str(len(self.simulator.cache)))
        return(currents)
    def runit(self,numsteps,initial_leads,initial_gates,final_leads,final_gates, stop='time', timeinterval=1000.0,maxsteps=1000, charge = ("0",100), debug=False):
    #returns list lists of currents to each lead
        listofcurrents=[]
        for i in range(numsteps+1):
            if (debug):
                print("Run "+str(i+1))
            #setting voltages
            self.set_lead_voltages(initial_leads*(1-i/numsteps)+final_leads*(i/numsteps))
            self.set_gate_voltages(initial_gates*(1-i/numsteps)+final_gates*(i/numsteps))
            curr=self.give_a_run(stop=stop, timeinterval=timeinterval,maxsteps=maxsteps, charge = charge,debug=debug)
            listofcurrents.append(curr)
        return(listofcurrents)
    
    
