# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 12:07:52 2023

@author: Yuli
This is a generic simulator class.
It works with a list of Actors, each can shot Events changing the State
At each simulation step, the simulator asks the Actors about the rates of their
Events and decides which Event happens. If cashing is on, it uses the cashed rates
It commands the corresponding Actor to change the state. 
During the simulation cycle, it keeps track of times spent in each state
and number of Events occured.
27-12-2023 change implementation of cashing
Version 0.0.1 Jan 31 2024
25-2-2024 Changed output of run: exprired, steps
Changing stopping function: time, steps, events, difference
So it becomes
Version 0.0.2 February 25
added movie_theater version 0.0.3 March 16
"""

"""Abstract Actor class"""

import numpy as np

class State:
    def __init__(self):
        self.id=None
    def setstate(self,arg):
        self.id=arg
    def prepare_to_compute(self):
        pass
    #access to the state by the variable state.id   
        

class Actor:
    def __init__(self,numevents,state):
        self.numevents=numevents
        self.state=state #access to info about the state
    def rates(self):
        return(np.ones[self.numevents])
    def work(self, eventid):
        pass

class Simulator:
    def __init__(self,actorlist,stateref):
        self.actors=actorlist
        self.stateref=stateref
        self.dim=sum([x.numevents for x in self.actors])
        self.allrates=np.zeros(self.dim)
        #now we need a correspondence between the index of this array and events
        self.eventdic={}
        self.indexlist=[]
        i=0
        for a in range(len(self.actors)):
            for e in range(self.actors[a].numevents):
                self.eventdic[(a,e)]=i
                self.indexlist.append((a,e))
                i+=1
        #this is to accumulate the results
        self.eventsoccured=np.zeros(self.dim,dtype=np.int64)   
        self.statetimes={}
        self.cache={}
        self.debug=False
        self.caching=True
        
    def _compute_rates(self):
        #here is the trick: some common operations with the state before computing the rates
        self.stateref.prepare_to_compute()
        i=0
        g=0
        arr=np.zeros(self.dim)
        for a in range(len(self.actors)):
            rates=self.actors[a].rates()
            for e in range(self.actors[a].numevents):
                g+=rates[e]
                arr[i]=g
                i+=1
        return(arr)
        
    def compute_rates(self):
        r=tuple(self.stateref.id)
        if(self.caching and (r in self.cache.keys())):
            return(self.cache[r])
        else:
            arr=self._compute_rates()
            if(self.caching):
                self.cache[r]=arr
        return(arr)
    
    def simstep(self):
        arr=self.compute_rates()
        i=np.searchsorted(arr,(arr[-1]*np.random.rand()))
        self.eventsoccured[i]+=1
        (a,e)=self.indexlist[i]
        timeincrement=-np.log(np.random.rand())/arr[-1]
        #collecting info about the state
        n=tuple(self.stateref.id)
        if (n in self.statetimes.keys()):
            self.statetimes[n]+=timeincrement
        else:
            self.statetimes[n]=timeincrement
        #let the actor change the state
        self.actors[a].work(e)
        if(self.debug):
            print("Actor:"+str(a)+" Event:"+str(e)+" Time:"+"{:.2f}".format(timeincrement)+ " New state:"+str(self.stateref.id))
        return(timeincrement)
    def movie_theater(self,fi,runs=25):
        movie=[]
        for i in range(runs):
            st=tuple(self.stateref.id)
            time=self.simstep()
            movie.append(str(st)+":duration:"+"{:.2f}".format(time)+"\n")
        f=open(fi,"w")
        f.writelines(movie)
        f.close()
        
    def simcycle(self, stop='time', timeinterval=1000,maxsteps=1000,events=(0,100), difference = (0,1,100), cleanup=True):
        #clean up if nessesary
        if(cleanup):
            self.eventsoccured[:]=0
            self.statetimes={}
            self.cache={}
        expired=0.0
        steps=0
        while{True}:
            expired+=self.simstep() 
            steps+=1
            #print(steps)
            #stopping it
            if( stop=='time' and expired>timeinterval):
                break
            if( stop=='steps' and steps>maxsteps):
                break
            if( stop=='events' and self.eventsoccured[events[0]]>events[1]): #tuple (actor id, maxevents)
                break
            #for difference: tuple (first id, second id, maxevents)
            if( stop=='difference' and np.abs(self.eventsoccured[difference[0]]-self.eventsoccured[difference[1]])>difference[2]):
                break
        return(expired,steps)
            
            
        
            
    