# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 19:19:47 2023

@author: Yuli
"""
import numpy as np
import CBSetup
from myutils import printall #for writting data to a file

#this is a probe file for Coulomb blockade simulator
#creates N-junction array with junction capacitances C, ground capacitances Cg
N=8
C=1
Cg=0.5
G=1

Setup= CBSetup.CBSetup()

#two leads
Setup.add_lead('Left')
Setup.add_lead('Right')

#creates N-2 islands and corresponding gates
for i in range(N-1):
    Setup.add_island('Island'+str(i).zfill(3))
    Setup.add_gate('Gate'+str(i).zfill(3),'Island'+str(i).zfill(3),Cg )
#junctions to the leads
Setup.add_junctionBIL("LJunction",'Island000','Left',G,C) 
Setup.add_junctionBIL('RJunction','Island'+str(N-2).zfill(3),'Right',G,C) 

#junctions in between 
for i in range(N-2):
    Setup.add_junctionBTI('Junction'+str(i).zfill(3),'Island'+str(i).zfill(3),'Island'+str(i+1).zfill(3),G,C)
#finishing
Setup.finish()
#Setup.print_info()
#print(Setup.islandstogates)
#setting voltages:
V=0.5 #under cb threshold
Setup.temperature=0.01
Setup.set_lead_voltages(np.array([V/2,-V/2])) #make a single run
currents=Setup.give_a_run(stop='steps',maxsteps=2000)
print(currents)
#let us take an I-Vcurve from 0 to 20 with 1000 steps, time interval 1000. Took about 15 min for N=8 T=0.27
zeros=np.zeros(N-1)
Vmax=20.0
Numsteps=100
voltages=np.linspace(0,Vmax, num=Numsteps+1)
# a long run of simulation cycles
listofcurrents=Setup.runit(Numsteps,np.array([0,0]),zeros,np.array([Vmax/2,-Vmax/2]),zeros,timeinterval=1000,debug=True)
#saving to disk
printall("CBSetupprobe3.dat",[voltages,listofcurrents],Numsteps+1)

