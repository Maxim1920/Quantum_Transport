# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 19:52:00 2026

Example file for the package supercurcuits
 
@author: dhyul
"""
#
# the calculation calls are commented out in this file
#

from supercircuits import SuperCircuit
import numpy as np
import time

############
#Minimum example: two terminals, one ballistic connector, no nodes
############
Minim=SuperCircuit()
#add terminals
Minim.addNormalTerminal("N")
Minim.addSuperTerminal("S")
#add connector 
Gb=2.0
Minim.addConnector("Ba", 'TT', 'Ballistic', 'N', 'S', conductance=Gb)
#let us compute the differential conductance versus voltage collecting the 
#results in the list
voltages=np.linspace(0,3.0,1000)
results=[]
for v in voltages:
    Minim.setNormalTerminal("N",voltage=v) #sets to voltage v
    #uncomment to run. Execution time 1.6 sec
    #results.append(Minim.fastDifferentialConductances("N"))
#print(results[3]) 
#expected output:
#[np.float64(-3.99960002376781), np.float64(3.99960002376781)] 
#the conductance approximately doubles at low voltage,

###end of the first result    

#########
#Chain of tunnel juncions between two superconducting terminals
########
Nd=3 #number of nodes.
Gt=1.0 #normal chain conductance
Eth=10.0 #for leakage connectors
Chain=SuperCircuit()
#add terminals
Chain.addSuperTerminal("Left", Gamma=0.01) #higher Gamma simplifies integration
Chain.addSuperTerminal("Right", Gamma=0.01) 
#add nodes
for i in range(Nd):
    Chain.addNode('N'+str(i))
#add connectors between the nodes
for i in range(Nd-1):
    Chain.addConnector("C"+str(i+1), "NN", 'Tunnel', 'N'+str(i), 'N'+str(i+1), conductance=Gt*(Nd+1))
#add connectors from end nodes to terminals
Chain.addConnector("C0", 'NT', 'Tunnel', 'N0', 'Left')
Chain.addConnector("C"+str(Nd), 'NT', 'Tunnel', 'N'+str(Nd-1), 'Right')
#add leakage connectors
for i in range(Nd):
    Chain.addLeakageConnector('L'+str(i),'N'+str(i),Eth*Nd)
#let us apply phase difference between the terminals
p=np.pi/2
Chain.setSuperTerminal('Left', phase=-p/2)
Chain.setSuperTerminal('Right', phase=p/2)
#Let us compute supercurrents to each terminal
#uncomment to execute. Ex time 20 sec for N=3
#print(Chain.FullCurrents())
#expected: [-0.68506052  0.68506052]

#increase temperature, compute supercurrent once again
#Chain.setGlobalTemperature(0.75)
#uncomment to execute. Ex time 22 sec for N=3
#print(Chain.FullCurrents())
#expected: [-0.30368223  0.30368223] significant decrease

#end of the second example

#################
#Threepole made of diffusive connectors: two superterminals, one normal.
#################
#conductances
G1,G2,G3=0.7,0.9,0.5
Threepole=SuperCircuit()
#add terminals
Threepole.addNormalTerminal('N',voltage=0.8) #set voltage right away.
Threepole.addSuperTerminal('S1')
Threepole.addSuperTerminal('S2')
#add node
Threepole.addNode('Node')
#add connectors
Threepole.addConnector('C3', 'NT','Diffusive','Node', 'N',conductance=G3)
Threepole.addConnector('C1', 'NT','Diffusive','Node', 'S1',conductance=G1)
Threepole.addConnector('C2', 'NT','Diffusive','Node', 'S2',conductance=G2)

#let us compute the currents for a set of phases
st=time.time()
Threepole.changeAdjust(0.2) #to remedy a possible overflow
for p in np.linspace(0,np.pi,10):
    Threepole.setSuperTerminal('S1', phase=-p/2)
    Threepole.setSuperTerminal('S2', phase=p/2)
    #uncomment to execute. Expected time: 56 sec
    #print("phase: ",p, "currents: ",Threepole.FullCurrents())
#expected output 
#phase:  0.0 currents:  [-0.33678964  0.14742982  0.18955263]
#phase:  0.3490658503988659 currents:  [-0.33517831  0.00822076  0.33129688]
#phase:  0.6981317007977318 currents:  [-0.34036577 -0.13763172  0.48056831]
#phase:  1.0471975511965976 currents:  [-0.3480389  -0.22678513  0.57593267]
#phase:  1.3962634015954636 currents:  [-0.36498948 -0.19489281  0.56007452]
#phase:  1.7453292519943295 currents:  [-0.37476808 -0.06357044  0.43910101]
#phase:  2.0943951023931953 currents:  [-0.36828313  0.01867815  0.34975382]
#phase:  2.443460952792061 currents:  [-0.35039601  0.06666046  0.28284445]
#phase:  2.792526803190927 currents:  [-0.32638595  0.09925233  0.2265193 ]
#phase:  3.141592653589793 currents:  [-0.31409758  0.13523843  0.17840148]

#Let us inspect the node at energy =0.8
Threepole.setEnergy(0.8)
Threepole.solve() #solves the circuit. takes 0.1 sec
#uncomment to execute
#print(Threepole.inspectNode('Node'))
#print(Threepole.currentsdensities(0.8))
#expected output (0.9942522841861563, (1.5707963267948954+8.881784197001256e-16j), 
#0.2035136929773898, 0.7042076111645286)
#[-0.20460402  0.08595361  0.11865041]

###end of the 3rd example#####


