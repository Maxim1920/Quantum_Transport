# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 14:23:19 2026

example file for the countingcircuit package

@author: dhyul
"""

from countingcircuit import CountingCircuit
import numpy as np
import pickle #to pickle long results to the disk

#we want to build H type circuit with four terminals (left/right, up/down)
#two nodes (left and right) and five connectors
#initializing the setup
H=CountingCircuit() #default energy and global temperature
#
#Terminals at default voltage, global temperature and counting fields
H.addTerminal('LeftUp')
H.addTerminal('RightUp')
H.addTerminal('LeftDown')
H.addTerminal('RightDown')
#nodes
H.addNode('LeftNode')
H.addNode('RightNode')
#now let us add connectors.
#Their conductances 
Glu, Gru, Gld, Grd = 1.0,1.0,1.0,1.0
Gb=0.001
#Diffusive connectors in the left branch
H.addConnector('LUconn','NT','Diffusive','LeftNode','LeftUp',conductance=Glu)
H.addConnector('LDconn','NT','Diffusive','LeftNode','LeftDown',conductance=Gld)
#Ballistic and diffusive connectors in the right branch
H.addConnector('RUconn','NT','Ballistic','RightNode','RightUp',conductance=Gru)
H.addConnector('RDconn','NT','Diffusive','RightNode','RightDown',conductance=Grd)
#Tunnel connector bridges the branches
H.addConnector('Bridge','NN','Tunnel','LeftNode','RightNode',conductance=Gb)
#
# The circuit is ready.
#Let's change something.
#more reasonable conductance Gb=2.0 for the bridge
H.setConductance('Bridge', 2.0)
#Elevated global temperature
H.setGlobalTemperature(2.0)
#Let us apply some voltages to the terminals
H.setTerminal('LeftUp',voltage=1.0) #could do this in very beginning
H.setTerminal('RightUp',voltage=-1.0) # opposite voltages
#Down terminals remain at default voltage, that is, are grounded
#Now computations. First the action, let us set some counting fields, 
#otherwise it is  zero
H.setTerminal('LeftUp',cf=0.3) #could do this in very beginning
H.setTerminal('RightUp',cf=0.5)

#let us compute the action. Uncomment to run it at your computer

#action per interval (the energy is set to default zero at the moment)
#a=H.actionPI()
#print(a)
#answer (-0.021614176016078103+0.13571170122975804j)
#Execution time =0.04

#Now integrated action assuming zero temperature: method Z
#a=H.actionZ()
#print(a)
#answer (-0.035403804503300905-0.16659979740155645j)
#Execution time 0.09

#Now integrate with temperature. method INT using default integration parameters
#    a=H.actionINT()
#    print(a)
#answer (-0.3396165803099158-0.15781177174829766j)
#Execution time 5.3  roughly 50 times longer than in the previous example

#Ok. Let us compute the average currents in all terminals. let us first define the list 
#of all terminals
allterminals=list(H.Terminals.keys())
#first at zero temperature method Z
#a = H.currents(allterminals,method='Z')
#print(a)
#answer [-0.83333333  0.83333333  0.16666667 -0.16666667]
#Execution time 0.18
#the currents are TO terminals, the first terminal is at positive voltage, 
#so the current TO it is negative

#now taking temperature into account. Let's improve integration parameters
H.integrationparameters=[2.0,0.0,400]
#a = H.currents(allterminals,method='INT')
#print(a)
#answer [-0.83333321  0.83333322  0.16666665 -0.16666664]
#Execution time 15.6
# much longer than in the previous example, the currents are almost the same 
#Indeed, they should not depend on temperature. The difference 
#indicates numerical error in the integration
#
#Let's address the noises and concentrate on down terminals only
#assuming zero temperature
#let's increase accuracy - good for comuting the derivaties
H.solvepars.goal=1e-12
#
#a = H.noises(['LeftDown','RightDown'],method='Z')
#print(a)
#answer [[ 0.30324074 -0.08487654]   4 values, 3 independent
# [-0.08487654  0.29012346]]
#Execution time 3.1
#taking temperature into account
#a = H.noises(['LeftDown','RightDown'],method='INT')
#print(a)
#answer =[[ 7.11010139 -1.83748588]
# [-1.83748588  6.87181596]]   Bigger noise due to temperature
# Execution time 399.3 !
#
#let us compute the noises per interval assuming f=1 in LeftUp and f=0 
#in all other terminals. We set the mask
H.putMask(['LeftUp'])
# and use the PI method
#a = H.noises(['LeftDown','RightDown'],method='PI')
#print(a)
#answer = [[ 0.18595679 -0.03858025]
# [-0.03858025  0.12037037]]
# Execution time 1.1
#
# third order cumulants - 8 values, 4 independent 
#a = H.thirdcumulants(['LeftDown','RightDown'],method='PI')
#print(a)
#answer = [[[ 0.05383802 -0.00492684]
#  [-0.00492684 -0.01754115]]
# [[-0.00492684 -0.01754115]
#  [-0.01754115  0.0613283 ]]]
# Execution time 22.670249399961904
#
#forth order cumulants 16 values, 5 independent
#a = H.fourthcumulants(['LeftDown','RightDown'],method='PI')
#print(a)
#
#answer = [[[[-0.01821533  0.00889538]
#   [ 0.00889538 -0.00762128]]
#  [[ 0.00889538 -0.00762128]
#   [-0.00762128  0.00241472]]]
# [[[ 0.00889538 -0.00762128]
#   [-0.00762128  0.00241472]]
#  [[-0.00762128  0.00241472]
#   [ 0.00241472  0.00294689]]]]
# Execution time 488.2
#do not need the mask any more
H.removeMask()
#
#Let us illustrate transmission distribution
# first we look at transmission matrix from LeftUp to RightDown
#values of transmissions to look at
ts=np.linspace(0.1,0.9,100)
#res = H.partialTransmissionDistribution(['LeftUp'],['RightDown'], ts) 
#pickle.dump([ts,res],open('transmission1.p','wb')) #pickles the result to the disk
#Execution time 47.8
#
# full transmission distibution: from LeftUp to all other terminals
#res = H.fullTransmissionDistribution(['LeftUp'], ts) 
#pickle.dump([ts,res],open('transmission2.p','wb')) #pickles the result to the disk
# Execution time 22.4
#
# gigantic fluctuations along a ray, method Z
# values of chi to check
chis=np.linspace(-3,3,200)
#ray: changing chi in LeftUp only
ray=np.array([1.0,0.0,0.0,0.0]) #4 terminals = 4 dim vector
#actions,currents = H.gfRay(chis,ray,method='Z') 
#pickle.dump([chis,actions,currents],open('gfRay.p','wb')) #pickles the result to the disk
#Execution time 40.67589179996867
#
# Gigantic fluctuation at given current settings
# dictionary of terminals and current settings
dic={'LeftUp': -0.6,'RightUp':0.2,'LeftDown':0.3}
#method 'Z'
#let us reduce the iteration agressivity: could be divergencies
H.changeAdjust(0.2)
#
#action,cfs,currents = H.gfIs(dic,method='Z') #default bound
#print('action: ',float(action),'cfs: ',cfs, "currents: ",currents)
# Output 
#action:  -0.9203351915138651 #log of probability
#cfs:  [-0.12824214  3.51962387  0.09312151] within bounds: reliable value
#currents:  [-0.60000136  0.20000053  0.30000008  0.10000074]  close to requested values
# Execution time: 27.4
#
#Finallly, let us try the node inspection
#Set energy
H.setEnergy(0.3)
#Set  counting fields conform to the previously obtained set
#of purely imaginary values
H.setTerminal('LeftUp',cf=1j*(-0.12824214))
H.setTerminal('LeftDown',cf=1j*(3.51962387))
H.setTerminal('RightUp', cf=1j*0.09312151)
H.setTerminal('RightDown',cf=0.0)
#solve
#H.solve()
#Execution time: 0.23
#inspect both nodes:
#for nName in ['LeftNode','RightNode']
#    f,chi=H.inspectNode(nName)
#    print(nName+': '+ 'filling factor:'+str(np.real(f))+' imaginary cfield:+'str(np.imag(chi)))
# we print real and imaginary part
#output:
#LeftNode: filling factor:0.4099530722134471 imaginary cfield:1.4746984513884913
#RightNode: filling factor:0.21402315697861624 imaginary cfield:0.89923892775182
#Execution time: 0.002
#

