# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 15:18:42 2026

checks spintronics

@author: dhyul
"""
#import numpy as np
from spintronics import SpintronicCircuit

Setup=SpintronicCircuit()
#parameters:
Gn=1.0 #for normal connectors
#for a connector between normal and ferro
GFN=[1.0,0.9,0.5,0.5] 
#members of the list: G - conductance P- polarization
# ReGm - real part of mixing conductance ImGm - imaginary part of mixing conductance

#adding normal terminals
Setup.addNormalTerminal("Left")
Setup.addNormalTerminal("Right")
#adding ferro terminals
Setup.addFerroTerminal("F1",magnetization=[1.0,0,0])
Setup.addFerroTerminal("F2",magnetization=[0,1.0,0])
Setup.addFerroTerminal("F3",magnetization=[0,0,1.0])
#adding normal nodes
for i in range(1,4):
    Setup.addNormalNode("N"+str(i))
#adding normal-normal connectors
Setup.addConductor("C0", 'NT', 'N1', 'Left', Gn)
Setup.addConductor("C3", 'NT', 'N3', 'Right', Gn)
Setup.addConductor("C1", 'NN', 'N1', 'N2', Gn)
Setup.addConductor("C2", 'NN', 'N2', 'N3', Gn)
#adding ferro-normal connectors
for i in range(1,4):
    Setup.addConductor('CF'+str(i), 'NT', 'N'+str(i), 'F'+str(i), GFN)
Setup.finalize()

#setting magnetization of the nodes
Setup.setMagnetization('T', 'F1', [-1.0,0,0])
Setup.setMagnetization('T', 'F2', [0,-1.0,0])
Setup.setMagnetization('T', 'F3', [0,0,-1.0])

Setup.compute()
#printing conductance matrix
print("Conductance matrix: ", Setup.ConductanceMatrix)
Setup.setVoltage(['Left','F1'], [0,1.0])
Setup.computeCurrents()
#printing current in F2
print("Current in F2", Setup.current('F2'))
#printing currents for two spin directions in F2
print("Currents in F2",Setup.gencurrent('F2'))
#printing node voltages
Setup.computeNodeVoltages()
for i in range(1,4):
    print("Matrix voltage in node N"+str(i)+": ", Setup.nodevoltage('N'+str(i)))



    


