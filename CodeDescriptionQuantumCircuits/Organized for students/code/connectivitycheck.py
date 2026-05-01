# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 13:05:30 2026

Connectivity report

@author: dhyul
"""
#We start with 3 dictionaries made from ids of elements
# T: terminal : list of connectors
# N: nodes: list of connectors
# C: connectors: [list of terminals, list of nodes]
def connectivityReport(T,N,C):
    #print('Terminals',T)
    #print('Nodes',N)
    #print('Connectors',C)
    #to start with, let us make dictionaries to put colors in
    Tcolor={ t: None for t in T.keys()}
    Ncolor={n: None for n in N.keys()}
    #intitial color
    currentcolor=0
    #report
    report=''
    def colorT(t):
        nextnodes,nextterminals=[],[]
        Tcolor[t]=currentcolor
        for c in T[t]:
            #check if really connected 
            assert (t in C[c][0]), "Error in dictionaries: terminal "+str(t)+" is not connected to the connector "+str(c)
            nextnodes+=list(filter(lambda x: Ncolor[x] == None, C[c][1]))
            nextterminals+=list(filter(lambda x: Tcolor[x] == None, C[c][0]))
        nextnodes=list(set(nextnodes))
        nextterminals=list(set(nextterminals))
        for ne in nextnodes:
            colorN(ne)
        for te in nextterminals:
            colorT(te)
    def colorN(n):
        nextnodes,nextterminals=[],[]
        Ncolor[n]=currentcolor
        for c in N[n]:
            assert (n in C[c][1]), "Error in dictionaries: node "+str(n)+" is not connected to the connector "+str(c)
            nextnodes+=list(filter(lambda x: Ncolor[x] == None, C[c][1]))
            nextterminals+=list(filter(lambda x: Tcolor[x] == None, C[c][0]))
        nextnodes=list(set(nextnodes))
        nextterminals=list(set(nextterminals))
        for ne in nextnodes:
            colorN(ne)
        for te in nextterminals:
            colorT(te)
    for t in Tcolor.keys():
        if(Tcolor[t] == None): #colors all terminals and nodes connected to this terminal to current 
            colorT(t)
            currentcolor+=1
    report+="All terminals colored. Number of groups: "+str(currentcolor)+"\n"
    termgroups=currentcolor
    for n in Ncolor.keys():
        if (Ncolor[n]==None):
            colorN(n)
            currentcolor+=1
    report+="All nodes colored. Extra groups revealed: "+str(currentcolor-termgroups)+"\n"
    if (currentcolor ==1):
        report+="Lucky you are. All your elements are neatly connected to each other\n"
    report+="Reporting the groups:\n"
    #sorting out the groups
    #dic={}
    for g in range(currentcolor):
        ter=list(filter(lambda x: Tcolor[x]==g,Tcolor.keys()))
        nod=list(filter(lambda x: Ncolor[x]==g,Ncolor.keys()))
        report+="Group #"+str(g+1)+": "+str(len(ter)+len(nod))+" elements.\n"
        if (len(ter)>0):
            report+=str(len(ter))+" terminals: "+str(ter)+"\n"
        if (len(nod)>0):   
            report+=str(len(nod))+" nodes: "+str(nod)+"\n"
        if (len(ter)+len(nod) == 0):
            report+="WARNING! Zero element group. Cannot be!\n"
        if (len(ter)+len(nod) == 1):
            report+="WARNING! Isolated group. Must be a mistake!\n"
        if (len(ter)==1):
            report+="WARNING! A single terminal in the group. Most likely a mistake!\n"
        if (len(ter)==0 and len(nod)>0):
            report+="WARNING! Isolated group of nodes. Most likely a mistake!\n"    
    return(report)

#lets run some checks
#Te={'1':['c1'],'2':['c2']}
#Ne={'n':['c1','c2'],'n1':[]}
#Ce={'c1':[['1'],['n']],'c2':[['2'],['n']]}

#print(connectivityReport(Te, Ne, Ce))