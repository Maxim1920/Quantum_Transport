# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 17:43:10 2022

@author: Yuli
"""
import numpy as np

#this function is to write the data to disk. has nothing to do with the code
def printall (filename,liste,N,chunk=0,extraline="\n"): 
    def stringrep(a):
        if type(a) is list:
            ll=[str(x) for x in a]
            s=" ".join(ll)
        else:
            s=str(a)
        return(s)
    allines=[]
    for i in range(N):
        if (i>0)and(chunk>0)and(i%chunk==0):
            allines.append(extraline)
        line=""
        for c in liste:
            line=line+" "+stringrep(c[i])
        allines.append(line+"\n")
    f=open(filename,"w")
    f.writelines(allines)
    f.close()
    
def twoblock(a,bvec):
    return np.array([[a+bvec[2],bvec[0]-1j*bvec[1]],[bvec[0]+1j*bvec[1],a-bvec[2]]])

def diablock(a):
    return twoblock(a,np.zeros(3))

def find_bracket(fu,a0,b0,c0): #a<b<c
    a,b,c=a0,b0,c0
    i=0
    while(i<20):
        A,B,C=fu(a),fu(b),fu(c)
        if((B<A)and(B<C)):
            break
        if(A>C):
            cn,an,bn=2*c-b,b,c
            a,b,c=an,bn,cn
        else:
            an,cn,bn=(2*a-b),b,a
            a,b,c=an,bn,cn
        i=i+1
    return a,b,c,i
            
def tracefunction(file,fun,linsp):
    a=np.apply_along_axis(fun, 0, linsp)  
    printall(file,[linsp,a],len(linsp))     
    
#tracefunction("sin.dat",np.sin,np.linspace(0,2*np.pi,num=1000))

def gnulines(fil, newf=None):
    #insert newlines required by gnu grid
    f=open(fil,'r')
    s=''
    #reads and appends first line
    line =f.readline()
    first=line.split()[0]
    s+=line
    while True:
        line=f.readline()
        if not line:
            break
        if not (first==line.split()[0]):
            s+="\n"
            first=line.split()[0]
        s+=line
    f.close()
    if (newf):
        out=newf
    else:
        out=fil+"gnu.dat"
    f=open(out,'w')
    f.write(s)
    f.close()
    
#gnulines('data\\484.dat','484.dat')
    
def filterit(fil, nfil, n,  upper, down): #assuming numbers everywhere
    f=open(fil,'r')
    lines=f.readlines()
    f.close()
    s=[]
    for lin in lines:
        first=float(lin.split()[n])
        if (first<upper)and(first>down):
            s.append(lin)
    f=open(nfil,'w')
    f.writelines(s)
    f.close()

def measuregrid(fil):
    f=open(fil,'r')
    s=f.read()
    f.close()
    first=s.split()[0]
    A=s.count(first)
    B=s.count('\n')
    return(int(B/A),A)
        
#filterit('data\\484.dat','484.dat',1,2.305,2.28)
#gnulines('484.dat','484f.dat')
#print(measuregrid('484.dat'))

def griddata(fil):
    A,B=measuregrid(fil)
    data=np.zeros((A,B,3),dtype=float)
    f=open(fil,'r')
    for i in range(A):
        for j in range(B):
            s=f.readline()
            li=s.split()
            for k in range(3):
                data[i][j][k]=float(li[k])
    return(data)

#D=griddata('484.dat')
#print(D[3][5])

        
