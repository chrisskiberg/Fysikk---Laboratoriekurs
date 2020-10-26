# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 10:24:50 2020

@author: Andre
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
import csv


# Horisontal avstand mellom festepunktene er 200 mm
h = 200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])
# Vi begrenser starthøyden (og samtidig den maksimale høyden) til
# å ligge mellom 250 og 300 mm
ymax = 300
yfast=np.asarray([251, 183, 205, 199, 140,  65, 144, 129])

# inttan: tabell med 7 verdier for (yfast[n+1]-yfast[n])/h (n=0..7); dvs
# banens stigningstall beregnet med utgangspunkt i de 8 festepunktene.
inttan = np.diff(yfast)/h
attempts=1

# Omregning fra mm til m:
xfast = xfast/1000
yfast = yfast/1000

cs = CubicSpline(xfast, yfast, bc_type='natural')
xmin = 0.000
xmax = 1.401
dx = 0.001
x = np.arange(xmin, xmax, dx)
Nx = len(x)
y = cs(x)
dy = cs(x,1)
b = np.arctan(dy);
d2y = cs(x,2)

g=9.81;
c= 2/5;

"""
npX =np.array(ekspX)
npY =np.array(ekspY)
dy = np.zeros(len(npY)-2)
for i in range(len(dy)):
    dy[i] = (npY[i+2]-npY[i])/(npX[i+2]-npX[i])
npX2 = np.delete(npX, 0)
npX2 = np.delete(npX2, -1)
d2y = np.zeros(len(dy)-2)
for i in range(len(d2y)):#skrive at det er lav numerisk presisjon om det ikke passer
    d2y[i] = (dy[i+2]-dy[i])/(npX2[i+2]-npX2[i])*/"""




#finner fart
v_x=[]
for i in range(1401):
    v_x.append(((2*g*(0.251-y[i]))/(1+c))**(1/2))
    

#horisontal fart
horzV = v_x*np.cos(b)    
#gjennomsnittlig horisontal fart
avrghorzV =  []
for i in range(1, 1401):
    avrghorzV.append(0.5*(horzV[i]+horzV[i-1]))
#endring i X verdier
deltT=[]
deltX = []
for i in range(1, 1401):
    deltX.append(x[i]-x[i-1])
    
#endring i tid for hver x
deltT.append(0);
for i in range(0, 1400):
    deltT.append(deltX[i]/avrghorzV[i])
deltT[0]=0.01
#nytt tidsutrykk
t=[]
t.append(0);
for i in range(1,1400):
    t.append(t[i-1]+deltT[i+1])
t.append(t[1399])
print(horzV)


#krumning:
krumning = (d2y)/((1+(dy)**2))**(3/2)


#regner ut normalkraften
sentripetal_aks = []
for i in range(1401):
    sentripetal_aks.append(((v_x[i])**2)*krumning[i])
normalkraft = (sentripetal_aks+g*np.cos(b))
#regner ut friksjonskraften
friksjon = (c*g*np.sin(b))/(1+c)
f_N = friksjon/normalkraft

#praktisk
header = []
data = []



filename = "fysikk_3.csv"
with open(filename) as csvfile:
    csvreader = csv.reader(csvfile)

    header = next(csvreader)

    for datapoint in csvreader:

        values = [float(value) for value in datapoint]
        data.append(values)
header[1] ="Simulert"
header[2] ="Eksperimentelt"
#header.append("Eksperimentelt")



time = [p[0]  for p in data]
ekspX = [p[1]+0.025 for p in data]
ekspY = [p[2]-0.02 for p in data]




npX =np.array(ekspX)
npY =np.array(ekspY)
dy = np.zeros(len(npY)-2)
for i in range(len(dy)):
    dy[i] = (npY[i+2]-npY[i])/(npX[i+2]-npX[i])
npX2 = np.delete(npX, 0)
npX2 = np.delete(npX2, -1)
d2y = np.zeros(len(dy)-2)
for i in range(len(d2y)):#skrive at det er lav numerisk presisjon om det ikke passer
    d2y[i] = (dy[i+2]-dy[i])/(npX2[i+2]-npX2[i])

helnig= np.arctan(dy)

v_x_eksp=[]
for i in range(len(npY)):
    v_x_eksp.append(((2*g*(0.251-npY[i]))/(1+c))**(1/2))

average = []
for i in range(len(helnig)-2):
    average.append((helnig[i+2]+helnig[i+1]+helnig[i])/3)
    
averagex = []
for i in range(len(npX2)-2):
    averagex.append((npX2[i+2] + npX2[i+1]+npX2[i])/3)
    
print(len(d2y))
dy =np.delete(dy,0)
dy = np.delete(dy,-1)
krumning2 = (d2y)/((1+(dy)**2))**(3/2)


#Plotting

baneform = plt.figure('y(x)',figsize=(12,3))
plt.plot(x,y, ekspX, ekspY,xfast,yfast,'*')
plt.legend(header[1:],  bbox_transform=plt.gcf().transFigure)
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$y$(m)',fontsize=20)
plt.ylim(0,0.350)
plt.grid()
plt.show()



helningsvinkel = plt.figure('b(x)',figsize=(12,3))
plt.plot(x,b, npX2, helnig)
plt.title('Banens helningsvinkel')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$\u03B2$ (radianer)',fontsize=20)
plt.ylim(-1,1)
plt.grid()
plt.show()

helningsvinkel = plt.figure('b(x)',figsize=(12,3))
plt.plot(x,b, averagex, average)
plt.legend(header[1:],  bbox_transform=plt.gcf().transFigure)
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$\u03B2$ (radianer)',fontsize=20)
plt.ylim(-1,1)
plt.grid()
plt.show()


npX2 = np.delete(npX2,0)
npX2 = np.delete(npX2,-1)
banekrumning = plt.figure('k(x)',figsize=(12,3))
plt.plot(x,krumning, npX2,  krumning2)
plt.title('Banens krummningsradius')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$k$(1/m)',fontsize=20)
plt.grid()
plt.show()




fart = plt.figure('v(x)',figsize=(12,3))
plt.plot(x,v_x, npX, v_x_eksp)
plt.legend(header[1:],  bbox_transform=plt.gcf().transFigure)
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$v(m/s)$',fontsize=20)
plt.ylim(0,2)
plt.grid()
plt.show()

print(normalkraft)
normalkraften = plt.figure('N(N)',figsize=(12,3))
plt.plot(x,normalkraft)
plt.title('Normalkraft')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$N(N)$',fontsize=20)
plt.grid()
plt.show()

friksjonskraft = plt.figure('|f(N)|',figsize=(12,3))
plt.plot(x,abs(f_N))
plt.title('Friksjonskraft')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$|f(N)|$',fontsize=20)
plt.grid()
plt.show()

distansetid = plt.figure('x(m)',figsize=(12,3))
plt.plot(t,x,time, ekspX)
plt.legend(header[1:],  bbox_transform=plt.gcf().transFigure)
plt.xlabel('$t$ (s)',fontsize=20)
plt.ylabel('$x(m)$',fontsize=20)
plt.grid()
plt.show()

vt =[]
vt.append(0)
for i in range(1400):
    vt.append(deltX[i]/deltT[i])

aks = np.zeros(len(vt)-2)

npT = np.array(t)
npT = np.delete(npT,0)
npT = np.delete(npT,0)
for i in range(len(d2y)):#skrive at det er lav numerisk presisjon om det ikke passer
    aks[i] = (vt[i+2]-vt[i])/(npT[i+2]-npT[i])
print(vt[1400])


print(aks)
farttid = plt.figure('v(m/s)',figsize=(12,3))
plt.plot(t,vt)
plt.legend(header[1:],  bbox_transform=plt.gcf().transFigure)#simulert og eksperimentell
plt.xlabel('$t$ (s)',fontsize=20)
plt.ylabel('$v(m/s)$',fontsize=20)
plt.grid()
plt.show()

akselerasjon = plt.figure('v(m/s)',figsize=(12,3))
plt.plot(npT,aks)
plt.legend(header[1:],  bbox_transform=plt.gcf().transFigure)#simulert og eksperimentell
plt.xlabel('$t$ (s)',fontsize=20)
plt.ylabel('$v(m/s)$',fontsize=20)
plt.grid()
plt.show()


baneform.savefig("baneform.pdf", bbox_inches='tight')
banekrumning.savefig("banekrumning.pdf", bbox_inches='tight')
helningsvinkel.savefig("helningsvinkel.pdf", bbox_inches='tight')
fart.savefig("fart.pdf", bbox_inches='tight')
normalkraften.savefig("normalkraf.pdf", bbox_inches='tight')
friksjonskraft.savefig("friksjonskraft.pdf", bbox_inches='tight')
distansetid.savefig("distansetid.pdf", bbox_inches='tight')
farttid.savefig("farttid.pdf", bbox_inches='tight')