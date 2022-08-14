"""
@author: DianaNtz
"""
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
#some initial values
t0=0
#tfinal=95#equal mass
tfinal=115#not equal mass
dt=0.005
m1=10**(25)#2*10**(26)
m2=2*10**(26)
G=6.67259*10**-11
steps=int((tfinal-t0)/dt)
t=np.empty(steps, dtype='double')
tn=t0
def fx1(x1,y1,z1,x2,y2,z2):
    return -m2*G*(x1-x2)*1/np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**3
def fy1(x1,y1,z1,x2,y2,z2):
    return -m2*G*(y1-y2)*1/np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**3    
def fz1(x1,y1,z1,x2,y2,z2):
    return -m2*G*(z1-z2)*1/np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**3    
def fx2(x1,y1,z1,x2,y2,z2):
    return m1*G*(x1-x2)*1/np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**3
def fy2(x1,y1,z1,x2,y2,z2):
    return m1*G*(y1-y2)*1/np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**3    
def fz2(x1,y1,z1,x2,y2,z2):
    return m1*G*(z1-z2)*1/np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**3 
x1=np.empty(steps+1, dtype='double')  
y1=np.empty(steps+1, dtype='double')  
z1=np.empty(steps+1, dtype='double') 
x2=np.empty(steps+1, dtype='double')  
y2=np.empty(steps+1, dtype='double')  
z2=np.empty(steps+1, dtype='double') 
x10=0
y10=0
z10=0
vx10=-7.5*1000
vy10=-20*1000
vz10=-15*1000
x20=3000*1000
y20=0
z20=0
vx20=7.5*1000*m1/m2
vy20=20*1000*m1/m2
vz20=15*1000*m1/m2
x1n=x10
y1n=y10
z1n=z10
vx1n=vx10
vy1n=vy10
vz1n=vz10
x2n=x20
y2n=y20
z2n=z20
vx2n=vx20
vy2n=vy20
vz2n=vz20
