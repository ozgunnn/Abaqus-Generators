import math
import numpy as np

axis='Strong'
e=15.0 #load eccentricity
ez=50.0 #rp distance from edge
tf=13.0 #flange thickness
tw=8.0 #web thickness
h=160.0 #section depth
t=10.0 #analysis time
ms=100 #mass scale
u=10.0 #assigned deformation
b=160.0 #section width
d=300.0 #concrete width/depth
cc=20.0  #clear cover
L=3000.0 #extrude length
#
lrd=20.0 #longitudinal rebar diameter
nr=12 #number of lrebars
std=8.0 #stirrup diameter
sts=100.0 #stirrup spacing
fs=600*1.1 #rebaryield
fy=460*1.1 #profileyield
fcm=50 #concrete compressive strength
#

steel_density=7.85e-9 #steel density
Esteel_profile=210000 #Esteelprofile
steel_poisson=0.3 #steel poisson

Esteel_rebar=200000 #Esteelrebar

fu=1.4*fy
esh=min(max(0.015,0.1*fy/fu-0.055),0.03)
eu=0.06
C1=(esh+0.25*(eu-esh))/eu
C2=(esh+0.4*(eu-esh))/eu
Esh=(fu-fy)/(C2*eu-esh)
Ey=210000
ey=fy/Ey
fc1eu=fy+(C1*eu-esh)*Esh

ec1=0.001*min(0.7*fcm**0.31,2.8)
ecu1=0.001*min(3.5,2.8+27*((98-fcm)/float(100))**4)
Ecm=1000*(22*(fcm/10)**0.3)
conc_denst=2.5e-9
conc_poisson=0.2

if fcm<58:
    fctm=0.3*(fcm-8)**(2.0/3)
else:
    fctm=2.12*math.log(1+fcm/10)
    
gfi=1.0/1000*73*fcm**(0.18)
    
ec=np.linspace(0,ecu1,num=20)

sigc=[]
k=1.05*Ecm*ec1/fcm
eta=ec/ec1
sigc=fcm*(k*eta-np.power(eta,2))/(1+(k-2)*eta)
sigc=np.append(sigc,fcm/10)
ec=np.append(ec,ecu1+0.0005)
sigc.tolist

cctable=np.zeros((len(ec),2))

for i in range(len(ec)):
    cctable[i][0]=sigc[i]
    cctable[i][1]=ec[i]-cctable[i][0]/Ecm


i=0
while i in range(len(cctable[:,1])):
    if cctable[i][0]>=0.4*fcm and cctable[i][1]>=0:
        break
    else:
        i=i+1

    
cuttab=cctable[i:]

i=0
while i in range(len(cuttab[:,1])):
    if cuttab[i][0]>fcm:
        cuttab=np.delete(cuttab,i,axis=0)
    else:
        i=i+1
        
sigma_eps=np.zeros((5,2))
sigma_eps[0,:]=(0,0)
sigma_eps[1,:]=(ey-ey,fy)
sigma_eps[2,:]=(esh-ey,fy)
sigma_eps[3,:]=(C1*eu-ey,fc1eu)
sigma_eps[4,:]=(eu-ey,fu)

table_st_eps=np.zeros((4,2))
table_st_eps[:,0]=sigma_eps[1:,1]
table_st_eps[:,1]=sigma_eps[1:,0]

table_st_fin=tuple(map(tuple, table_st_eps))

i=0
while i in range(len(cuttab[:,1])):
    cuttab[i][1]=(cuttab[i][1]-cuttab[0][1])
    i=i+1
    
cctab=tuple(map(tuple, cuttab))

