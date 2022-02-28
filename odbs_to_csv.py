from abaqusConstants import *
from abaqus import *
import csv
from odbAccess import *
from caeModules import *

lis=['A11','A13']
l = []

for i in range(182):
    l.append([])

for j in lis: 
    l[0].append(j+'U')
    l[0].append(j+'F')

for j in lis: 
    od = openOdb(path='C:/temp/'+'Job_'+j+'_geoimp.odb')
    o7 = session.odbs['C:/temp/'+'Job_'+j+'_geoimp.odb']
    session.viewports['Viewport: 1'].setValues(displayedObject=o7)

    if j=='A11':
        session.xyDataListFromField(odb=od, outputPosition=NODAL, variable=(('U', 
        NODAL, ((COMPONENT, 'U1'), )), ), nodeLabels=(('Profile Instance', ('2018', 
        )), ))
    else:
        session.xyDataListFromField(odb=od, outputPosition=NODAL, variable=(('U', 
        NODAL, ((COMPONENT, 'U2'), )), ), nodeLabels=(('Profile Instance', ('2018', 
        )), ))        

    session.xyDataListFromField(odb=od, outputPosition=NODAL, variable=(('RF', 
    NODAL, ((COMPONENT, 'RF3'), )), ), nodeLabels=(('ASSEMBLY', ('1', )), ))

    if j=='A11':
        xy1 = session.xyDataObjects['U:U1 PI: Profile Instance N: 2018']
    else:
        xy1 = session.xyDataObjects['U:U2 PI: Profile Instance N: 2018']

    xy2 = session.xyDataObjects['RF:RF3 PI: ASSEMBLY N: 1']

    for i in range(len(xy1)):
        l[i+1].append(abs(xy1[i][1]))
        l[i+1].append(abs(xy2[i][1])/1000)

    session.viewports['Viewport: 1'].setValues(displayedObject=o7)
    del session.xyDataObjects['RF:RF3 PI: ASSEMBLY N: 1']

    if j=='A11':
        del session.xyDataObjects['U:U1 PI: Profile Instance N: 2018']
    else:
        del session.xyDataObjects['U:U2 PI: Profile Instance N: 2018']

    session.odbs['C:/temp/'+'Job_'+j+'_geoimp.odb'].close()
    
with open('test.csv','wb') as f:
    writer = csv.writer(f)
    for i in range(len(l)):
         writer.writerow(l[i])
