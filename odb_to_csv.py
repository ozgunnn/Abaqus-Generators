from abaqusConstants import *
import csv
odb = session.odbs['C:/temp/Job_A13_geoimp.odb']

#session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
#    NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=("SET-1", ))

session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
    NODAL, ((COMPONENT, 'U1'), )), ), nodeLabels=(('Profile Instance', ('2018', 
    )), ))

#session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
#    NODAL, ((INVARIANT, 'Magnitude'), )), ), nodeSets=(
#    "REFERENCE_POINT_        2", ))

session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
    NODAL, ((COMPONENT, 'RF3'), )), ), nodeLabels=(('ASSEMBLY', ('1', )), ))

#session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
#    NODAL, ((INVARIANT, 'Magnitude'), )), ), nodeLabels=(('ASSEMBLY', ('2', )), 
#    ))
xy1 = session.xyDataObjects['U:U1 PI: Profile Instance N: 2018']
xy2 = session.xyDataObjects['RF:RF3 PI: ASSEMBLY N: 1']

with open('test.csv','wb') as f:
    writer = csv.writer(f)
    for i in range(len(xy1)):
        writer.writerow([abs(xy1[i][1]),abs(xy2[i][1])])