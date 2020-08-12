from abaqus import *
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *

odb = session.openOdb(name='C:/temp/Job-v122-e200.odb')

rpcoordlb=odb.rootAssembly.nodes[0].label

a=[0.0,0.0,1500.0]

distin=a[2]
for i in range(len(odb.rootAssembly.instances['PROFILESHELLPART-1'].nodes)):
    x1=odb.rootAssembly.instances['PROFILESHELLPART-1'].nodes[i].coordinates[0]
    x2=odb.rootAssembly.instances['PROFILESHELLPART-1'].nodes[i].coordinates[1]
    x3=odb.rootAssembly.instances['PROFILESHELLPART-1'].nodes[i].coordinates[2]
    distdum=sqrt((x1-a[0])**2+(x2-a[1])**2+(x3-a[2])**2)
    if distdum<distin:
        distin=distdum
        label=odb.rootAssembly.instances['PROFILESHELLPART-1'].nodes[i].label

session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', NODAL, ((COMPONENT, 'RF3'), )), ), nodeLabels=(('ASSEMBLY', ('1', )), ))
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', NODAL, ((COMPONENT, 'U2'), )), ), nodeLabels=(('PROFILESHELLPART-1', (str(label), )), ))
#

#xy1 = session.xyDataObjects['U:U2 PI: ASSEMBLY N: '+str(myNode1)+'']
#xy2 = session.xyDataObjects['RF:RF3 PI: ASSEMBLY N: '1']
#xy3 = combine(abs(xy1), abs(xy2))
#xyp = session.xyPlots['XYPlot-1']
#chartName = xyp.charts.keys()[0]
#chart = xyp.charts[chartName]
#c1 = session.Curve(xyData=xy3)
#chart.setValues(curvesToPlot=(c1, ), )