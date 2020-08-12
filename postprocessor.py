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

odb = session.openOdb(name='C:/temp/Job-try.odb')

allNodes = mdb.models["Model-1"].rootAssembly.nodes

delta = 1.0
(x,y,z)=(0.0,0.0,0.0)

xmin, ymin, zmin = x-delta, y-delta, z-delta
xmax, ymax, zmax = x+delta, y+delta, z+delta

myNode1 = allNodes.getByBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)

(x,y,z)=(0.0,0.0,120.0)
xmin, ymin, zmin = x-delta, y-delta, z-delta
xmax, ymax, zmax = x+delta, y+delta, z+delta

myNode2 = allNodes.getByBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)

##session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', NODAL, ((COMPONENT, 'RF2'), )), ), nodePick=(('ASSEMBLY', 1, ('[#2 ]', )),  ), )
#
#session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
#    NODAL, ((COMPONENT, 'RF2'), )), ), nodeLabels=(('ASSEMBLY', ('1', )), ))

#session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
#    NODAL, ((COMPONENT, 'U2'), )), ),nodeLabels=(('ASSEMBLY', (str(myNode1), )), ))

#session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
#    NODAL, ((COMPONENT, 'RF2'), )), ),nodeLabels=(('ASSEMBLY', (str(myNode2), )), ))

#xy1 = session.xyDataObjects['U:U2 PI: ASSEMBLY N: '+str(myNode1)+'']
#xy2 = session.xyDataObjects['RF:RF2 PI: ASSEMBLY N: '+str(myNode2)+'']
#xy3 = combine(abs(xy1), abs(xy2))
#xyp = session.xyPlots['XYPlot-1']
#chartName = xyp.charts.keys()[0]
#chart = xyp.charts[chartName]
#c1 = session.Curve(xyData=xy3)
#chart.setValues(curvesToPlot=(c1, ), )