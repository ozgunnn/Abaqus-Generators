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
        print[i]