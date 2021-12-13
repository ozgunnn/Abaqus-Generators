from odbAccess import *

odb = openOdb(path='C:/temp/Job_B14.odb')

mycenter = odb.rootAssembly.nodeSets['REFERENCE_POINT_        1']

max=0
#for i in range(len(odb.steps['Step-1'].frames)):
for i in range(100):
    #u = odb.steps['Step-1'].frames[i].fieldOutputs['U'].getSubset(region=mycenter).values[0].data[2]
    f = abs(odb.steps['Step-1'].frames[i].fieldOutputs['RF'].getSubset(region=mycenter).values[0].data[2])
    if f > max:
        max = f

f