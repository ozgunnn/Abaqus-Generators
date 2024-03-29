from abaqus import *
from abaqusConstants import *
import regionToolset
import math
import numpy as np
import copy
import sketch
import part
import section
import material
import assembly
import mesh 
import regionToolset
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *

name='B14'  #comment when default name is used , line 136
perfectstep=1 #only one of these can be 1
bucklestep=0 #only one of these can be 1
geoimpstep=0 #only one of these can be 1, imp values defined at the end of script
resstrstep=0 #switch, 1 or 0 regardless of others
axis='Strong'  #Strong or Weak
shape='Rect' #Circular or Rect
e=15.0 #load eccentricity
ez=25.0 #rp distance from edge

t=10.0 #analysis time
ms=100 #mass scale
u=10.0 #assigned deformation

b=160.0 #section width
h=160.0 #section depth
tf=13.0 #flange thickness
tw=8.0 #web thickness

d=350.0 #concrete diameter (if circular)

bc=300.0 #concrete width (if rectangle)
hc=300.0 #concrete depth (if rectangle)

cc=20.0  #clear cover

L=6000.0 #extrude length (full length)
nocores=6
#
lrd=20.0 #longitudinal rebar diameter
nr=12 #number of lrebars
std=10.0 #stirrup diameter
sts=75.0 #stirrup spacing
fs=500.0 #rebaryield 
fy=525.0 #profileyield 
fcm=92.7 #concrete compressive strength
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

meshsize=30.0

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

session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])

i=0
a=copy.deepcopy(cuttab[0][1])
while i in range(len(cuttab[:,1])):
    cuttab[i][1]=(cuttab[i][1]-a)
    i=i+1
    
cctab=tuple(map(tuple, cuttab))

session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])

#column_model = mdb.models['Model-1']
#name='_'.join([shape[:3],'d',str(int(d)),'fcm',str(int(fcm)),'h',str(int(h)),'fy',str(int(fy)),'sts',str(int(sts)),axis])

if perfectstep==1:
    column_model=mdb.Model(name=name, modelType=STANDARD_EXPLICIT)
elif bucklestep==1:
    column_model=mdb.Model(name=name+'_buckle', modelType=STANDARD_EXPLICIT)
elif geoimpstep==1:
    column_model=mdb.Model(name=name+'_geoimp', modelType=STANDARD_EXPLICIT)

import material

concmaterial=column_model.Material(name='mat_concrete')
concmaterial.ConcreteDamagedPlasticity(table=((35,0.1,1.16,0.67,0),))
concmaterial.concreteDamagedPlasticity.ConcreteCompressionHardening(table=cctab)
concmaterial.concreteDamagedPlasticity.ConcreteTensionStiffening(table=((fctm, gfi), ), type=GFI)
concmaterial.Density(table=((conc_denst, ), ))
concmaterial.Elastic(table=((Ecm, conc_poisson), ))

column_model.Material(name='mat_profile')
column_model.materials['mat_profile'].Density(table=((steel_density, ), ))
column_model.materials['mat_profile'].Elastic(table=((Esteel_profile, steel_poisson), ))
column_model.materials['mat_profile'].Plastic(table=table_st_fin)

column_model.Material(name='mat_rebar')
column_model.materials['mat_rebar'].Density(table=((steel_density, ), ))
column_model.materials['mat_rebar'].Elastic(table=((Esteel_rebar, steel_poisson), ))
column_model.materials['mat_rebar'].Plastic(table=((fs,0),))

import sketch
import part

concrete_sketch=column_model.ConstrainedSketch(name='beam',sheetSize=160)

concrete_sketch.Line(point1=(b/2,h/2),point2=(b/2,h/2-tf))
concrete_sketch.Line(point1=(b/2,h/2-tf),point2=(tw/2,h/2-tf))
concrete_sketch.Line(point1=(tw/2,h/2-tf),point2=(tw/2,-(h/2-tf)))
concrete_sketch.Line(point1=(tw/2,-(h/2-tf)),point2=(b/2,-(h/2-tf)))
concrete_sketch.Line(point1=(b/2,-(h/2-tf)),point2=(b/2,-h/2))
concrete_sketch.Line(point1=(b/2,-h/2),point2=(-b/2,-h/2))
concrete_sketch.Line(point1=(-b/2,-h/2),point2=(-b/2,-(h/2-tf)))
concrete_sketch.Line(point1=(-b/2,-(h/2-tf)),point2=(-tw/2,-(h/2-tf)))
concrete_sketch.Line(point1=(-tw/2,-(h/2-tf)),point2=(-tw/2,(h/2-tf)))
concrete_sketch.Line(point1=(-tw/2,(h/2-tf)),point2=(-b/2,(h/2-tf)))
concrete_sketch.Line(point1=(-b/2,(h/2-tf)),point2=(-b/2,h/2))
concrete_sketch.Line(point1=(-b/2,h/2),point2=(b/2,h/2))

if shape == "Rect":
    concrete_sketch.rectangle(point1=(bc/2, hc/2), point2=(-bc/2, -hc/2))
elif shape == "Circular":
    concrete_sketch.CircleByCenterPerimeter((0,0),(0,d/2))

concrete_part=column_model.Part(name='concrete',dimensionality=THREE_D,type=DEFORMABLE_BODY)
concrete_part.BaseSolidExtrude(sketch=concrete_sketch,depth=L)
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(b/2,h/2,6),point2=(b/2,-h/2,5),point3=(b/2,-h/2,4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,h/2,6),point2=(-b/2,-h/2,5),point3=(-b/2,-h/2,4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,h/2,6),point2=(b/2,h/2,5),point3=(-b/2,h/2,4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,-h/2,6),point2=(b/2,-h/2,5),point3=(-b/2,-h/2,4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,-(h/2-tf),6),point2=(b/2,-(h/2-tf),5),point3=(-b/2,-(h/2-tf),4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,(h/2-tf),6),point2=(b/2,(h/2-tf),5),point3=(-b/2,(h/2-tf),4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(0,(h/2-tf),6),point2=(0,-(h/2-tf),5),point3=(0,(h/2-tf),4))

beam_sketch=column_model.ConstrainedSketch(name='beam',sheetSize=160)

beam_sketch.Line(point1=(-b/2,(h-tf)/2),point2=(0,(h-tf)/2))
beam_sketch.Line(point1=(0,(h-tf)/2),point2=(b/2,(h-tf)/2))
beam_sketch.Line(point1=(-b/2,-(h-tf)/2),point2=(0,-(h-tf)/2))
beam_sketch.Line(point1=(0,-(h-tf)/2),point2=(b/2,-(h-tf)/2))
beam_sketch.Line(point1=(0,-1*(h-tf)/2),point2=(0,(h-tf)/2))

beam_part=column_model.Part(name='Beam',dimensionality=THREE_D,type=DEFORMABLE_BODY)
beam_part.BaseShellExtrude(sketch=beam_sketch,depth=L)

rebar_sketch=column_model.ConstrainedSketch(name='rebar',sheetSize=160)
rebar_sketch.Line(point1=(0,0),point2=(L,0))
rebar_part=column_model.Part(name='Rebar',dimensionality=THREE_D,type=DEFORMABLE_BODY)
rebar_part.BaseWire(sketch=rebar_sketch)

stirrup_sketch=column_model.ConstrainedSketch(name='rebar',sheetSize=160)
if shape=="Rect":
    stirrup_sketch.rectangle((bc/2-cc-std/2,hc/2-cc-std/2),(-(bc/2-cc-std/2),-(hc/2-cc-std/2)))
elif shape=="Circular":
    stirrup_sketch.CircleByCenterPerimeter((0,0),(0,d/2-cc-std/2))
stirrup_part=column_model.Part(name='Stirrup',dimensionality=THREE_D,type=DEFORMABLE_BODY)
stirrup_part.BaseWire(sketch=stirrup_sketch)

import section

column_model.HomogeneousShellSection(name='sec_flange', preIntegrate=OFF, material='mat_profile', thicknessType=UNIFORM, thickness=tf, thicknessField='', nodalThicknessField='', idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, useDensity=OFF,integrationRule=SIMPSON, numIntPts=5)
column_model.HomogeneousShellSection(name='sec_web', preIntegrate=OFF, material='mat_profile', thicknessType=UNIFORM, thickness=tw, thicknessField='', nodalThicknessField='', idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, useDensity=OFF,integrationRule=SIMPSON, numIntPts=5)

f=beam_part.faces
if shape=="Rect":
    faces_fl=f.findAt(((b/4.0,h/2.0-tf/2.0,L/4),),((-b/4.0,h/2.0-tf/2.0,L/4),),((-b/4.0,-(h/2.0-tf/2.0),L/4),),((b/4.0,-(h/2.0-tf/2.0),L/4),))
elif shape=="Circular":
    faces_fl=f.findAt(((b/4,h/2-tf/2,L/4),),((-b/4,h/2-tf/2,L/4),),((-b/4,-(h/2-tf/2),L/4),),((b/4,-(h/2-tf/2),L/4),))

faces_web=f.findAt((0,0,L/4),)

region_fl = (faces_fl,)
region_web = (faces_web,)

beam_part.SectionAssignment(region=region_fl, sectionName='sec_flange', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)
beam_part.SectionAssignment(region=region_web, sectionName='sec_web', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)

conc_sec=column_model.HomogeneousSolidSection(name='Concrete_Section',material='mat_concrete')
conc_region=(concrete_part.cells,)
concrete_part.SectionAssignment(region=conc_region,sectionName='Concrete_Section')

column_model.CircularProfile(name='longrebar_profile', r=lrd/2)
column_model.BeamSection(name='Lrebar_Section', integration=DURING_ANALYSIS, poissonRatio=0.0, profile='longrebar_profile', material='mat_rebar', temperatureVar=LINEAR, consistentMassMatrix=False)
rebar_edges=(rebar_part.edges,)
rebar_part.SectionAssignment(region=rebar_edges,sectionName='Lrebar_Section')

column_model.CircularProfile(name='stirrup_profile', r=std/2)
column_model.BeamSection(name='Stirrup_Section', integration=DURING_ANALYSIS, poissonRatio=0.0, profile='stirrup_profile', material='mat_rebar', temperatureVar=LINEAR, consistentMassMatrix=False)
stirrup_edges=(stirrup_part.edges,)
stirrup_part.SectionAssignment(region=stirrup_edges,sectionName='Stirrup_Section')

import assembly

columnAssembly=column_model.rootAssembly
concreteInstance=columnAssembly.Instance(name='Concrete Instance',part=concrete_part,dependent=ON)
profileInstance=columnAssembly.Instance(name='Profile Instance',part=beam_part,dependent=ON)

if shape=="Rect":
    lrebarInstance=columnAssembly.Instance(name='Lrebar Instance',part=rebar_part,dependent=ON)
    lrebarInstance=columnAssembly.Instance(name='Lrebar Instance1',part=rebar_part,dependent=ON)
    lrebarInstance=columnAssembly.Instance(name='Lrebar Instance2',part=rebar_part,dependent=ON)
    lrebarInstance=columnAssembly.Instance(name='Lrebar Instance3',part=rebar_part,dependent=ON)
elif shape=="Circular":
    lrebarInstance=columnAssembly.Instance(name='Lrebar Instance',part=rebar_part,dependent=ON)

stirrupInstance=columnAssembly.Instance(name='Stirrup Instance',part=stirrup_part,dependent=ON)
columnAssembly.regenerate()

if shape=="Rect":
    columnAssembly.rotate(instanceList=('Lrebar Instance', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, -1.0, 0.0), angle=90.0)
    columnAssembly.rotate(instanceList=('Lrebar Instance1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, -1.0, 0.0), angle=90.0)
    columnAssembly.rotate(instanceList=('Lrebar Instance2', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, -1.0, 0.0), angle=90.0)
    columnAssembly.rotate(instanceList=('Lrebar Instance3', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, -1.0, 0.0), angle=90.0)

    columnAssembly.translate(instanceList=('Lrebar Instance', ), vector=(bc/2-cc-std-lrd/2, hc/2-cc-std-lrd/2, 0.0))
    columnAssembly.translate(instanceList=('Lrebar Instance1', ), vector=(bc/2-cc-std-lrd/2, -(hc/2-cc-std-lrd/2), 0.0))
    columnAssembly.translate(instanceList=('Lrebar Instance2', ), vector=(-bc/2+cc+std+lrd/2, -hc/2+cc+std+lrd/2, 0.0))
    columnAssembly.translate(instanceList=('Lrebar Instance3', ), vector=(-bc/2+cc+std+lrd/2, hc/2-cc-std-lrd/2, 0.0))

    columnAssembly.translate(instanceList=('Stirrup Instance', ), vector=(0.0, 0.0, sts/2))

    columnAssembly.LinearInstancePattern(instanceList=('Lrebar Instance', ), direction1=(0.0, -1.0, 0.0), direction2=(0.0, 1.0, 0.0), number1=int(nr/4), number2=1, spacing1=(hc-2*cc-2*std-lrd)/(nr/4), spacing2=sts)
    columnAssembly.LinearInstancePattern(instanceList=('Lrebar Instance1', ), direction1=(-1.0, 0.0, 0.0), direction2=(0.0, 1.0, 0.0), number1=int(nr/4), number2=1, spacing1=(bc-2*cc-2*std-lrd)/(nr/4), spacing2=sts)
    columnAssembly.LinearInstancePattern(instanceList=('Lrebar Instance2', ), direction1=(0.0, 1.0, 0.0), direction2=(0.0, 1.0, 0.0), number1=int(nr/4), number2=1, spacing1=(hc-2*cc-2*std-lrd)/(nr/4), spacing2=sts)
    columnAssembly.LinearInstancePattern(instanceList=('Lrebar Instance3', ), direction1=(1.0, 0.0, 0.0), direction2=(0.0, 1.0, 0.0), number1=int(nr/4), number2=1, spacing1=(bc-2*cc-2*std-lrd)/(nr/4), spacing2=sts)

    columnAssembly.LinearInstancePattern(instanceList=('Stirrup Instance', ), direction1=(0.0, 0.0, 1.0), direction2=(0.0, 1.0, 0.0), number1=int((L-sts/2)/sts)+1, number2=1, spacing1=sts, spacing2=sts)

    #stirrupsset=columnAssembly.edges.getByBoundingCylinder(center1=(0,d/2,-5) ,center2=(0,d/2,L+5),radius=500)
    #stirrupsset=columnAssembly.instances[findAt(0.0, d/2-cc-std/2, sts/2)]
    #columnAssembly.InstanceFromBooleanMerge(name='Stirrups', instances=(stirrupsset,), originalInstances=DELETE, domain=GEOMETRY)
elif shape=="Circular":
    columnAssembly.rotate(instanceList=('Lrebar Instance', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, -1.0, 0.0), angle=90.0)
    columnAssembly.translate(instanceList=('Lrebar Instance', ), vector=(0.0, d/2-cc-std-lrd/2, 0.0))
    columnAssembly.translate(instanceList=('Stirrup Instance', ), vector=(0.0, 0.0, sts/2))
    columnAssembly.LinearInstancePattern(instanceList=('Stirrup Instance', ), direction1=(0.0, 0.0, 1.0), direction2=(0.0, 1.0, 0.0), number1=int((L-sts/2)/sts)+1, number2=1, spacing1=sts, spacing2=sts)
    columnAssembly.RadialInstancePattern(instanceList=('Lrebar Instance', ), point=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0), number=nr, totalAngle=360.0)

    #stirrupsset=columnAssembly.edges.getByBoundingCylinder(center1=(0,d/2,-5) ,center2=(0,d/2,L+5),radius=500)
    #stirrupsset=columnAssembly.instances[findAt(0.0, d/2-cc-std/2, sts/2)]
    #columnAssembly.InstanceFromBooleanMerge(name='Stirrups', instances=(stirrupsset,), originalInstances=DELETE, domain=GEOMETRY)
    
AllInstances_List = columnAssembly.instances.keys()
stirrupinstances=[]

for i in AllInstances_List:
    if 'Stirrup' in i:
        stirrupinstances.append(i)
    else:
        continue
        
lrebarinstances=[]

for i in AllInstances_List:
    if 'Lrebar' in i:
        lrebarinstances.append(i)
    else:
        continue

    
columnAssembly.InstanceFromBooleanMerge(name='Stirrups', instances=([columnAssembly.instances[stirrupinstances[i]]
    for i in range(len(stirrupinstances))] ), originalInstances=DELETE, domain=GEOMETRY)
    
columnAssembly.InstanceFromBooleanMerge(name='Lrebars', instances=([columnAssembly.instances[lrebarinstances[i]]
    for i in range(len(lrebarinstances))] ), originalInstances=DELETE, domain=GEOMETRY)

p = column_model.parts['Stirrups']
region=regionToolset.Region(edges=p.edges)
p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(1.0, 0.0, -1.0))

column_model.parts['Stirrups'].SectionAssignment(region=region,sectionName='Stirrup_Section')

p = column_model.parts['Lrebars']
region=regionToolset.Region(edges=p.edges)
p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 1.0, -1.0))

column_model.parts['Lrebars'].SectionAssignment(region=region,sectionName='Lrebar_Section')

AllInstances_List = columnAssembly.instances.keys()
stirrupsinstances=[]
lrebarsinstances=[]

for i in AllInstances_List:
    if 'Stirrups' in i:
        stirrupsinstances.append(i)
    else:
        continue
        
for i in AllInstances_List:
    if 'Lrebars' in i:
        lrebarsinstances.append(i)
    else:
        continue     

e1 = columnAssembly.instances[stirrupsinstances[0]].edges
e2 = columnAssembly.instances[lrebarsinstances[0]].edges
edges = (e1,e2,)

column_model.EmbeddedRegion(name='Con_emb_reg', embeddedRegion=edges, hostRegion=None, weightFactorTolerance=1e-06, absoluteTolerance=0.0, fractionalTolerance=0.05, toleranceMethod=BOTH)

if axis=='Strong':
    columnAssembly.ReferencePoint(point=(0, e, -ez))
    refcoord=(0, e, -ez)
    columnAssembly.ReferencePoint(point=(0, e, L+ez))
    refcoord1=(0, e, L+ez)
    
elif axis=='Weak':
    columnAssembly.ReferencePoint(point=(e, 0, -ez))
    refcoord=(e, 0, -ez)
    columnAssembly.ReferencePoint(point=(e, 0, L+ez))
    refcoord1=(e, 0, L+ez)

f1 = columnAssembly.instances['Concrete Instance'].faces
faces5 = f1.getByBoundingBox(zMin=-1,zMax=0)
r1 = columnAssembly.referencePoints.findAt(refcoord)
e1 = columnAssembly.instances['Profile Instance'].edges
edges5 = e1.getByBoundingBox(zMin=-1,zMax=0)

column_model.RigidBody(name='rigid_body', refPointRegion=(r1,), tieRegion=(faces5,edges5,))

f11 = columnAssembly.instances['Concrete Instance'].faces
faces51 = f11.getByBoundingBox(zMin=L,zMax=L+1)
r11 = columnAssembly.referencePoints.findAt(refcoord1)
e11 = columnAssembly.instances['Profile Instance'].edges
edges51 = e11.getByBoundingBox(zMin=L,zMax=L+1)

column_model.RigidBody(name='rigid_body1', refPointRegion=(r11,), tieRegion=(faces51,edges51,))

column_model.ContactProperty('Int_Fric')
column_model.interactionProperties['Int_Fric'].TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((0.5, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005, elasticSlipStiffness=None)
column_model.interactionProperties['Int_Fric'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)

if bucklestep==1:
    column_model.ContactStd(createStepName='Initial', name='Inter')
else:
    column_model.ContactExp(name='Inter', createStepName='Initial')
column_model.interactions['Inter'].includedPairs.setValuesInStep(stepName='Initial', useAllstar=ON)
column_model.interactions['Inter'].contactPropertyAssignments.appendInStep(stepName='Initial', assignments=((GLOBAL, SELF, 'Int_Fric'), ))

#f1 = columnAssembly.instances['Concrete Instance'].faces
#faces1 = f1.getByBoundingBox(zMin=L,zMax=L+1)
#e1 = columnAssembly.instances['Profile Instance'].edges
#edges1 = e1.getByBoundingBox(zMin=L,zMax=L+1)
#column_model.ZsymmBC(name='BC-sim', createStepName='Initial', region=(faces1,edges1,), localCsys=None)

if bucklestep==1:
    column_model.BuckleStep(name='Step-1', numEigen=2, previous='Initial', vectors=4)
else:
    column_model.SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((0.0, 0.0), (t, 1.0)))
    column_model.ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=t, massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, ms, 0.0, None, 0, 0, 0.0, 0.0, 0, None), ), improvedDtMethod=ON)

r1 = columnAssembly.referencePoints.findAt(refcoord)

if bucklestep==1:
    column_model.DisplacementBC(name='BC-2', createStepName='Step-1', region=(r1,), u1=0, u2=0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=0, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    column_model.DisplacementBC(name='BC-3', createStepName='Step-1', region=(r11,), u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=0, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    column_model.ConcentratedForce(cf3=1.0, createStepName='Step-1', distributionType=UNIFORM, field='', localCsys=None , name='Load-1', region=(r1,), )
else:
    column_model.DisplacementBC(name='BC-2', createStepName='Step-1', region=(r1,), u1=0, u2=0, u3=u, ur1=UNSET, ur2=UNSET, ur3=0, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    column_model.DisplacementBC(name='BC-3', createStepName='Step-1', region=(r11,), u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=0, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    column_model.fieldOutputRequests['F-Output-1'].setValues(numIntervals=500)

import mesh 

region_profile_mesh=(f,)
elem_type_profile=mesh.ElemType(elemCode=S4R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=DEFAULT)
beam_part.setElementType(regions=region_profile_mesh, elemTypes=(elem_type_profile, ))
beam_part.seedPart(size=meshsize, deviationFactor=0.1, minSizeFactor=0.1)
beam_part.generateMesh()

region_concrete_mesh=(concrete_part.cells,)
elem_type_concrete=mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, hourglassControl=DEFAULT, distortionControl=DEFAULT)
concrete_part.setElementType(regions=region_concrete_mesh, elemTypes=(elem_type_concrete, ))
concrete_part.seedPart(size=meshsize, deviationFactor=0.1, minSizeFactor=0.1)
concrete_part.generateMesh()

region_stirrups_mesh=(column_model.parts['Stirrups'].edges,)
elem_type_bars=mesh.ElemType(elemCode=B31, elemLibrary=EXPLICIT)
column_model.parts['Stirrups'].setElementType(regions=region_stirrups_mesh, elemTypes=(elem_type_bars, ))
column_model.parts['Stirrups'].seedPart(size=meshsize, deviationFactor=0.1, minSizeFactor=0.1)
column_model.parts['Stirrups'].generateMesh()

region_lrebars_mesh=(column_model.parts['Lrebars'].edges,)
elem_type_bars=mesh.ElemType(elemCode=B31, elemLibrary=EXPLICIT)
column_model.parts['Lrebars'].setElementType(regions=region_lrebars_mesh, elemTypes=(elem_type_bars, ))
column_model.parts['Lrebars'].seedPart(size=meshsize, deviationFactor=0.1, minSizeFactor=0.1)
column_model.parts['Lrebars'].generateMesh()

session.viewports['Viewport: 1'].setValues(displayedObject=columnAssembly)

fy=235.0 #comment when defining the residual stresses based on strengths
if resstrstep==1:
    #fy=100.0
    #el = column_model.parts['Beam'].elements.getByBoundingCylinder(center1=(b/2-15,(h-tf)/2,0),center2=(b/2-15,(h-tf)/2,40),radius =20)
    el = column_model.parts['Beam'].elements.getByBoundingBox(xMin=0, xMax=b/2, yMin=(h-tf)/2-1, yMax=(h-tf)/2+1, zMin=0, zMax=40)

    fl_mesh_a=meshsize
    for i in range(4):
	    #c1=el[0].getElemEdges()[i].getNodes()[0].coordinates
	    #c2=el[0].getElemEdges()[i].getNodes()[1].coordinates
	    #d=sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2)
	    #if d != fl_mesh_a:
	    #	fl_mesh_a = d
	    #	break

        z1=el[0].getElemEdges()[i].getNodes()[0].coordinates[2]
        z2=el[0].getElemEdges()[i].getNodes()[1].coordinates[2]
        x1=el[0].getElemEdges()[i].getNodes()[0].coordinates[0]
        x2=el[0].getElemEdges()[i].getNodes()[1].coordinates[0]
        
        if z1==z2:
            fl_mesh_a= abs(x1-x2)
            break

    #el = column_model.parts['Beam'].elements.getByBoundingCylinder(center1=(0,0,0),center2=(0,0,40),radius =40)
    el = column_model.parts['Beam'].elements.getByBoundingBox(xMin=0, xMax=0, yMin=-h/2, yMax=h/2, zMin=0, zMax=40)

    wb_mesh_a=meshsize
    for i in range(4):
        #c1=el[0].getElemEdges()[i].getNodes()[0].coordinates
        #c2=el[0].getElemEdges()[i].getNodes()[1].coordinates
        #d=sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2)
        #if d != wb_mesh_a:
        #	wb_mesh_a = d
        #	break

        z1=el[0].getElemEdges()[i].getNodes()[0].coordinates[2]
        z2=el[0].getElemEdges()[i].getNodes()[1].coordinates[2]
        y1=el[0].getElemEdges()[i].getNodes()[0].coordinates[1]
        y2=el[0].getElemEdges()[i].getNodes()[1].coordinates[1]
        
        if z1==z2:
            wb_mesh_a= abs(y1-y2)
            break

    nef=round(b/fl_mesh_a) #no of elements in the width of flange
    ned=round((h-tf)/wb_mesh_a) #no of elements in the depth of web

    import regionToolset

    column_model.parts['Beam'].MaterialOrientation( additionalRotationField='', additionalRotationType=ROTATION_ANGLE, angle=90.0, axis=AXIS_1, fieldName='', localCsys=None, orientationType=SYSTEM, region=regionToolset.Region(faces=column_model.parts['Beam'].faces))

    for i in range(int(ned)):
        y=-1*((h-tf)/2-wb_mesh_a/2)+i*wb_mesh_a
        el = column_model.parts['Beam'].elements.getByBoundingCylinder(center1=(0,y,0),center2=(0,y,L),radius =wb_mesh_a/2+1)
        column_model.parts['Beam'].Set(elements= el, name='W_Set_y='+str(int(y)))
        column_model.Stress( distributionType=UNIFORM, name='Predefined Field-'+'W_Set_y='+str(int(y)), region=column_model.rootAssembly.instances['Profile Instance'].sets['W_Set_y='+str(int(y))], sigma11=-0.5*fy+abs(y)*fy/((h-tf)/2), sigma12=0.0, sigma13=None, sigma22=0.0, sigma23=None, sigma33=None)

    for i in range(int(nef)):
        x=-1*(b/2-fl_mesh_a/2)+i*fl_mesh_a
        el1 = column_model.parts['Beam'].elements.getByBoundingCylinder(center1=(x,(h-tf)/2,0),center2=(x,(h-tf)/2,L),radius =fl_mesh_a/2+1)
        el2 = column_model.parts['Beam'].elements.getByBoundingCylinder(center1=(x,-1*(h-tf)/2,0),center2=(x,-1*(h-tf)/2,L),radius =fl_mesh_a/2+1)
        column_model.parts['Beam'].Set(elements= el1+el2, name='FL_Set_x='+str(int(x)))
        column_model.Stress( distributionType=UNIFORM, name='Predefined Field-'+'FL_Set_x='+str(int(x)), region=column_model.rootAssembly.instances['Profile Instance'].sets['FL_Set_x='+str(int(x))], sigma11=0.5*fy-abs(x)*fy/(b/2), sigma12=0.0, sigma13=None, sigma22=0.0, sigma23=None, sigma33=None)
	

if bucklestep==1:
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model=name+'_buckle', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Job_'+name+'_buckle', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

    column_model.keywordBlock.synchVersions(storeNodesAndElements=False)
    column_model.keywordBlock.setValues()
    line_num = 0
    for n, line in enumerate(column_model.keywordBlock.sieBlocks):
        if line == '*Restart, write, frequency=0':
            line_num = n
            break
    if line_num:
        column_model.keywordBlock.insert(position=line_num,text='*NODE FILE \nU')
    else:
        e = ("Error: Part '{}' was not found".format(partname),
            "in the Model KeywordBlock.")
        raise Exception(" ".join(e))

    
elif perfectstep==1:
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model=name, modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Job_'+name, nodalOutputPrecision=SINGLE, 
    numCpus=nocores, numDomains=nocores, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

elif geoimpstep==1:
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model=name+'_geoimp', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Job_'+name+'_geoimp', nodalOutputPrecision=SINGLE, 
    numCpus=nocores, numDomains=nocores, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

    column_model.keywordBlock.synchVersions(storeNodesAndElements=False)
    column_model.keywordBlock.setValues()
    line_num = 0
    for n, line in enumerate(column_model.keywordBlock.sieBlocks):
        if line == '** ----------------------------------------------------------------\n** \n** STEP: Step-1\n** ':
            line_num = n
            break
    if line_num:
        column_model.keywordBlock.insert(position=line_num,text='*IMPERFECTION,FILE=Job_'+name+'_buckle,STEP=1 \n1,'+str(-L/500)+'\n2,'+str(-L/500)+'')
    else:
        e = ("Error: Part '{}' was not found".format(partname),
            "in the Model KeywordBlock.")
        raise Exception(" ".join(e))