import sketch
import part
from abaqus import *
from abaqusConstants import *
import regionToolset
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
nocores=6
#
lrd=20.0 #longitudinal rebar diameter
nr=12 #number of lrebars
std=8.0 #stirrup diameter
sts=100.0 #stirrup spacing
fs=500.0*1.1 #rebaryield
fy=550.0*1.1 #profileyield
fcm=100.0*1.3 #concrete compressive strength
#

name='_'.join(['rect','d',str(int(d)),'fcm',str(int(fcm)),'h',str(int(h)),'fy',str(int(fy)),'sts',str(int(sts)),axis])
column_model=mdb.Model(name=name, modelType=STANDARD_EXPLICIT)

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
concrete_sketch.rectangle(point1=(d/2, d/2), point2=(-d/2, -d/2))

concrete_part=column_model.Part(name='concrete',dimensionality=THREE_D,type=DEFORMABLE_BODY)
concrete_part.BaseSolidExtrude(sketch=concrete_sketch,depth=L)
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(b/2,h/2,6),point2=(b/2,-h/2,5),point3=(b/2,-h/2,4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,h/2,6),point2=(-b/2,-h/2,5),point3=(-b/2,-h/2,4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,h/2,6),point2=(b/2,h/2,5),point3=(-b/2,h/2,4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,-h/2,6),point2=(b/2,-h/2,5),point3=(-b/2,-h/2,4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,-(h/2-tf),6),point2=(b/2,-(h/2-tf),5),point3=(-b/2,-(h/2-tf),4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(-b/2,(h/2-tf),6),point2=(b/2,(h/2-tf),5),point3=(-b/2,(h/2-tf),4))
concrete_part.PartitionCellByPlaneThreePoints(cells=concrete_part.cells,point1=(0,(h/2-tf),6),point2=(0,-(h/2-tf),5),point3=(0,(h/2-tf),4))
