#fy=100.0
#el = column_model.parts['Beam'].elements.getByBoundingCylinder(center1=(b/2-15,(h-tf)/2,0),center2=(b/2-15,(h-tf)/2,40),radius =20)
el = column_model.parts['Beam'].elements.getByBoundingBox(xMin=0, xMax=b/2, yMin=(h-tf)/2, yMax=(h-tf)/2, zMin=0, zMax=40)

fl_mesh_a=30
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

wb_mesh_a=30
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
	