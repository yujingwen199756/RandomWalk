import os
import vtk
import math
import glob
import random
import numpy as np

# Diffusion Coefficient
diffusionCoefficient = 0.05
# Time internval
time = 3
# Size of the data
timePoint = 2000
# Working Directory
workdir = '/Users/yujingwen199756/Desktop/ER-5-300-DEF'
# Could be this for Tatsu's computer: '/Users/Tatsuhisa/Desktop/080516_ATP3_2_R3_33_G4_33_300nm_x100_5_Pos0'

#-------------------------------------------------------------------------
#:: AUXILIAR FUNCTION FOR CALCULATING THE DISTANCES
#-------------------------------------------------------------------------


def CalculateDistance(mito_surface_name,cell_inner_surface_name, cell_outer_surface_name):

    #Open the surface of Mito
    SurfaceMito = LegacyVTKReader(FileNames=[mito_surface_name])
    SurfaceMito = GetActiveSource()
    SurfaceMito = servermanager.Fetch(SurfaceMito) 
  
    #Open the outer surface of Cell
    OuterSurfaceCell = LegacyVTKReader(FileNames=[cell_outer_surface_name])
    OuterSurfaceCell = GetActiveSource()
    OuterSurfaceCell = servermanager.Fetch(OuterSurfaceCell) 
    geometryFilterOuterCell = vtk.vtkGeometryFilter()
    geometryFilterOuterCell.SetInputData(OuterSurfaceCell)
    geometryFilterOuterCell.Update()
    polydataOuterCell = geometryFilterOuterCell.GetOutput()
    #print polydataOuterCell.GetNumberOfPoints()

    #Open the inner surface of Cell
    InnerSurfaceCell = LegacyVTKReader(FileNames=[cell_inner_surface_name])
    InnerSurfaceCell = GetActiveSource()
    InnerSurfaceCell = servermanager.Fetch(InnerSurfaceCell) 
    geometryFilterInnerCell = vtk.vtkGeometryFilter()
    geometryFilterInnerCell.SetInputData(InnerSurfaceCell)
    geometryFilterInnerCell.Update()
    polydataInnerCell = geometryFilterInnerCell.GetOutput()
    
    
    #Get the bounds of the cell(xmin,xmax,ymin,ymax,zmin,zmax)
    bounds = [0]*6
    OuterSurfaceCell.GetBounds(bounds)
   
    #Creating the point locator 
    LocatorMito = vtk.vtkPointLocator()
    LocatorMito.SetDataSet(SurfaceMito)
    LocatorMito.BuildLocator()
   
    #Vector to store the distance from foci to mito
    DistanceToMito = []
    DistanceMoved = []

 
    selectEnclosedPointsOuterCell = vtk.vtkSelectEnclosedPoints()
    selectEnclosedPointsOuterCell.Initialize(polydataOuterCell)

    selectEnclosedPointsInnerCell = vtk.vtkSelectEnclosedPoints()
    selectEnclosedPointsInnerCell.Initialize(polydataInnerCell)

    selectEnclosedPointsMito = vtk.vtkSelectEnclosedPoints()
    selectEnclosedPointsMito.Initialize(SurfaceMito)    

    insideOuterCell = 0
    insideInnerCell = 1 
    
    renderWindow = vtk.vtkRenderWindow()
    renderer = vtk.vtkRenderer()
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    mapper = vtk.vtkPolyDataMapper()
    mapperMito = vtk.vtkPolyDataMapper()
    mapperMito.SetInputData(SurfaceMito)
    actorMito = vtk.vtkActor()
    actorMito.SetMapper(mapperMito)
    mapper.SetInputData(polydataOuterCell)#need to change
    mapperInner = vtk.vtkPolyDataMapper()
    mapperInner.SetInputData(polydataInnerCell)
    actor =vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(0.5)
    actorinner =vtk.vtkActor()
    actorinner.SetMapper(mapperInner)
    actorinner.GetProperty().SetOpacity(0.5)
    renderWindow.SetSize(3000,2000)
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderer.AddActor(actorMito)
    renderer.AddActor(actor)
    renderer.AddActor(actorinner)
    renderWindow.Render()
    
    while insideOuterCell == 0 and insideInnerCell == 1:
        x = random.uniform(bounds[0], bounds[1])
        y = random.uniform(bounds[2], bounds[3])
        z = random.uniform(bounds[4], bounds[5])
       
        insideOuterCell = selectEnclosedPointsOuterCell.IsInsideSurface(x, y, z)
        insideInnerCell = selectEnclosedPointsInnerCell.IsInsideSurface(x, y, z)
        #Check to see if the random foci is inside the cell
        if insideOuterCell ==1 and insideInnerCell==0:
            insideMito = selectEnclosedPointsMito.IsInsideSurface(x,y,z)
            #Check to see if the random foci is inside the mitochroniral
            if insideMito:
                insideOuterCell = 0
                insideInnerCell = 1
                continue
            else:  
                #Calculate the distance between the foci (x,y,z) and the Surface                
                
                foci = vtk.vtkSphereSource()
                foci.SetCenter(x,y,z)
                foci.SetRadius(0.05)
                mapper1 = vtk.vtkPolyDataMapper()
                mapper1.SetInputConnection(foci.GetOutputPort())
                actor1 =vtk.vtkActor()
                actor1.GetProperty().SetColor(1,0,0)
                actor1.SetMapper(mapper1)  
                renderWindow.AddRenderer(renderer) 
                renderWindowInteractor.SetRenderWindow(renderWindow)
                renderer.AddActor(actor1)
                renderWindow.Render()  
               
                r = [x, y, z]
                ptId = LocatorMito.FindClosestPoint(r)
                u = SurfaceMito.GetPoints().GetPoint(ptId)
                distance = math.sqrt((r[0]-u[0])**2+(r[1]-u[1])**2+(r[2]-u[2])**2)
                DistanceToMito.append(distance)
                DistanceMoved.append(0)
        else:
            insideOuterCell = 0
            insideInnerCell = 1
            continue
    
    
    for randomNumber in range(timePoint):   
        # mean of x,y,z
        muX = x; 
        muY = y; 
        muZ = z; 
        # standard deviation
        sigma = math.sqrt(2 * diffusionCoefficient * time)  
        #Use
        xi = np.random.normal(muX, sigma, 1) 
        yi = np.random.normal(muY, sigma, 1) 
        zi = np.random.normal(muZ, sigma, 1) 

        radius = math.sqrt((xi-x)**2+(yi-y)**2+(zi-z)**2)
        insideOuterCell = selectEnclosedPointsOuterCell.IsInsideSurface(xi, yi, zi)
        insideInnerCell = selectEnclosedPointsInnerCell.IsInsideSurface(xi, yi, zi)
        if insideOuterCell == 1 and insideInnerCell == 0:
            insideMito = selectEnclosedPointsMito.IsInsideSurface(xi,yi,zi)
            #Check to see if the random foci is inside the mitochroniral
            if insideMito:
                continue
            else:  
                #Calculate the distance between the foci (x,y,z) and the Surface
                
                focii = vtk.vtkSphereSource()
                focii.SetCenter(xi,yi,zi)
                focii.SetRadius(0.05)                   
                mapper2 = vtk.vtkPolyDataMapper()
                mapper2.SetInputConnection(focii.GetOutputPort())
                actor2 =vtk.vtkActor()
                actor2.SetMapper(mapper2)
                renderer.AddActor(actor2)
                path = vtk.vtkLineSource()
                path.SetPoint1(x,y,z)
                path.SetPoint2(xi,yi,zi)
                mapper3 = vtk.vtkPolyDataMapper()
                mapper3.SetInputConnection(path.GetOutputPort())
                actor3 =vtk.vtkActor()
                actor3.SetMapper(mapper3)
                renderer.AddActor(actor3)
                renderWindow.Render()
                
                r = [xi, yi, zi]
                ptId = LocatorMito.FindClosestPoint(r)
                u = SurfaceMito.GetPoints().GetPoint(ptId)
                distance = math.sqrt((r[0]-u[0])**2+(r[1]-u[1])**2+(r[2]-u[2])**2)
                DistanceToMito.append(distance)
                DistanceMoved.append(radius)
                x = xi
                y = yi
                z = zi
        else:
            continue
    
    renderWindow.Render()
    renderWindowInteractor.Start()
    
    Delete(GetActiveSource())

    del SurfaceMito
    del OuterSurfaceCell
    del InnerSurfaceCell
    del LocatorMito

    return DistanceToMito,DistanceMoved

#-------------------------------------------------------------------------
#:: MAIN FUNCTION STARTS HERE
#-------------------------------------------------------------------------


#Main directory
#workdir = '/Users/yujingwen199756/Desktop/ER-5-300-DEF'
#workdir = '/Users/Tatsuhisa/Desktop/080516_ATP3_2_R3_33_G4_33_300nm_x100_5_Pos0'

 
#List the main directory content
for item in os.listdir(workdir):
 
    #If the item corresponds to a subfolder
    if os.path.isdir(os.path.join(workdir,item)):
 
        subdir = os.path.join(workdir,item)
 
        print subdir

        #File where the result are going to be written down
        fsave = open("%s/results.txt"%(subdir), 'w')

        fsave.write("Folder\t\tSurface\t\t\tFoci\tDistanceToMito\tDistanceMoved\n")
 
        #Vector to store surfaces name
        SurfaceNamesMito = glob.glob(os.path.join(workdir,item,'*00_surface.vtk'))
        SurfaceNamesInnerCell = glob.glob(os.path.join(workdir,item,'InnerCell*.vtk'))
        SurfaceNamesOuterCell = glob.glob(os.path.join(workdir,item,'OuterCell*.vtk'))


        #Read the cell file surface_OM*.vtk and the mito file surface_IM*.vtk
        for s in range(len(SurfaceNamesMito)):
            DistanceToMito, DistanceMoved = CalculateDistance(SurfaceNamesMito[s], SurfaceNamesInnerCell[s], SurfaceNamesOuterCell[s])
            
            for d in range(len(DistanceToMito)):
                fsave.write("%s\t%s\t%d\t%1.3f\t%1.3f\n" %(os.path.split(subdir)[-1:][0],os.path.basename(SurfaceNamesMito[s]),d,DistanceToMito[d],DistanceMoved[d]))
fsave.close()

print "Analysis complete."