import pygalmesh
import vtk
from cmath import sqrt, acos
from collections import OrderedDict
from vmtk import pypes
from vtkmodules.all import (
    vtkAppendPolyData,
    vtkCellArray,
    vtkCleanPolyData,
    vtkIdList,
    vtkPoints,
    vtkPolyData,
    vtkSTLWriter,
    vtkXMLPolyDataReader,
    vtkXMLPolyDataWriter,
    vtkXMLUnstructuredGridReader)
import numpy as np
from core import getDistance

def getNodeIdOnTheBorder(polydata):
    """
    Returns a list of indexes of points lying on the instream and outstream curves.\n
    - polydata - vtkPolyData object, vessel lumen.
    """
    node2cell = OrderedDict()
    for i in range(polydata.GetNumberOfCells()):
        for j in range(3):
            if (polydata.GetCell(i).GetPointId(j) in node2cell):
                node2cell[polydata.GetCell(i).GetPointId(j)] = node2cell[polydata.GetCell(i).GetPointId(j)] + 1
            else:
                node2cell[polydata.GetCell(i).GetPointId(j)] = 1
    nodeId = list()

    # check which points are on the border
    for key, value in node2cell.items():
        if (value == 3):
            nodeId.append(key)
    return nodeId

def sort(polydata, ptsId):
    """
    Sort nodes on border to create cells.\n
    Returns an ordered list of indexes.\n
    - polydata - vtkPolyData object, vessel lumen,
    - ptsId - idxs points on border.
    """
    listLen = list()
    listLen.append(0.0)
    pt0 = polydata.GetPoint(ptsId[0])
    for i in range(1, len(ptsId)):
        pt1 = polydata.GetPoint(ptsId[i])
        listLen.append(getDistance(pt0, pt1))

    newlistId = ptsId
    sort_flg = False
    while (not sort_flg):
        sort_flg = True
        for i in range(len(listLen) - 1):
            if listLen[i] > listLen[i + 1]:
                listLen[i], listLen[i + 1] = listLen[i + 1], listLen[i]
                newlistId[i], newlistId[i + 1] = newlistId[i + 1], newlistId[i]
                sort_flg = False

    listAngle = list()
    listAngle.append(0.0)
    listAngle.append(0.0)
    pt0 = polydata.GetPoint(newlistId[0])
    pt1 = polydata.GetPoint(newlistId[1])
    v01 = [pt1[0] - pt0[0], pt1[1] - pt0[1], pt1[2] - pt0[2]]
    for i in range(2, len(ptsId)):
        pt2 = polydata.GetPoint(newlistId[i])
        v02 = [pt2[0] - pt0[0], pt2[1] - pt0[1], pt2[2] - pt0[2]]
        listAngle.append(abs(acos(np.dot(v02, v01) / sqrt(np.dot(v01, v01) * np.dot(v02, v02)))))

    sort_flg = False
    while (not sort_flg):
        sort_flg = True
        for i in range(1, len(listAngle) - 1):
            if listAngle[i] > listAngle[i + 1]:
                listAngle[i], listAngle[i + 1] = listAngle[i + 1], listAngle[i]
                newlistId[i], newlistId[i + 1] = newlistId[i + 1], newlistId[i]
                sort_flg = False
    return newlistId

def closeSurface(polyLumen, polyWall):
    """
    Returns the volumetric closed surface of the vessel.\n
    (Since the size of lumen and wall do not match, the outer surface 
    is created according to the average thickness of the vessel)
    - polyLumen - vtkPolyData, lumen vessel
    - polyWall - vtkPolyData, wall vessel 
    """
    appendFilter = vtkAppendPolyData()
    appendFilter.AddInputData(polyLumen)
    appendFilter.AddInputData(polyWall)

    cleanFilter = vtkCleanPolyData()
    cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
    cleanFilter.Update()

    polydata = cleanFilter.GetOutput()
    ids = getNodeIdOnTheBorder(polydata)
    idsL1 = ids[0:int(len(ids) / 4)]
    idsL2 = ids[int(len(ids) / 4):int(len(ids) / 2)]
    idsW1 = sort(polydata, ids[int(len(ids) / 2):int(len(ids) * 3 / 4)])
    idsW2 = sort(polydata, ids[int(len(ids) * 3 / 4):int(len(ids))])

    for i in range(len(idsW1)):
        dists1 = list()
        dists2 = list()
        for j in range(len(idsL1)):
            dists1.append(getDistance(polydata.GetPoint(idsW1[i]), polydata.GetPoint(idsL1[j])))
            dists2.append(getDistance(polydata.GetPoint(idsW2[i]), polydata.GetPoint(idsL2[j])))
        idx1 = 0
        idx2 = 0
        for j in range(1, len(dists1)):
            if dists1[idx1] > dists1[j]:
                idx1 = j
            if dists2[idx2] > dists2[j]:
                idx2 = j
        idsL1[i], idsL1[idx1] = idsL1[idx1], idsL1[i]
        idsL2[i], idsL2[idx2] = idsL2[idx2], idsL2[i]

    idList = vtkIdList()
    idList.SetNumberOfIds(3)
    for i in range(len(idsL1) - 1):
        idList.SetId(0, idsL1[i])
        idList.SetId(1, idsL1[i + 1])
        idList.SetId(2, idsW1[i])
        polydata.InsertNextCell(vtk.VTK_TRIANGLE, idList)

        idList.SetId(0, idsW1[i])
        idList.SetId(1, idsW1[i + 1])
        idList.SetId(2, idsL1[i + 1])
        polydata.InsertNextCell(vtk.VTK_TRIANGLE, idList)

        idList.SetId(0, idsL2[i])
        idList.SetId(1, idsL2[i + 1])
        idList.SetId(2, idsW2[i])
        polydata.InsertNextCell(vtk.VTK_TRIANGLE, idList)

        idList.SetId(0, idsW2[i])
        idList.SetId(1, idsW2[i + 1])
        idList.SetId(2, idsL2[i + 1])
        polydata.InsertNextCell(vtk.VTK_TRIANGLE, idList)

    idList.SetId(0, idsL1[len(idsL1) - 1])
    idList.SetId(1, idsL1[0])
    idList.SetId(2, idsW1[len(idsL1) - 1])
    polydata.InsertNextCell(vtk.VTK_TRIANGLE, idList)

    idList.SetId(0, idsW1[len(idsL1) - 1])
    idList.SetId(1, idsW1[0])
    idList.SetId(2, idsL1[0])
    polydata.InsertNextCell(vtk.VTK_TRIANGLE, idList)

    idList.SetId(0, idsL2[len(idsL2) - 1])
    idList.SetId(1, idsL2[0])
    idList.SetId(2, idsW2[len(idsW2) - 1])
    polydata.InsertNextCell(vtk.VTK_TRIANGLE, idList)

    idList.SetId(0, idsW2[len(idsW2) - 1])
    idList.SetId(1, idsW2[0])
    idList.SetId(2, idsL2[0])
    polydata.InsertNextCell(vtk.VTK_TRIANGLE, idList)

    return polydata

def calcCenterLine(polydata, path):
    """
    Writing polydata in *.vtp for VMTK scripts.\n
    - polydata - vtkPolyData, lumen vessel,
    - path - output path
    """
    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(path + 'lumen.vtp')
    writer.SetInputData(polydata)
    writer.Write()

    myArguments = 'vmtkcenterlines -seedselector profileidlist -sourceids 0 -ifile ' + path + 'lumen.vtp -ofile ' + path + 'centerline.vtp'
    pypes.PypeRun(myArguments)

def readCenterLineVTP(path):
    """
    Read 'centerline.vtp' from path.
    """
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(path + 'centerline.vtp')
    reader.Update()
    return reader.GetOutput()

def generateWall(polyLumen, polyWall, centerLine):
    """
    Returns the outer surface of the vessel, plotted by the average radius.\n
    - polyLumen - vtkPolyData, lumen vessel,
    - polyWall - vtkPolyData, wall vessel,
    - centerLine - vtkPolyData, centerline.
    """
    tempWall = vtkPolyData()
    tempWall.DeepCopy(polyLumen)
    scale = getAverageRadius(polyWall, centerLine) / getAverageRadius(polyLumen, centerLine)
    nPtsW = tempWall.GetNumberOfPoints()
    nPtsCL = centerLine.GetNumberOfPoints()

    for i in range(nPtsW):
        pt0 = tempWall.GetPoint(i)
        pt1 = centerLine.GetPoint(0)
        rad = getDistance(pt0, pt1)
        for j in range(1, nPtsCL):
            newRad = getDistance(pt0, centerLine.GetPoint(j))
            if (newRad < rad):
                rad = newRad
                pt1 = centerLine.GetPoint(j)
        pt3 = list()
        for k in range(3):
            pt3.append(scale * pt0[k] + pt1[k] * (1 - scale))
        tempWall.GetPoints().SetPoint(i, pt3)

    idNodeWB = getNodeIdOnTheBorder(tempWall)

    pts = vtkPoints()
    cells = vtkCellArray()
    newWall = vtkPolyData()

    for i in range(tempWall.GetNumberOfPoints()):
        pts.InsertNextPoint(tempWall.GetPoint(i))
    for i in range(tempWall.GetNumberOfCells()):
        idlist = vtkIdList()
        tempWall.GetCellPoints(i, idlist)
        if (len(set(idNodeWB).intersection([idlist.GetId(0), idlist.GetId(1), idlist.GetId(2)])) == 0):
            cells.InsertNextCell(idlist)
    newWall.SetPoints(pts)
    newWall.SetPolys(cells)

    clean = vtkCleanPolyData()
    clean.SetInputData(newWall)
    clean.Update()
    return clean.GetOutput()

def getAverageRadius(polydata, centerLine):
    """
    Returns the average distance between the surface and the line.\n
    - polydata - vtkPolyData, surface
    - centerLine - vtkPolyData, line
    """
    nPtsPD = polydata.GetNumberOfPoints()
    nPtsCL = centerLine.GetNumberOfPoints()
    radiuses = list()

    for i in range(nPtsPD):
        pt0 = polydata.GetPoint(i)
        dists = list()
        for j in range(nPtsCL):
            pt1 = centerLine.GetPoint(j)
            dists.append(getDistance(pt0, pt1))
        radiuses.append(min(dists))
        dists.clear()
    return sum(radiuses) / len(radiuses)

def vtu2csv(filename, outpath):
    """
    Convert *.vtu to *.csv file.
    """
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(outpath + filename)
    reader.Update()
    ugrid = reader.GetOutput()

    nodesfile = open(outpath + 'nodesMesh.csv', 'w')
    for i in range(ugrid.GetNumberOfPoints()):
        pt = ugrid.GetPoint(i)
        nodesfile.write(f'{i:d}' + '\t' + f'{pt[0]:.6f}' + '\t' + f'{pt[1]:.6f}' + '\t' + f'{pt[2]:.6f}' + '\n')
    nodesfile.close()

    elemsfile = open(outpath + 'elemsMesh.csv', 'w')
    idx = 0
    for i in range(ugrid.GetNumberOfCells()):
        cell = ugrid.GetCell(i)
        if (cell.GetCellType() == vtk.VTK_TETRA):
            elemsfile.write(f'{idx:d}' + '\t' + f'{cell.GetPointId(0):d}' + '\t' + f'{cell.GetPointId(1):d}' + '\t'f'{cell.GetPointId(2):d}' + '\t'f'{cell.GetPointId(3):d}' + '\n')
            idx += 1
    elemsfile.close()

def runMesher(lumenStl, wallStl, outpath):
    """
    Generate a volume mesh.
    """
    print("\nMesher: Center line VMTK calculation...")
    calcCenterLine(lumenStl, outpath)
    centerLine = readCenterLineVTP(outpath)

    print("Mesher: Vessel generation...")
    newWall = generateWall(lumenStl, wallStl, centerLine)
    vesselPolyData = closeSurface(lumenStl, newWall)

    writerVTP = vtkSTLWriter()
    writerVTP.SetFileName(outpath + 'vessel.stl')
    writerVTP.SetInputData(vesselPolyData)
    writerVTP.Write()

    # Triangulate:
    print("Mesher: Volume mesh generation....")
    mesh = pygalmesh.generate_volume_mesh_from_surface_mesh(
        outpath + 'vessel.stl',
        lloyd=True,
        min_facet_angle=25.0,
        max_radius_surface_delaunay_ball=1.5,
        max_facet_distance=0.03,
        max_circumradius_edge_ratio=3.0,
        reorient=True
    )

    print("Mesher: Writing results to a file.")
    mesh.write(outpath + "volumeMesh.vtu", 'vtu')
    vtu2csv("volumeMesh.vtu", outpath)

    print("Mesher:done!")
