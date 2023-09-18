import sys
import numpy as np
from vtkmodules.all import (
    vtkCellArray,
    vtkPoints,
    vtkPolyData)
from writer import writeDisplacementsCSV,writePolyDataAsSTL

#-------------------------------------#
#some operations on vectors
def getFaceNormal(segm: list) -> list:
    """
    Returns the normal to the surface using three points [0, 1, 2].
    """
    dir1 = np.diff([segm[1], segm[0]], axis=0)[0]
    dir2 = np.diff([segm[1], segm[2]], axis=0)[0]
    segm_normal = np.cross(dir1, dir2)
    return segm_normal/np.linalg.norm(segm_normal)

def getDistance(pt1: list, pt2: list) -> float:
    vec = np.diff([pt2, pt1], axis=0)[0]
    return np.sqrt(np.dot(vec, vec))

def normalize(pt: list) -> list:
    return np.divide(pt, np.linalg.norm(pt))

#-------------------------------------#
#preparation the data

def createCenterline(segms: list) -> list:
    centerline = []
    for pts in segms:
        newPt = [0., 0., 0.]
        for pt in pts:
            newPt = np.add.reduce([newPt, pt])
        newPt = np.divide(newPt, len(pts))
        centerline.append(newPt)  
    return centerline  

def calcOffset(bdsSegments: list, data_list: list, functionName: str) -> list:
    offset = []
    N = bdsSegments[3] - bdsSegments[2]
    for i in range(0, N):
        offset.append(funcOffset(i, bdsSegments, data_list, functionName))    
    return offset
   
def funcOffset(k, bdsSegments: list, data_list: list, functionName: str) -> float:
    match functionName:
        case "sin":
            # f(k) = (delta_min_diam/2)*sin(pi*k/N_s)
            return (np.abs(data_list[5] - data_list[6]) *
                    (np.sin(np.pi * k / (bdsSegments[3] - bdsSegments[2]))))
        case "s-curve":
            print("Warring: \'s-curve\' has not been implemented yet!")
            return 0
        case _:
            print("Warring: \'" + functionName + "\' not found, check \'init.ini\'!")
            sys.exit()

def createDisplacementWall(lumenStl, segms3DLumen, bdsSegments, centerline, offset):
    cellArray= vtkCellArray()
    cellArray.DeepCopy(lumenStl.GetCells())
    pts = vtkPoints()
    pts.DeepCopy(lumenStl.GetPoints())

    listPts = [] #idx pts for preparation

    segm_init = segms3DLumen[bdsSegments[2]]
    segm_init_normal = getFaceNormal(segm_init)
    segm_fin = segms3DLumen[bdsSegments[3]]
    segm_fin_normal = getFaceNormal(segm_fin)
        
    for i in range(0, pts.GetNumberOfPoints()):
        flg = True
        pt = pts.GetPoint(i)

        p2f_init = np.diff([pt, segm_init[0]], axis=0)[0]
        d_init = np.dot(p2f_init, segm_init_normal)
        flg *= d_init < 0
                
        p2f_fin = np.diff([pt, segm_fin[0]], axis=0)[0]
        d_fin = np.dot(p2f_fin, segm_fin_normal)
        flg *= d_fin > 0

        if flg:
            listPts.append(i)

    offset2vec = []
    for i in listPts:
        pt = pts.GetPoint(i)
        listdist = []
        for j in range(bdsSegments[2], bdsSegments[3]):
            listdist.append(getDistance(pt, centerline[j]))
        minIdx = np.argmin(listdist)
        
        norm_vec = normalize(np.diff([pt, centerline[int(minIdx + bdsSegments[2])]], axis=0)[0])
        offset2vec.append(np.multiply(norm_vec, offset[minIdx]))

    for j in range(0, 10):
        tempPts = vtkPoints()
        tempPts.DeepCopy(pts)
        for i in range(0, len(listPts)):
            pt = tempPts.GetPoint(listPts[i])
            tempPts.SetPoint(listPts[i], np.add.reduce([pt, np.multiply(offset2vec[i], j/10)]))
        tempPts.Modified()
        newPD = vtkPolyData()
        newPD.SetPoints(tempPts)
        newPD.SetPolys(cellArray)

        writePolyDataAsSTL('data/output/offset_' + str(j) + '.stl', newPD)  #TODO: весь вывод должен быть в writer.py
    
    writeDisplacementsCSV('data/output/disp.csv', pts, [[0.0]*3] * pts.GetNumberOfPoints())
