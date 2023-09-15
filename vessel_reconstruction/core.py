import sys
import numpy as np
from vtkmodules.all import (
    vtkCellArray,
    vtkPoints,
    vtkPolyData)
from writer import writePolyDataAsSTL

#-------------------------------------#
#some operations on vectors
def sub(pt1: list, pt2: list) -> list:
    """
    Vector difference.
    """
    return [a - b for a, b in zip(pt1, pt2)]

def add(pt1: list, pt2: list) -> list:
    """
    Sum of vectors.
    """
    return [a + b for a, b in zip(pt1, pt2)]

def truediv(pt: list, num: int) -> list:
    """
    Dividing a vector by a number.
    """
    return [a / num for a in pt]

def multy(pt: list, num: float) -> list:
    return [a * num for a in pt]

def getNormal(segm: list) -> list:
    """
    Returns the normal to the surface using three points [0, 1, 2].
    """
    dir1 = sub(segm[0], segm[1])
    dir2 = sub(segm[2], segm[1])
    segm_normal = np.cross(dir1, dir2)
    return segm_normal/np.linalg.norm(segm_normal)

def getDistance(pt1: list, pt2: list) -> float:
    vec = sub(pt1, pt2)
    return np.sqrt(np.dot(vec, vec))

def normalize(pt: list) -> list:
    return truediv(pt, np.linalg.norm(pt))

#-------------------------------------#
#preparation the data

def createCenterline(segms: list) -> list:
    centerline = []
    for pts in segms:
        newPt = [0., 0., 0.]
        for pt in pts:
            newPt = add(newPt, pt)
        newPt = truediv(newPt, len(pts))
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


def tempSTL(lumenStl, segms3DLumen, bdsSegments, centerline, offset):
    cellArray= vtkCellArray()
    cellArray.DeepCopy(lumenStl.GetPolys())
    pts = vtkPoints()
    pts.DeepCopy(lumenStl.GetPoints())

    listPts = [] #idx pts por preparation

    segm_init = segms3DLumen[bdsSegments[2]]
    segm_init_normal = getNormal(segm_init)
    segm_fin = segms3DLumen[bdsSegments[3]]
    segm_fin_normal = getNormal(segm_fin)
        
    for i in range(0, pts.GetNumberOfPoints()):
        flg = True
        pt = pts.GetPoint(i)

        p2f_init = sub(segm_init[0], pt)
        d_init = np.dot(p2f_init, segm_init_normal)
        flg *= d_init < 0
                
        p2f_fin = sub(segm_fin[0], pt)
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

        offset2vec.append(multy(normalize(sub(centerline[int(minIdx + bdsSegments[2])], pt)), offset[minIdx]))

    for j in range(1, 11):
        tempPts = vtkPoints()
        tempPts.DeepCopy(pts)
        for i in range(0, len(listPts)):
            pt = tempPts.GetPoint(listPts[i])
            tempPts.SetPoint(listPts[i], add(pt, truediv(offset2vec[i], j)))
        tempPts.Modified()
        newPD = vtkPolyData()
        newPD.SetPoints(tempPts)
        newPD.SetPolys(cellArray)

        writePolyDataAsSTL('data/result/offset_' + str(j) + '.stl', newPD)