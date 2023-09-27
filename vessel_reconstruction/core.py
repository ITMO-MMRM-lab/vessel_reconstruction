import sys
import numpy as np
from vtkmodules.all import (
    vtkCellArray,
    vtkPoints,
    vtkPolyData,
    vtkPlane,
    vtkCutter)
from scipy.interpolate import lagrange

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

class DataAlgorithms(object):
    def __init__(self, volumeMesh, lumenStl, segms3DLumen, bdsSegments, data_list, funcName="sin", steps=10):
        self.cline = self.createCenterline(segms3DLumen)
        self.funcName = funcName
        self.offsets = self.calcOffsets(bdsSegments, data_list, self.funcName)
        [self.displs, self.trajs] = self.createDisplacementWall(volumeMesh, segms3DLumen, bdsSegments, self.cline, self.offsets, steps)
        self.prelumen = self.createPreLumen(lumenStl, segms3DLumen, bdsSegments, self.cline, self.offsets)
        self.cur_diams = self.getDiameters(self.prelumen, self.cline, bdsSegments, data_list)
    
    def createCenterline(self, segms: list) -> list:
        centerline = []
        for pts in segms:
            newPt = [0., 0., 0.]
            for pt in pts:
                newPt = np.add.reduce([newPt, pt])
            newPt = np.divide(newPt, len(pts))
            centerline.append(newPt)  
        return centerline  

    def calcOffsets(self, bdsSegments: list, data_list: list, functionName: str) -> list:
        """
        The function calculates scalar offsets for each segment using a given function.
        """
        offsets = []
        N = bdsSegments[3] - bdsSegments[2]
        for i in range(0, N):
            offsets.append(self.funcOffset(i, bdsSegments, data_list, functionName))    
        return offsets
    
    def funcOffset(self, k, bdsSegments: list, data_list: list, functionName: str) -> float:
        N_s = bdsSegments[3] - bdsSegments[2]
        deltaMinR = np.abs(data_list[5] - data_list[6])/2.
        match functionName:
            case "sin":
                # f(k) = R * sin(pi * k / N_s)
                return deltaMinR* (np.sin(np.pi * k / (N_s)))
            
            case "arctan":
                # f(k) = ± (R / 2) * (2 * arctan(k-n1) / pi) ± 1)
                n1 = N_s / 6.
                n2 = 5. * n1
                if (k <= N_s / 2.):
                    return  (deltaMinR / 2.) * ((np.arctan(k - n1) / (np.pi / 2.)) + 1)
                else:
                    return -(deltaMinR / 2.) * ((np.arctan(k - n2) / (np.pi / 2.)) - 1)

            case "polynomial":
                if not hasattr(self, 'poly_coef'):
                    x = [0., N_s/2., N_s]
                    y = [0., deltaMinR, 0.]
                    self.poly_coef = lagrange(x, y).coef
                return  self.poly_coef[0]*k*k + self.poly_coef[1]*k + self.poly_coef[2]
            
            case _:
                print("Warring: \'" + functionName + "\' not found, check \'init.ini\'!")
                sys.exit()

    def createDisplacementWall(self, volumeMesh, segms3DLumen, bdsSegments, cline, offsets, steps=10):
        cellArray= vtkCellArray()
        cellArray.DeepCopy(volumeMesh.GetCells())
        pts = vtkPoints()
        pts.DeepCopy(volumeMesh.GetPoints())

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
                listdist.append(getDistance(pt, cline[j]))
            minIdx = np.argmin(listdist)
            
            norm_vec = normalize(np.diff([pt, cline[int(minIdx + bdsSegments[2])]], axis=0)[0])
            offset2vec.append(np.multiply(norm_vec, offsets[minIdx]))

        visTraj = []
        for j in range(0, steps):
            tempPts = vtkPoints()
            tempPts.DeepCopy(pts)
            for i in range(0, len(listPts)):
                pt = tempPts.GetPoint(listPts[i])
                tempPts.SetPoint(listPts[i], np.add.reduce([pt, np.multiply(offset2vec[i], j/10)]))
            newPD = vtkPolyData()
            newPD.SetPoints(tempPts)
            newPD.SetPolys(cellArray)
            visTraj.append(newPD)

        allOffsets = [[0.0]*3] * pts.GetNumberOfPoints()
        for i in range(0, len(listPts)):
            allOffsets[listPts[i]] = offset2vec[i]

        return [allOffsets, visTraj]
    
    def createPreLumen(self, lumenStl, segms3DLumen, bdsSegments, cline, offsets):
        cellArray= vtkCellArray()
        cellArray.DeepCopy(lumenStl.GetPolys())
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
                listdist.append(getDistance(pt, cline[j]))
            minIdx = np.argmin(listdist)
            
            norm_vec = normalize(np.diff([pt, cline[int(minIdx + bdsSegments[2])]], axis=0)[0])
            offset2vec.append(np.multiply(norm_vec, offsets[minIdx]))

        newPts = vtkPoints()
        newPts.DeepCopy(pts)
        for i in range(0, len(listPts)):
            pt = newPts.GetPoint(listPts[i])
            newPts.SetPoint(listPts[i], np.add.reduce([pt, offset2vec[i]]))

        newLumen = vtkPolyData()
        newLumen.SetPoints(newPts)
        newLumen.SetPolys(cellArray)
        return newLumen

    def getDiameters(self, prelumen, cline, bdsSegms, data_list):
        '''
        The function returns the diameters of the vessel.
        - prelumen: lumen surface, genetared polydata 
        - cline: centerLine from createCenterline(...), ordered array of points
        - bdsSegms: segment numbers [initVessel, finVessel, initStent, finStent]
        '''   
        diams = []
        for i in range(0, len(cline) - 1):
            plane = vtkPlane()
            plane.SetOrigin(cline[i])
            plane.SetNormal(normalize(np.diff([cline[i+1], cline[i]], axis=0)[0]))

            cutter = vtkCutter()
            cutter.SetInputData(prelumen)
            cutter.SetCutFunction(plane)
            cutter.Update()

            pts = cutter.GetOutput().GetPoints()
            npts = pts.GetNumberOfPoints()

            sumRad = 0
            for j in range(0, npts):
                sumRad += getDistance(cline[i], pts.GetPoint(j))
            diams.append(sumRad*2./npts)
        return diams