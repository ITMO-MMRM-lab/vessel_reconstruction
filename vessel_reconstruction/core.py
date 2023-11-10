import sys

import numpy as np
import vtk
from progress.bar import FillingCirclesBar
from scipy.interpolate import lagrange
from vtkmodules.all import (
    vtkCellArray, 
    vtkCutter, 
    vtkPlane, 
    vtkPoints,
    vtkPointSet, 
    vtkPolyData, 
    vtkSmoothPolyDataFilter,
    vtkTransform, 
    vtkTransformFilter,
    vtkUnstructuredGrid)


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

    def tranformSegmentOLD(self, pointSet, p1, p2):
        # 1 - moving the segment to the origin
        center = pointSet.GetCenter()
        transform1 = vtkTransform()
        transform1.Translate(np.multiply(pointSet.GetCenter(), -1))
        transformFilter1  = vtkTransformFilter()
        transformFilter1.SetInputData(pointSet)
        transformFilter1.SetTransform(transform1)
        transformFilter1.Update()
        pointSet = transformFilter1.GetOutput()

        # 2 - moving the segment to the vessel
        midPt = np.divide(np.sum([p1, p2], axis=0), 2.)
        # midPt = p1

        v1 = np.diff([p1, p2], axis=0)[0]
        v2 = [0., 0., 1.]
        axe = normalize(np.cross(v1, v2))

        v1_u = v1 / np.linalg.norm(v1)
        v2_u = v2 / np.linalg.norm(v2)

        alpha = -np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

        transform = vtkTransform()
        transform.Translate(np.diff([pointSet.GetCenter(), midPt], axis=0)[0])
        # transform.Translate(center)
        transform.RotateWXYZ(alpha, axe)
        
        transformFilter = vtkTransformFilter()
        transformFilter.SetInputData(pointSet)
        transformFilter.SetTransform(transform)
        transformFilter.Update()
        return [transformFilter.GetOutput(), axe]
    

    def tranformSegment(self, segmsList, clinePts, rightId):
        axes = []
        betta = 0.0
        for i in range(0, len(segmsList)):
            [segmsList[i], axeI] = self.tranformSegmentOLD(segmsList[i], clinePts.GetPoint(rightId + i), clinePts.GetPoint(rightId + i + 1))
            axes.append(axeI)
        
        # for i in range(0, len(segmsList)):
        #     v1 = axes[i]
        #     v2 = [0., 1., 0.]
        #     v1_u = v1 / np.linalg.norm(v1)
        #     v2_u = v2 / np.linalg.norm(v2)
        #     betta = -np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
        #     center = segmsList[i].GetCenter()

        #     transform  = vtkTransform()
        #     transform.Translate(np.multiply(center, 1))
        #     transform.RotateWXYZ(betta, normalize(np.diff([clinePts.GetPoint(rightId + i), clinePts.GetPoint(rightId + i + 1)], axis=0)[0]))
        #     transform.Translate(np.multiply(center, -1))

        #     transformFilter = vtkTransformFilter()
        #     transformFilter.SetInputData(segmsList[i])
        #     transformFilter.SetTransform(transform)
        #     transformFilter.Update()

        #     segmsList[i] = transformFilter.GetOutput()
        return segmsList




    def tranformStent(self, stent, bdsSegments, centerLine):
        bdsStent = stent.GetBounds()
        lenStent = bdsStent[5] - bdsStent[4]
        midleSegm = int((bdsSegments[2] + bdsSegments[3]) / 2.) #TODO: see in the Issue 11.

        clinePts = centerLine.GetPoints()
        dist = np.inf
        midleId = 0
        for i in range(0, clinePts.GetNumberOfPoints()):
            curdist = getDistance(self.cline[midleSegm], clinePts.GetPoint(i))
            if(curdist < dist):
                dist = curdist
                midleId = i
        midleSegm = midleId
        
        leftcline = []
        rightcline = []
        sum = 0
        i = 0
        while sum < lenStent/2.:
            sum += getDistance(clinePts.GetPoint(midleSegm - i), clinePts.GetPoint(midleSegm - (i + 1)))
            leftcline.append([0., 0., -sum])
            i += 1
        leftId = midleSegm - i

        sum = 0
        i = 0
        while sum < lenStent/2.:
            sum += getDistance(clinePts.GetPoint(midleSegm + i), clinePts.GetPoint(midleSegm + (i + 1)))       
            rightcline.append([0., 0., sum])     
            i += 1
        rightId = midleSegm + i

        clineStent = list(reversed(leftcline)) + [[0., 0., 0.]] + rightcline

        newPts = vtkPoints()
        newPts.Allocate(stent.GetNumberOfPoints())
        newPts.SetNumberOfPoints(stent.GetNumberOfPoints())
        pts = stent.GetPoints()

        suffix = '%(percent)d%% [%(elapsed_td)s / %(eta_td)s]'
        progress_bar = FillingCirclesBar('Stent transformation: ', suffix=suffix, max = len(clineStent)-1)

        segmsList = []
        ListsIDS = []
        # for i in range(0, len(clineStent)- 1):    
        #     progress_bar.next()        
        #     tempListIds = []
        #     segmPts = vtkPoints()
        #     for j in range(0, stent.GetNumberOfPoints()):
        #         pt = pts.GetPoint(j)
        #         flg = True

        #         p2f_init = np.diff([pt, clineStent[i]], axis=0)[0]
        #         d_init = np.dot(p2f_init, [0., 0., 1.])
        #         flg *= d_init < 0
                        
        #         p2f_fin = np.diff([pt, clineStent[i+1]], axis=0)[0]
        #         d_fin = np.dot(p2f_fin, [0., 0., 1.])
        #         flg *= d_fin > 0

        #         if flg:
        #             tempListIds.append(j)
        #             segmPts.InsertNextPoint(pt)
        #     pointSet = vtkPointSet()
        #     pointSet.SetPoints(segmPts)

        #     segmsList.append(pointSet)
        #     ListsIDS.append(tempListIds)

            # trasformPts = self.tranformSegment(pointSet, clinePts.GetPoint(rightId + i), clinePts.GetPoint(rightId + i + 1))

            # for j in range(0, trasformPts.GetNumberOfPoints()):
            #     newPts.SetPoint(tempListIds[j], trasformPts.GetPoint(j))


        pointSet = vtkPointSet()
        pointSet.SetPoints(pts)
        segmsList.append(pointSet)
        segmsList = self.tranformSegment(segmsList, clinePts, rightId)
        for i in range(0, len(segmsList)):
            for j in range(0, segmsList[i].GetNumberOfPoints()):
                newPts.SetPoint(j, segmsList[i].GetPoint(j))


        print('\n')
        transformStent = vtkUnstructuredGrid()
        transformStent.SetPoints(newPts)
        transformStent.SetCells(vtk.VTK_HEXAHEDRON, stent.GetCells())
        return transformStent
    
    def tranformStentOLD(self, stent, bdsSegments, centerLine):
        bdsStent = stent.GetBounds()
        lenStent = bdsStent[5] - bdsStent[4]
        midleSegm = int((bdsSegments[2] + bdsSegments[3]) / 2.) #TODO: see in the Issue 11.

        clinePts = centerLine.GetPoints()
        dist = np.inf
        midleId = 0
        for i in range(0, clinePts.GetNumberOfPoints()):
            curdist = getDistance(self.cline[midleSegm], clinePts.GetPoint(i))
            if(curdist < dist):
                dist = curdist
                midleId = i
        midleSegm = midleId
        
        leftcline = []
        rightcline = []
        sum = 0
        i = 0
        while sum < lenStent/2.:
            sum += getDistance(clinePts.GetPoint(midleSegm - i), clinePts.GetPoint(midleSegm - (i + 1)))
            leftcline.append([0., 0., -sum])
            i += 1
        leftId = midleSegm - i

        sum = 0
        i = 0
        while sum < lenStent/2.:
            sum += getDistance(clinePts.GetPoint(midleSegm + i), clinePts.GetPoint(midleSegm + (i + 1)))       
            rightcline.append([0., 0., sum])     
            i += 1
        rightId = midleSegm + i

        clineStent = list(reversed(leftcline)) + [[0., 0., 0.]] + rightcline

        newPts = vtkPoints()
        newPts.Allocate(stent.GetNumberOfPoints())
        newPts.SetNumberOfPoints(stent.GetNumberOfPoints())
        pts = stent.GetPoints()

        suffix = '%(percent)d%% [%(elapsed_td)s / %(eta_td)s]'
        progress_bar = FillingCirclesBar('Stent transformation: ', suffix=suffix, max = len(clineStent)-1)

        for i in range(0, len(clineStent)- 1):    
            progress_bar.next()        
            tempListIds = []
            segmPts = vtkPoints()
            for j in range(0, stent.GetNumberOfPoints()):
                pt = pts.GetPoint(j)
                flg = True

                p2f_init = np.diff([pt, clineStent[i]], axis=0)[0]
                d_init = np.dot(p2f_init, [0., 0., 1.])
                flg *= d_init < 0
                        
                p2f_fin = np.diff([pt, clineStent[i+1]], axis=0)[0]
                d_fin = np.dot(p2f_fin, [0., 0., 1.])
                flg *= d_fin > 0

                if flg:
                    tempListIds.append(j)
                    segmPts.InsertNextPoint(pt)
            pointSet = vtkPointSet()
            pointSet.SetPoints(segmPts)

            trasformPts = self.tranformSegment(pointSet, clinePts.GetPoint(rightId + i), clinePts.GetPoint(rightId + i + 1))

            for j in range(0, trasformPts.GetNumberOfPoints()):
                newPts.SetPoint(tempListIds[j], trasformPts.GetPoint(j))
        print('\n')
        transformStent = vtkUnstructuredGrid()
        transformStent.SetPoints(newPts)
        transformStent.SetCells(vtk.VTK_HEXAHEDRON, stent.GetCells())
        return transformStent

    def smoothCenterline(self, centerline):
        smooth = vtkSmoothPolyDataFilter()
        smooth.SetInputData(centerline)
        smooth.SetNumberOfIterations(150)
        smooth.SetRelaxationFactor(0.15)
        smooth.FeatureEdgeSmoothingOn()
        smooth.BoundarySmoothingOn()
        smooth.Update()
        return smooth.GetOutput()

    def smoothCenterline2(self, centerline, refinements=10):
        pts = []
        vtkPts = centerline.GetPoints()
        for i in range(vtkPts.GetNumberOfPoints()):
            pts.append(vtkPts.GetPoint(i))

        # Chaikin’s corner cutting algorithm:
        pts = np.array(pts)

        for _ in range(refinements):
            L = pts.repeat(2, axis=0)
            R = np.empty_like(L)
            R[0] = L[0]
            R[2::2] = L[1:-1:2]
            R[1:-1:2] = L[2::2]
            R[-1] = L[-1]
            pts = L * 0.75 + R * 0.25

        # Densification:
        min_dist = 30.0 * 0.036 # Length (max) of FE of stent 

        newPts = []
        i = 0
        count = 0
        while i < len(pts):
            if count == 0:
                newPts.append(pts[i])
                count += 1
            dist = getDistance(newPts[-1], pts[i + count])
            if dist >= min_dist:
                i = i + count
                count = 0
            else:
                count += 1        
            if i + count >= len(pts):
                newPts.append(pts[-1])
                break

        newVtkPts = vtkPoints()
        for pt in newPts:
            newVtkPts.InsertNextPoint(pt)
        newCenterline = vtkPolyData()
        newCenterline.SetPoints(newVtkPts)
        
        return newCenterline


        
        
