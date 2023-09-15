import configparser
import pandas
import numpy as np
from vtkmodules.all import (
    vtkCellArray,
    vtkPoints,
    vtkPolyData,
    vtkSTLReader,
    vtkSTLWriter)


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

def normalize(pt: list) ->list:
    return truediv(pt, np.linalg.norm(pt))


class Reader(object):
    def __init__(self, path_cfg):
        self.config = configparser.ConfigParser()
        if not self.config.read(path_cfg):
            print(f"E: File \"{path_cfg}\" not found.")

        #TODO: read param from file? 
        self.param = [
            'OBSTRUCTION: DIAMETER STENOSIS PRE (%)',
            'STENT: DIAMETER STENOSIS POST (%)',
            'VESSEL: DIAMETER STENOSIS POST (%)',
            'VESSEL: MEAN LUMEN DIAMETER PRE (mm)',
            'VESSEL: MEAN LUMEN DIAMETER POST (mm)',
            'OBSTRUCTION: MINIMUM LUMEN DIAMETER PRE (mm)',
            'IN-SEGMENT: MINIMUM LUMEN DIAMETER POST (mm)']
        # parameter values
        self.data_list = []
        # segment numbers [initVessel, finVessel, initStent, finStent]
        self.bdsSegments = [0, 0, 0, 0]
        # list of segments3d [[[point3,..],...],[[point3,..],...],...]
        self.segms3DLumen = []
        self.segms3DWall = []
        self.isPrint = True
        self.centerline = []

    def update(self):
        self.readContours2D()
        self.readXLSX()
        self.readContours3D()

    def correctingOrderOfSegments(self):
        """
        In the vessel lumen data, part of the segments is duplicated. 
        This function replaces the extra segments from the middle 
        with segments from the end of the list.
        """
        if self.isPrint:
            print('\nCorrecting order of segments...') 

        newInitVessel = 0
        npts = len(self.segms3DLumen[0])
        for i in range(1, len(self.segms3DLumen)):
            if len(self.segms3DLumen[i]) != npts:
                newInitVessel = i
                break

        newFinVessel = len(self.segms3DLumen) - 1

        segm_init = self.segms3DLumen[newInitVessel]
        segm_init_normal = getNormal(segm_init)
        segm_fin = self.segms3DLumen[newFinVessel]
        segm_fin_normal = getNormal(segm_fin)

        part1 = []
        part2 = []
        part3 = []
        for i in range(0, newInitVessel):
            segment = self.segms3DLumen[i]
            flg_init = True
            flg_fin = True
            for pt in segment:
                p2f_init = sub(segm_init[0], pt)
                d_init = np.dot(p2f_init, segm_init_normal)
                flg_init *= d_init < 0
                
                p2f_fin = sub(segm_fin[0], pt)
                d_fin = np.dot(p2f_fin, segm_fin_normal)
                flg_fin *= d_fin > 0

            if flg_init:
                part1.append(segment)
            if flg_fin:
                part3.append(segment)
        part2 = self.segms3DLumen[newInitVessel: newFinVessel]
        for segms in part2:
            segms.reverse()
        
        self.segms3DLumen = part1 + part2 + part3
        
        if self.isPrint:
            print('New number of lumen contours: ', len(self.segms3DLumen))
    
    def checkIntersections(self):
        """
        Checking for intersections between adjacent sections.
        """
        if self.isPrint:
            print('\nCheck intersections...') 

        badlist = []

        for i in range(1, len(self.segms3DLumen)-1, 2):
            flg = True
            normal = getNormal(self.segms3DLumen[i])
            for pt1 in self.segms3DLumen[i]:
                for pt2 in self.segms3DLumen[i - 1]:
                    p2f = sub(pt1, pt2)
                    d = np.dot(p2f, normal)
                    flg *= (d > 0)

                for pt2 in self.segms3DLumen[i + 1]:
                    p2f = sub(pt1, pt2)
                    d = np.dot(p2f, normal)
                    flg *= (d < 0)

            if not flg:
                badlist.append(i)

        for i in sorted(badlist, reverse=True):
            del self.segms3DLumen[i]         
        if self.isPrint:
            print('New number of lumen contours: ', len(self.segms3DLumen))   

    def simplifySegments(self, dist=0.5, n=1, maxPts = 40):
        """
        Delete some points and segments.

        dist - required distance between segments,\n
        n - number of iterations,\n
        maxPts - required number of points in the segment.
        """
        if self.isPrint:
            print('\nSimplify segments...')  
        
        npts = []
        for segm in self.segms3DLumen:
            npts.append(len(segm))
        if (np.min(npts) > maxPts):
            maxPts = np.min(npts)


        for segm in self.segms3DLumen:
            k = int(len(segm)/maxPts)
            for i in range(len(segm)-1, 0, -1):
                if i%k !=0:
                    del segm[i]
        
        for j in range(0, n):
            badlist = []
            for i in range(0, len(self.segms3DLumen)-1, 2):
                flg = True
                for pt1 in self.segms3DLumen[i]:
                    for pt2 in self.segms3DLumen[i - 1]:
                        p2f = sub(pt1, pt2)
                        flg *= np.linalg.norm(p2f) > dist

                    for pt2 in self.segms3DLumen[i + 1]:
                        p2f = sub(pt1, pt2)
                        flg *= np.linalg.norm(p2f) > dist
                if not flg:
                    badlist.append(i)

            for i in sorted(badlist, reverse=True):
                del self.segms3DLumen[i]           
        if self.isPrint:
            print('New number of lumen contours: ', len(self.segms3DLumen))       

    def readXLSX(self):
        """
        Search for parameter values in the table.
        """
        table = pandas.read_excel(self.config["DATA"]["PathDataSet"])
        maskID = table['Subject ID'].str.contains(self.config['DATA']['ID'], na=True)
        for i in self.param:
            self.data_list.append(float(table.loc[maskID, i].iat[0]))
        if self.isPrint:
            print('\nParameters:')
            print(self.data_list)

    def readContours2D(self):
        """
        Read .ctr file.
        Search for the initial and final numbers segments of the stented vessel
        param.
        """
        with open(self.config["DATA"]["PathLumenContours"]) as file:
            lines = file.readlines()
            self.bdsSegments[0] = int(lines[2].split(' ')[0])
            self.bdsSegments[1] = int(lines[len(lines)-1].split(' ')[0])

        with open(self.config["DATA"]["PathStentContours"]) as file:
            lines = file.readlines()
            self.bdsSegments[2] = int(lines[2].split(' ')[0])
            self.bdsSegments[3] = int(lines[len(lines)-1].split(' ')[0])

        if self.isPrint:
            print('segment numbers [initVessel, finVessel, initStent, finStent]:')
            print(self.bdsSegments)
    
    def readContours3D(self):
        """
        Read .ctr file.
        Creates a field of points based on data from the file.
        """
        self.__readContours3D('PathLumen3DContours', self.segms3DLumen)
        self.__readContours3D('PathWall3DContours', self.segms3DWall)   
        if self.isPrint:
            print('Number of lumen contours: ', len(self.segms3DLumen))
            print('Number of wall contours: ', len(self.segms3DWall))

    def __readContours3D(self, pathName, segms3D):
        with open(self.config['DATA'][pathName]) as file:
            lines = file.readlines()
            npts = -1
            points3 = []
            for i in range(1, len(lines)):
                line = lines[i].split(' ')
                if len(line) < 3: # skip an empty line
                    continue

                if line[1] == 'Contour':
                    segms3D.append(points3.copy())
                    points3.clear()
                    continue

                if line[1] == 'Number':
                    npts = int(line[5])
                    continue
                points3.append([float(line[0]), float(line[1]), float(line[2])])
            segms3D.append(points3.copy())

    def createCenterline(self):
        for pts in self.segms3DLumen:
            newPt = [0., 0., 0.]
            for pt in pts:
                newPt = add(newPt, pt)
            newPt = truediv(newPt, len(pts))
            self.centerline.append(newPt)

    def readStl(self):
        reader = vtkSTLReader()
        reader.SetFileName(self.config["DATA"]["PathLumenSTL"])
        reader.Update()
        self.lumenStl = reader.GetOutput()
    
    def writeStl(self, polydata, filename):
        writer = vtkSTLWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polydata)
        writer.Write()
       
    def calcOffset(self):
        self.offset = []
        N = self.bdsSegments[3]-self.bdsSegments[2]
        for i in range(0, N):
            self.offset.append(self.funcSin(i))
    
    def funcSin(self, k):
        return (np.abs(self.data_list[5] - self.data_list[6]) *
                (np.sin(np.pi*k/(self.bdsSegments[3]-self.bdsSegments[2]))))
    
    def tempSTL(self):
        cellArray= vtkCellArray()
        cellArray.DeepCopy(self.lumenStl.GetPolys())
        pts = vtkPoints()
        pts.DeepCopy(self.lumenStl.GetPoints())

        listPts = [] #idx pts por preparation

        segm_init = self.segms3DLumen[self.bdsSegments[2]]
        segm_init_normal = getNormal(segm_init)
        segm_fin = self.segms3DLumen[self.bdsSegments[3]]
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
            for j in range(self.bdsSegments[2], self.bdsSegments[3]):
                listdist.append(getDistance(pt, self.centerline[j]))
            minIdx = np.argmin(listdist)

            offset2vec.append(multy(normalize(sub(self.centerline[int(minIdx+self.bdsSegments[2])], pt)), self.offset[minIdx]))

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

            self.writeStl(newPD ,'data/result/offset_' + str(j) + '.stl')
