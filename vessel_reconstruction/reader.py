from core import *
import configparser
import pandas
import numpy as np
import vtk
from vtkmodules.all import (
    vtkSTLReader, 
    vtkXMLUnstructuredGridReader, 
    vtkPoints, 
    vtkUnstructuredGrid,
    vtkIdList)
from progress.bar import IncrementalBar


class Reader(object):
    def __init__(self, path_cfg):
        self.config = configparser.ConfigParser()
        if not self.config.read(path_cfg):
            print(f"E: File \"{path_cfg}\" not found.")

        self.outpath = self.config["CONFIG"]["PathOutput"]
        # column names in the table:
        self.param = [
            'OBSTRUCTION: DIAMETER STENOSIS PRE (%)',           #0
            'STENT: DIAMETER STENOSIS POST (%)',                #1
            'VESSEL: DIAMETER STENOSIS POST (%)',               #2
            'VESSEL: MEAN LUMEN DIAMETER PRE (mm)',             #3
            'VESSEL: MEAN LUMEN DIAMETER POST (mm)',            #4
            'OBSTRUCTION: MINIMUM LUMEN DIAMETER PRE (mm)',     #5
            'IN-SEGMENT: MINIMUM LUMEN DIAMETER POST (mm)',     #6
            'VESSEL: MINIMUM LUMEN DIAMETER POST (mm)',         #7
            'VESSEL: MAXIMUM LUMEN DIAMETER PRE (mm)',          #8
            'VESSEL: MAXIMUM LUMEN DIAMETER POST (mm)',         #9
            'DIST: MAXIMUM LUMEN DIAMETER POST (mm)',           #10
            'PROX: MAXIMUM LUMEN DIAMETER POST (mm)',           #11
            'STENT: MAXIMUM LUMEN DIAMETER POST (mm)',          #12
            'IN-SEGMENT: MAXIMUM LUMEN DIAMETER POST (mm)'      #13
            ]
        
        # parameter values
        self.data_list = []

        # segment numbers [initVessel, finVessel, initStent, finStent]
        self.bdsSegments = [0, 0, 0, 0]

        # list of segments3d [[[point3,..],...],[[point3,..],...],...]
        self.segms3DLumen = []
        self.segms3DWall = []

        self.isPrint = True

        # list of centers segments
        self.centerline = []

        # name of the offset function
        self.funcOffset = self.config["CONFIG"]["functionOffset"]

        # data .stl
        self.lumenStl = self.readStl('data/' + self.config["CONFIG"]["ID"] + '/' + self.config["CONFIG"]["ID"] + '-BL/LUMEN/Lumen_mesh.stl')
        self.wallStl = self.readStl('data/' + self.config["CONFIG"]["ID"] + '/' + self.config["CONFIG"]["ID"] + '-BL/VESSEL WALL/VesselWall_mesh.stl')
        self.stent = self.readStent(self.config["CONFIG"]["PathStentNodes"], self.config["CONFIG"]["PathStentElements"])

    def update(self):
        """
        Reading all data and corrects the order of the segments.
        """
        self.readContours2D()
        self.readXLSX()
        self.readContours3D()
        self.correctsOrderOfSegments()

    def correctsOrderOfSegments(self):
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
        segm_init_normal = getFaceNormal(segm_init)
        segm_fin = self.segms3DLumen[newFinVessel]
        segm_fin_normal = getFaceNormal(segm_fin)

        part1 = []
        part2 = []
        part3 = []
        for i in range(0, newInitVessel):
            segment = self.segms3DLumen[i]
            flg_init = True
            flg_fin = True
            for pt in segment:
                p2f_init = np.diff([pt, segm_init[0]], axis=0)[0]
                d_init = np.dot(p2f_init, segm_init_normal)
                flg_init *= d_init < 0
                
                p2f_fin = np.diff([pt, segm_fin[0]], axis=0)[0]
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
            print('New number of lumen contours: ', len(self.segms3DLumen),'\n')
    
    def readXLSX(self):
        """
        Search for parameter values in the table.
        """
        table = pandas.read_excel(self.config["CONFIG"]["PathDataSet"])
        maskID = table['Subject ID'].str.contains(self.config['CONFIG']['ID'], na=True)
        for i in self.param:
            self.data_list.append(float(table.loc[maskID, i].iat[0]))
        if self.isPrint:
            print('\nParameters:')
            print(["%.2f" % var for var in self.data_list])

    def readContours2D(self):
        """
        Read .ctr file.
        Search for the initial and final numbers segments of the stented vessel
        param.
        """
        with open('data/' + self.config["CONFIG"]["ID"] + '/' + self.config["CONFIG"]["ID"] + '-BL/CONTOURS/LumenExport.ctr') as file:
            lines = file.readlines()
            self.bdsSegments[0] = int(lines[2].split(' ')[0])
            self.bdsSegments[1] = int(lines[len(lines)-1].split(' ')[0])

        with open('data/' + self.config["CONFIG"]["ID"] + '/' + self.config["CONFIG"]["ID"] + '-BL/CONTOURS/StentExport.ctr') as file:
            lines = file.readlines()
            self.bdsSegments[2] = int(lines[2].split(' ')[0])
            self.bdsSegments[3] = int(lines[len(lines)-1].split(' ')[0])

        if self.isPrint:
            print('Segment numbers:\n[initVessel, finVessel, initStent, finStent] = ',self.bdsSegments)

    def readContours3D(self):
        """
        Read .ctr file.
        Creates a field of points based on data from the file.
        """
        self.__readContours3D('data/' + self.config["CONFIG"]["ID"] + '/' + self.config["CONFIG"]["ID"] + '-BL/LUMEN/Lumen_contours_3d.txt', self.segms3DLumen)
        self.__readContours3D('data/' + self.config["CONFIG"]["ID"] + '/' + self.config["CONFIG"]["ID"] + '-BL/VESSEL WALL/VesselWall_contours_3d.txt', self.segms3DWall)   
        if self.isPrint:
            print('Number of lumen contours: ', len(self.segms3DLumen))
            print('Number of wall contours: ', len(self.segms3DWall))

    def __readContours3D(self, pathName, segms3D):
        with open(pathName) as file:
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
    
    def readStl(self, filename):
        reader = vtkSTLReader()
        reader.SetFileName(filename)
        reader.Update()
        return reader.GetOutput()
    
    def readUnstructuredGrid(self, filename):
        reader = vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()
        self.volumeMesh = reader.GetOutput()

    def readStent(self, fileNodes, fileElements):
        ptIds = []
        pts = vtkPoints()
        stent = vtkUnstructuredGrid()

        with open(fileNodes) as fNodes:
            lines = fNodes.readlines()
            for i in range(0, len(lines)):
                line = lines[i].split('\t')
                ptIds.append(int(line[0]))
                pts.InsertNextPoint([float(line[1]), float(line[2]), float(line[3])])   
        stent.SetPoints(pts)

        with open(fileElements) as fElems:
            lines = fElems.readlines()
            stent.Allocate(len(lines))
            print(len(lines))
            print(pts.GetNumberOfPoints())
            progress_bar = IncrementalBar('Read node_elements.csv: ', max = len(lines))
            for i in range(0, len(lines)):
                progress_bar.next()
                line = lines[i].split('\t')[1:]
                idList = vtkIdList()
                idList.SetNumberOfIds(len(line))
                for j in range(0, len(line)):
                    idList.InsertId(j, ptIds.index(int(line[j])))
                stent.InsertNextCell(vtk.VTK_HEXAHEDRON, idList)
        return stent
