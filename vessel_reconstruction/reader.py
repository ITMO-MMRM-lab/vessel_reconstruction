import configparser
import pandas
import numpy as np


def sub(pt1: list, pt2: list) -> list:
    return ([pt1[0] - pt2[0], 
             pt1[1] - pt2[1], 
             pt1[2] - pt2[2]])

def getNormal(segm: list) -> list:
    dir1 = sub(segm[0], segm[5])
    dir2 = sub(segm[10], segm[5])
    segm_normal = np.cross(dir1, dir2)
    return segm_normal/np.linalg.norm(segm_normal)

class Reader(object):
    def __init__(self, path_cfg):
        self.config = configparser.ConfigParser()
        if not self.config.read(path_cfg):
            print(f"E: File \"{path_cfg}\" not found.")

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
        
        self.segms3DLumen = part1 + part2 + part3
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