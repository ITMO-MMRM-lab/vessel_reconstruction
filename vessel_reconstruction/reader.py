#!/usr/bin/env python
import configparser
import os
import pandas
import numpy as np


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
        # list of segments3d
        self.segms3DLumen = []
        self.segms3DWall = []

    def update(self):
        self.readContours2D()
        self.readXLSX()
        self.readContours3D()

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
    
    def readContours3D(self):
        """
        Read .ctr file.
        Creates a field of points based on data from the file.
        """
        self.__readContours3D('PathLumen3DContours', self.segms3DLumen)
        self.__readContours3D('PathWall3DContours', self.segms3DWall)    

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

    def readXLSX(self):
        """
        Search for parameter values in the table.
        """
        param = [
            'OBSTRUCTION: DIAMETER STENOSIS PRE (%)',
            'STENT: DIAMETER STENOSIS POST (%)',
            'VESSEL: DIAMETER STENOSIS POST (%)',
            'VESSEL: MEAN LUMEN DIAMETER PRE (mm)',
            'VESSEL: MEAN LUMEN DIAMETER POST (mm)',
            'OBSTRUCTION: MINIMUM LUMEN DIAMETER PRE (mm)',
            'IN-SEGMENT: MINIMUM LUMEN DIAMETER POST (mm)']

        table = pandas.read_excel(self.config["DATA"]["PathDataSet"])
        maskID = table['Subject ID'].str.contains(self.config['DATA']['ID'], na=True)
        data1 = float(table.loc[maskID, 'VESSEL: DIAMETER STENOSIS POST (%)'].iat[0])

        for i in param:
            self.data_list.append(float(table.loc[maskID, i].iat[0]))


def main():
    reader = Reader(os.getcwd() + '/config/init.ini')
    reader.update()


if __name__ == '__main__':
    main()