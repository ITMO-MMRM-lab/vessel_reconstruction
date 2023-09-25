# !/usr/bin/env python
from reader import Reader
from core import DataAlgorithms
from mesher import runMesher
import os


def main():
    reader = Reader(os.getcwd() + '/config/init.ini')
    reader.update()

    dataAlg = DataAlgorithms()

    centerline = dataAlg.createCenterline(reader.segms3DLumen)

    dataAlg.analysisOfVessel(reader.lumenStl, centerline, reader.bdsSegments, reader.data_list)
    
    offset = dataAlg.calcOffset(reader.bdsSegments, reader.data_list, reader.funcOffset)
  
    if not os.path.isfile(reader.outpath + 'volumeMesh.vtu'):
        runMesher(reader)

    reader.readUnstructuredGrid(reader.outpath + 'volumeMesh.vtu')

    dataAlg.createDisplacementWall(reader.volumeMesh, reader.segms3DLumen, reader.bdsSegments, centerline, offset)


if __name__ == '__main__':
    main()