#!/usr/bin/env python
from reader import Reader
from core import *
from mesher import runMesher
import os


def main():
    reader = Reader(os.getcwd() + '/config/init.ini')
    reader.update()
    
    centerline = createCenterline(reader.segms3DLumen)
    offset = calcOffset(reader.bdsSegments, reader.data_list, reader.funcOffset)
    
    if not os.path.isfile('data/result/volumeMesh.vtu'):
        runMesher(reader)

    reader.readUnstructuredGrid('data/result/volumeMesh.vtu')
    createDisplacementWall(reader.volumeMesh, reader.segms3DLumen, reader.bdsSegments, centerline, offset)

if __name__ == '__main__':
    main()