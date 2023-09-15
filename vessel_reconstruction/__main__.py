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
    tempSTL(reader.lumenStl, reader.segms3DLumen, reader.bdsSegments, centerline, offset)
    runMesher(reader)
    reader.readUnstructuredGrid('data/result/volumeMesh.vtu')


if __name__ == '__main__':
    main()