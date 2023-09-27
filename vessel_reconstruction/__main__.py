# !/usr/bin/env python
from reader import Reader
from core import DataAlgorithms
from mesher import runMesher
from writer import writeDisplacementsCSV, writePolyDataAsSTL, writeComparisonMeasurements
import os


def main():
    reader = Reader(os.getcwd() + '/config/init.ini')
    reader.update()
  
    if not os.path.isfile(reader.outpath + 'volumeMesh.vtu'):
        runMesher(reader.lumenStl, reader.wallStl, reader.outpath)

    reader.readUnstructuredGrid(reader.outpath + 'volumeMesh.vtu')

    dataAlg = DataAlgorithms(reader.volumeMesh, reader.lumenStl, reader.segms3DLumen, reader.bdsSegments, reader.data_list, reader.funcOffset)

    writeComparisonMeasurements(reader.data_list, dataAlg.cur_diams)

    writeDisplacementsCSV('data/output/vessel_disp.csv', reader.volumeMesh.GetPoints(), dataAlg.displs, 10)

    for i in range(0, len(dataAlg.trajs)):
        writePolyDataAsSTL('data/output/offsetVessel_' + str(i) + '.stl', dataAlg.trajs[i])


if __name__ == '__main__':
    main()