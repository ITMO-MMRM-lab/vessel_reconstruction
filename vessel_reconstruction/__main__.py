# !/usr/bin/env python
import os

from core import DataAlgorithms
from mesher import runMesher
from reader import Reader
from writer import (
    printComparisonMeasurements,
    writeDisplacementsCSV,
    writePolyDataAsSTL,
    writePolyDataAsVTP,
    writeUnstructuredGrid)


def main():
    reader = Reader(os.getcwd() + '/config/init.ini')
    reader.update()
  
    if not os.path.isfile(reader.outpath + 'volumeMesh.vtu'):
        runMesher(reader.lumenStl, reader.wallStl, reader.outpath)

    reader.readVolumeMesh(reader.outpath + 'volumeMesh.vtu')

    dataAlg = DataAlgorithms(reader.volumeMesh, reader.lumenStl, reader.segms3DLumen, reader.bdsSegments, reader.data_list, reader.funcOffset)

    printComparisonMeasurements(reader.data_list, dataAlg.cur_diams)

    writeDisplacementsCSV('data/output/vessel_disp.csv', reader.volumeMesh.GetPoints(), dataAlg.displs, 10)

    for i in range(0, len(dataAlg.trajs)):
        writePolyDataAsSTL('data/output/offsetVessel_' + str(i) + '.stl', dataAlg.trajs[i])

    reader.readStent()
    if not os.path.isfile(reader.config["CONFIG"]["PathStentVTU"]):
        writeUnstructuredGrid(reader.config["CONFIG"]["PathStentVTU"], reader.stent)
    centerline = reader.readPolyData(reader.outpath + 'centerline.vtp')

    cline_smooth = dataAlg.smoothCenterline2(centerline)

    writePolyDataAsVTP(reader.outpath + 'centerline_smooth.vtp', cline_smooth)
    writeUnstructuredGrid(reader.outpath + 'tranformStent.vtu', dataAlg.tranformStent(reader.stent, reader.bdsSegments, cline_smooth))


if __name__ == '__main__':
    main()