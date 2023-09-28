from vtkmodules.all import (
    vtkPoints,
    vtkPolyData,
    vtkUnstructuredGrid,
    vtkXMLPolyDataWriter,
    vtkXMLUnstructuredGridWriter,
    vtkSTLWriter)
import numpy as np


def writeSegmentsCSV(filename, segments:list):
    with open(filename, 'w', encoding='utf-8') as file:
        for segment in segments:
            for pt3 in segment:
                file.write(str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + '\n')

def writePointsCSV(filename, points: list):
    with open(filename, 'w', encoding='utf-8') as file:
        for pt3 in points:
            file.write(str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + '\n')

def writePolyDataAsSTL(filename, polydata):
    writer = vtkSTLWriter()
    writer.SetInputData(polydata)
    writer.SetFileName(filename)
    writer.Write()

def writePolyDataAsVTP(filename, polydata: vtkPolyData):
    writer = vtkXMLPolyDataWriter()
    writer.SetInputData(polydata)
    writer.SetFileName(filename)
    writer.Write()

def writeUnstructuredGrid(filename, ugrid: vtkUnstructuredGrid):
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(filename)
    writer.Write()

def writeDisplacementsCSV(filename, points: vtkPoints, offsets: list, steps=10):
    '''
    filename - path to file\n
    poitns - is an array of vx-vy-vz triplets accessible by (point or cell) id\n
    offset - is an array of vx-vy-vz triplets, characterizes the displacement of a point\n
    steps - responsible for the number of displacements (linearly) when installing the stent\n
    ! the structure is similar to exporting displacements from ABAQUS !
    '''
    with open(filename, 'w', encoding='utf-8') as file:
        file.write('Initial coordinates\n')

        for i in range(0, points.GetNumberOfPoints()):
            pt = points.GetPoint(i)
            file.write(str(i+1) + '\t' + str(pt[0]) + '\t' +  str(pt[1]) + '\t' +  str(pt[2]) + '\t')
        file.write('\nDisplacements\n')

        #TODO: думаю, стоит пересмотреть траекторию смещений. В данном случае все линейно.
        for i in range(0, steps):
            line = []
            tempPts = vtkPoints()
            tempPts.DeepCopy(points)
            for j in range(0, tempPts.GetNumberOfPoints()):
                pt = np.add.reduce([tempPts.GetPoint(j), np.multiply(offsets[j], i/10.0)], axis=0)
                diff = np.diff([pt, tempPts.GetPoint(j)], axis=0)[0]
                file.write(str(j+1) + '\t' + str(diff[0]) + '\t' +  str(diff[1]) + '\t' +  str(diff[2]) + '\t')
            file.write('\n')
        
def printComparisonMeasurements(data_list, cur_diams):
    print('VESSEL: Δ MAXIMUM LUMEN DIAMETER: | ', '%0.2f' % data_list[8], ' -  ', '%0.2f' % np.max(cur_diams), '| = ' , '%0.2f' % abs(data_list[8] - np.max(cur_diams)))
    print('VESSEL: Δ MEAN LUMEN DIAMETER:    | ', '%0.2f' % data_list[3], ' -  ', '%0.2f' % np.mean(cur_diams), '| = ', '%0.2f' % abs(data_list[3] - np.mean(cur_diams)))
    print('VESSEL: Δ MINIMUM LUMEN DIAMETER: | ', '%0.2f' % data_list[5], ' -  ', '%0.2f' % np.min(cur_diams), '| = ',  '%0.2f' % abs(data_list[5] - np.min(cur_diams)))

    # https://www.nejm.org/doi/10.1056/NEJM199108153250701
    # https://radcalculators.org/carotid-artery-stenosis-nascet-and-ecst-calculator/
    pr_stenosis = (1 - np.min(cur_diams)/np.mean(cur_diams))*100
    print('STENT:  Δ DIAMETER STENOSIS (%):  |', '%0.2f' % data_list[0], ' - ', '%0.2f' % pr_stenosis,  '| = ', '%0.2f' % abs(data_list[0] - pr_stenosis))