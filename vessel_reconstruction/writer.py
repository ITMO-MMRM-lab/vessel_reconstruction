from vtkmodules.all import (
    vtkPoints,
    vtkPolyData,
    vtkXMLPolyDataWriter,
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
        

