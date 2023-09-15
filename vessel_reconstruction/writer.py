from vtkmodules.all import (
    vtkPoints,
    vtkXMLPolyDataWriter,
    vtkSTLWriter)
from core import *

def writeSegmentsCSV(filename, segments):
    with open(filename, 'w', encoding='utf-8') as file:
        for segment in segments:
            for pt3 in segment:
                file.write(str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + '\n')

def writePointsCSV(filename, points):
    with open('output2.txt', 'w', encoding='utf-8') as file:
        for pt3 in points:
            file.write(str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + '\n')

def writeCenterLine(filename, cline):
    with open('data/result/centerline.txt', 'w', encoding='utf-8') as file:
        for pt3 in cline:
            file.write(str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + '\n')

def writePolyDataAsSTL(filename, polydata):
    writer = vtkSTLWriter()
    writer.SetInputData(polydata)
    writer.SetFileName(filename)
    writer.Write()

def writePolyDataAsVTP(filename, polydata):
    writer = vtkXMLPolyDataWriter()
    writer.SetInputData(polydata)
    writer.SetFileName(filename)
    writer.Write()

def writeDisplacementsCSV(filename, points: vtkPoints, offsets):
    '''
    filename - path to file\n
    poitns - is an array of vx-vy-vz triplets accessible by (point or cell) id\n
    offset - is an array of vx-vy-vz triplets, characterizes the displacement of a point
    '''
    with open(filename, 'w', encoding='utf-8') as file:
        file.write('Initial coordinates')
        line = ''
        for i in range(0, points.GetNumberOfPoints()):
            pt = points.GetPoint(i)
            #TODO: FIX IT! ->
            line += str(i+1) + '\t' + str(pt[0]) + '\t' +  str(pt[1]) + '\t' +  str(pt[2]) + '\t'
        file.write(line)
        file.write('Displacements')

        for i in range(0, 10):
            line = []
            tempPts = vtkPoints()
            tempPts.DeepCopy(points)
            for j in range(0, tempPts.GetNumberOfPoints()):
                pt = add(tempPts.GetPoint(j), multy(offsets[j], i/10))
                diff = sub(tempPts.GetPoint(j), pt)
                line += str(j+1) + '\t' + str(pt[0]) + '\t' +  str(pt[1]) + '\t' +  str(pt[2]) + '\t'


