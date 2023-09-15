from vtkmodules.all import (
    vtkXMLPolyDataWriter,
    vtkSTLWriter)

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