#!/usr/bin/env python
from reader import Reader
import os

def writeCSV(segments, filename):
    with open(filename, 'w', encoding='utf-8') as file:
        for segment in segments:
            for pt3 in segment:
                file.write(str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + '\n')

def writeCSV_2(segment):
    with open('output2.txt', 'w', encoding='utf-8') as file:
        for pt3 in segment:
            file.write(str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + '\n')

def writeCenterLine(cline):
    with open('data/result/centerline.txt', 'w', encoding='utf-8') as file:
        for pt3 in cline:
            file.write(str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + '\n')

def main():
    reader = Reader(os.getcwd() + '/config/init.ini')
    reader.update()
    reader.correctingOrderOfSegments()
    reader.createCenterline()
    writeCSV(reader.segms3DLumen, 'data/result/output.csv')
    writeCenterLine(reader.centerline)
    reader.readStl()
    reader.calcOffset()
    reader.tempSTL()


if __name__ == '__main__':
    main()