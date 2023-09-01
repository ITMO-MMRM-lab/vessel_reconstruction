#!/usr/bin/env python
from reader import Reader
import preparation
import os

def writeCSV(segments):
    with open('output.txt', 'w', encoding='utf-8') as file:
        for segment in segments:
            for pt3 in segment:
                file.write(str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + '\n')


def main():
    reader = Reader(os.getcwd() + '/config/init.ini')
    reader.update()
    reader.correctingOrderOfSegments()
    writeCSV(reader.segms3DLumen)


if __name__ == '__main__':
    main()