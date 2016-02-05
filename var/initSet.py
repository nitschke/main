#!/usr/bin/python

import csv
import argparse

parser = argparse.ArgumentParser(description='Set values in an AMDiS init file')

parser.add_argument('-f', '--file', required=True, help='init file')
parser.add_argument('-k', '--key', required=True, help='key, e.g., userParameter->K0')
parser.add_argument('-v', '--value', required=True, help='value')

args = parser.parse_args()

file1 = open(args.file, 'rb')
reader = csv.reader(file1, delimiter=':', skipinitialspace=True)
rowList = []
for row in reader:
    rowList.append(row)
file1.close()

maxKeyLength = 0
for row in rowList:
    if row != [] : 
        kl = len(row[0]) 
    else: 
        kl = 0
    if maxKeyLength < kl: maxKeyLength = kl

for row in rowList:
    if row != [] :
        if ('%' not in row[0]) and (args.key in row[0]):
            if len(row) > 1:
                keyAndCom = row[1].split('%')
                keyAndCom[0] = args.value
                row[1] = keyAndCom.pop(0)
                for s in keyAndCom:
                    row[1] += " %" + s
            else:
                row.append(args.value)
        kl = len(row[0])
        whiteSpaceLen = maxKeyLength - kl
        row[0] += whiteSpaceLen*" "
        if len(row) > 1: row[1] = " " + row[1]

file2 = open(args.file, 'wb')
writer = csv.writer(file2, delimiter=':', lineterminator='\n')
writer.writerows(rowList)
file2.close()
