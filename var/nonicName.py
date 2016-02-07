#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser(description='Map the NonicPressed ID to name, e.g., 6 -> nonic095r125c04375b')

parser.add_argument('--id', type=int, required=True, help='the NonicPressed ID')
parser.add_argument('--getPress','-gb', help='get PRESS parameter as float', action="store_true")
parser.add_argument('--getStretch','-gc', help='get STRETCH parameter as float', action="store_true")

args = parser.parse_args()

ch = 0.25
bh = 0.0875

c = -ch
b = -bh
depth = 0
for at in str(args.id):
    atomid = float(at)
    c += atomid * ch / 10**depth
    b += atomid * bh / 10**depth
    depth += 1

cstr = str(c).replace('.','').rstrip('0')
bstr = str(b).replace('.','').rstrip('0')

if cstr == '':
    cstr = '0'
if bstr == '':
    bstr = '0'
    
if args.getStretch:
    print c
elif args.getPress:
    print b
else:
    print 'nonic095r' + cstr + 'c' + bstr + 'b'
