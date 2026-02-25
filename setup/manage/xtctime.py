#!/usr/bin/env python3

from struct import unpack
import traceback
import sys

for xtc in sys.argv[1:]:
    try:
        f = open(xtc, 'rb')
 
        tag = f.read(8)                           # Tag: magic number and number of atoms
        n = 92 + unpack(">l",f.read(84)[-4:])[0]  # Size of frame in bytes
 
        f.seek(int(-5*n/4), 2)                    # This should contain a complete frame
        frame = f.read()                          # Read the remaining part in
        n = frame.index(tag)                      # Find the tag
 
        print(xtc, 0.001*unpack(">f",frame[n+12:n+16])[0])
    except Exception as excpt:
        print(xtc, "-1")
        #raise excpt
