#!/usr/bin/python

import sys
import shutil
print 'Input:', str(sys.argv)

if ( len(sys.argv) < 3 ) :
	sys.exit()

fromDir = sys.argv[ 1 ]
toDir = sys.argv[ 2 ]

print "copying : ", fromDir + "params.dat", " to ", toDir
shutil.copyfile( fromDir+"params.dat", toDir+"vpd_4DB.dat" )

print "copying : ", fromDir + "rQA.pdf", " to ", toDir
shutil.copyfile( fromDir+"rQA.pdf", toDir+"qa.pdf" )

print "copying : ", fromDir + "log.txt", " to ", toDir
shutil.copyfile( fromDir+"log.txt", toDir+"log.txt" )
