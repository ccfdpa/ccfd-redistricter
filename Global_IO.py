# Global Settings
import socket
import Global_Settings as G

# Key folder navigation for all files

thisRoot=G.thisRoot
thisPath=G.thisPath
this_IO_File=G.this_IO_File

# Key files/variables
# Input County Data and Shapes
# Data structural design
global rawDemoData
rawDemoData=['Size']

# Output
global filetoQGIS
filetoQGIS='Initial District Assignments.csv'   # This file has district name and color directions appended to the csv infile above and saves it.
global districtSummary
districtSummary='Initial District Summary.csv'  # this file is a csv file with one row per district and columns of miscellaneous district info.
global districtShapefile
districtShapefile='PA State Senate Districts.shp'  # Output: An ESRI shapefile with district information and color directions

# The part below is saved in individual txt files
# Note: the ID is the first element in the variable list

# Input State Data and Shape
global StateShapeFile, StateVariableDict, StateCoder

# Input County Data and Shapes
global CountyShapeFile, CountyVariableDict,CountyCoder

# Input Municipal Data and Shapes
global MunicipalShapefile, MunicipalVariableDict,MunicipalCoder

# Input Detailed Attributes as CSV (at VTD level)
global VTDShapeFile, VTDVariableDict, VTDCoder

# Input scale adjustments
global VTDinputAreaScaleFactor
global MunicipalinputAreaScaleFactor

# Input initial forced inclusions (usually missing neighbors)
global inputNeighborForced
inputNeighborForced=[]

def setGlobal(var,value):
    globals()[var]=deepcopy(value)
    return

def getGlobal(var):
    return globals()[var]

