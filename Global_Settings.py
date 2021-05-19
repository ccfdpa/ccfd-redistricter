# Global Settings
import socket
from shapely.geometry import Polygon
from copy import deepcopy


# Special
global g
# Key design criteria
global populationTotal
populationTotal=12702379  # This is used as a control total.
global numDistricts
numDistricts=18
global districtPrefix
districtPrefix='C'              # The convention I'm using is 'C': US congressional; 'S': PA Senate; 'H': PA House
global allowedPopVariance
allowedPopVariance=0.05

# Key folder navigation for all files
global thisPath,dataPath,thisRoot
if socket.gethostname()=='HessVAIO':
    thisRoot='L:\\Hess Work\\Districting\\'
elif socket.gethostname()=='StudyDesktop':
    thisRoot='L:\\HessViao\\Hess Work\\Districting\\'
else:
    thisPath='C:\\'

thisPath=thisRoot+'MemberList Assemblers Solutions\\'
dataPath=thisRoot+'MemberList Assemblers Raw Data\\'
dataInitID='2010'

# IO mapping (small special purpose py file)
global this_IO_File
this_IO_File='Turzai_2010_IO'

# Precedence relationship/codelength for political subdivisions (Steps 1 and 2)
# This is created as a list of tuples because the order in the list sets the precedence
# List the precedence ID and the start character (regionMemberList variables at this precedence will start with this character
global regionPrecedence
regionPrecedence=[('state','S',1),('county','C',4),('municipal','M',9),('VTD','V',13)]


#**************************** Do Not Modify Below This Line *****************************
#****************************************************************************************

# Some data sources require special adjustments. Enter these in the special IO file.
global regionExclusions
regionExclusions=[]
global neighborExclusions   # Used to remove manually neighbors that share extremely short common boundaries
neighborExclusions=[]
global neighborForced # Connects two regions as neighbors in a graph representation for special reasons. They might become a group.
neighborForced=[]

# Output
global filetoQGIS
filetoQGIS='Initial District Assignments.csv'   # This file has district name and color directions appended to the csv infile above and saves it.
global districtSummary
districtSummary='Initial District Summary.csv'  # this file is a csv file with one row per district and columns of miscellaneous district info.
global districtShapefile
districtShapefile='PA State Senate Districts.shp'  # Output: An ESRI shapefile with district information and color directions

# Operating assumptions
global verbose
verbose=False
global showPlots
showPlots=False
global calculateGeoData
calculateGeoData=False
global populationCheckFlag
populationCheckFlag=True


# Data structural design and initial values
global regionDemogData
regionDemogData=['TAPersons']
# Note elements of regionData not in regionDemoData are considered geometry data
# They can be recomputed with function CCFD_shapes.calcRegionGeoData(inRegionSet)
global regionData
global dataDict
regionData=regionDemogData+['ID','Perimeter','Area','North','South','East','West','Latitude','Longitude','BdgSqrIndex']
dataDict=dict.fromkeys(regionData)
# Initial values (not too sensitive)
for x in regionDemogData:
    dataDict[x]=0
for x in ['Perimeter','Area','BdgSqrIndex']:
    dataDict[x]=0.0
dataDict['North']=0.0
dataDict['South']=90.0
dataDict['East']=0.0
dataDict['West']=-180.0

global regionMemberSet
regionMemberSet={'ID':'','Name':'','Group':set(),'Districts':{districtPrefix:None},'Neighbors':set(),'Data':dataDict,'geometry':Polygon([])}
global districtMemberSet
districtMemberSet={'ID':'','Regions':set(),'Neighbors':set(),'distNeighbors':set(),'Data':dataDict,'geometry':Polygon([])}

# Initial storage locations for use during a run
# Note there should be a global xxxMemberList for each of the above.
global stateMemberList
stateMemberList={}
# global countyMemberList
# countyMemberList={}
global municipalMemberList
municipalMemberList={}
global VTDMemberList
VTDMemberList={}

# Intermediate storage locations for use during a run
# Global IO Flags
global memberListOpenFlags
memberListOpenFlags={
            'boundaryDict':False,
            'county':False,           # Keeps track whether countyMemberList is Open
            'municipal':False,        # Keeps track: municipalMemberList
            'VTD':False,              # Keeps track: VTDMemberList
            'state':False             # Keeps track: stateMemberList
}
global regionMemberList
global holdMunis
global districtMemberList
districtMemberList={}
regionMemberList={}
global groupMemberList
groupMemberList={}
global boundaryLengthDict
boundaryLengthDict=None

# Working constants computed here as necessary
global popTarget
popTarget=int(populationTotal/numDistricts)
global sigDigits
sigDigits=8   # significant digits in Lat and Long values
global crs_cea_proj_str,projType,projUnits
# projection settings
projType='+proj=cea '
projUnits='+units=km '
crs_cea_proj_str='epsg:3395'
global perimeterToMiles
perimeterToMiles=0.621371  # This adjustment must be consistent with the units choice in projection 'projType' above
global areaToSqMiles
areaToSqMiles=perimeterToMiles**2  # Approximate lat x long units in area to square miles in PA
 
# Special global variable set functions
def setGlobal(var,value):
    globals()[var]=deepcopy(value)
    return

def getGlobal(var):
    return globals()[var]
