import math
from statistics import mean,pstdev,median
from numpy import isnan
from random import sample
from pandas import read_json, DataFrame
from copy import deepcopy,copy
from shapely.geometry import Polygon
import Connectedness as Cnct

import Global_Settings as G
from File_Handlers import LoadList
from MemberList_Creators import calcDistrictBoundaries

import Region_Checking_Functions as CCFD_checks
import Shape_Handler as CCFD_shapes
import Goals_Tests as CCFD_tests

# Unit conversion functions
# Start: Source: John Cook Consulting (https://www.johndcook.com/blog/2009/04/27/converting-miles-to-degrees-longitude-or-latitude/)
earth_radius = 3960.0
degrees_to_radians = math.pi/180.0
# Converts a latitude change to its equivalent in miles
def deltaLAT_to_miles(deltaLAT):
    return (deltaLAT*degrees_to_radians)*earth_radius
# Converts a longitude change to its equivalent in miles
#    Requires a latitude as an argument
def deltaLONG_to_miles(deltaLONG,latitude):
    r=earth_radius*math.cos(latitude*degrees_to_radians)
    return (deltaLONG*degrees_to_radians)*r
# Converts a distance in miles to its latitude change equivalent
def change_in_latitude(miles):
    "Given a distance north, return the change in latitude."
    return (miles/earth_radius)*radians_to_degrees
# converts a distance and a latitude to its longitude change equivalent
def change_in_longitude(latitude, miles):
    "Given a latitude and a distance west, return the change in longitude."
    # Find the radius of a circle around the earth at given latitude.
    r = earth_radius*math.cos(latitude*degrees_to_radians)
    return (miles/r)*radians_to_degrees
# End: Source: John Cook Consulting

# Formula for the binomial coefficient (number of combinations of taking k items from a population of n, ignoring ordering)
def binom(n,k):
    if not 0<=k<=n: return 0
    b=1
    for t in range(min(k,n-k)):
        b*=n
        b/=t+1
        n-=1
    return int(b)
# End: Source: keithbriggs  https://gist.github.com/rougier/ebe734dcc6f4ff450abf

# inputs a string of items separated by a delimiter, returns a set of items
def neighborsList(nstring,delimiter='+'):
# No longer used. Employs pysal code instead.
# This function takes the string of neighbors from the QGIS input file and creates a set from them
# The input string has the form *nnn*+*nnn*+*nnn*... where * are spaces
# The input is a string with plus sign as separator
    xlist=set()
    i=0
    xstring=nstring+delimiter
    xstring=xstring.replace(' ','')
    len_in=len(xstring)
    while i<len_in:
        j=xstring.find(delimiter,i,len(xstring))
        xlist.add(xstring[i:j])
        i=j+1
    xlist.discard('')
    return xlist

# inputs a set or list, returns a single string with elements separated by the specified delimiter
def neighborsString(inList,delimiter=','):
# This function inserts a delimiter between elements of a string list.
# It returns a single string
    xstr=delimiter.join(inList)
    return xstr

# Creates an ID from some component pieces:
    # a prefix string
    # an input integer
    # the number of character places reserved in the ID for the integer. It is left-filled with zeros as necessary
        # if not specified, no left filling takes place
    # a suffix. If not specified, no suffix is added
# Typical area ID is, for example, Mcccmmmmm.0
def areaID(prefix,inNumber,numPlaces=None,suffix=None):
# This function creates an ID of the form Xnnn where X is a prefix
# The inputs are a character and a number
    if suffix is None: suffix=''
    if numPlaces is None: numPlaces=len(str(inNumber))
    xstr='0'*numPlaces+str(inNumber)
    xstr=prefix+xstr[len(xstr)-numPlaces:len(xstr)]+suffix
    return xstr

# the coreID is the ID before the decimal point.
def coreID(inRegion):
    xseq=inRegion.find('.')
    if xseq>0:
        return inRegion[:xseq]
    else:
        return inRegion

# returns a member list at the requested boundary level
def getInitialMemberList(inchar):
    if inchar=='S':
        memberList=G.stateMemberList

    elif inchar=='C':
        memberList=G.countyMemberList

    elif inchar=='M':
        memberList=G.municipalMemberList

    elif inchar=='V':
        if len(G.VTDMemberList)==0:
            results,G.VTDMemberList=LoadList('VTD')
        memberList=G.VTDMemberList
    else:
        return None

    return memberList

# returns the region code after the decimal point
def subRegSeq(inRegion):
    xseq=inRegion.find('.')
    if xseq>0:
        return inRegion[xseq+1:]
    else:
        return ''


# Builds a dictionary of all area IDs associated with their initial letters
def BuildExclusions(inList):
    # Exclusion precedence is 'C','M','W','V'     (County, Municipality, Ward, VTD)
    # The first character of the exclusion applies to that level of precedence and all subsequent ones.
    # (e.g., a 'C' applies to all regions; a 'M' applies to 'M','W','V' regions)
    # Note: the 'W' and 'V' regions are not yet implemented.
    precedence=['C','M','W','V']
    lenPrec=len(precedence)
    
    def getSubstr(invar,precList):
            if invar[0] in precedence:
                return [x+invar[1:] for x in precedence[precedence.index(invar[0]):]]
            else:
                return[invar]

    # Build the exclusion set
    exclSet={'Len':None,'NbrSet':set()}
    exclDict={}
    for xitem in inList:
        for xseq in range(0,2):
            yseq=(1,0)[xseq]

            for xvar in getSubstr(xitem[xseq],precedence):
                if xvar not in exclDict.keys(): exclDict[xvar]=deepcopy(exclSet)
                exclDict[xvar]['Len']=len(xvar)
                for yvar in getSubstr(xitem[yseq],precedence):
                    exclDict[xvar]['NbrSet'].update({yvar})

    return exclDict

# Removes any specifically identified regions from inDict, and reallocate any relevant demo data
# inTuple: (excluded Region, region receiving data)
def forcedRegionExclusions(inTupleList,inDict):
    inRegionList=deepcopy(inDict)
    for excl in inTupleList:
        if excl[0] in inRegionList.keys():
            if excl[1] not in inRegionList.keys():
                print('Error: receiving region not found:',excl[1])
                return None
            print(' ',excl[0],'. Moving demog data to ',excl[1])
            exclRegion=inRegionList[excl[0]]
            # Move the demog data to the receiving region
            for demogData in G.regionDemogData:
                inRegionList[excl[1]]['Data'][demogData]+=exclRegion['Data'][demogData]
            # Changes to references to the excluded region
            for nbr in exclRegion['Neighbors']-{'Outside'}:
                # The excluded region will become part of 'Outside' (may be a hole)
                inRegionList[nbr]['Neighbors'].discard(excl[0])
                inRegionList[nbr]['Neighbors'].update({'Outside'})
                # The boundary between the neighbor and the excluded region will now be associated with 'Outside'  
                if 'Outside' not in G.boundaryLengthDict[nbr].keys():
                    G.boundaryLengthDict[nbr]['Outside']=G.boundaryLengthDict[nbr][excl[0]]
                else:
                    G.boundaryLengthDict[nbr]['Outside']+=G.boundaryLengthDict[excl[0]]['Outside']
#                del G.boundaryLengthDict[nbr][excl[0]]
            del inRegionList[excl[0]]

            # Finally, if there is just one remaining subregion (alpha prior to decimal point), remove it
            if coreID(excl[1])[-1] in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                testCore=coreID(excl[1])[:-1]
                subregions=[x for x in inRegionList.keys() if x[:len(testCore)]==testCore]
                if len(subregions)==1:
                    newCode=testCore+'.'+str(subRegSeq(excl[1]))
                    print('  One subregion remains. Rename',excl[1],' to ',newCode)
                    # Rename the subregion code with out subregion character:
                    inRegionList[newCode]=deepcopy(inRegionList[excl[1]])
                    del inRegionList[excl[1]]
                    inRegionList[newCode]['ID']=newCode
                    inRegionList[newCode]['Data']['ID']=newCode
                    inRegionList[newCode]['Group']={newCode}
                    # Rename it in the surrounding regions' Neighbors sets.
                    for nbr in inRegionList[newCode]['Neighbors']-{'Outside'}:
                        inRegionList[nbr]['Neighbors'].discard(excl[1])
                        inRegionList[nbr]['Neighbors'].update({newCode})
                    # Fix the boundary length file. Involves using the new name in the list and its transposes
                    G.boundaryLengthDict[newCode]=deepcopy(G.boundaryLengthDict.pop(excl[1]))
                    for btranspose in set(G.boundaryLengthDict[newCode].keys())-{'Outside'}:
                        G.boundaryLengthDict[btranspose][newCode]=G.boundaryLengthDict[newCode][btranspose]

        else:
            if G.verbose: print('  Region not found. No exclusion: ',excl[0])

    return inRegionList

# Remove any specifically identified neighbor exclusions
# intupe contains tuples of regions for which we state are not neighbors of each other
def forcedExclusions(inExclusions,inDict):

    for xreg in inExclusions.keys():

        if xreg in ['xxx']:
            print('...')

        xlen=inExclusions[xreg]['Len']
        nbrSet=inExclusions[xreg]['NbrSet']
        regionSet={x for x in inDict.keys() if x[:xlen]==xreg}
        for thisreg in regionSet:
            for thisnbr in copy(inDict[thisreg]['Neighbors']):
                if thisnbr[:xlen] in nbrSet:
                    print('  Disconnecting:',thisreg, thisnbr)
                    inDict[thisreg]['Neighbors'].discard(thisnbr)
                    inDict[thisnbr]['Neighbors'].discard(thisreg)
                
    return inDict

# Adds specifically required neighbors
# input contains tuples of regions for which we make neighbors of each other
def forcedInclusions(intuple,inRegionList):
    for incl in intuple:
        if incl[0] in inRegionList.keys() and incl[1] in inRegionList:
            if incl[0] not in inRegionList[incl[1]]['Neighbors'] or incl[1] not in inRegionList[incl[0]]['Neighbors']:
                print('  Connection created between:',incl[0],'and',incl[1])
                inRegionList[incl[0]]['Neighbors'].update({incl[1]})
                inRegionList[incl[1]]['Neighbors'].update({incl[0]})
            else:
                print('  Requested connection already exists. No action:',incl[0],incl[1])
    return inRegionList


# Constructs a districtMemberList from a geodataframe
def build_districtMemberList(inFrame):
    # Rebuild districtMemberList from a saved dataFrame

    regionData=G.regionData
    districtMemberSet=G.districtMemberSet
    outList={}

    for i in range(0,len(inFrame['RegionIDs'])):

        if i in {'S043'} and G.verbose:
            print('Pt build_districtMemberList 1')

        thisDistrict=inFrame['RegionIDs'][i]
        outList.update({thisDistrict:deepcopy(districtMemberSet)})

        outList[thisDistrict]['ID']=inFrame['ID'][i]
        outList[thisDistrict]['geometry']=inFrame['geometry'][i]

        for xd in regionData:
            outList[thisDistrict]['Data'][xd]=inFrame[xd][i]
            if xd=='TAPersons' and not isnan(inFrame[xd][i]):
                outList[thisDistrict]['Data'][xd]=int(inFrame[xd][i])

        xNeighbors={name for name in neighborsList(inFrame['Neighbors'][i])}
        outList[thisDistrict]['Neighbors'].update(xNeighbors)

        xNeighbors={name for name in neighborsList(inFrame['distNeighbors'][i])}
        outList[thisDistrict]['distNeighbors'].update(xNeighbors)
        

        xNeighbors={name for name in neighborsList(inFrame['Regions'][i])}
        outList[thisDistrict]['Regions'].update(xNeighbors)
    return outList

# iterates to the next region in sequence as determined by the integer following the decimal point
def region_and_Number(inRegion):
    return inRegion[:inRegion.find('.')],int(inRegion[inRegion.find('.')+1:])

# Specialized function -- creates a new regionSet based on the name and ID of a parent set
# Iterates the suffix digit of the parent ID by one, adds modifier to the parent name
def MakeSubRegionSet(inRegion,regionMemberList=None):
    if regionMemberList is None: regionMemberList=G.regionMemberList

    xdot=inRegion.find('.')
    # To add a new subregion, we need to find the highest existing subregion number for this region
    xnum=1+max([int(x[xdot+1:]) for x in regionMemberList.keys() if x[:xdot]==inRegion[:xdot]])
    # Add the new region
    newRegion=inRegion[:xdot+1]+str(xnum)
    regionSet=deepcopy(G.regionMemberSet)
    regionSet['ID']=newRegion
    regionSet['Name']=regionMemberList[inRegion]['Name']+' Part'+str(xnum)
    regionSet['Districts'][G.districtPrefix]=copy(regionMemberList[inRegion]['Districts'][G.districtPrefix])
    print('  New region and set created: ',newRegion)

    return regionSet

# Routine to assign colors to geographical regions so that no neighbor has the same color.
def colorize(inList,neighborStr,colorCount):

    colorNumbers=[x for x in range(0,colorCount)]

    # Make a list of the districts
    colorDict={x:None for x in inList.keys()}

    loops=0
    # Run through the districts, assigning colors randomly from among the colors not used for neighbors
    while True:
        for xdist in colorDict.keys():

            reDo=False    
            if colorDict[xdist]==None:
                usedColors=set([colorDict[x] for x in inList[xdist][neighborStr]])
                availableColors=set(colorNumbers)-usedColors
                if len(availableColors)>0:
                    xcolor=sample(availableColors,1)[0]
                    colorDict[xdist]=xcolor
                else:
                    # We've encountered a member whose neighbors already have colors. Set its color and repeat process.
                    # This should not happen very often.
                    colorDict={x:None for x in inList.keys()}  # Set this for debugging purposes.
                    xcolor=sample(colorNumbers,1)[0]
                    colorDict[xdist]=xcolor
                    reDo=True
                    break # Breaks from the for loop
        loops=loops+1
        if loops>5 or not reDo:
            break # Breaks from the while loop

    if reDo and G.verbose: print('Coloring problem for: ',xdist)

    return colorDict

def colorize2(inDict,numColors=6):
    colorNumbers={x for x in range(numColors)}

    # Make a list of the districts
    colorDict={x:None for x in inDict.keys()}

    # Run through the districts, assigning colors randomly from among the colors not used for neighbors
    # This way, we always use all colors (or color all districts differently)
    unAssigned=set(inDict.keys())

    # Run through the first numColors of inDict members, assigning a number to each
    xnum=0
    sortDict=list(inDict.keys())
    sortDict.sort()
    for x in sortDict[:numColors]:
        colorDict[x]=xnum
        unAssigned-={x}
        xnum+=1

    # Now assign colors to the rest of the districts
    while len(unAssigned)>0:
        thisRegion=unAssigned.pop()
        if colorDict[thisRegion] is None:
            usedColors={colorDict[x] for x in inDict[thisRegion] if x is not None}
            availableColors=colorNumbers-usedColors
            while len(availableColors)==0:
                # release the color of one of the neighbors
                pickOne=sample(inDict[thisRegion]-{thisRegion},1)[0]
                colorDict[pickOne]=None
                unAssigned.update({thisRegion,pickOne})
                usedColors={colorDict[x] for x in inDict[thisRegion]}
                availableColors=colorNumbers-usedColors

            xcolor=sample(availableColors,1)[0]
            colorDict[thisRegion]=xcolor

    return colorDict

# Creates a sublist of regions in inMemberList for whom the 'Outside' region is a neighbor
def createOutsideList(inMemberList):
    outSideList=[]
    for xreg in inMemberList.keys():
        if 'Outside' in inMemberList[xreg]['Neighbors']:
            outSideList.append(xreg)
    return outSideList

# Returns a dictionary of inlist with the 'Outside' region stripped from all 'Neighbors' sets
def removeOutsideList(inList):
    inMemberList=deepcopy(inList)
    for x in inMemberList.keys():
        inMemberList[x]['Neighbors'].discard('Outside')
    if 'Outside' in inMemberList.keys():
        del inMemberList['Outside']
    
    return inMemberList

# Aggregate inregionDataSet data into indistrictDataSet data
# Some data values are summed, others involve various mathematical operations
# the 'Perimeter' data value involves the use of dictionary G.boundaryLengthDict
def aggregate_Data(indistrictDataSet,regionDataSet,varlist=None):

    districtDataSet=deepcopy(indistrictDataSet)

    def First(var,xin,yin):
        return xin[var]

    def Sum(var,xin,yin):
        return xin[var]+yin[var]

    def Max(var,xin,yin):
        return max(xin[var],yin[var])

    def Min(var,xin,yin):
        return min(xin[var],yin[var])

    def Perimeter(var,xin=0.0,yin=0.0,shared=None):
        if shared is None:
            try:
                shared=G.boundaryLengthDict[yin['ID']][xin['ID']]*G.perimeterToMiles   # boundaryLengths measured in natural lat-long units 
            except:
                try:
                    shared=G.boundaryLengthDict[xin['ID']][yin['ID']]*G.perimeterToMiles
                except:
                    shared=0.0
        return xin[var]+yin[var]-2.0*shared

    def Coordinate(var,xin,yin):
        return (xin[var]*xin['Area']+yin[var]*yin['Area'])/(xin['Area']+yin['Area'])

    def BdgSqrIndex(var,xin,yin):
        return min(1.0,(xin['Area']+yin['Area'])/bounding_square(max(xin['North'],yin['North']),min(xin['South'],yin['South']),max(xin['East'],yin['East']),min(xin['West'],yin['West'])))

    aggr_Type={
        'ID':First,
        'TAPersons':Sum,
        'Area':Sum,
        'Perimeter':Perimeter,
        'North':Max,
        'South':Min,
        'East':Max,
        'West':Min,
        'Latitude':Coordinate,
        'Longitude':Coordinate,
        'BdgSqrIndex':BdgSqrIndex
        }

    if varlist is None:
        varlist=G.regionData
    for xdata in varlist:
        func=aggr_Type.get(xdata,'nothing')
        if districtDataSet[xdata] is None:
            districtDataSet[xdata]=regionDataSet[xdata]
        elif xdata not in aggr_Type:
            print('No aggregate_data action on ',xdata,'Unknown data ID.')
        else:
            if districtDataSet[xdata] is None:
               districtDataSet[xdata]=0.0
               print('Missing data for District set to zero',xdata)
            if regionDataSet[xdata] is None:
                regionDataSet[xdata]=0.0
                print('Missing data for region set to zero',xdata)
            districtDataSet[xdata]=func(xdata,indistrictDataSet,regionDataSet)
        
    return districtDataSet

# Creates a mew districtMemberSet
def makeNewDistrict(distID):
    thisMemberSet=deepcopy(G.districtMemberSet)
    thisMemberSet['ID']=distID
    thisMemberSet['Data']['ID']=thisMemberSet['ID']
    thisMemberSet['geometry']=None

    G.boundaryLengthDict[distID]={}
    return thisMemberSet

# Returns a basic member list. Loads it if necessary
def confirmListAtLevel(thisLevel):
    precedenceTuple=G.regionPrecedence[thisLevel]
    thisPrecedence=precedenceTuple[0]
    if not G.memberListOpenFlags[thisPrecedence]:
    # Load the start member list
        results,memberList=LoadList(thisPrecedence)
        if not results:
            print('Terminating. Bad member list for',thisPrecedence)
            sys.exit(0)
        G.memberListOpenFlags[thisPrecedence]=True
        G.setGlobal(thisPrecedence+'MemberList',memberList)
    else:
        memberList=G.getGlobal(thisPrecedence+'MemberList')
    return memberList

# Constructs district data from a list of regions
# if the regionMemberList includes 'geometry' elements, uses Shapely as needed for district data entries
def build_Data(districtMemberSet,regionMemberList,varlist=None,useGeometry=True):
    if varlist is None:
        varlist=districtMemberSet['Data'].keys()

    reglist=list(districtMemberSet['Regions'])

    # Initialize the dataSet

    districtMemberSet['Data']={x:0.0 for x in G.regionData}
    for x in G.regionDemogData:
        districtMemberSet['Data'][x]=0
    districtMemberSet['Data']['ID']=districtMemberSet['ID']

    if len(reglist)>0:
        for xb in reglist:
            districtMemberSet['Data']=aggregate_Data(districtMemberSet['Data'],regionMemberList[xb]['Data'],varlist)

        if useGeometry:
            districtMemberSet['geometry']=CCFD_shapes.tryUnion(reglist,regionMemberList)
            districtMemberSet=CCFD_shapes.calcRegionGeoData(districtMemberSet)
    else:
        if G.verbose: print('Warning: district has no regional members:',districtMemberSet['ID'])

    return districtMemberSet['Data']

# Constructs a full or partial list of district Neighbors from a regionMemberList
def build_distNeighbors(districtMemberSet,regionMemberList):
    distNeighbors=set()
    districtPrefix=districtMemberSet['ID'][0]
    for x in districtMemberSet['Regions']:
        for y in regionMemberList[x]['Neighbors']:
            if y!='Outside':
                if districtPrefix in regionMemberList[y]['Districts'].keys():
                    if regionMemberList[y]['Districts'][districtPrefix] is not None:
                        distNeighbors.update({regionMemberList[y]['Districts'][districtPrefix]})
            else:
                distNeighbors.update({'Outside'})
    distNeighbors.discard(districtMemberSet['ID'])

    return distNeighbors

# Takes as inputs 4 extremum points on the compass as lat/long, and computes the appproximate area of a minimum bounding square.
def bounding_square(north,south,east,west):
    if north is None or south is None or east is None or west is None:
        print('Error: bounds value(s): N',north,' S',south,' E',east,' W',west)
        return 0.0

    # Compute E-W and N-S 
    del_lat = north-south
    del_lng = east-west
    r=7917.5/2.0  # radius of the earth in miles (12742/2 km)

    # Convert these angles to radians
    north,south=map(math.radians,[north,south])
    del_lat,del_lng=map(math.radians,[del_lat,del_lng])

    # N-S distance (with d_lng = 0)
    NSdist=r*del_lat

    # E-W distance (d_lat = 0, more complex because it's not on a great circle
    # Using a formula for small fractions of the great circle.
    # Using the southern extremum as latitude.
    EWdist=r*math.asin(math.cos(south)*math.sin(del_lng))
    

    # Use the larger to compute the area
    return max(abs(NSdist),abs(EWdist))**2.0

# Returns distance in miles between two points given by lat,long tuples
def get_distance(point1,point2): 
    d_lat = point2[0] - point1[0]
    d_lng = point2[1] - point[1]
    r=7917.5/2.0  # radius of the earth in miles (12742/2 km)
  
    temp = (    
         math.sin(d_lat / 2) ** 2 
       + math.cos(lat_1) 
       * math.cos(lat_2) 
       * math.sin(d_lng / 2) ** 2
    )

    return r * (2 * math.atan2(math.sqrt(temp), math.sqrt(1 - temp)))

# Prints two population stats for a memberList with a 'Data' dictionary: a max-min spread and RMS spread wrt the global population target
def populationSpread(inMemberList):
    popList=[inMemberList[x]['Data']['TAPersons'] for x in inMemberList.keys()]
    minPop=min(popList)
    maxPop=max(popList)
    popSpread=math.sqrt(sum((x/G.popTarget)**2 for x in popList))
    print('District population pct vs Target: max-min, RMS spread:',round((maxPop-minPop)/G.popTarget*100,4),round(popSpread*100,6))
    return

# Returns multiple stats for an input memberList with a 'Data' dictionary
def indexSpread(metric,memberList,targetDict=None):
    if targetDict is None:
        targetDict={x:0.0 for x in districtMemberList.keys()}
    indexValues=[CCFD_tests.metricCalc(metric,memberList[x]['Data'],targetDict[x]) for x in memberList.keys()]
    minIndex=min(indexValues)
    maxIndex=max(indexValues)
    medIndex=median(indexValues)
    meanIndex=mean(indexValues)
    pstdIndex=pstdev(indexValues)
    print('Stats for ',metric,' index mean:',round(meanIndex,4),', max-min from target:',round((maxIndex-minIndex)/meanIndex,4),', RMS Spread to target:',round(pstdIndex/meanIndex,4))
    return {'min':minIndex,'max':maxIndex,'med':medIndex,'mean':meanIndex,'pstd':pstdIndex}


# This function moves a region "inreg" from its current district to a new district "newDist"
# It makes all necessary changes to lists "regionMemberList" and "districtMemberList"
# It returns "True" if all the changes succeed, "False" if not.
# It prints some notes associated with several specific reasons for failure.
# The variable "singleRegionOverride" allows the function to remove the last region from a district list and delays computation of the empty district data
    # In context, this is valid only just before another invocation of this function that moves a new region into it.
def MoveRegionTo(inregSet,newDist,singleRegionOverride,districtMemberList,regionMemberList,useGeometry=True):

    # Create a set of regions that will be moved as a group.
    # If input is a string, assume it is a single region ID
    if isinstance(inregSet,str):
        inreg={inregSet}
    else:
        inreg=inregSet

    # Initial working info
    districtPrefix=G.districtPrefix
    startList=deepcopy(districtMemberList)

    # Check to make sure that each of the regions has data
    if not inreg.issubset(set(regionMemberList.keys())):
        missingIDs-=set(regionMemberList.keys())
        print('Incorrect region IDs: ',missingIDs)
        return False,startList

    # Check to make sure the district ID is valid.                    
    if newDist not in districtMemberList.keys():
        print('Incorrect district ID: ',newDist)
        return False,startList

    # For test purposes compute the total population for the subset of districts in the current collection
    if len(districtMemberList.keys())>0:
        testPop=[x for x in districtMemberList.keys() if districtMemberList[x]['Data']['TAPersons'] is not None]
        if len(testPop)==0:
            print('Stop: no population entries.')
        controlPopulation=sum([districtMemberList[x]['Data']['TAPersons'] for x in testPop])
    else:
        controlPopulation=0

    # Determine beginning region membership
    # The input regions must all reside in the same district
    oldDist=set()
    for xreg in inreg:
        tryDist=regionMemberList[xreg]['Districts'][districtPrefix]
        if tryDist is None:
            print('Error: no district in record for: ',xreg)  # Interrupt for debugging.
            return False,startList
        oldDist.update({tryDist})
    if len(oldDist)>1:
        print('Error: input regions reside in multiple districts',{x:regionMemberList[x]['Districts']['districtPrefix'] for x in inreg})
        return False,startList
    elif len(oldDist)==0:
        print('Error: no districts identified in regions',{x:regionMemberList[x]['Districts']['districtPrefix'] for x in inreg})
    else:
        oldDist=list(oldDist)[0]


    # check to make sure the regions are currently in the old District regions
    if  not inreg.issubset(districtMemberList[oldDist]['Regions']):
        print('Error: some of region set ',inreg,' not in district ',oldDist,':',districtMemberList[oldDist]['Regions'])
        print()
        return False,startList

    # Under normal conditions, emptying a district is disallowed. There is an override capability, however.
    if len(districtMemberList[oldDist]['Regions'])==len(inreg):
        if singleRegionOverride:
            if G.verbose: print('   ',oldDist,'; ',districtMemberList[oldDist]['Regions'],' (Override: ',singleRegionOverride,')')
        else:
            print('Error: Cannot remove all members of a district.')
            return False, startList

    # Update the region assignments
    for xreg in inreg:
        regionMemberList[xreg]['Districts'][districtPrefix]=newDist

    # Update the district memberships
    districtMemberList[newDist]['Regions'].update(inreg)
    districtMemberList[oldDist]['Regions']-=inreg

    # Prepare to recalculate the district neighborhood and its Data values
    # If we're allowing non-compliant move, only address the receiving district
    if singleRegionOverride and districtMemberList[oldDist]['Regions']==inreg:
        # We're creating an empty district. Reset its info.
        districtMemberList[oldDist]=makeNewDistrict(oldDist)
        G.boundaryLengthDict[oldDist]={}
        fixList=[newDist]
    else:
        fixList=[newDist,oldDist]

    # Update each district's info to include the updated membership
    for bx in fixList:

        # Regional neighbors
        districtMemberList[bx]['Neighbors']=set()
        for reg in districtMemberList[bx]['Regions']:
            districtMemberList[bx]['Neighbors'].update(regionMemberList[reg]['Neighbors'])
        districtMemberList[bx]['Neighbors']-=districtMemberList[bx]['Regions']

        # Meighboring districts
        districtMemberList[bx]['distNeighbors']=set()
        for thisRegion in districtMemberList[bx]['Neighbors']:
            if thisRegion=='Outside':
                districtMemberList[bx]['distNeighbors'].update({'Outside'})
            else:
                districtMemberList[bx]['distNeighbors'].update({regionMemberList[thisRegion]['Districts'][districtPrefix]})
        # Remove self-referenced district
        districtMemberList[bx]['distNeighbors'].discard(bx)

        # Compute the data
        districtMemberList[bx]['Data']=build_Data(districtMemberList[bx],regionMemberList,varlist=G.regionData)

        # recalculate the district/region boundary lengths for the newly modified district
        results=calcDistrictBoundaries(districtMemberList[bx],regionMemberList)
    
    # Perform post-move checks.

    # test overall population. should not change.
    if not CCFD_checks.CheckPopTotal(districtMemberList,controlPopulation,verbose=True):
        print('    Population test failed.')
        results=CCFD_checks.CheckRegionDistrictConsistency(districtMemberList,regionMemberList)        
        return False,startList
    # Test Connectedness
    if not CCFD_checks.CheckConnectedness({newDist:districtMemberList[newDist]},regionMemberList,verbose=True):
        print('  Connectedness issue:',newDist,districtMemberList[newDist]['Regions'])
    # check district boundaries
    results=CCFD_checks.checkRegionBoundaries({newDist:districtMemberList[newDist]},verbose=False)
    if not results:
        CCFD_checks.boundaryDecomposer(districtMemberList[newDist],regionMemberList)
        print()

    return True,districtMemberList


# This function returns lists of connected elements
# Inputs: inRegionList -- a set of region codes
#         regionNeighborsDict: a dictionary of region:{neighbors}
# Returns: a list of lists of connected regions
def getConnected(inRegionList,regionNeighborsDict):

    def numConverters(inList):
        outDict={}
        backDict={}
        seq=0
        for x in inList:
            if x not in outDict.keys():
                outDict[x]=seq
                backDict[seq]=x
                seq+=1
        return outDict,backDict
    def numberToName(inList,inDict):
        outList=[]
        for x in inList:
            outList.append(inDict[x])
        return outList

    # Alternative graph code
    nameDict,numDict=numConverters(inRegionList)
    nodes=len(nameDict)
    numNeighborsDict={}
    numList=[x for x in range(nodes)]
    for x in inRegionList:
        thisSeq=nameDict[x]
        numNeighborsDict[thisSeq]=[]
        if x in regionNeighborsDict.keys():
            for y in regionNeighborsDict[x]:
                if y in inRegionList:
                    numNeighborsDict[thisSeq].append(nameDict[y]) 

    g=Cnct.Graph(nodes)
    for x in range(nodes):
        for y in range(x+1,nodes):
            if y in numNeighborsDict[x]:
                g.addEdge(x,y)
    cc=g.connectedComponents()

    outLists=[]
    for x in cc:
        outLists.append(numberToName(x,numDict))

    return outLists


# This function returns True if "inDistrict" is connected. False otherwise.
# Empty inGroup is defined as connected.
def testConnectedness(inGroup,regionMemberList):
    if len(inGroup)==0:
        return True
    xneighbors={x:regionMemberList[x]['Neighbors'] for x in regionMemberList.keys()}
    cc=getConnected(inGroup,xneighbors)
    return len(cc)==1

