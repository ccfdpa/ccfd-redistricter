# File handling functions
import Global_Settings as G
import CCFD_Functions as CCFD
import File_Handlers as CCFD_files
import Shape_Handler as CCFD_shapes
import subRegion_Handler as CCFD_subregions
import Region_Checking_Functions as CCFD_checks

import sys
import pandas as pd
import geopandas as gp
from numpy import isnan,nan
from random import sample
from pathlib import Path

from shapely import geometry
from shapely.geometry import Point,Polygon
from copy import deepcopy, copy

global districtMemberList
global regionMemberList

def check_file_open(file_path):
# From stackoverflow contributor Razbi: March 12, 2021

    path = Path(file_path)
    
    if not path.exists():
        return False
    
    try:
        path.rename(path)
    except PermissionError:
        print('File is in use.')
        return True
    else:
        return False

def build_regionMemberList(inFrame):
    # Rebuild regionMemberList

    outList={}
    vars=list(inFrame.columns)
    rcount=0

    for i in range(len(inFrame['RegionIDs'])):

        regSet=inFrame.iloc[i,:]
        thisRegion=regSet.loc['RegionIDs']
        outList[thisRegion]={'Data':{}}

        for var in vars:

            if var in ['ID','Name','geometry']:
                # Var is a single entity
                outList[thisRegion][var]=regSet.loc[var]

            elif var in ['Neighbors','Group','Regions','distNeighbors']:
                # Var is a set
                outList[thisRegion][var]=CCFD.neighborsList(regSet.loc[var])

            elif var[:8]=='District':
                outList[thisRegion]['Districts']={}
                xpos=var.find('_')+1
                if xpos<0:
                    # Legacy districts form
                    # Convert to new form
                    distSet=CCFD.neighborsList(regSet.loc[var])
                    outList[thisRegion][var]={x[0]:x for x in distSet}
                else:
                    # Var is a dictionary with the key following an underline character in var
                    xprefix=var[xpos:]
                    outList[thisRegion]['Districts'].update({xprefix:regSet.loc[var]})

            elif var in G.regionData:
                value=regSet.loc[var]
                if isnan(value): value=0.0
                if var in G.regionDemogData:
                    # For integer values
                    outList[thisRegion]['Data'][var]=int(value)
                else:
                    # For floating point values
                    outList[thisRegion]['Data'][var]=value

        outList[thisRegion]['Data']['ID']=outList[thisRegion]['ID']
        rcount+=1
        if rcount>500:
            print('.',end='')
            rcount=0
    print('.')
    return outList

def JsonToList(thisPath,returnFrame=False):

    print('Loading from JSON: ',thisPath,' .',end='')
    
    try:
        xFrame=gp.read_file(thisPath+'.json')
        if returnFrame:
            return xFrame
    except:
        print('  Error: file not found.')
        return None

    xlist=build_regionMemberList(xFrame)

    return xlist

def fileInLists(thisPath,runID=0):
    # Read in the info for regionMemberList and districtMemberList
    regionMemberList=JsonToList('region',thisPath,runID)
    if regionMemberList is not None: print('  regionMemberList reconstructed.')

    # Rebuild districtMemberList
    districtMemberList=JsonToList('district',thisPath,runID)
    if districtMemberList is not None: print('  districtMemberList reconstructed.')

    # We're finished with the input phase.
    print('Input phase completed.')
    return True

# Create a new dataframe of the regionMemberList with key information for transfer purposes
def build_regionFrame(inList):
    regionList=list(inList.keys())

    if not isinstance(inList[regionList[0]],dict):                               # the input list is not a memberList
        regionFrame=pd.DataFrame({'ID':list(inList.keys()),'Value':[inList[x] for x in inList.keys()]})
        return regionFrame

    dataList=list(inList[regionList[0]].keys())
    regionFrame=gp.GeoDataFrame({'RegionIDs':regionList})

    for xdata in dataList:

        if xdata in ['Neighbors','distNeighbors','Regions','Group']:      # to create a string of set members
            if isinstance(inList[regionList[0]][xdata],dict):
                toList={reg:[x[0]+':'+inList[reg][xdata][x] for x in inList[reg][xdata].keys()]for reg in regionList}
            else:
                toList={x:list(inList[x][xdata]) for x in regionList}
            regionFrame[xdata]=[CCFD.neighborsString(toList[reg],'+') for reg in regionList]
        elif xdata in ['Districts']:                                       # to create a districts dictionary
            for xd in inList[regionList[0]][xdata].keys():
                regionFrame['District_'+xd]=[inList[val][xdata][xd] for val in regionList]                
        elif xdata in ['Data']:                                             # to create a data dictionary
            for xd in inList[regionList[0]][xdata].keys():          # to split subdictionary
                regionFrame[xd]=[inList[val][xdata][xd] for val in regionList]
                regionFrame[xd].fillna(nan,inplace=True)

        elif xdata=='geometry':                                             # to create the geometry element
            regionFrame['geometry']=[inList[val]['geometry'] for val in regionList]


        else:                                                              # every other element in a memberList
            regionFrame[xdata]=[inList[val][xdata] for val in regionList]


    print('')
    
    return regionFrame

def fileOutCSV(inList,thisPath=None,runID='QGIS CSV Transfer'):
    if thisPath is None: thisPath=G.thisPath

#    try:
    regionFrame=build_regionFrame(inList)
    fileCSV=thisPath+runID+'.csv'
    regionFrame.to_csv(fileCSV,index=False)
    print(fileCSV+' created.')
    return True
#    except:
#        print('Failed: save to CSV')

#    return False

def fileOutSHP(inList,thisPath=None,runFileName='QGIS SHP Transfer',distColors=True,districtPrefix=None,regionMemberList=G.regionMemberList):
    def DistrictsFromRegions(inList):
        districts={}
        for region in inList.keys():
            thisDistrict=inList[region]['Districts'][districtPrefix]
            if thisDistrict not in districts.keys():
                districts[thisDistrict]=set()
            for trynbr in inList[region]['Neighbors']:
                if trynbr in inList.keys():
                    tryDistrict=inList[trynbr]['Districts'][districtPrefix]
                    if tryDistrict!=thisDistrict:
                        districts[thisDistrict].update({tryDistrict})
        return districts

    if districtPrefix is None: districtPrefix=G.districtPrefix
    if thisPath is None: thisPath=G.thisPath
    results=True

    # Geometries need to be constructed if this List consists of districts
    if 'Regions' in inList.items()[0].keys():
        for thisReg in inList.keys():
            inList[thisReg]['geometry']=unary_union([regionMemberList[xreg]['geometry'] for xreg in inList[thisReg]['Regions']])
            if not checkPolygon(inList[thisReg]['geometry']):
                print('  Warning: geometry not a polygon:',thisReg)
    try:
        regionFrame=build_regionFrame(inList)
    except:
        failStep='regionFrame creation failed.'
        results=False

    if distColors and results:
        try:
            districtList=DistrictsFromRegions(inList)

            distColors=CCFD.colorize2(districtList,6)
            regionFrame['Colorize']=[distColors[inList[val]['Districts'][districtPrefix]] for val in regionFrame['RegionIDs']]
        except:
            failStep='colorcode setting failed.'
            results=False

    if results:
        region_ID=thisPath+runFileName+".shp"
        regionFrame.to_file(region_ID,driver='ESRI Shapefile')
        print('Shapefile: ',region_ID, 'created.')
    else:
        print('Shapefile creation failed...',failStep)

    return results

def dictToJson(inDict,thisPath,runID='',inType=None):

    try:
        regionFrame=build_regionFrame(inDict)
    except:
        print('GeoDataFrame not created.')
        return False
    try:
        Export=regionFrame.to_file(thisPath+str(runID)+'.json',driver='GeoJSON')
        print('Filed: ',thisPath+runID+'.json')        
        return True
    except:
        print('Error saving',thisPath+str(runID)+'.json')
        check_file_open(thisPath+str(runID)+'.json')
    return False

def fileOutResults(thisPath,runID,toCSV=False):
    print('Starting from File_Handlers fileOutResults ...')

    regionMemberList=G.regionMemberList
    districtMemberList=G.districtMemberList

    regionFrame_ID='regionFrame '+str(runID)
    print('Saving from regionMemberList')
    results1=dictToJson(regionMemberList,thisPath,regionFrame_ID)

    districtFrame_ID='districtFrame '+str(runID)
    print('Saving from districtMemberList')
    results2=dictToJson(districtMemberList,thisPath,districtFrame_ID)

    # Create a csv file for transfer region info back to QGIS
    if toCSV:
        distColors=CCFD.colorize(districtMemberList,'distNeighbors',6)

        # For now, any region broken into multiple subregions gets a color of 6 (grey)
        xtry={val:distColors[list(regionMemberList[val]['Districts'])[0]] for val in regionFrame['Regions']}
        for val in xtry.keys():
            if int(val[5])>0:
                xtry[val]=6
                xtry[val[0:5]+'0']=6

        regionFrame['Colorize']=[xtry[val] for val in regionFrame['RegionIDs']]

        regionCSV_ID='Senate Region Assignments '+str(runID)
        regionFrame.to_csv(thisPath+regionCSV_ID+'.csv',index=False)
        print(thisPath+regionCSV_ID+'.csv'+' created.')

        districtCSV_ID='Senate District Assignments '+str(runID)
        districtFrame.to_csv(thisPath+districtCSV_ID+'.csv',index=False)
        print(thisPath+districtCSV_ID+'.csv'+' created.')

    return results1 and results2

def build_statsFrame():
    statsFrame=pd.DataFrame({'Stats':[]})
    statsFrame.set_index('Stats',inplace=True)
    return statsFrame

def add_statsToFrame(resultsFrame,inStats,runID):

    if type(runID)!=str: runID=str(runID)
    runNUM='Run'+runID
    if runNUM not in resultsFrame.columns: resultsFrame[runNUM]=None

    for x in list(inStats.keys()):
        resultsFrame.loc[x,runNUM]=inStats[x]
    
    return resultsFrame

def fixIslands(inMemberList):
            # Special fix for municipality neighbors using inefficient method
    islandList=[x for x in inMemberList.keys() if len(inMemberList[x]['Neighbors'])==0]
    for id in islandList:
        print('Testing: ',id)
        thisGeo=inMemberList[id]['geometry']
        thisCounty=inMemberList[id]['ID'][1:4]
        testDict={x:inMemberList[x]['geometry'] for x in inMemberList if inMemberList[x]['ID'][1:4]==thisCounty}
        for cty in inMemberList['C'+thisCounty+'.0']['Neighbors']:
            xcode=cty[1:4]
            testDict.update({x:inMemberList[x]['geometry'] for x in inMemberList if G.municipalMemberList[x]['ID'][1:4]==xcode})        
        testDict.pop(id,None)
        for x in testDict.keys():
            # Special test for specific muni here. For debugging only.
            if id=='Mcccmmmmm.0':   # Coding is 3 digit county FIPS followed by 5 digit municipal FIPS
                print('Checking: ',id,x)
            if IsNeighbor(thisGeo,testDict[x]):
                print('Adding neighbors ',id,x)
                inMemberList[id]['Neighbors'].update({x})
                inMemberList[x]['Neighbors'].update({id})
    
    return inMemberList

def fixMunicipalNeighbors(inCountyFips,countyMemberList=None,municipalMemberList=None):
    if countyMemberList is None: countyMemberList=G.countyMemberList
    if municipalMemberList is None: municipalMemberList=G.municipalMemberList

    # This function assumes that the neighbors for the counties is correct.
    # Fixes two possible issues with the municipal list:
    # 1) The initial list must be fully connected.
    # 2) It must have a full set of municipal neighbors in its own county and the surrounding counties


    thisMuniList={x:municipalMemberList[x]['Neighbors'] for x in municipalMemberList.keys() if x[:4]=='M'+inCountyFips}

    # Special fix for '003'
    if inCountyFips=='003':
        for x in thisMuniList.keys():
            municipalMemberList[x]['Neighbors'].discard('Outside')

    # check to make sure the county's own munis are fully connected.
    this_Connected=CCFD.getConnected(list(thisMuniList.keys()),thisMuniList)
    if len(this_Connected)>1:
        # Use fixNeighbors to try to resolve them.
        xlist=list(thisMuniList.keys())
        counter=(len(xlist)*len(xlist)/2+len(xlist))/10
        xcount=0
        xdecile=10
        for x in range(0,len(xlist)):
            for y in range(x+1,len(xlist)):
                if CCFD_shapes.IsNeighbor(municipalMemberList[xlist[x]]['geometry'],municipalMemberList[xlist[y]]['geometry']):
                    if xlist[y] not in municipalMemberList[xlist[x]]['Neighbors']:
                        municipalMemberList[xlist[x]]['Neighbors'].update({xlist[y]})
                        municipalMemberList[xlist[y]]['Neighbors'].update({xlist[x]})
                        print('Adding Neighbors:',xlist[x],',',xlist[y])
                xcount+=1
                if xcount>counter:
                    print(xdecile,' pct completed.')
                    xcount=0
                    xdecile+=10
    # Does this resolve internal connectedness?
    this_Connected=CCFD.getConnected(list(thisMuniList.keys()),thisMuniList)
    if len(this_Connected)>1:
        print('Error: connectedness issues remain.')
    else:
        print('Connectedness confirmed.')

    # Create a list of the boundary munis for this county and for the surrounding counties
    thisMuniList={x:municipalMemberList[x]['Neighbors'] for x in municipalMemberList.keys() if x[:4]=='M'+inCountyFips}
    boundarySet=CCFD_subregions.BoundaryMembers(thisMuniList)
    for cty in countyMemberList['C'+inCountyFips+'.0']['Neighbors']:
        boundarySet.update(CCFD_subregions.BoundaryMembers({x:municipalMemberList[x]['Neighbors'] for x in municipalMemberList.keys() if x[:4]=='M'+cty[1:4]}))

    xlist=list(boundarySet)
    counter=(len(xlist)*len(xlist)/2+len(xlist))/10
    xcount=0
    xdecile=10
    for x in range(0,len(xlist)):
        for y in range(x+1,len(xlist)):
            if CCFD_shapes.IsNeighbor(municipalMemberList[xlist[x]]['geometry'],municipalMemberList[xlist[y]]['geometry']):
                if xlist[y] not in municipalMemberList[xlist[x]]['Neighbors']:
                    municipalMemberList[xlist[x]]['Neighbors'].update({xlist[y]})
                    municipalMemberList[xlist[y]]['Neighbors'].update({xlist[x]})
                    print('Adding Neighbors:',xlist[x],',',xlist[y])
            xcount+=1
            if xcount>counter:
                print(xdecile,' pct completed.')
                xcount=0
                xdecile+=10


    # No resolve the boundary neighbors

    return municipalMemberList

def LoadList(regionLevel,applyExclusions=True,dropOutside=True): #Valid inputs are County,Municipal,VTD
    print()
    print('Loading Initial Member Lists for ',regionLevel)

    # Build the exclusions dictionary
    if applyExclusions:
        exclusionDict=CCFD.BuildExclusions(G.neighborExclusions)

    # Load the list
    try:
        outList=CCFD_files.JsonToList(G.thisPath+regionLevel+'MemberList_Initial'+G.dataInitID)
    except:
        print('Failed to load:',G.thisPath+regionLevel+'MemberList_Initial'+G.dataInitID)
        return False,None

    # The 'Outside' region should appear only in neighbors lists
    if 'Outside' in outList.keys() and dropOutside:
        del outList['Outside']

    if applyExclusions:
        # Make any specifically identified neighbor exclusions and inclusions
        print('Disconnecting edges between specific regions:',)
        outList=CCFD.forcedExclusions(exclusionDict,outList)

        # Force removal of specifically identified regions (special discontinuities with zero demographics)
        print('Removing specific entire regions:')
        outList=CCFD.forcedRegionExclusions(G.regionExclusions,outList)

    # Add neighbor connections not identified in original build.
    print('Adding edges between specific regions')
    outList=CCFD.forcedInclusions(G.neighborForced,outList)

    if regionLevel=='municipal':
        print('Municipal miscellaneous initial repairs ...')
  
        print('  Checking for small geometry irregularities.')
        for xgroup in ['M041','M045','M061']:
            print('    Checking ',xgroup)
            CCFD_shapes.ScrubPolygons(xgroup,outList)

        # Temporary entry here while its existence is ambiguous
        special='000000000.0'
        if special in outList.keys():
            print('  Removing',special)
            del outList[special]

    # Initial checks
    print('Initial input checks:',end='')
    results=True
    results=results and CCFD_checks.CheckGeometries(outList)
    results=results and CCFD_checks.CheckRegionNeighbors(outList)
    results=results and CCFD_checks.CheckNeighborConsistency(outList)
    results=results and CCFD_checks.CheckPopTotal(outList)

    if not results:
        print(' failed.')

        holdVerbose=copy(G.verbose)
        G.verbose=True
        GEOresults=CCFD_checks.CheckGeometries(outList)
        NBRresults=CCFD_checks.CheckRegionNeighbors(outList)
        NBR2results=CCFD_checks.CheckNeighborConsistency(outList)
        POPresults=CCFD_checks.CheckPopTotal(outList)
        G.verbose=copy(holdVerbose)

        results=GEOresults and NBRresults and NBR2results
        results=results and (POPresults or not G.populationCheckFlag)
    else:
        print(' OK')

    # Check for termination
    if not results:
        print('Exiting ...')
        sys.exit(0)
    else:
        print('Initial loading and checks complete for',regionLevel)
        print()

    return results,outList
