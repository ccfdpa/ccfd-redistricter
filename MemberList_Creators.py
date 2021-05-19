import os
import pandas as pd
import geopandas as gp
import time
import sys
import json
import fiona

from math import pi
from copy import deepcopy, copy
from operator import itemgetter
from shapely.geometry import Polygon
from shapely.ops import unary_union,linemerge,transform
from pysal.lib import weights # Needed to constract contiguity matrix

import File_Handlers as CCFD_files
import Shape_Handler as CCFD_shapes
import Region_Checking_Functions as CCFD_checks
import Global_IO as F
import Global_Settings as G
import CCFD_Functions as CCFD


def inCSVFile(inFile,inKeepList=None):        
    rawData=pd.read_csv(inFile)
        
    if inKeepList is None:
        return rawData
    # Keep only what we need
    else:
        return  rawData[[x for x in inKeepList.keys()]]


# Load the existing boundaryLengthDict
def tryOpenBoundaryDict(newBoundary=False):

    if newBoundary:
        G.boundaryLengthDict={}
        print('Creating new boundaryLengthDict file.')
        return

    try:
        print('Loading boundary length list:',end='')
        G.boundaryLengthDict=json.load(open(G.thisPath+'boundaryLengthDict'+G.dataInitID+'.json','r'))
        results=CCFD_checks.checkRegionBoundaries(G.municipalMemberList)
        if not results:
            print()
            print(' Terminating. Bad boundary dictionary.')
            sys.exit(0)
        else:
            G.memberListOpenFlags['boundaryDict']=True
            print(' OK')
    except:
        print('Error: boundaryLengthDict not found at: ',G.thisPath+'boundaryLengthDict'+G.dataInitID+'.json')
        print('Terminating ...')
        sys.exit(0)
    return

# Controller for computing the boundaries
# regionCode: one of ['state','county','municipal','VTD']
# ctyList: None or a List of 3 digit county codes. The included counties must be contiguous.
# inBoundaries: If True, loads an existing boundary dictionary
# inBad: if True, loads an existing dictionary of regions whose boundary assembly failed.
def boundaryRepair(regionCode,ctyList=None,inBoundaries=True,inBad=False):

    # Open an existing bad length dictionary
    if inBad:
        try:
            inBadFrame=inCSVFile(G.thisPath+'Bad Boundaries'+G.dataInitID+'.csv')
            allBadDict={inBadFrame.loc[x,'ID']:inBadFrame.loc[x,'Value'] for x in inBadFrame.index}
            print('File opened: allBadDict:',len(allBadDict),' bad length sets.')
        except:
            print('Creating new allBadDict file.')
            allBadDict={}
    else:
        allBadDict=None

    if inBoundaries:
        try:
            boundaries=json.load(open(G.thisPath+regionCode+'BoundaryLengths'+G.dataInitID+'.json','r'))
            print('File opened: Boundary Length Dictinary.')
        except:
            print('Failed: Boundary load for ',regionCode)
            print('  :',G.thisPath+regionCode+'BoundaryLengths'+G.dataInitID+'.json')
            boundaries=None
    else:
        boundaries=None

    if ctyList is not None:
        focusRegionList=['42'+x for x in ctyList]
    else:
        focusRegionList=None

    levelBoundaryLengthDict=makeBoundary(regionCode,
                                          fixRegionList=ctyList,
                                          focusRegionList=focusRegionList,
                                          allBadDict=allBadDict,
                                          boundaries=boundaries,
                                          saveJSON=True)

    # Update the full boundary dictionary
    thisFile=G.thisPath+'boundaryLengthDict'+G.dataInitID+'.json'
    if not G.memberListOpenFlags['boundaryDict']:
        try:
            tryOpenBoundaryDict()
        except:
            print('Error: boundaryLengthDict not opened.')
            print('Terminating ...')
            sys.exit(0)
    try:
        G.boundaryLengthDict.update(levelBoundaryLengthDict)
        json.dump(G.boundaryLengthDict,open(thisFile,'w'))
        print('JSON filed:',thisFile)
    except:
        print('Failed: Update of:',thisFile)

    return

# Computes district-level boundary length given a set of correct lengths in the memberList
def calcDistrictBoundaries(inSet,memberList=None):

    if memberList is None:
        memberList=G.regionMemberList

    # Make computations:
    thisDist=inSet['ID']
    G.boundaryLengthDict[thisDist]={}
    for bndyReg in inSet['Neighbors']:
        G.boundaryLengthDict[thisDist][bndyReg]=0.0
        if bndyReg!='Outside':
            # Care must be taken here to include in the calculations from excluded neighbors when needed
            if bndyReg not in G.boundaryLengthDict.keys():
                print('Warning: No boundaries exist for ',bndyReg)
                G.boundaryLengthDict[bndyReg]={}
            nbrSet=set(G.boundaryLengthDict[bndyReg].keys()) & inSet['Regions']
            for xreg in nbrSet:
                if xreg in G.boundaryLengthDict[bndyReg].keys():
                    G.boundaryLengthDict[thisDist][bndyReg]+=G.boundaryLengthDict[bndyReg][xreg]
                else:
                    if G.verbose: print('Pair already in boundary dictionary:',bndyReg,xreg)
            # for normal boundaries, add the transpose length as well
            G.boundaryLengthDict[bndyReg][thisDist]=G.boundaryLengthDict[thisDist][bndyReg]
        else:
            for xreg in inSet['Regions']:
                if 'Outside' in G.boundaryLengthDict[xreg]:
                    G.boundaryLengthDict[thisDist]['Outside']+=G.boundaryLengthDict[xreg]['Outside']

    return True

def calcSetBoundaries(inSet,memberList=None):
    holdVerbose=G.verbose
    if 'Group' in inSet.keys():
        if 'Mcccmmmmm.0' in inSet['Group']:
            G.verbose=True

    if memberList is None:
        memberList=G.regionMemberList

    if not G.memberListOpenFlags['boundaryDict']:
        tryOpenBoundaryDict()

    # Initialize
    levelsDict={G.regionPrecedence[x][1]:x for x in range(len(G.regionPrecedence))}

    # Three possible cases:
    # 1) xreg is a new region, use its group members;
    # 2) xreg is an existing group, deal with it in a special way
    # 3) xreg is a district, use its regions:
    setReg=inSet['ID']
    if 'Group' in inSet.keys():
        if len(inSet['Group'])==0:
            print('No elements in Group for ',inSet['ID'])

    def getSubList(inSet,baseList):
        # Breaks the set membership to a set of individual subegions (no groups or districts)

        # Initialize
        outSet=set()
        checkSet={inSet['ID']}

        # Loop until we get a full set of regions from the raw member set
        while len(checkSet)>0:
            checkReg=checkSet.pop()

            # Find the proper List element to use
            if checkReg==inSet['ID']:
                trySet=copy(inSet)             # Always does this in first loop

            elif checkReg in baseList.keys():
                trySet=baseList[checkReg]
            elif checkReg[0] in levelsDict.keys():
                memberList=CCFD.confirmListAtLevel(levelsDict[checkReg[0]])
                trySet=memberList[checkReg]
            else:
                print('Problem in getSubList',inSet['ID'])
                return None

            if 'Group' in trySet.keys():
                trySub=trySet['Group']
            elif 'Regions' in trySet.keys():
                trySub=trySet['Regions']
            else:
                print('Unknown List structure in',inSet['ID'])
                return None

            if len(trySub)>1:
                checkSet.update(trySub)
            else:
                outSet.update(trySub)


        return outSet

    # Further break down members to their group constituents if necessary
#    for thisreg in copy(xregMembers):
#        if thisreg[4]=='G' and thisreg not in G.boundaryLengthDict.keys():
#            xregMembers.update(memberList[thisreg]['Group'])
#            xregMembers.discard(thisreg)

    # Special settings to turn on verbose
    # G.verbose=False
    # if len(setReg)>4:
    #    if setReg[4]=='G':
    #       for x in inSet['Neighbors']:
    #          if len(x)>4:
    #              if x[4]=='G':
    #                G.verbose=True
    # if setReg[0]=='_':
    #    G.verbose=True

    # Check the boundaries we already have.
    if setReg in G.boundaryLengthDict.keys():
        knownBoundaries=set(G.boundaryLengthDict[setReg].keys()) & inSet['Neighbors']
    else:
        knownBoundaries=set()
        G.boundaryLengthDict[setReg]={}
    newBoundaries=inSet['Neighbors']-knownBoundaries
    # These are the absolute new boundaries. # There are also partial new boundaries
    #   associated with the regions these boundaries touch (i.e. some known boundaries are currently partial.

    fullMemSet=getSubList(inSet,memberList)
    touchingRegions=set()
    for xreg in newBoundaries-{'Outside'}:
        touchingRegions.update(fullMemSet & memberList[xreg]['Neighbors'])
    touchingBoundaries=set()
    for xreg in touchingRegions:
        touchingBoundaries.update(memberList[xreg]['Neighbors'] & inSet['Neighbors'])
    newBoundaries.union(touchingBoundaries)
    updateSetDict={xnbr:0.0 for xnbr in newBoundaries}
    updateSetDict.update(G.boundaryLengthDict[setReg])

    # Check to see if this set region has defined boundaries with its neighbors
    Spaces=' '*4
    print(Spaces,'Creating boundary dictionary for ',setReg)
    if G.verbose:
        print(Spaces,'  Region collection:',fullMemSet)
        print(Spaces,'  nbrs:',inSet['Neighbors'])
        print(Spaces,'  known:',knownBoundaries)
        print(Spaces,'  touching:',touchingBoundaries)
        print(Spaces,'  new:',newBoundaries)

    # The 'Outside' region is handled in a special way. It is not in boundaryLengthDict.keys()
    #     but it is in some subkeys (it is not a region but it is a neighbor)

    if 'Outside' in newBoundaries:
        updateSetDict['Outside']=sum([G.boundaryLengthDict[xreg]['Outside'] for xreg in fullMemSet if 'Outside' in G.boundaryLengthDict[xreg].keys()])
        print('       Outside bound:',round(updateSetDict['Outside'],3))
        
    for thisNeighbor in newBoundaries-{'Outside'}:

        if 'Group' in memberList[thisNeighbor]:
            if 'Mcccmmmmm.0' in memberList[thisNeighbor]['Group']:
                G.verbose=True

        if G.verbose: print(Spaces,'  thisNeighbor:',thisNeighbor)

        keepRegFlag=True
        keepNbrFlag=True        
        holdNbrs=set()
        if thisNeighbor not in G.boundaryLengthDict.keys():
            fullNbrMembers=getSubList(memberList[thisNeighbor],memberList)
        elif setReg not in G.boundaryLengthDict[thisNeighbor].keys():
            fullNbrMembers=getSubList(memberList[thisNeighbor],memberList)
        else:
            fullNbrMembers={thisNeighbor}

        nbrMembers=copy(fullNbrMembers)

        # Sum up all the boundaries between the members of thisNeighbor
        #    and the members if this region
        if G.verbose==True:       
            print(Spaces,'  Boundary for:',setReg,end='')
            print(fullMemSet,end='')
            print(' with ',thisNeighbor,end='')
            print(nbrMembers,end='')
            print()

        while len(nbrMembers)>0:
            xnbr=nbrMembers.pop()
            if G.verbose: print('        Processing:',xnbr)
            if xnbr in G.boundaryLengthDict.keys():
                memSet=copy(fullMemSet)
                while len(memSet)>0:
                    xmem=memSet.pop()
                    if xmem in G.boundaryLengthDict[xnbr].keys():
                        if G.verbose: print('          add type 1 using ',xmem,',',xnbr,':',round(G.boundaryLengthDict[xmem][xnbr]*G.perimeterToMiles,4))
                        updateSetDict[thisNeighbor]+=G.boundaryLengthDict[xnbr][xmem]
                    elif xnbr in G.boundaryLengthDict[xmem].keys():
                        if G.verbose: print('        add type 2 using ',xnbr,',',xmem,':',round(G.boundaryLengthDict[xnbr][xmem]*G.perimeterToMiles,4))
                        updateSetDict[thisNeighbor]+=G.boundaryLengthDict[xmem][xnbr]
                    else:
                        if G.verbose: print('        No boundary pair for ',xmem,xnbr)
                        # Is this pair really separated? If so, one must be enclosed.
#                        memSet.update(getSubList(memberList[xmem]))
#                        memSet.discard(xmem)
            else:
                if G.verbose: print('        No boundary values for ',xnbr)


#                if G.verbose: print('        Retrying with neighbor members',fullMemSet,nbrMembers,'(',holdNbrs,')')

        if updateSetDict[thisNeighbor]<0.000001:
            print('  Missing boundary distance:',setReg,thisNeighbor)
            print('    Used member Set:',fullMemSet)
            print('         neighbor set:',fullNbrMembers)
            # What is the base list for these
            thisRegionCode=inSet['ID'][0]
            thisLevel=levelsDict[thisRegionCode]
            baseList=CCFD.confirmListAtLevel(thisLevel)

            if baseList is not None:
                plotSet=set()
                for x in getSubList(inSet,memberList):
                    if x[4]=='G':
                        plotSet.update(getSubList(memberList[x],memberList))
                    else:
                        plotSet.update(baseList[x]['Group'])
                CCFD_shapes.GroupPlot({setReg:{x:baseList[x]['geometry'] for x in plotSet},
                                                thisNeighbor:{x:baseList[x]['geometry'] for x in memberList[thisNeighbor]['Group']}},
                                        showPlots=True)

    if setReg=='Outside':
        print(' Should not occur.')
    G.boundaryLengthDict[setReg]=updateSetDict
    # Set the paired length
    for xnbr in set(updateSetDict.keys())-{'Outside'}:
        if xnbr not in G.boundaryLengthDict.keys():
            G.boundaryLengthDict[xnbr]={}
        if G.verbose: print('       Setting reverse boundary for',setReg,'with',xnbr)
        G.boundaryLengthDict[xnbr][setReg]=updateSetDict[xnbr]

    if G.verbose: print(Spaces,'  ',end='')
    results=CCFD_checks.checkRegionBoundaries({setReg:inSet})
    if not results:
        print('Error: inconsistent boundary for',setReg)
 
    G.verbose=holdVerbose
    return results

# Assemble boundaries dictionary only
# Argument regionCode below can be one of ['VTD','municipal','county', or 'state']
# Argument 'fixRegionList' (optional) recalculates boundary lengths for specific counties only by countycode (e.g. '003')
# Argument 'focusRegionList' (optional) creates a small memberList before computing boundaries. Can be faster than use of fixRegionList
# Argument 'allBadDict' (optional) works only on regions with failing checkes. Uses a specific dictionary.
# Argument 'boundaries' inputs and existing boundary Length dictionary, which will be modified then returned. If 'None', creates a new one.
# Argument 'saveJSON' resaves the regional boundaryDict before return
def makeBoundary(regionCode,fixRegionList=None,focusRegionList=None,allBadDict={},boundaries=None,saveJSON=True):

    if boundaries is None:
        boundaries={}

    # Sometimes its quicker to do a quick assembly of a smalll member list rather than load all of an existing one
    if focusRegionList is not None:
        memberList=CreateList(regionCode,Exit=False,rawSubList=focusRegionList)
    else:
        results,memberList=CCFD_files.LoadList(regionCode,applyExclusions=False,dropOutside=False)

    newBoundaries=makeBoundaryLengthDict(memberList,boundaries,calcRegionList=fixRegionList,allBadDict=allBadDict)
    boundaries.update(newBoundaries)
    if saveJSON:
        thisFile=G.thisPath+regionCode+'boundaryLengthDict'+G.dataInitID+'.json'
        try:
            json.dump(boundaries,open(thisFile,'w'))
            print('JSON filed:',thisFile)
        except:
            print('Failed: JSON file:',thisFile)
    return boundaries

def makeBoundaryLengthDict(inList,inDict=None,IDSnippetStart=None,calcRegionList=None,allBadDict={}):

    if inDict is None and G.boundaryLengthDict is not None:
        inDict=G.boundaryLengthDict
    if allBadDict is None:
        allBadDict={}

    def prtIndent(indent):
        return ' '+'  '*(indent-1)

    indent=0
    prtSpaces=prtIndent(indent)

#    G.verbose=True
    print()
    print('Creating boundary lengths ')
    if len(allBadDict)>0:
        print('  Starting with bad list of ',len(allBadDict),'regions.')

    indent+=1
    prtSpaces=prtIndent(indent)

    allowedBoundarySumVariance=0.01 # In km, roughly the length in km of the narrow dimension of a typical urban lot (10 meters, 33 feet)
    bufferSize=0.000004  # Somewhat arbitrary buffer length for line-polygon intersections (set by experimentation)

    def calcPreciseLength(inGEO):

        # Set to a cartesian projection from the approximate centroid
        gShape=gp.GeoSeries(inGEO)
        centroids=gShape.centroid
        gShape=CCFD_shapes.changeCRS(gShape,'cartesian',zeroLongLatTuple=(centroids.x,centroids.y))

        thisGeometry=gShape.geometry
        thisExterior=thisGeometry.boundary[0]
        return thisExterior.length


    # Set up to automatically use a shortened ID for selected ID lengths
    # The snippet is only used for output.
    # Note: Boundaries created county by county.
    if IDSnippetStart is None:
        IDSnippedStart=0
        # Pre-selected snippets are based on length
        IDLen=max([len(x) for x in list(inList.keys())[:10]])
        # County snippets start on character zero
        if IDLen<8:
            IDSnippetStart=0
        # Municipal snippets start on character 4
        elif IDLen<13:
            IDSnippetStart=4
        # VTD snippets start on character 9 
        elif IDLen<19:
            IDSnippetStart=9

    regList=list(inList.keys())

    # Run through the input list and create a list of county codes
    # Note: this step needs to be generalized
    print('Computing boundaries using county groups:')
    countyList=sorted(list({x[1:4] for x in inList.keys()}))
    if calcRegionList is not None:
        countyList=sorted([x for x in countyList if x in calcRegionList])

    outDict={}
    countyCount=0
    allBounds=0
    seqList=countyList.copy()
    if 'Outside'[1:4] in countyList:
        seqList.remove('Outside'[1:4])

    reTries=0
    while reTries<2:
        if reTries>0:
            print()
            print('Rerunning the boundary length creation. Retry',reTries)
            indent=1
            prtSpaces=prtIndent(indent)
        reTries+=1
        for thisCounty in seqList:

            # Create a selection list. Include all elements of the county and all out-of-county neighbors
            totalBounds=0
            countyCount+=1
            badDict={}
            print(prtSpaces,'County',thisCounty)

            indent+=1
            prtSpaces=prtIndent(indent)

            #Pare down the calculation frame:
            # Include all regions in this county:
            regionList={x for x in inList.keys() if x[1:4]==thisCounty}
            # Also include all of the neighbors of this list. The easy way is just recreate the set
            # Also, include any neighbors deleted using the exlusion list
            nbrList=set()
            for thisRegion in regionList:
                nbrList.update(inList[thisRegion]['Neighbors'])
            IDList=list(nbrList|regionList)

            # Create a dataframe wih IDs and geometries from the input memberList
            calcFrame=gp.GeoDataFrame()

            # Add the geometry and set the CRS
            calcFrame['geometry']=[inList[x]['geometry'] for x in IDList]
            calcFrame=CCFD_shapes.changeCRS(calcFrame,'default')

            # Input the location information and the input areas and lengths
            calcFrame['InArea']=[inList[x]['Data']['Area'] for x in IDList]
            calcFrame['InPerimeter']=[inList[x]['Data']['Perimeter'] for x in IDList]
            calcFrame['Latitude']=[inList[x]['Data']['Latitude'] for x in IDList]
            calcFrame['Longitude']=[inList[x]['Data']['Longitude'] for x in IDList]
            calcFrame['ID']=IDList
            calcFrame.set_index('ID',inplace=True)

            # add Outside zeros. These values are never used. Only 'geometry' is used.
            calcFrame.loc['Outside','InArea']=0.0
            calcFrame.loc['Outside','InPerimeter']=0.0

            # Boundary calcs are performed in KM and sq KM.
            calcFrame['InArea']=calcFrame['InArea']/G.areaToSqMiles
            calcFrame['InPerimeter']=calcFrame['InPerimeter']/G.perimeterToMiles
        
            # Loop over all but 'Outside'
            for row in calcFrame.loc[list(regionList)].itertuples():
            
                # thisReg='C117.0'
                # Reset outDict for this region
                thisReg=row.Index
                thisRegSnippet=thisReg[IDSnippetStart:]
    #            calcFrame.to_crs(G.crs_cea_proj_str,inplace=True)

                targetBoundary=row.InPerimeter

                if thisReg not in outDict.keys():
                    outDict[thisReg]={}
                nbrDict={}
                thisCalc='Unknown'

                thisNbrList=list(inList[thisReg]['Neighbors'])
                if len(set(thisNbrList)-set(calcFrame.index))>0:
                    print('Error')
                xList=copy(thisNbrList)
                xList.append(thisReg)
                nbrFrame=calcFrame.loc[xList,['geometry','InArea','InPerimeter']]

                nbrFrame=CCFD_shapes.changeCRS(nbrFrame,'cartesian',zeroLongLatTuple=(row.Longitude,row.Latitude))

    #            nbrFrame=CCFD_shapes.CalcPreciseGeoms(nbrFrame)
    #            nbrFrame.rename(columns={'Area':'AreaCRS','Perimeter':'AreaPerimeter'},inplace=True)
    #            print(nbrFrame.loc[xList,['InArea','AreaCRS','InPerimeter','AreaPerimeter']])

                IntersectFrame=nbrFrame.drop(index=[thisReg])
                IntersectFrame['geometry']=IntersectFrame['geometry'].buffer(bufferSize,join_style=1)
                thisExterior=nbrFrame.loc[thisReg,'geometry'].exterior


                IntersectFrame['geometry']=IntersectFrame['geometry'].intersection(thisExterior)
                # Note: each end of the intersected line is too long by the size of the buffer addition
                IntersectFrame['Length']=IntersectFrame['geometry'].length-2.0*bufferSize
                IntersectFrame['geom_type']=IntersectFrame.geom_type

                # Case: the region encloses a neighbor (intersection length with target is 0 and target is the only neighbor)
                enclosed=[x for x in IntersectFrame.index if IntersectFrame.loc[x,'Length']<bufferSize and len(inList[x]['Neighbors'])==1]
                # If so, fix
                for xReg in enclosed:
                    xLength=calcFrame.loc[xReg,'InPerimeter']
                    IntersectFrame.loc[xReg,'Length']=xLength
                    IntersectFrame.loc[xReg,'geom_type']='Enclosed.'
                    IntersectFrame.loc[xReg,'geometry']=calcFrame.loc[xReg,'geometry']

                nbrSum=sum([IntersectFrame.loc[x,'Length'] for x in IntersectFrame.index])
                currentVariance=nbrSum-targetBoundary
    #            print('      ',thisReg,'Targ:',round(targetBoundary,3),'; var:',round(currentVariance,3))
                if G.verbose:

                    prtSpaces=prtIndent(indent+1)

                    print(prtSpaces,'Preiminary lengths with',thisReg[IDSnippetStart:])
                    print(prtSpaces,IntersectFrame[['Length','geom_type']].copy())
                    print(prtSpaces,'Region:',round(targetBoundary,3),'; Sum of nbr boundaries:',round(nbrSum,3),'; diff:',round(currentVariance,3))

                    prtSpaces=prtIndent(indent)

                # Note Any problems:
                indent+=1
                prtSpaces=prtIndent(indent)
                if abs(currentVariance)>allowedBoundarySumVariance:

                    print(prtSpaces,'Attempting to resolve large variance for ',thisReg,':',currentVariance)

                    for thisNbr in IntersectFrame.index:
                        print(prtSpaces,'  ',thisNbr,round(IntersectFrame.loc[thisNbr,'Length'],3),IntersectFrame.loc[thisNbr,'geom_type'])
                    print(prtSpaces,'Region:',round(targetBoundary,3),'; Sum of nbr boundaries:',round(nbrSum,3),'; diff:',round(currentVariance,3))

    # Run through occasional remedies here.

                # Look to see if this issue was resolved in a previous makeBoundaryList run.
                if inDict is not None and abs(currentVariance)>allowedBoundarySumVariance:
                    # Sum up the lengths of the neighber intersections
                    if thisReg in inDict.keys():
                        boundaryNbrs=set(inDict[thisReg].keys())
                        listNbrs=inList[thisReg]['Neighbors']
                        if len(listNbrs-boundaryNbrs)==0:
                            # all of the neighbors are in the boundaryDict
                            checkNbrs=sum(inDict[thisReg][xnbr] for xnbr in listNbrs)
                            if abs(checkNbrs-targetBoundary)<allowedBoundarySumVariance:
                                # Just substitute the earlier lengths:
                                IntersectFrame['Length']=[inDict[thisReg][xNbr] for xNbr in IntersectFrame.index]
                                nbrSum=sum([IntersectFrame.loc[x,'Length'] for x in IntersectFrame.index])
                                currentVariance=nbrSum-targetBoundary
                                print(prtSpaces,'Repaired with substitution from existing boundaryDict')

                # Try eliminating spurious lines from target
                if abs(currentVariance)>allowedBoundarySumVariance: 
                    # Try Scrubbing the target and recalculate its boundary
                    fixedRegion=CCFD_shapes.ScrubPolygons(CCFD.coreID(thisReg),inList)
                    if len(fixedRegion)>0:
                        newBoundary=calcPreciseLength(fixedRegion[thisReg]['geometry'])
                        if newBoundary<targetBoundary:
                            # Seems to be a fix
                            currentVariance-=targetBoundary-newBoundary
                            targetBoundary=copy(newBoundary)
                            # No boundary recalcs are needed.
                            print(prtSpaces,'Repairing by scrubbing out Target interior lines:',round(currentVariance,0))

                # Try using the reverse calc if it has already been calculated.
                if abs(currentVariance)>allowedBoundarySumVariance: 
                        # This calculation is redundant in [thisReg][thisNbr] and [thisNbr][thisReg]
                        #    unless the spatial intersections differ. Note when this is the case
                    for thisNbr in IntersectFrame.index:
                        if thisNbr in outDict.keys():
                            if thisReg in outDict[thisNbr].keys():
                                boundaryVariance=outDict[thisNbr][thisReg]-IntersectFrame.loc[thisNbr,'Length']
                                # If this variance reduces the nbrSum (currentVariance) use the new one
                                if abs(currentVariance+boundaryVariance)+0.15*allowedBoundarySumVariance<abs(currentVariance):
                                    outDict[thisReg][thisNbr]=outDict[thisNbr][thisReg]
                                    currentVariance+=boundaryVariance
                                    print(prtSpaces,'Repairing by substituting the reverse intersection length from',thisNbr,'; new Variance:',round(currentVariance,3))

                # Try using the area intersection (snake method).
                if abs(currentVariance)>allowedBoundarySumVariance: 
                    # Create an intersection area between the two with a small buffer.
                    # Then assume that the area of the intersection is a long rectangle.
                    holdVerbose=copy(G.verbose)
    #                G.verbose=True
                    snakeBuffer=0.005
                    xFrame=nbrFrame.drop(index=[thisReg])
                    xFrame['ID']=xFrame.index
                    xFrame['geometry']=xFrame['geometry'].buffer(snakeBuffer,join_style=1)
                    xFrame.reset_index
                    targetFrame=nbrFrame.loc[[thisReg]]
                    targetFrame['geometry']=targetFrame['geometry'].buffer(snakeBuffer,join_style=1)
                    snakes=gp.overlay(xFrame,targetFrame,how='intersection')
                    # Divide by the twice the width of the snake (the buffer size) yields a length
                    #   less 0.5 * pi * buffer to account for extra area at ends.
                    snakes['snakeLength']=snakes.area/2.0/snakeBuffer - 0.5*pi*snakeBuffer
                    snakes['geom_type']=snakes.geom_type
                    snakes.set_index('ID',inplace=True)
                    if G.verbose:
                        print(prtSpaces,'Intersection approach:')
                        xgeom={thisReg:targetFrame.loc[thisReg,'geometry']}
                        xgeom.update({x:xFrame.loc[x,'geometry'] for x in xFrame.index})
                        CCFD_shapes.QuickPlot(xgeom,showPlots=True)
                    G.verbose=copy(holdVerbose)
                    if abs(sum(snakes.snakeLength)-targetBoundary)+0.15*allowedBoundarySumVariance<abs(currentVariance):
                        # This approach worked, at least partially
                        print(prtSpaces,'Applying intersection approach estimates:')
                        for thisNbr in IntersectFrame.index:
                            if abs(currentVariance+snakes.loc[thisNbr,'snakeLength']-IntersectFrame.loc[thisNbr,'Length'])<abs(currentVariance):
                                newValue=snakes.loc[thisNbr,'snakeLength']
                                currentVariance+=newValue-IntersectFrame.loc[thisNbr,'Length']
                                IntersectFrame.loc[thisNbr,'Length']=copy(newValue)
                                print(prtSpaces,'  ','Substituting new value for ',thisNbr,'; new Variance:',round(currentVariance,3))
                        nbrSum=sum([IntersectFrame.loc[x,'Length'] for x in IntersectFrame.index])
                        currentVariance=nbrSum-targetBoundary
                        if G.verbose:
                            for thisNbr in IntersectFrame.index:
                                print(prtSpaces,' ',thisNbr,round(IntersectFrame.loc[thisNbr,'Length'],3),IntersectFrame.loc[thisNbr,'geom_type'])
                        print(prtSpaces,'Region:',round(targetBoundary,3),'; Sum of nbr boundaries:',round(nbrSum,3),'; diff:',round(currentVariance,3))
 
                # This region is entirely enclosed by one neighbor. Then, the boundary with the neighbor is its exterior length
                # Probably only occurs in the case of 'S.0' region and 'Outside'
                if abs(currentVariance)>allowedBoundarySumVariance: 
                    if len(thisNbrList)==1:
                        thisNbr=list(thisNbrList)[0]
                        print(prtSpaces,'This region is entirely enclosed by ',thisNbr,' Use its length as the boundary.')
                        IntersectFrame.loc[thisNbr,'Length']=targetBoundary
                        nbrSum=targetBoundary
                        currentVariance=0.0

                # If all the above have failed, keep track of failed region. (There should be only a few.)
                if abs(currentVariance)>allowedBoundarySumVariance:
                    print(prtSpaces,'Boundary length calculation failure for',thisReg)                    
                    if G.verbose:
                        plotDict={thisReg:thisExterior}
                        plotDict.update({x:IntersectFrame.loc[x,'geometry'] for x in IntersectFrame.index})
                        CCFD_shapes.QuickPlot(plotDict,showPlots=True)
                    badDict[thisReg]=currentVariance

                # Enter the computed and corrected data for this row.
                if thisReg not in outDict.keys():
                    outDict[thisReg]={}
                outDict[thisReg].update({x:IntersectFrame.loc[x,'Length'] for x in IntersectFrame.index})
                # Note: the full outDict could include some 'Outside' entries associated with regions on county borders
                #          the 'Outside' length will be the length along the county border.
                totalBounds+=len(IntersectFrame.index)

                indent-=1
                prtSpaces=prtIndent(indent)

            indent-=1
            prtSpaces=prtIndent(indent)
                    
            badList=list(badDict.keys())
            numBad=len(badList)
            if numBad>0:
                print(prtSpaces,'Attempting to match existing decompositions from full dictionary:',numBad)
            # Keep only combinations already in outDict

            indent+=1
            prtSpaces=prtIndent(indent)
            badList.sort()
            for i in range(0,numBad):
                thisReg=badList[i]
                print(prtSpaces,'Adjusting ',thisReg,'; Variance:',round(badDict[thisReg],3))
                currentVariance=badDict[thisReg]
                # Go through its neighbors to see if they differ
                for thisNbr in outDict[thisReg].keys():
                    if thisNbr in outDict.keys():
                        if thisReg in outDict[thisNbr].keys():
                            lenDiff=outDict[thisNbr][thisReg]-outDict[thisReg][thisNbr]
                            if abs(currentVariance+lenDiff)+0.15*allowedBoundarySumVariance<abs(currentVariance):
                                # The reverse length improves the result
                                outDict[thisReg][thisNbr]=outDict[thisNbr][thisReg]
                                totalBounds+=1
                                currentVariance+=lenDiff
                                print(prtSpaces,' ','Repairing with ',thisNbr,':',round(currentVariance,3))
                # Recalculate difference and evaluate:
                nbrSum=sum([outDict[thisReg][x] for x in outDict[thisReg].keys()])
                targetBoundary=calcFrame.loc[thisReg,'InPerimeter']
                currentVariance=nbrSum-targetBoundary
                if abs(currentVariance)<allowedBoundarySumVariance:
                    print(prtSpaces,' ','Variance resolved using existing neighbor length for',thisReg,'; Variance:',round(currentVariance,3))
                    del badDict[thisReg]
                    if thisReg in allBadDict.keys():
                        del allBadDict[thisReg]
                else:
                    badDict.update({thisReg:currentVariance})

            indent-=1
            prtSpaces=prtIndent(indent)

            # Adjust allBadDict for rerunning this county by deleting all references for this county
            #   then updating with this county's badDict
            for thisBad in list(allBadDict.keys()):
                if thisBad[1:4]==thisCounty:
                    del allBadDict[thisBad]

            if len(badDict)>0:
                print(prtSpaces,'  Unresolved boundary lengths:',len(badDict),'of',totalBounds)
                allBadDict.update(badDict)
                print(prtSpaces,'  Bad:',list(badDict.keys()))
                seqList=sorted(list({x[1:4] for x in allBadDict.keys()}))
            else:
                print(prtSpaces,'  Total bounds (no errors):',totalBounds)
                reTries=1000

            allBounds+=totalBounds
            indent-=1
            prtSpaces=prtIndent(indent)


        # First, save the regions with bad boundaries.
        if len(allBadDict)>0:
            try:
                CCFD_files.fileOutCSV(allBadDict,runID='Bad Boundaries'+G.dataInitID)
                print(prtSpaces,'Bounds created for:',allBounds,'regions; bad Count:',len(allBadDict))
                for x in allBadDict.keys():
                    print(prtSpaces,'  ',x,round(allBadDict[x],3))
            except:
                print(' Fail: bad boundary csv file not created.')

        # Provide summary stat    
        print('Boundary Dictionary entries completed for ',countyCount,'counties.')
        print('  Bounds:',allBounds,'; Bad:',len(allBadDict))

    # Returns a dictionary of dictionaries: outDict[region1 ID][region2 ID]=boundarylength in miles
    return outDict

def Initial_Data_Input(xfile):

    alphalist='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    splits=0
    # Fix any multipolygon issues
    badList=[x for x in xfile.index if not CCFD_shapes.checkPolygon(xfile.loc[x,'geometry'],verbose=False)]
    if len(badList)>0:
        print('  Splitting ',str(len(badList)),'instances of non-polygonal geometries')
    for i in badList:
        # Create a small dictionary of the polygons in badList
        badPolys={x:xfile.loc[i,'geometry'][x] for x in range(len(xfile.loc[i,'geometry']))}
        # Sort them into a list, largest area to smallest area
        sortList=[x[0] for x in sorted({x:badPolys[x].area for x in badPolys.keys()}.items(),key=lambda item: item[1],reverse=True)]

        # Now add each of these polygons to the region List, distributing the population
        splits+=len(badPolys)-1

        baseID=xfile.loc[i,'regionID']
        basename=xfile.loc[i,'Name']
        totpop=xfile.loc[i,'Size']
        totarea=xfile.loc[i,'geometry'].area
        addbad=pd.DataFrame([[baseID+alphalist[x],basename+' Pt '+alphalist[x],
                              int(badPolys[sortList[x]].area/totarea*(float(totpop)+0.45)),
                              badPolys[sortList[x]]] 
                              for x in range(len(sortList))],columns=['regionID','Name','Size','geometry'])
        # Make sure split populations still sum
        popvar=totpop-addbad['Size'].sum()
        while popvar!=0:
            # Add or subtract one from each split in order from largest to smallest until equal
            for seq in range(popvar):
                addbad.loc[sortList[seq],'Size']+=1
            popvar=totpop-addbad['Size'].sum()

        
        addbad.set_index('regionID',drop=False,inplace=True)
        xfile=xfile.append(addbad,sort=True)
        xfile.drop([i],inplace=True)

    xfile.set_index('regionID',drop=False,inplace=True)
    print('  Increased base geometries by ',splits)

    if xfile['Size'].sum()!=G.populationTotal:
        print('Warning: Population sum variance:',xfile['Size'].sum(),' (Actual=',G.populationTotal,')')
    else:
        print('  Population total after split additions is confirmed.')

    print('  Computing member data ...')
    xfile=CCFD_shapes.CalcPreciseGeoms(xfile)
    xbounds=xfile.bounds
    xfile['North']=xbounds['maxy']
    xfile['South']=xbounds['miny']
    xfile['East']=xbounds['maxx']
    xfile['West']=xbounds['minx']

    # Add the "outside" region
    outsideID='000000000'
    outside=CCFD_shapes.AddOutsideShape(list(xfile['geometry']))
    xfile=xfile.append({'regionID':outsideID,'geometry':outside,'Size':0},ignore_index=True)
    xfile.set_index('regionID',drop=False,inplace=True)


    # print(xfile[['regionID','geometry']])
    # Use pysal to create neighbors list (defines network)
    print('  Creating neighbor sets ...')
    rW = weights.contiguity.Rook.from_dataframe(xfile, idVariable = 'regionID')

    xlist=[set(rW[x].keys()) for x in list(xfile['regionID'])]

    xfile['Neighbors']=xlist
    # rW is a dictionary of county codes whose values are lists of neighbor county codes

    print('  Initial number of regions (including Outside): ',len(xfile))
    print(' ')

    return xfile

def inSHPFile(inFile,invarDict=None):

    print('  Reading geometry file: ',inFile)
    # Read in transfer file into a GeoDataFrame
    try:
        region_shp = os.path.join(F.dataPath, inFile)
        xfile = gp.read_file(region_shp)
        print('    Attribute file import complete:',len(xfile),'regions.')
    except:
        print('File input error:')
        print('path: ',F.dataPath)
        print('filename: ',inFile)
        return None

    # First, lets keep only what we need.
    tryVars=set(invarDict.values())-set(xfile.columns)
    if len(tryVars)>0:
        print('Error on input from',inFile)
        print('  Variables not found:',tryVars)
    xfile=xfile.loc[:,[invarDict[x] for x in invarDict.keys()]]
    # Rename to used column Names
    xfile.rename(columns={invarDict[x]:x for x in invarDict.keys()},inplace=True)

    # If we dont have some demographics or area scalars yet, set them to a selected value.
    valueSpecialSet={'Size':1,'ALAND':0.0,'AWATER':0.0,'CountyRawCode':'ccc','MuniRawCode':'mmmmm','VTDRawCode':'vvvv'}
    for thisVar in valueSpecialSet.keys():
        if thisVar not in xfile.columns:
            print('    Setting ',thisVar,' value to ', valueSpecialSet[thisVar])
            xfile[thisVar]=valueSpecialSet[thisVar]

    return xfile

def DF_to_memberList(xfile,codePrefix,scaleParameter=None):

    print('Creating member list ...')

    if scaleParameter is None:
        scaleParameter=G.perimeterToMiles
    # Set up the region dictionary.
    # The county key has form Cnnn.x
    # The municipal key has form Mcccmmmmm.x
    # The VTD key has form Vcccmmmmmvvvvv.x


    numRegions=len(xfile)
    regionMemberSet=G.regionMemberSet
    regionMemberList={}
    outsideCode=codePrefix+'000000000.0'

    # Create an empty list of region records
    xfile.set_index('regionID',inplace=True)
    for xcode in xfile.index:

        if xcode=='00143056A':
            print('')
        # Base the region code on its FIPS code
        codestr=CCFD.areaID(codePrefix,xcode)  # We'll use the suffix "C" for counties, followed by the 3 digit county FIPS code (state FIPS excluded)

        # Build lists of regions and their neighbors
        #    Region Name consists of a name plus '.n' in cases with multiple districts in one region
        thisRegion=codestr+'.0'
        regionMemberList.update({thisRegion:deepcopy(G.regionMemberSet)})

        regionMemberList[thisRegion]['ID']=thisRegion
        # Initial group is itself
        regionMemberList[thisRegion]['Group'].update({thisRegion})

        regionMemberList[thisRegion]['Name']=xfile.loc[xcode,'Name']

        regionMemberList[thisRegion]['Data']['TAPersons']=xfile.loc[xcode,'Size']
        regionMemberList[thisRegion]['Data']['Area']=xfile.loc[xcode,'Area']*scaleParameter*scaleParameter # geoPandas ESPG 4269 to square miles adjustment in vicinity of PA
        regionMemberList[thisRegion]['Data']['Perimeter']=xfile.loc[xcode,'Perimeter']*scaleParameter # geoPandas ESPG 4269 to miles adjustment in vicinity of PA
        regionMemberList[thisRegion]['Data']['North']=xfile.loc[xcode,'North']
        regionMemberList[thisRegion]['Data']['South']=xfile.loc[xcode,'South']
        regionMemberList[thisRegion]['Data']['East']=xfile.loc[xcode,'East']
        regionMemberList[thisRegion]['Data']['West']=xfile.loc[xcode,'West']
        regionMemberList[thisRegion]['Data']['Longitude']=xfile.loc[xcode,'Longitude']
        regionMemberList[thisRegion]['Data']['Latitude']=xfile.loc[xcode,'Latitude']
        regionMemberList[thisRegion]['geometry']=xfile.loc[xcode,'geometry']

        # Change neighbor set to list of regionIDs
        regionMemberList[thisRegion]['Neighbors']={CCFD.areaID(codePrefix,x)+'.0' for x in list(xfile.loc[xcode,'Neighbors'])}

        # Compute the bounding-square index from this data
        regionMemberList[thisRegion]['Data']['BdgSqrIndex']=min(1.0,regionMemberList[thisRegion]['Data']['Area']/
            CCFD.bounding_square(regionMemberList[thisRegion]['Data']['North'],
                                    regionMemberList[thisRegion]['Data']['South'],
                                    regionMemberList[thisRegion]['Data']['East'],
                                    regionMemberList[thisRegion]['Data']['West'])                                                   
                                )
        if regionMemberList[thisRegion]['Data']['BdgSqrIndex']>1.0:
            print('  Sqr index issue:',thisRegion)

        # For understandability, rename the '000000000' region to 'Outside'
        if  outsideCode in regionMemberList[thisRegion]['Neighbors']:
            # Remove this from the list and add 'Outside'
            regionMemberList[thisRegion]['Neighbors'].remove(outsideCode)
            regionMemberList[thisRegion]['Neighbors'].update({'Outside'})

    # Rename the outside region in the Dictionary
    xreg=regionMemberList[outsideCode]
    xreg['ID']='Outside'
    del regionMemberList[outsideCode]
    regionMemberList['Outside']=xreg

    # Apply any initial forced.
    if len(F.inputNeighborForced)>0:
        regionMemberList=CCFD.forcedInclusions(F.inputNeighborForced,regionMemberList)

    return regionMemberList

def CreateStateList():

    regionCode='S.0'
    # The stateMemberList is created from the countyMemberList
    try:
        if not G.memberListOpenFlags['boundaryDict']:
            G.boundaryLengthDict=json.load(open(G.thisPath+'boundaryLengthDict.json','r'))
            G.memberListOpenFlags['boundaryDict']=True
    except:
        print('Error: boundaryDict not opened.')
        return False

    try:
        if not G.memberListOpenFlags['county']:
            results,G.countyMemberList=CCFD_files.LoadList('County')
            if results:
                G.memberListOpenFlags['county']=True
    except:
        print('Error: countyMemberList not opened .')
        return False

    countyList=list(G.countyMemberList.keys())
    if 'Outside' in countyList:
        countyList.remove('Outside')

    # Create the stateMemberList dictionary
    stateMemberList={regionCode:deepcopy(G.regionMemberSet)}

    # The state geometry is a polygon of the union of the state geometries
    stateMemberList[regionCode]['geometry']=unary_union([G.countyMemberList[x]['geometry'] for x in countyList])
    results=CCFD_shapes.checkPolygon(stateMemberList[regionCode]['geometry'])
    if not results:
        print('Error: geometry union for state total failed.')
        return False

    # Complete the creation of the data for this list
    # Geometry and data vars related to it
    stateMemberList[regionCode]=CCFD_shapes.calcRegionGeoData(stateMemberList[regionCode])
    # Demography variables
    for thisVar in G.regionDemogData:
        stateMemberList[regionCode]['Data'][thisVar]=sum([G.countyMemberList[x]['Data'][thisVar] for x in countyList])
    
    # Create outside region
    stateMemberList['Outside']=deepcopy(G.regionMemberSet)
    stateMemberList['Outside']['geometry']=CCFD_shapes.AddOutsideShape(stateMemberList[regionCode]['geometry'])
    stateMemberList['Outside']['ID']='Outside'
    stateMemberList['Outside']['Neighbors']={regionCode}
    # The Neighbor list is simple
    stateMemberList[regionCode]['Neighbors']={'Outside'}
    # Initial District Assignment:
    stateMemberList[regionCode]['Districts'][G.districtPrefix]=None
    # Initial ID
    stateMemberList[regionCode]['ID']=regionCode
    stateMemberList[regionCode]['Name']='STATE'
    stateMemberList[regionCode]['Group']={regionCode}

    results1=CCFD_files.dictToJson(stateMemberList,G.thisPath,'stateMemberList_Initial'+G.dataInitID)
    if not results1: print('Error: JSON file for stateMemberList_Initial not created.')

    return True,stateMemberList

def CreateCountyList(inCountyfile,inDataDict):

    # We input the area county file only to get the Census Geocode and name
    inFrame=inSHPFile(inCountyfile,inDataDict)
    # Create a small dictionary with the raw county code and Name
    xDict={inFrame.at[x,'regionID']:inFrame.at[x,'Name'] for x in inFrame.index}

    # Countyies are built up from VTD regions
    try:
        VTDFrame=CCFD_files.JsonToList(G.thisPath+'VTDMemberList_Initial'+G.dataInitID,returnFrame=True)
        VTDMemberList=CCFD_files.build_regionMemberList(VTDFrame)
    except:
        print('Error: need VTD geodataframe to create initial county member List')
        print('Exiting ....')
        sys.exit(0)

    # Remove the Outside region prior to aggregating
    VTDFrame=VTDFrame[VTDFrame['ID']!='Outside']

    # Add county code columns to the VTD
    VTDFrame['countyCode']=VTDFrame['ID'].str.slice(1,4)

    # Keep a minimal set of variables -- only those that sum plus the geometry
    keepList=['countyCode','geometry','Area','Perimeter','Neighbors','ID']
    keepList=keepList+G.regionDemogData
    keepList=[x for x in keepList if x in VTDFrame.columns]
    VTDFrame=VTDFrame[keepList]

    # Dissolve to county areas
    countyFrame=VTDFrame.dissolve(by='countyCode',aggfunc='sum')
    countyFrame['Name']=[xDict[x] for x in countyFrame.index]
    countyFrame['regionID']=list(countyFrame.index)
    countyFrame.rename(columns={'TAPersons':'Size','Area':'area'},inplace=True)


    # Special perimeter calculations
    print('  Loading VTD boundary length list:')
    G.boundaryLengthDict=json.load(open(G.thisPath+'VTDBoundaryLengths'+G.dataInitID+'.json','r'))
    results=CCFD_checks.checkRegionBoundaries(VTDMemberList)
    if not results:
        print('  Bad boundary dictionary')
        print()

    # The notion is, you add up VTD perimeters then substract the boundary value of each in-county neighbor
    # Has the effect of subtracting twice the common internal boundaries

    countyFrame['perimeter']=0.0
    for countyCode in countyFrame.index:
        smallFrame=VTDFrame.loc[VTDFrame['countyCode']==countyCode,['ID','Perimeter','Neighbors']]
        boundarySum1=smallFrame['Perimeter'].sum()
        boundarySum2=0.0
        for smallRow in smallFrame.index:
            thisVTD=smallFrame.at[smallRow,'ID']
            neighborSet=CCFD.neighborsList(smallFrame.at[smallRow,'Neighbors'])
            for thisNeighbor in neighborSet:
                if thisNeighbor[1:4]==countyCode:
                    if thisNeighbor in G.boundaryLengthDict[thisVTD]:
                        boundarySum1-=2.0*G.boundaryLengthDict[thisVTD][thisNeighbor]
                else:
                    boundarySum2+=G.boundaryLengthDict[thisVTD][thisNeighbor]
        print(countyCode,'Netted internal boundaries:',boundarySum1,'External boundaries: ',boundarySum2)
        countyFrame.at[countyCode,'Perimeter']=boundarySum1

    # Create the memberList
    G.countyMemberList=DF_to_memberList(countyFrame,'C')

    CountyBoundaries=makeBoundaryLengthDict(G.countyMemberList)
    G.boundaryLengthDict.update(CountyBoundaries)
    json.dump(CountyBoundaries,open(G.thisPath+'countyBoundaryLengths'+G.dataInitID+'.json','w'))
    json.dump(G.boundaryLengthDict,open(G.thisPath+'boundaryLengthDict'+G.dataInitID+'.json','w'))

    return G.countyMemberList

def CreateMuniList(inMuniFile,inDataDict,countySubList=None):

    # Raw VTD input typically comes from Census shapefiles
    try:
        inFrame=inSHPFile(inMuniFile,inDataDict)
        # For debugging, reduce file membership
        if countySubList is not None:
            inFrame=inFrame[inFrame['CountyRawCode'].str[1:4].isin(countySubList)]
            print('  For debugging',len(inFrame),'records retained.')
    except:
        print('Muni shapefile input failed.')
        print('Error exit ...')
        sys.exit(0)

    # Create a muni ID. Call it "regionID"
    muniCode=['M'+('000'+str(inFrame.loc[xcode,'CountyRawCode']))[-3:]+('00000'+str(inFrame.loc[xcode,'MuniRawCode']))[-5:] for xcode in inFrame.index]
    inFrame['regionID']=muniCode
    
    # Special: Define Philadelphia Wards as the Muni Level withn Philadelphia
    specialIndexes=[xcode for xcode in inFrame.index if str(inFrame.loc[xcode,'MuniRawCode'])[:5]=='60000']
    for i in specialIndexes:
        wardNo=inFrame.loc[i,'VTDName'][16:]
        wardNo=wardNo[:wardNo.find('PCT')].strip()
        inFrame.loc[i,'regionID']='M101WD'+('0000'+wardNo)[-3:]
        inFrame.loc[i,'MuniName']='Philadelphia Ward '+wardNo
    nameDict={inFrame.loc[i,'regionID']:inFrame.loc[i,'MuniName'] for i in inFrame.index}
    
    # Special: Define Pittsburgh Wards as the Muni Level within Pittsburgh
    specialIndexes=[xcode for xcode in inFrame.index if str(inFrame.loc[xcode,'MuniRawCode'])[:5]=='61000']
    for i in specialIndexes:
        wardNo=inFrame.loc[i,'VTDName'][14:]
        wardNo=wardNo[:wardNo.find('DIST')].strip()
        inFrame.loc[i,'regionID']='M003WD'+('0000'+wardNo)[-3:]
        inFrame.loc[i,'MuniName']='Pittsburgh Ward '+wardNo
    nameDict={inFrame.loc[i,'regionID']:inFrame.loc[i,'MuniName'] for i in inFrame.index}

    # Aggregate geometries using regionID
    print('  Aggregating VTD to Muni geometries ...')
    inMunifile=inFrame[['regionID','Size','geometry']].dissolve(by='regionID',aggfunc='sum')

    print('Initial Muni files after aggregation:',len(inMunifile))

    if inMunifile['Size'].sum()!=G.populationTotal:
        print('Warning: Population sum variance:',inMunifile['Size'].sum(),' (Actual=',G.populationTotal,')')
    else:
        print('  Population total after dissolve is confirmed.')

    inMunifile['Name']=[nameDict[x] for x in inMunifile.index]
    inMunifile.reset_index(inplace=True)

    # Input initial VTD-level geometry and data
    inFrame=Initial_Data_Input(inMunifile)

    if inFrame['Size'].sum()!=G.populationTotal:
        print('Warning: Population sum variance:',inFrame['Size'].sum(),' (Actual=',G.populationTotal,')')
    else:
        print('  Population total for inFrame is confirmed.')

    G.municipalMemberList=DF_to_memberList(inFrame,'')


    # Special Fix for the "outside" region
    outSide='Outside'
    # Special for Allegheny county. None of regions is a neighbor of outSide.
    for xreg in {x for x in G.municipalMemberList.keys() if x[:4]=='M003'}:
        G.municipalMemberList[xreg]['Neighbors'].discard(outSide)

 
    # Special check of neighbors for municipalities
    def FixNeighbors(inList):

        seqList=list(inList)
        seqList.sort()

        for xseq in range(len(seqList)):
            id=seqList[xseq]

            if G.verbose: print('  Testing: ',id)
            thisGeo=G.municipalMemberList[id]['geometry']

            # Remove neighbors if incorrect
            for x in {xreg for xreg in G.municipalMemberList[id]['Neighbors'] if xreg!=outSide}:
                if not CCFD_shapes.IsNeighbor(thisGeo,G.municipalMemberList[x]['geometry']):
                    print('    Removing neighbors ',id,x)
                    G.municipalMemberList[id]['Neighbors'].remove(x)
                    G.municipalMemberList[x]['Neighbors'].discard(id)

            # Now, add neighbors if they are missing. This requires a larger comparison set.
            # To reduce the comparisons, limit neighbors checks to munis in this county and the county's neighbors
            # Also we eliminate in-county comparisons we've already made (uses if a is a neighbor of b then b is a neighbor of a)
            thisCounty=G.municipalMemberList[id]['ID'][1:4]
            testDict={x:G.municipalMemberList[x]['geometry'] for x in G.municipalMemberList if G.municipalMemberList[x]['ID'][1:4]==thisCounty and x in seqList[xseq+1:]}
            for cty in G.countyMemberList['C'+thisCounty+'.0']['Neighbors']:
                xcode=cty[1:4]
                testDict.update({x:G.municipalMemberList[x]['geometry'] for x in G.municipalMemberList if G.municipalMemberList[x]['ID'][1:4]==xcode})        

            # Add neighbors if missing
            for x in testDict.keys():
#                if (id=='M003CD005.0' and x=='M00335424.0') or (x=='M003CD005.0' and id=='M00335424.0'):
#                    print('Checking: ',id,x)
                if CCFD_shapes.IsNeighbor(thisGeo,testDict[x]):
                    if x not in G.municipalMemberList[id]['Neighbors']:
                        print('    Adding neighbor ',id,x)
                        G.municipalMemberList[id]['Neighbors'].update({x})
                        G.municipalMemberList[x]['Neighbors'].update({id})
                
    # Special fix for municipality neighbors using inefficient method
    islandList=[x for x in G.municipalMemberList.keys() if len(G.municipalMemberList[x]['Neighbors'])==0]    
    print('  Neighbors fix for islands.')
    FixNeighbors(islandList)
    # Special fix for subregions in Philadelphia
    print('  Neighbors fix for Philadelphia')
    FixNeighbors([x for x in G.municipalMemberList.keys() if G.municipalMemberList[x]['ID'][1:4]=='101'])
    # Special fix for subregions in Allegheny county
    print('  Neighbors fix for Allegheny county')
    FixNeighbors([x for x in G.municipalMemberList.keys() if G.municipalMemberList[x]['ID'][1:4]=='003'])

    
    MunicipalBoundaries=makeBoundaryLengthDict(G.municipalMemberList)
    G.boundaryLengthDict.update(MunicipalBoundaries)
    json.dump(MunicipalBoundaries,open(G.thisPath+'MunicipalBoundaryLengths'+G.dataInitID+'.json','w'))
    json.dump(G.boundaryLengthDict,open(G.thisPath+'boundaryLengthDict'+G.dataInitID+'.json','w'))

    return G.municipalMemberList

def CreateVTDList(inVTDFile,inDataDict,regionSubList=None):


    # Raw VTD input typically comes from Census shapefiles
    try:
        rawFrame=inSHPFile(inVTDFile,inDataDict)
        # For debugging, reduce file membership
        if regionSubList is not None:
            subLen=len(regionSubList[0])
            rawFrame=rawFrame[rawFrame['rawRegionCode'].str[:subLen].isin(regionSubList)]
            print('  For debugging',len(rawFrame),'records retained.')
            print('    To modify, remove or change argument *rawSubList* in the CreateList line.')
    except:
        print('VTD shapefile input failed.')
        print('Error exit ...')
        sys.exit(0)

    rawFrame.set_index('rawRegionCode',inplace=True)
    # Create key geometricc characteristics for comparison with Census

    # Create a memberList for VTDs
    # At this stage, we don't know the Municipal codes yet'
    rawFrame['regionID']=rawFrame['CountyRawCode']+rawFrame['MuniRawCode']+('00000'+str(inFrame.loc[xcode,'VTDRawCode']))[-4:]

    # At this point, the 'regionID' variables have to be unique. Test this.
    IDCheck=sorted(list(rawFrame.regionID))
    badID=set()
    for seq in range(1,len(IDCheck)):
        if IDCheck[seq]==IDCheck[seq-1]:
            badID.update({IDCheck[seq]})
    if len(badID)>0:
        print('Error: duplicate IDs found')
        print(badID)
        print('Terminating ...')
        sys.exit(0)

    # From census, total area  in square kilometers
    rawFrame['AreaCensus']=(rawFrame['ALAND']+rawFrame['AWATER'])*F.VTDinputAreaScaleFactor


    #Create memberList
    # Input initial adjacency and data for the regions
    rawFrame=Initial_Data_Input(rawFrame)

    # Document area variance: computed vs Census
    # Use this as evidence for meaningful perimeter calculation
    rawFrame['area_Variance']=rawFrame['Area']-rawFrame['AreaCensus']

    xFrame=rawFrame[['Name','AreaCensus','area_Variance']]
    yFrame='V'+rawFrame.regionID
    outFrame=pd.concat([xFrame,yFrame],axis=1)
    outFrame.set_index('regionID',inplace=True)
    try:
        outFrame.to_csv(G.thisPath+'VTD_areas.csv',mode='w')
        print('  Area variance file created: ',G.thisPath+'VTD_areas.csv')
    except:
        print('  Area variance file action failed.')

    memberList=DF_to_memberList(rawFrame,'V')


    return memberList

def CreateRegionList(inRegionFile,inDataDict,regionSubList=None,regionLevel='Region',inputAreaScaleFactor=1.0):


    try:
        rawFrame=inSHPFile(inRegionFile,inDataDict)
        if rawFrame is None:
            sys.exit(0)
        # For debugging, reduce file membership
        if regionSubList is not None:
            subLen=len(regionSubList[0])
            rawFrame=rawFrame[rawFrame['rawRegionCode'].str[:subLen].isin(regionSubList)]
            print('  For debugging',len(rawFrame),'records retained.')
            print('    To modify, remove or change argument *rawSubList* in the CreateList line.')
    except:
        print(regionLevel,'shapefile input failed.')
        print('Error exit ...')
        sys.exit(0)

    # At this stage, we don't know the region codes yet'
    if regionLevel=='VTD':
        rawFrame['regionID']=[('000'+str(rawFrame.loc[xcode,'CountyRawCode']))[-3:]+('00000'+str(rawFrame.loc[xcode,'MuniRawCode']))[-5:]+('00000'+str(rawFrame.loc[xcode,'VTDRawCode']))[-4:] for xcode in rawFrame.index]
        regionCode='V'
    elif regionLevel=='municipal':
        rawFrame['regionID']=[('000'+str(rawFrame.loc[xcode,'CountyRawCode']))[-3:]+('00000'+str(rawFrame.loc[xcode,'MuniRawCode']))[-5:] for xcode in rawFrame.index]
        regionCode='M'
    elif regionLevel=='county':
        rawFrame['regionID']=[('000'+str(rawFrame.loc[xcode,'CountyRawCode']))[-3:] for xcode in rawFrame.index]
        regionCode='C'
    else:
        print('Error: invalid region level:',regionLevel)
        print('Exiting...')
        sys.exit(0)

    # At this point, the 'regionID' variables have to be unique. Test this.
    IDCheck=sorted(list(rawFrame.regionID))
    badID=set()
    for seq in range(1,len(IDCheck)):
        if IDCheck[seq]==IDCheck[seq-1]:
            badID.update({IDCheck[seq]})
    if len(badID)>0:
        print('Error: duplicate IDs found')
        rawFrame.reset_index(inplace=True)
        badRows=[x for x in range(len(rawFrame)) if rawFrame.loc[x,'regionID'] in badID]
        badFrame=rawFrame.loc[badRows,['regionID','rawRegionCode','CountyRawCode','MuniRawCode','VTDRawCode','Size']]
        print(badFrame)
        print('Terminating ...')
        sys.exit(0)

    # From census, total area  in square kilometers
    rawFrame['AreaCensus']=(rawFrame['ALAND']+rawFrame['AWATER'])*inputAreaScaleFactor


    #Create memberList
    # Input initial adjacency and data for the regions
    rawFrame=Initial_Data_Input(rawFrame)

    # Document area variance: computed vs Census
    # Use this as evidence for meaningful perimeter calculation
    rawFrame['area_Variance']=rawFrame['Area']-rawFrame['AreaCensus']

    xFrame=rawFrame[['Name','AreaCensus','area_Variance']]
    yFrame=regionCode+rawFrame.regionID
    outFrame=pd.concat([xFrame,yFrame],axis=1)
    outFrame.set_index('regionID',inplace=True)
    try:
        outFrame.to_csv(G.thisPath+'VTD_areas.csv',mode='w')
        print('  Area variance file created: ',G.thisPath+regionLevel+'_areas.csv')
    except:
        print('  Area variance file action failed.')

    memberList=DF_to_memberList(rawFrame,regionCode)


    return True,memberList

def CreateList(thisLevel,rawSubList=None,Exit=True):

    print()
    results=True

    if thisLevel=='state':
        # Note: in this implementation, the stateMemberList is created from the countyMemberList
        results,outList=CreateStateList()

    if thisLevel=='county':
        results,outList=CreateRegionList(F.CountyShapeFile,F.CountyVariableDict,rawSubList,regionLevel=thisLevel,
                                    inputAreaScaleFactor=F.CountyinputAreaScaleFactor)

    if thisLevel=='municipal':
        results,outList=CreateRegionList(F.MunicipalShapeFile,F.MunicipalVariableDict,rawSubList,regionLevel=thisLevel,
                                    inputAreaScaleFactor=F.MunicipalinputAreaScaleFactor)
    
    if thisLevel=='VTD':
        results,outList=CreateRegionList(F.VTDShapeFile,F.VTDVariableDict,rawSubList,regionLevel=thisLevel,
                                 inputAreaScaleFactor=F.VTDinputAreaScaleFactor)

    if not results:
        print('Error: creation of',thisLevel,'MemberList failed.')
        print('Exiting ...')
        sys.exit(0)

    # Initial checks
    print('Initial input checks: ',end='')
    results=True
    results=results and CCFD_checks.CheckGeometries(outList)
    results=results and CCFD_checks.CheckRegionNeighbors(outList)
    results=results and CCFD_checks.CheckNeighborConsistency(outList)
    if rawSubList is None:
        results=results and CCFD_checks.CheckPopTotal(outList)

    if not results:
        print('Failed.')
        print('Error: Assembly of ',thisLevel,'failed at least one check.')

        holdVerbose=copy(G.verbose)
        G.verbose=True
        GEOresults=CCFD_checks.CheckGeometries(outList)
        NBRresults=CCFD_checks.CheckRegionNeighbors(outList)
        NBR2results=CCFD_checks.CheckNeighborConsistency(outList)
        if rawSubList is None:
            POPresults=CCFD_checks.CheckPopTotal(outList)
        G.verbose=copy(holdVerbose)

        results=GEOresults and NBRresults and NBR2results
        results=results and (POPresults or not G.populationCheckFlag)
    else:
        print('OK')

    # Check for termination
    if results:
        JSONresults=CCFD_files.dictToJson(outList,G.thisPath,thisLevel+'MemberList_Initial'+G.dataInitID)        
        if not JSONresults:
            print('Error: JSON file for'+regionLevel+MemberList_Initial+G.dataInitID+' not created.')
 
    results=results and JSONresults
    if Exit or not results:
        print('Exiting ...')
        sys.exit(0)

    return outList

def CreateBoundaryDict(inSelectionList=['Municipal'],regionSubList=None,checkOverride=False):

    if 'Municipal' in inSelectionList:
        if len(G.municipalMemberList)==0:
            try:
                # County list is needed to create municipal list
                fullMemberList=CCFD_files.JsonToList(G.thisPath+'municipalMemberList_Initial'+G.dataInitID)
                if regionSubList is not None:
                    subLen=len(regionSubList[0])
                    thisMemberList={x:fullMemberList[x] for x in fullMemberList.keys() if x[:subLen] in regionsSubList}
                    # Add all of the Neighbors
                    addNbrs=set()
                    for thisReg in thisMemberList.keys():
                        addNbrs.update(thisMemberList[thisReg]['Neighbors'])
                    thisMemberList.update({x:fullMemberList[x] for x in addNbrs})
                    thisMemberList['Outside']=deepcopy(G.regionMemberSet)
                    thisMemberList['Outside']['geometry']=CCFD_shapes.AddOutsideShape([thisMemberList[x]['geometry'] for x in thisMemberList.keys()])
                else:
                    thisMemberList=fullMemberList
            except:
                print('Municipal MemberList needed and not found.')
                print('Error exit ...')
                sys.exit(0)

            theseBoundaries=makeBoundaryLengthDict(thisMemberList,calcRegionList=countySubList)
            G.boundaryLengthDict.update(theseBoundaries)

        # Check the Boundaries

        results=CCFD_checks.checkRegionBoundaries(thisMemberList,G.boundaryLengthDict)
        if not results and not checkOverride:
            print('Boundary check failed.')
        else:
            print('  Boundary checks completed.')
            try:
                saveFile=G.thisPath+'selectedBoundaryLengths'+G.dataInitID
                json.dump(theseBoundaries,open(saveFile+'.json','w'))
                print('Json file: '+saveFile+' created.')
            except:
                print('Failed to save:municipalBoundarLengths')
                results=False
            try:
                saveFile=G.thisPath+'boundaryLengthDict'+G.dataInitID
                json.dump(G.boundaryLengthDict,open(saveFile+'.json','w'))
                print('Json file: '+saveFile+' created.')
            except:
                print('Failed to save: BoundarLengths')
                results=False

    return results