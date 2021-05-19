# SubRegion Handler
import operator
import numpy as np
from math import atan2, pi
from statistics import mean
from copy import deepcopy, copy
import random as Rnd

from shapely import geometry
from shapely.geometry import Polygon, Point, LineString
from shapely.ops import unary_union

import Global_Settings as G
import Shape_Handler as CCFD_shapes
import CCFD_Functions as CCFD
import MemberList_Creators as CCFD_inputs
import Region_Checking_Functions as CCFD_checks


# From a collection of regions, returns the set of regions with neighbors not in the collection
# Inputs: dictionary of region:its neighbors
def BoundaryMembers(inDict):
    # inDict is a dictionary regionID:its neighbors
    # A Boundary of a region has at least one neighbor not a member of itself
    members=set(inDict.keys())
    return {x for x in members if len(inDict[x]-members)>0}

# The fake region "Outside" exists only as a neighbor of regions on the state-border
# This routine strips off those neighbors and makes a list of regions which had them
def StripOutside(inDict):

    outSideList=[]
    for x in inDict.keys():
        if 'Outside' in inDict[x]['Neighbors']:
            outSideList+=[x]
            inDict[x]['Neighbors'].discard('Outside')

    return inDict,outSideList

def getElement(inSet):
    # Implemented to allow the Districts set in regionMemberList to contain more than one district type
    # Currently, this function assumes inSet has at most one element. If not, it returns the first one.
    if inSet is None:
        return None
    elif len(inSet)==0:
        return set()
    else:
        return next(iter(inSet))
 
def OrderedEdges(inDict,first=None,endMember=None):
# inDict is a regionMemberList
# Creates a list of edge members starting with member "start"

    # Clockwise direction in range 0 to 360 (where 0 the positive y direction)
    def ClockwiseAngle(delta_x,delta_y):
        totradians=2.0*pi
        inradians=(atan2(delta_x,delta_y)+ totradians) % totradians
        return inradians*(360.0/totradians)

    fullEdgeSet=BoundaryMembers({x:inDict[x]['Neighbors'] for x in inDict.keys()})
    edgeSet=set(fullEdgeSet)
    

    # Start with a member of the edgelist
    if first is None:
        first=list(inDict.keys())[0]
    if endMember is None: endMember=first

    if first not in inDict.keys() or endMember not in inDict.keys():
        print('Error: region  not found in input dictionary: ',first,',',endMember)
        return None,None

    start=first
    if G.verbose: print('Initial region: ',start)
    outList=[]

    # Find the neighboring district farthest to the right of the start district
    # Start with a vector pointed from the start muni using its neighbors
    # If the start muni is on a convex area, the vector will point inward and the function will list the edges in clockwise order.
    # If the start muni is on a concave area, the vector will point outward and the function wil list the edges in counter-clockwise order.
    avgLongitude=mean([inDict[x]['Data']['Longitude'] for x in inDict.keys() if x in inDict[first]['Neighbors'] & fullEdgeSet])
    avgLatitude=mean([inDict[x]['Data']['Latitude'] for x in inDict.keys() if x in inDict[first]['Neighbors'] & fullEdgeSet])

    toCenter=ClockwiseAngle((avgLongitude-inDict[start]['Data']['Longitude']),(avgLatitude-inDict[start]['Data']['Latitude']))
    testdirs={x:ClockwiseAngle((inDict[x]['Data']['Longitude']-inDict[start]['Data']['Longitude']),(inDict[x]['Data']['Latitude']-inDict[start]['Data']['Latitude'])) for x in inDict[start]['Neighbors'] & edgeSet}
    testangles={x:(testdirs[x]-toCenter+360.0) % 360.0 for x in testdirs.keys()}

    lastMember=[x[0] for x in sorted(testangles.items(),key=lambda y:y[1])][0]  
    
    thisMember=''
    enclosedList=[]
    
    # Find the next member in clockwise direction
    # It is a neighbor of the start region.
    # The first clockwise neighbor is the one with the sharpest left turn relative to the ray pointing inward from the current start region
    while thisMember!=endMember:
        thisMember=''

        nextgroup=inDict[start]['Neighbors'] & edgeSet
        nextgroup=nextgroup - set(outList)
        nextgroup.discard(lastMember)
        nextgroup=list(nextgroup)
        # This should alway be valid. If not, invalid: start was a neighbor of lastMember but lastMember is not a neighbor of start.

        if len(nextgroup)==1:    # This is the typical case
            thisMember=nextgroup[0]
            outList.append(thisMember)
            lastMember=start
            start=thisMember

        elif len(nextgroup)>1:
            nextset=set(nextgroup)
            preselection={x:inDict[x]['Neighbors'] & edgeSet for x in nextgroup}
            # Find any enclosed.

            iter=list(preselection.keys())
            for x in iter:
                if len(preselection[x])==1:
                    # This region is enclosed by thisMember. Add it to the enclosed list.
                    enclosedList.append(x)
                    edgeSet.remove(x)
                    # Then delete it from preselection
                    del preselection[x]

            selection={x:len(preselection[x] - nextset) for x in preselection.keys()}

            if len(selection)==0:
                print('Logic error: No selections')
                return None,None

            thisMember=[x[0] for x in sorted(selection.items(),key=lambda y:y[1])][0]  
            outList.append(thisMember)
            lastMember=start
            start=thisMember

        else:
            # This will happen if the first region is enclosed.
            # Use its enclosing region as the first region
            if len(outList)==0:
                if G.verbose: print('  Initial region reset to: ',lastMember)
                start=lastMember
                first=lastMember
                endMember=first

                avgLongitude=mean([inDict[x]['Data']['Longitude'] for x in inDict.keys() if x in inDict[first]['Neighbors'] & fullEdgeSet])
                avgLatitude=mean([inDict[x]['Data']['Latitude'] for x in inDict.keys() if x in inDict[first]['Neighbors'] & fullEdgeSet])

                toCenter=ClockwiseAngle((avgLongitude-inDict[start]['Data']['Longitude']),(avgLatitude-inDict[start]['Data']['Latitude']))

                testdirs={x:ClockwiseAngle((inDict[x]['Data']['Longitude']-inDict[start]['Data']['Longitude']),(inDict[x]['Data']['Latitude']-inDict[start]['Data']['Latitude'])) for x in inDict[start]['Neighbors'] & edgeSet}
                testangles={x:(testdirs[x]-toCenter+360.0) % 360.0 for x in testdirs.keys()}
                # The below selects for clockwise ordering at a convex point, and counterclockwise ordering at a concave point
                lastMember=[x[0] for x in sorted(testangles.items(),key=lambda y:y[1])][0]  
    
                outList=[]
                thisMember=''
                edgeSet.discard(start)
                start=first
            else:
                # The current region is a dead end.
                # Add the current region to the enclosed list, and delete it from edgeSet
                # The back up one region
                if G.verbose: print('  Enclosed: ',start)
                enclosedList.append(start)
                edgeSet.discard(start)
                outList.pop()
                start=outList[-1]
                lastMember=''

    outList.insert(0,first)
    # Make the enclosed list all in edgeSet not in outList
    skippedList=list(fullEdgeSet-set(outList))

    return outList,skippedList

def SubregionAssign(inLongLatTuple,geometryDictionary):
    # This routine places a point tuple (Longitude, Latitude) in the selection of geometries in the input dictionary {(regionID,geometry),...}
    # It returns the first regionID for which the point lies within its geometry

    try:
        regionList=[x for x in geometryDictionary.keys() if geometryDictionary[x].contains(Point(inLongLatTuple))]
    except:
        return None
    if len(regionList)==0:
        return None

    return regionList[0]

def MakeEnclosedRegions(regionPrefix, inDict=None):

    if G.verbose: print('Enclosure Substitutions in: ',regionPrefix)

    if inDict is None:
        if regionPrefix[0]=='M':
            inDict=G.municipalMemberList
        else:
            print('Exception: Unknown region prefix:',regionPrefix)

    else:
        # shortDict is inDict
        shortDict={x:inDict[x] for x in inDict.keys() if x[:len(regionPrefix)]==regionPrefix}

    if len(shortDict)==0:
        print('Error: no regions found with prefix:',regionPrefix)
        return None


#        nbrDict={x:inDict[x] for xreg in shortDict.keys() for x in shortDict[xreg]['Neighbors'] if x!='Outside'}
#        shortDict.update(nbrDict)


    # At the subregion level, there may be multiple entirely enclosed subregions.
    # For those fully embedded in a larger one, we must carry these forward in the analysis together with their enclosers
    # Substitute a group for any subregion that fully encloses any subregions.
    # The group will consist of the enclosing subregion and all of those it encloses.
    # Test 1: Most common
    # The enclosing subregion is any subregion for which it is the only neighbor of another subregion
    groupEnclosers1={list(shortDict[x]['Neighbors'])[0] for x in shortDict.keys() if len(shortDict[x]['Neighbors'])==1 and x[:4]==regionPrefix}
    # Test 2: Enclosing region has interior holes
    groupEnclosers={x for x in shortDict.keys() if len(shortDict[x]['geometry'].interiors)>0 and x[:4]==regionPrefix}

    if not groupEnclosers1.issubset(groupEnclosers):
        if G.verbose: print('  Warning: encloser tests yield inconsistent results:')
        xproblem=groupEnclosers1-groupEnclosers
        if G.verbose: print('  Regions',xproblem,' do not have holes but are enclosers. May be boundary enclosers.')
        groupEnclosers.union(groupEnclosers1)

    # Build a list group ID dictionary from this set, with this encloser as its first member
    groupSeq=0
    groupList={}
    enclosureDict={}
    encloserList=list(groupEnclosers)
    adjoiningEnclosers={x:set() for x in encloserList}

#    encloserList=['M01754192A.0','M01754184.0']
#    adjoiningEnclosers['M01754192A.0'].update({'M01754184.0'})
#    Above for debugging
    while len(encloserList)>0:
        while True:
            thisRegion=encloserList.pop(0)
            #Capture error ...
            if thisRegion[:4]!=regionPrefix:
                print('Error: external region:',thisRegion)
                print('')

            
        # Run through all the munis and assign the appropriate enclosed ones to this group.
            # First, collect all of the in-county neighbors of thisRegion

            if thisRegion in ['Mcccmmmmm.0']:
# Special: we need to process these in a particular order and the logic below sometimes doesn't get us there

                print('Focus:',thisRegion)
            xNeighbors={x for x in shortDict.keys() if thisRegion in shortDict[x]['Neighbors'] if x[:4]==regionPrefix}
            # the order of processing enclosers matters.
            # If any of these regions are enclosers, we need to group them before we group their enclosers.
            # We only proceed if the regions in-county neighbors have no holes.
            # Side-by-side regions with holes frustrate simple hole-counting. The dictionary adjoiningEnclosers deals with this issue.
            #   Nested hole-containing regions must be monitored carefully.
            enclNbrs={x for x in xNeighbors if len(shortDict[x]['geometry'].interiors)>0 and x not in adjoiningEnclosers[thisRegion]}
            if len(enclNbrs)>0:
                if G.verbose: print('    Reordering: adjoining or nested hole regions detected',thisRegion,enclNbrs)

                encloserList.append(thisRegion)
                for x in enclNbrs:
                    adjoiningEnclosers[x].update({thisRegion})
                    # Put these near the front of the list
                    if x in encloserList:
                        encloserList.remove(x)
                    encloserList.insert(0,x)


            else:
                # We can process this encloser now.
                # Members must remain in this region (eg. county). We exclude those already flagged as adjoining
                # Group members are included if they are in the same region as thisRegion
                newMembersx={x for x in xNeighbors-adjoiningEnclosers[thisRegion] if x[:4]==regionPrefix}
                # newMembersx may include neighbors of thisRegion that are external to thisRegion
                # Step 1 of eliminating them is to restrict members to those show share neighbors with each other.
                # This is an if condition but not an only if condition. Neighbors of thisRegion enveloped by other neighbors of thisRegion also pass this test.
                newMembers1={x for x in newMembersx if shortDict[x]['Neighbors']-newMembersx=={thisRegion}}
                # newMembers1 has reduced the set somewhat, but still may include external neighbors not yet identified as enclosers themselves.
                # Internal members also may be enclosers themselves. So, include all neighboring enclosers, whether internal or external, here.
                newMembers2={x for x in newMembersx if len(shortDict[x]['geometry'].interiors)>0}
                newMembers=newMembers1|newMembers2

                # The set may still include some regions that are not truly enclosed. We need to impose the "only if" condition
                tryGroup=newMembers|{thisRegion}
                newMembers={x for x in copy(newMembers) if len(shortDict[x]['Neighbors']-tryGroup)==0}

                if len(newMembers)>0:
                    break
                else:
                    # This can happen if the region contains interior holes that do not correspond to municipalities.
                    # This can happen with nested hole-containing regions.
                    # Rememdy at this point by reversing the adjoiningEnclosers assignments
                    if thisRegion not in adjoiningEnclosers.keys():
                        print('  Problem. Stop.')
                    xNested=list(adjoiningEnclosers[thisRegion])[0]
                    if G.verbose:
                        print('    Identified adjoining enclosing regions:',xNested,'and',thisRegion)
                    if xNested not in adjoiningEnclosers.keys(): adjoiningEnclosers[xNested]=set()
                    adjoiningEnclosers[xNested].update({thisRegion})
                    # Check to see if any enclosers in this region are really totally enclosed in this region.
                    # If so, they are not adjoining
                    for x in copy(adjoiningEnclosers[thisRegion]):
                        if shortDict[x]['Neighbors']=={thisRegion}:
                            adjoiningEnclosers[thisRegion].discard(xNested)
                    if len(encloserList)==0:
                        print('  Problem. Stop.')

        groupSeq+=1

        # Create a new regionSet for this. Here, has form MmmmGgggg.0
        regionID=regionPrefix+'G'+str(groupSeq).zfill(4)+'.0'
        groupList[regionID]=deepcopy(shortDict[thisRegion])
        enclosureDict[thisRegion]=regionID
                    

        if thisRegion in ['Mcccmmmmm.0']:
            print('Focus:',thisRegion)
        # If this is successful, we should be able verify that the areas of the interior holes is equal to the area of the enclosed regions
        # This comparison is complicated a little bit because some enclosed regions sit on the region boundary or are neighbors of  boundary region
        # Create the set of bounding Members
        # The set may still include some regions that are not truly enclosed. We need to impose the "only if" condition
        tryGroup=newMembers|{thisRegion}
        newBoundingMembers=BoundaryMembers({x:shortDict[x]['Neighbors'] for x in tryGroup})
        newBoundingMembers.discard(thisRegion)
        testBounds=copy(newBoundingMembers)
        testMembers=list(newMembers)
        while len(testMembers)>0:
            thisMember=testMembers.pop()
            # Neighbor of Simple Bounding
            if thisMember in shortDict.keys():
                if len(shortDict[thisMember]['Neighbors'] & testBounds)>0:
                    newBoundingMembers.update({thisMember})
            elif thisMember[4]!='G':
                # State Bounding
                if len({x for x in shortDict[thisMember]['Neighbors'] if x=='Outside'})>0:
                    newBoundingMembers.update({thisMember})
            else:
                # State-Bounding bounding Group members.
                testMembers+=shortDict[thisMember]['Group']

        holeArea=0.0
        for thishole in shortDict[thisRegion]['geometry'].interiors:
            holeArea+=Polygon(thishole).area
        for thishole in newBoundingMembers:
            holeArea+=shortDict[thishole]['geometry'].area
        regionArea=0.0
        for thisMember in newMembers:
            regionArea+=shortDict[thisMember]['geometry'].area
        if abs(regionArea-holeArea)*2324800.0>0.1:   # Area difference larger than about one-tenth acre
            print('  Hole/Region area inconsistency',round((regionArea-holeArea)*2324800.0,2),thisRegion,newMembers)

    # We've exited the while loop.
    # Update info and Aggregate the graph info, the geometry, and the data for each group
        # the graph info
        x=groupList[regionID]
        x['Group'].update(newMembers)
        x['ID']=regionID
        x['Data']['ID']=x['ID']
        x['Neighbors']-=newMembers
        x['Name']=x['Name']+' and Enclosed'
        # the geometry
        if len(x['Group'])!=len({xtry for xtry in x['Group'] if xtry in shortDict.keys()}):
            print('Error: Missing info for entry in group:',regionID,x['Group'])
        x['geometry']=unary_union([shortDict[xgrp]['geometry'] for xgrp in x['Group']])            
        # the data
        x=CCFD_shapes.calcRegionGeoData(x)
        x['Data']['TAPersons']=int(sum([shortDict[xreg]['Data']['TAPersons'] for xreg in x['Group']]))



        # Add this new group to the working dictionary
        shortDict[regionID]=x
        if G.verbose: print('  Group',regionID,' encloser created:',x['Group'])

        # Fix the neighbors lists for all affected munis (substitute the group for the encloser)
        # In this case, we are just dropping the original enclosed muni ID and substituting the Group ID
        # One complication. Some neighbors may not be in the target county. To remedy this, we add any non-County neighbors to munis for future substitution into regionMemberList
        # Fix the list of neighbors in the neighboring regions
        for thisMember in x['Group']:
            
            for nbr in shortDict[thisMember]['Neighbors']:
                # We only need to consider external neighbors
                if nbr not in x['Group'] and nbr!='Outside':
                    if nbr not in shortDict.keys():
                        if G.verbose: print('    Adding region',nbr,' to shortDict')
                        shortDict[nbr]=deepcopy(inDict[nbr])
                    # Add group to this region's neighbors
                    shortDict[nbr]['Neighbors'].update({regionID})
                    # Delete thisMember as appropriate
                    shortDict[nbr]['Neighbors'].discard(thisMember)

                    # Add this region to the neighbors of x
                    # x['Neighbors'].update({nbr})    ?

        # Remove munis that are now parts of the group
        dropRegions= {xreg for xreg in x['Group'] if xreg in shortDict.keys()}
        while len(dropRegions)>0:
            xreg=dropRegions.pop()
            # if xreg happens to be a group, move its members to the new group list before deleting
            if xreg[:5]==regionPrefix+'G':
                x['Group'].update(shortDict[xreg]['Group'])
                dropRegions.update(shortDict[xreg]['Group'])
            if xreg in shortDict.keys():
                del shortDict[xreg]

        # After creating a group, need to check all remaining munis to see if it is enclosed itself.
        # First remove its elements from adjoiningEnclosers. They no longer apply.
        for x in shortDict[regionID]['Group']:
            if x in adjoiningEnclosers.keys():
                if G.verbose: print('    Removing from adjoiningEnclosers:',x)
                del adjoiningEnclosers[x]

        #    (Note: it won't be an encloser because of the nature of enclosing groups.)
        # If so, need to add it to the list of enclosers so it will be checked again.
        # By the nature of groups, this group should have no holes
        if len(shortDict[regionID]['geometry'].interiors)>0:
            print('  Error: group has holes',regionID)
        # Check all of the neighbors
        inCountyNeighbors={x for x in shortDict[regionID]['Neighbors'] if x[1:4]==regionPrefix}
        # Is this group and encloser?
        if len({x for x in inCountyNeighbors if len(shortDict[x]['Neighbors'])==1})>0:
            if G.verbose: print('    Group',regionID,'encloses',list(shortDict[regionID]['Neighbors'])[0],'. Add to encloserList.')
            print('  Problem: enclosing groups should not be enclosers:')
            encloserList.append(regionID)

        # Is this group enclosed?
        if len(shortDict[regionID]['Neighbors'])==1:
            xencloser=list(shortDict[regionID]['Neighbors'])[0]
            if G.verbose: print('    Region',xencloser,'encloses',regionID,end='')
            # Make sure this encloser is in the encloser List -- group regionID will eventually be added.
            if xencloser in encloserList:
                if G.verbose: print(' is already in encloserList.')
            else:
                encloserList.append(xencloser)
                if G.verbose: print(' added to encloserList.')

        # Any members of this group in encloserList?

        # Finally, we need to update adjoiningEnclosers, both the keys and the values
        for thisReg in shortDict[regionID]['Group']:
            for thisAdj in adjoiningEnclosers.keys():
                if thisReg in adjoiningEnclosers[thisAdj]:
                    adjoiningEnclosers[thisAdj].remove(thisReg)
                    adjoiningEnclosers[thisAdj].update({regionID})
        for thisAdj in list(adjoiningEnclosers.keys()):
            if thisAdj in shortDict[regionID]['Group']:
                adjoiningEnclosers[regionID]=copy(adjoiningEnclosers[thisAdj])
                del adjoiningEnclosers[thisAdj]

    return shortDict

def MakeConnectedRegions(regionPrefix,dataDict=None):


    if G.verbose: print('Connector Substitutions in: ',regionPrefix)

    shortDict={}
    if dataDict is None:
        # Use the regionID info to build this
        if regionPrefix[0]=='M':
            dataDict=G.municipalMemberList
            # shortDict must include all regions with prefix "regionPrefix" and all of the regions listed as neighbors of these regions
        shortDict={x:dataDict[x] for x in dataDict.keys() if dataDict[x]['ID'][:len(regionPrefix)]==regionPrefix}
        for xreg in shortDict.keys():
            for nbr in shortDict[xreg]['Neighbors']:
                shortDict[nbr]=dataDict[nbr]
    else:
       shortDict=deepcopy(dataDict)

    if shortDict is None or len(shortDict)==0:
        print('Error: no subregions entered:',regionPrefix)
        return None

    # Upon creation of the initial database, split subregions are given an alpha sequence just to the left of the decimal point.
    splittedSubregions={x for x in shortDict.keys() if '0123456789'.find(x[x.find('.')-1])==-1 and x[:4]==regionPrefix}

    # There may be some splitted munis pieces that have already been grouped. Include their groups as well.
    xgroups={x:xgrp for xgrp in shortDict.keys() if xgrp[:5]==regionPrefix+'G' for x in shortDict[xgrp]['Group']}

    # for reference purposes, remember the grouped members and their groups
    
    
    splittedSubregions.update({x for x in xgroups.keys() if '0123456789'.find(x[x.find('.')-1])==-1})

    # For iteration purposes, we need to know what group sequence we're at at this point. Find it by looking through the input dictionary shortDict
    grpList=[int(x[5:9]) for x in shortDict.keys() if x[:5]==regionPrefix+'G']
    if len(grpList)>0:
        groupSeq=max(grpList)
    else:
        groupSeq=0

    # Collect these into the desired groups.
    splitDict={}
    rootNameDict={}
    for thisSubregion in splittedSubregions:
        xroot=thisSubregion[:thisSubregion.find('.')-1]

        if xroot not in splitDict.keys():
            splitDict[xroot]=set()

        # if the muni x is already grouped, use the group instead of the muni
        if thisSubregion in xgroups.keys():
            splitDict[xroot].update({xgroups[thisSubregion]})
            rootNameDict[xroot]=shortDict[xgroups[thisSubregion]]['Name']
        else:
            splitDict[xroot].update({thisSubregion})
            rootNameDict[xroot]=shortDict[thisSubregion]['Name'][:-5]

    # Create a new group for each root:
    for thisRoot in splitDict.keys():

        # Add for debugging purpose here.
        # watchSet={'M00356384','M04518160'}
        watchSet={'M091G0001'}
        if thisRoot in watchSet:
            G.verbose=True

        # Create a new regionSet for this. Here, has form MmmmGgggg.0
        groupSeq+=1
        regionID=thisRoot[:4]+'G'+str(groupSeq).zfill(4)+'.0'
        if G.verbose: print('Making ',regionID,splitDict[thisRoot])
        # Start the region info for this group. It doesn't matter which element to start with
        xgrp=deepcopy(G.regionMemberSet)
        # Set up initial info
        xgrp['ID']=regionID
        xgrp['Data']['ID']=xgrp['ID']

        xgrp['Group'].update(splitDict[thisRoot])
        # It is possible that a member of this group has already been consolidated into another group.
        # This would have happened in an earlier step within the regions that share "regionPrefix"
        # Make the requisite substitutions, the member's group instead of the member
        # All substitutions of any kind.
        for thisReg in copy(xgrp['Group']):
            if thisReg not in shortDict.keys():
                # This member has already been consolidated. Use the consolidated region instead.
                if thisReg not in xgroups.keys():
                    print('  Error: region not in shortDict or xgroups',thisReg)
                # Make the substitution
                xgrp['Group'].update({xgroups[thisReg]})
                xgrp['Group'].remove(thisReg)

                # Add to xgroups
                xgroups[thisReg]=regionID
                if thisReg[4]=='G':
                    # This is a group. We need to fix the reference of its members in xgroups
                    for yreg in xgroups.keys():
                        if xgroups[yreg]==thisReg: xgroups[yreg]=regionID


        # Create the neighbors for xgrp. All members should now be in shortDict
        # They are simply all the neighbors of the group members less any regions that are members of xgrp
        # Remove any references to the group members, substituting the group ID instead.
        for thisReg in xgrp['Group']:
            xgrp['Neighbors'].update(shortDict[thisReg]['Neighbors'])
        xgrp['Neighbors']-=xgrp['Group']

        # Select a name based on the root
        if len(rootNameDict[thisRoot])>0:
            xgrp['Name']=rootNameDict[thisRoot]+' and Connectors'
        else:
            print('  Problem with name.')


        # Find common neighbors for a dict of element:set of neighbors
        def CommonNeighbors(inDict):
            xlist=list(inDict.keys())
            neighbors=set(inDict[xlist[0]])
            for i in range(1,len(xlist)):
                neighbors=neighbors & set(inDict[xlist[i]])
            return neighbors

               
        # Add subregions until we get a fully connected group
        while True:
            # Check the connectivity of this set.

            # tryConnected is a list of lists of connected subregions.
            # We add neighboring subregions to this group until this falls to one, at which point all of the split segments are connected.
            theseConnectors={}
            for xreg in xgrp['Group']:
                if xreg in shortDict.keys():
                    theseConnectors[xreg]=shortDict[xreg]['Neighbors'] & xgrp['Group']
                elif xreg in xgroups.keys():
                    if xgroups[xreg] in shortDict.keys():
                        theseConnectors[xreg]=shortDict[xgroups[xreg]]['Neighbors'] & xgrp['Group']
                else:
                    print('  Unknown group member:',xreg,' of group',regionID)
            tryConnected=CCFD.getConnected(xgrp['Group'],theseConnectors)
            # tryConnected is a list of unconnected lists of connected subregions.

            chkRegions={'Mcccmmmmm.0','Mcccmmmmm.0'}
            if len(chkRegions & xgrp['Group'])>0:
                print('Checking:',chkRegions)
            if G.verbose:
                print('tryConnected',tryConnected)
                plotSet={}
                for x in range(len(tryConnected)):
                    plotSet['Piece'+str(x)]={y:shortDict[y]['geometry'] for y in tryConnected[x]}
                CCFD_shapes.GroupPlot(plotSet)
            if len(tryConnected)==1:
                if G.verbose: print('  Group',regionID,' connector created:',xgrp['Group'])
                break

            # Create a little dictionary with a temporary set name for each internally-connected list in the bigger disconnected list.
            tempDict={}
            for seq in range(len(tryConnected)):
                # For each connected group, make a list of in-county neighbors (to remain compliant with Rule 2)
                tempDict[seq]={z for xreg in tryConnected[seq] for z in shortDict[xreg]['Neighbors'] if z[:4]==regionPrefix}

            # Connect the pieces with one or more connecting neighbors
            xRemaining=set(tempDict.keys())
            thesePieces=copy(xRemaining)
            tryCombos=1
            counter=1
            pcounter=1
            alreadyTried={}
            keepNumber=len(xRemaining)

            # Below we look for a subregion that has at least two of the split regions as a neighbor.
            # Such a common neighbor would connect them if it is included in the group.
            # We start by looking for a common neighbor of all of them. If none exists, try smaller combinations of the split subregions
            # (For example, if there are 4 split subregions of this root, try all 4, then try each of the combinations of 3 of them, then combinations of two, etc., until a common region is found.)
            while True:
                if G.verbose: print('keepNumber',keepNumber)
                if keepNumber<2:
                    if G.verbose: print('No in-region connector found for:',xgrp['Group'])

                    groupList=list(xgrp['Group'])
                    # Remedy attempt. Add neighbor of one piece that a line between the centroids intersects.
                    thisRay=LineString([shortDict[groupList[0]]['geometry'].centroid,shortDict[groupList[1]]['geometry'].centroid])
                    countyNeighbors={xreg for xreg in xgrp['Neighbors'] if xreg[:4]==regionPrefix}
                    tryNeighbors=set()
                    for xreg in countyNeighbors:
                        if thisRay.intersects(shortDict[xreg]['geometry']):
                            tryNeighbors.update({xreg})

                    # tryNeighbors={xreg for xreg in countyNeighbors if thisRay.intersects(shortDict[xreg]['geometry'])}

                    # if tryNeighbors is empty, we've failed.
                    if len(tryNeighbors)==0:
                        if G.verbose: print('Deleting smaller member as outside the current region. Proceeding with larger only.')
                        # CCFD_shapes.QuickPlot({xreg:shortDict[xreg]['geometry'] for xreg in xgrp['Group']})
                        xDict={xreg:shortDict[xreg]['Data']['Area'] for xreg in xgrp['Group']}
                        tryNeighbors={max(xDict,key=xDict.get)}
                        for dropped in {xreg for xreg in xgrp['Group'] if xreg not in tryNeighbors}:
                            del shortDict[dropped]
                            xgrp['Group'].remove(dropped)
                        tryNeighbors=set()
                        break

                    else:
                        # Pick one of these Neighbors and add it to the group
                        pickNeighbor=Rnd.sample(tryNeighbors,1)[0]
                        if G.verbose: print('  Trying the addition of interveneing region: ',pickNeighbor)
                        tryNeighbors=[pickNeighbor]
                else:
 
                # look for a subregion that some of the split subregions have in common. If so, tryneighbors will have length>0.
                    while True:

                        tryNeighbors=CommonNeighbors({thisCollection:tempDict[thisCollection] for thisCollection in tempDict.keys() if thisCollection in thesePieces})
                        alreadyTried[pcounter]=thesePieces
                        pcounter+=1
                        if len(tryNeighbors)>0:
                            break
                        else:
                            # No common neighbor from any of these combos (xRemaining taking keepNumber without replacement)
                            # Keep track of how many combinations (must have at least two pieces to continue)
                            if len(alreadyTried)>=tryCombos:
                                keepNumber-=1
                                if keepNumber<2:
                                    break
                                tryCombos=CCFD.binom(len(xRemaining),keepNumber)
                            # We need to keep track of the combos we've already tried but only at the same tuple level.
                                counter=0
                                pcounter=1
                                alreadyTried={}
                                if G.verbose: print('    tried:',alreadyTried)
                            while True:
                                # Find an untried set at this combo level. We're choosing randomly, so we might often pick on already tried
                                counter+=1
                                thesePieces=set(Rnd.sample(xRemaining,keepNumber))
                                testPiece=False
                                for pSeq in range(1,pcounter):
                                    testPiece=testPiece or thesePieces==alreadyTried[pSeq]
                                if not testPiece:
                                    break
                                if counter>100:
                                    print('Counter limit',counter,' exceeded running MakeConnectedRegions.')

                # We've found a list of munis that connect some of the pieces, or the process has failed.
                # Pick the one with the smallest feature. We're using area here.
                if len(tryNeighbors)>0:

                    xDict={}
                    for xreg in tryNeighbors:
                        xDict[xreg]=shortDict[xreg]['Data']['Area']
                    pickNeighbor=min(xDict,key=xDict.get)    

                    # If necessary, add to shortList the neighbors of this new member
                    # If it's a group, its neighbors are already in shortDict
                    if pickNeighbor[4]!='G':
                        for x in dataDict[pickNeighbor]['Neighbors']:
                            if x not in shortDict.keys():
                                shortDict.update({x:dataDict[x]})
                    else:
                        if pickNeighbor not in shortDict.keys():
                            print('Error: Missing group in shortDict')

                    # Update the membership and neighbors for this group
                    # the graph
                    xgrp['Group'].update({pickNeighbor})
                    # The neighbors List
                    xgrp['Neighbors'].update(shortDict[pickNeighbor]['Neighbors'])
                    xgrp['Neighbors']-=xgrp['Group']         # Neighbors of group members do not include each other
                    break

        # We've got a fully connected group. Complete creating the MemberSet info for it.
        # the geometry
        if len(xgrp['Group'])<2:
            # This is not actually a group. No additional action. We don't need a connector.
            if G.verbose: print('  Not considered a group:',regionID,xgrp['Group'])
            # del shortDict[regionID]
            groupSeq-=1
        else:
            # Connecting success ...

            # Check for inadvertantly enclosed some regions inside this group
            if 'M01172824.0' in xgrp['Neighbors']:
                print()
            # Create a set of all neighbors of the newly assembled connected group (test: if group is a neighbor)
            testnbrs=set()
            for x in xgrp['Group']:
                testnbrs.update(shortDict[x]['Neighbors'])
            # If all neighbors of one of these neighbors are part of the group, it is enclosed. (test: if group is the only neighbor)
            xenclosed=set()
            for x in testnbrs:
                # Restrict to neighbors within the consolidating region (no boundary splits)
                if x[:4]==regionPrefix:
                    trial=shortDict[x]['Neighbors'] - xgrp['Group']
                    if len(trial)==0:
                        xenclosed.update({x})
            # Add any regions found to the group
            xgrp['Group'].update(xenclosed)
            # Adjust the neighbors list
            xgrp['Neighbors']-=xenclosed

            # Assemble the geomtry               
            unionList=[]
            for x in xgrp['Group']:
                unionList.append(shortDict[x]['geometry'])
            xgrp['geometry']=unary_union(unionList)
            if not CCFD_shapes.checkPolygon(xgrp['geometry']):
                print('Improper geometry:',xgrp['ID'],'Group:',xgrp['Group'])
                print({y:(shortDict[y]['Neighbors']&xgrp['Group']) for y in xgrp['Group']})
                if G.verbose:
                    CCFD_shapes.QuickPlot({xgrp['ID']:xgrp['geometry']})
            
            # Compute the data, many from the geometry
            xgrp=CCFD_shapes.calcRegionGeoData(xgrp)           
            xgrp['Data']['TAPersons']=sum([shortDict[xreg]['Data']['TAPersons'] for xreg in xgrp['Group']])

            # With all info collected, add this group to munis
            shortDict[regionID]=xgrp

            # Enter the new information and modified information in xgroups
            for thisReg in xgrp['Group']:
                xgroups[thisReg]=regionID
                # if thisReg is a group, then replace the individual region references to the new group from this one
                if thisReg[4]=='G':
                    for thisMember in xgroups.keys():
                        if xgroups[thisMember]==thisReg:
                            xgroups[thisMember]=regionID

            # in the neighboring regions
            # Any reference to a member of this group in neighbors' lists should be replaced with this group name
            xcleanup={nbr for x in xgrp['Group'] for nbr in shortDict[x]['Neighbors'] if nbr!='Outside'}
            # All of these regions should have the groupname as a neighbor and should have all group members eliminated as members
            for x in xcleanup:
                if x not in shortDict.keys() and x not in dataDict.keys():
                    print('Error: ',x,' not found in working dictionary.')
                    return None
                elif x not in shortDict.keys():
                        print('Why adding this so late?',x)
                        shortDict[x]=deepcopy(dataDict[x])

                # Members of the group are not Neighbors anymore, the group ID is
                # Fix neighbors lists
                for thisReg in xgrp['Group']:
                    if thisReg in shortDict[x]['Neighbors']:
                        shortDict[x]['Neighbors'].remove(thisReg)
 
                # Remove any members of this group from this neighbor
                shortDict[x]['Neighbors']-=xgrp['Group']
                # Add the group ID as a neighbor
                shortDict[x]['Neighbors'].update({regionID})

            # Remove the members of this group from shortDict and references to it
            for x in copy(xgrp['Group']):
                if x in shortDict.keys():
                    # Just before deleting, if x is a group, replace its members in xgrp['Group']
                    if x[4]=='G':
                        shortDict[regionID]['Group'].update(shortDict[x]['Group'])
                        shortDict[regionID]['Group'].remove(x)
                    del shortDict[x]
                elif x in xgroups.keys():
                    # The region has already been dropped. If so, it should be in xgroups
                    # Make sure it's group replacement is listed.
                    if xgroups[x] not in shortDict[regionID]['Group']:
                        # Make the substitutions
                        shortDict[regionID]['Group'].update(shortDict[xgroups[x]]['Group'])
                        shortDict[regionID]['Group'].remove(x)
                else:
                    print('  Warning: group member not in shortDict:',x)


    return shortDict,xgroups

def SubregionDistrictAssigner(regDict,districtGeometry=None):

    watchSet={'Mo03G0004.0'}
    if len(set(regDict.keys()) & watchSet)>0:
        G.verbose=True

#    G.verbose=True
    regionPrefix=list(regDict.keys())[0][:4]
    if districtGeometry is None:
        districtGeometry={x:x['geometry'] for x in G.districtMemberList.keys()}

    # Build a dictionary of subregions with District IDs that comprise this subcounty
    # Assign the these munis:

    print('  District Assignments: ')
    regionSubstitutes={}

    for thisregion in regDict.keys():
        thissubcounty=SubregionAssign((regDict[thisregion]['Data']['Longitude'],regDict[thisregion]['Data']['Latitude']),districtGeometry)

        if thissubcounty is None:
            print('Error in region assignment: ',thisregion)
            CCFD_shapes.QuickPlot(districtGeometry)
            CCFD_shapes.QuickPlot({**districtGeometry,**{thisregion:regDict[thisregion]['geometry']}},thisregion)

        # Add this muni to list of substitutes
        regionSubstitutes[thisregion]=regDict[thisregion]
        # Assign its district
        if 'Districts' not in regionSubstitutes[thisregion].keys():
            regionSubstitutes[thisregion]['Districts']={G.districtPrefix:None}
        regionSubstitutes[thisregion]['Districts'][G.districtPrefix]=thissubcounty

    # Check to make sure all of the county munis are assigned to a district.
    for x in regionSubstitutes.keys():
        if regionSubstitutes[x]['Districts'][G.districtPrefix] is None:
            print('  Assignment error: ',x,' has no district assignment.')

    # Check to make sure every district has at least one muni assigned to it.
    distDict={}
    for xdist in districtGeometry.keys():
        if xdist not in distDict.keys():
            distDict[xdist]=set()
        distDict[xdist].update({xreg for xreg in regionSubstitutes.keys() if regionSubstitutes[xreg]['Districts'][G.districtPrefix]==xdist})
        if len(distDict[xdist])==0:
            # Assign a single region to this district. The one with the largest overlapping area
            # This process isn't very efficient but won't be used very often
            # To speed up this process, first limit the the candidate regions to the boundary regions
            pickList=BoundaryMembers({x:regionSubstitutes[x]['Neighbors'] for x in regionSubstitutes.keys()})
            # From this list make another small dictionary with the areas of their overlap with district
            tryDict={x:regionSubstitutes[x]['geometry'].intersection(districtGeometry[xdist]).area for x in pickList}
            pickRegion=max(tryDict,key=tryDict.get)
            # assign it
            oldDist=regionSubstitutes[pickRegion]['Districts'][G.districtPrefix]
            regionSubstitutes[pickRegion]['Districts'][G.districtPrefix]=xdist
            print('  Error No regions for:',xdist,'. Assigned',pickRegion,'from',oldDist)

    
    # The district boundary can be ragged, leading to discontinuities. Fix them.

    # function build_Set creates a set of all values of the dictionary not counting its keys.
    # useful for create a set of neighbors
    def build_Set(inDict):
        if len(inDict.keys())==0:
            print(' Boundary dictionary empty. Exiting as None')
            return None
        outSet=set()
        for x in inDict.keys():
            outSet.update(inDict[x])
        outSet-=set(inDict.keys())
        return outSet

    xtest=True
    checkDistricts={regionSubstitutes[xreg]['Districts'][G.districtPrefix] for xreg in regionSubstitutes.keys()}
    denyRepeatDict={}         # A dictionary of subregion:from/to tuples to prevent infinite loopes
    while len(checkDistricts)>0:
        thisDist=getElement(checkDistricts)
        while True:
            theseRegions=[x for x in regionSubstitutes.keys() if regionSubstitutes[x]['Districts'][G.districtPrefix]==thisDist]
            inCountyNeighbors={}
            for x in theseRegions:
                inCountyNeighbors[x]={nbr for nbr in regionSubstitutes[x]['Neighbors'] if nbr[:4]==regionPrefix}
            # We must establish continuity within the county, so only look at the muni members of the district's regions
            check_Connected=CCFD.getConnected(theseRegions,inCountyNeighbors)

            icount=len(check_Connected)
            if icount==1:
                if G.verbose: print('  District ',thisDist,' is connected.')
                checkDistricts.remove(thisDist)
                break
            # We must fix unconnected subregions.
            # Two ways to fix:
            #   1) if lots of potential neighbors, find one in common with another connected subset.
            #   2) if unconnected group is small, remove them from the list
            # Generally, we will attempt the former rather than risk creating orphan subregions.
            if G.verbose:
                print('  Fixing unconnected:',thisDist)
                grp1={x:regionSubstitutes[x]['geometry'] for x in check_Connected[0]}
                grp2={x:regionSubstitutes[x]['geometry'] for x in theseRegions if x not in check_Connected[0]}
                CCFD_shapes.GroupPlot({'Set 1':grp1,'Set 2':grp2})

            # build_Set creates a set of all neighbors of the set 
            # Neighbors of the first connected set
            testSet=build_Set({x:regionSubstitutes[x]['Neighbors'] for x in check_Connected[0]})
            # Neighbors of the second connected set
            checkSet=build_Set({x:regionSubstitutes[x]['Neighbors'] for x in check_Connected[1]})

            # We only allow establishing connectedness within the county in STAGE 1
            testSet={x for x in testSet if x[:4]==regionPrefix}
            checkSet={x for x in checkSet if x[:4]==regionPrefix}

            # Create null entries in denyRepeat Dict if necessary to avoid error in the upcoming test
            # denyRepeatDict has a region key and a set of tuples of (oldDistrict,newDistrict) of reassignments tried for that region
            for x in testSet&checkSet:
                if x not in denyRepeatDict.keys():
                    denyRepeatDict[x]=set()

            # Build set of common Neighbors not tried yet
            commonSet=set()
            for x in checkSet & testSet:
                if (thisDist,getElement(regionSubstitutes[x]['Districts'][G.districtPrefix])) not in denyRepeatDict[x]:
                    commonSet.update({x})
                    
            if len(commonSet)>0:
                if G.verbose: print('  Connection options:',len(commonSet),commonSet)
                # We've found a common neighbor. Adding it to this district's regions will increased connectedness.
                # Select a region that accomplishes this
                thisRegion=list(commonSet)[0]
 
                # Change the district membership for this region
                oldDistrict=regionSubstitutes[thisRegion]['Districts'][G.districtPrefix]
                if thisRegion in denyRepeatDict.keys():
                    print('  Checking:',thisRegion,denyRepeatDict[thisRegion])
                    if (thisDist,oldDistrict) in denyRepeatDict[thisRegion] and len(commonSet)==1:
                        print('Problem.')
                        # Try combining oldDistrict and thisDist and splitting 50/50 some other way
                        tempSet={x for x in regionSubstitutes.keys() if List(regionSubstitutes[x]['Districts'])[0] in {thisDist,oldDistrict}}

                checkDistricts.update({oldDistrict})                                # Make sure we retest the giving regions for connectedness
                regionSubstitutes[thisRegion]['Districts'][G.districtPrefix]=thisDist

                if G.verbose: print('    Added to district',thisDist,': region',thisRegion,'(removed from: ',oldDistrict,')')

                # Forbid taking this action again on this region during this loop
                if thisRegion not in denyRepeatDict.keys():
                    denyRepeatDict[thisRegion]=set()
                denyRepeatDict[thisRegion].update({(oldDistrict,thisDist)})
                
            else:
                # There is no common neighbor.
                # As an alternative, find the neighbors of the smaller of the two comparison sets and change the district to that of one member of the selected neighbor
                if len(check_Connected[0])>len(check_Connected[1]):
                    fixIt=1
                    chgSet=checkSet
                else:
                    fixIt=0
                    chgSet=testSet
                chgSet.discard('Outside')  # Clears uncomputable neighbor
                chgDict={xreg:regionSubstitutes[xreg]['Districts'][G.districtPrefix] for xreg in chgSet if xreg[:4]==regionPrefix}

                # Remove any unacceptable candidates
                distSet=set()
                for xreg in chgDict.keys():
                    if xreg[:4]==regionPrefix:                      # Subsets must remain within this region
                        if chgDict[xreg]!=thisDist:                 # Must change the district
                            if xreg not in denyRepeatDict.keys():   # Not already tried
                                distSet.update({chgDict[xreg]})
                            elif (thisDist,chgDict[xreg]) not in denyRepeatDict[xreg]:   
                                distSet.update({chgDict[xreg]})

                if len(distSet)>0:
                    # Pick one of the district candidates randomly
                    newDistrict=Rnd.sample(distSet,1)[0]

                    if newDistrict==thisDist:
                        print('Error: moving from current District to itself.',xreg,',',thisDist)

                    # Now, simply change the district identifier for the chgSet regions to the new District
                    for xreg in check_Connected[fixIt]:
                        if G.verbose: print('  Added to district',newDistrict,': region',xreg,'(removed from: ',thisDist,')')
                        regionSubstitutes[xreg]['Districts'][G.districtPrefix]=newDistrict

                        # And prevent this from being tried again.
                        if xreg not in denyRepeatDict.keys():
                            denyRepeatDict[xreg]=set()
                        denyRepeatDict[xreg].update({(thisDist,newDistrict)})

                else:
                    #No options remaining. Need to reverse a previous step
                    for xreg in chgDict.keys():
                        if xreg in denyRepeatDict.keys():
                            for xtuple in denyRepeatDict[xreg]:
                                if xtuple[0]==thisDist:
                                    distSet.update({(xreg,xtuple[0],xtuple[1])})
                                    if G.verbose: print('  Need to reverse. Candidate = ',xreg,(thisDist,xtuple[1]))

    
                    # Don't allow repeat of this step
                    for xset in copy(distSet):
                        if (xset[2],xset[1]) in denyRepeatDict[xset[0]]:
                            distSet.discard(xset)

                    if len(distSet)==0:
                        if G.verbose: print('  No combinations available to select. Select a random neighbor ...')
                        for xreg in check_Connected[fixIt]:
                            if thisDist in regionSubstitutes[xreg]['Districts'][G.districtPrefix]:
                                for tryreg in regionSubstitutes[xreg]['Neighbors']:
                                    if tryreg[:4]==regionPrefix:
                                        if thisDist not in regionSubstitutes[tryreg]['Districts'][G.districtPrefix]:
                                            distSet.update({(xreg,regionSubstitutes[tryreg]['Districts'][G.districtPrefix],thisDist)})
                    else:
                        check_Connected[fixIt]=[xreg]


                    # Pick one of the district candidates randomly
                    pickTuple=Rnd.sample(distSet,1)[0]
                    xreg=pickTuple[0]
                    newDistrict=pickTuple[1]
                    thisDist=pickTuple[2]

                    #Make the change
                    for xreg in check_Connected[fixIt]:
                        if G.verbose: print('    Assigned to district',newDistrict,': region',xreg,'(removed from: ',thisDist,')')
                        regionSubstitutes[xreg]['Districts'][G.districtPrefix]=newDistrict

                    # And prevent this from being tried again.
                    if xreg not in denyRepeatDict.keys():
                        denyRepeatDict[xreg]=set()
                    denyRepeatDict[xreg].update({(thisDist,newDistrict)})

                    # Finally add newDistrict back to checkDistricts for analysis
                    checkDistricts.update({newDistrict})
                    checkDistricts.update({thisDist})

        # Return to test connectedness of remaining checkDistricts                    

    # Plot the final district assignments
    if G.verbose:
        plotDict={}
        for xdict in districtGeometry.keys():
            plotDict[xdict]={xreg:regionSubstitutes[xreg]['geometry'] for xreg in regionSubstitutes.keys() if regionSubstitutes[xreg]['Districts'][G.districtPrefix]==xdict}
        CCFD_shapes.GroupPlot(plotDict)
    print('Subregion assignments complete for ',regionPrefix)

    return regionSubstitutes