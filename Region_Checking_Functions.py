# This module contains a number of checking functions to verify that the regionMemberList is appropriately populated.
import Global_Settings as G
from Shape_Handler import QuickPlot, checkPolygon
import CCFD_Functions as CCFD

from copy import copy

# Checks to make sure ID values are the same for each record
def checkIDConsistency(inList):
    for xID in inList.keys():
        if xID!=inList[xID]['ID']: 
            print('ID Test failed:',xID,inList[xID]['ID'],inList[xID]['Data']['ID'])
            return False
        if xID!=inList[xID]['Data']['ID']:
            print('ID Test failed:',xID,inList[xID]['ID'],inList[xID]['Data']['ID'])
            return False
    return True

# Checks whether the computed perimeter of a region equals the sum of the boundary lengths between the region and its neighbors
# Returns 'True' or 'False'. If in vebose mode, True also prints 'OK"; False prints a list of regions potentially causing the problem.
def checkRegionBoundaries(inList,boundaryDict=None,allowedBoundarySumVariance=175.0,verbose=None,indent=''):

    if boundaryDict is None:
        boundaryDict=G.boundaryLengthDict

    if verbose is None:
        verbose=G.verbose

    # allowedBoundarySumVariance=175.0 # Roughly 5 times the length in feet of the narrow dimension of a typical urban lot, or 1/2 city block

    if verbose: print(indent,'Boundary checks for',list(inList.keys()),'... ',end='')
    badreg=[]
    listToCheck=set(inList.keys())
    listToCheck.discard('Outside')
    Totresults=True
    for xreg in listToCheck:
        results=True
        sumbounds=0.0
        if xreg not in boundaryDict.keys():
            badreg.append(xreg)
            results=False
        else:
#            if inList[xreg]['Data']['Perimeter'] is None:
#                print('Error')
            exteriorLength=inList[xreg]['Data']['Perimeter']
            for xnbr in inList[xreg]['Neighbors']:
                if xnbr not in boundaryDict[xreg].keys():
                    badreg.append(xreg+'/'+xnbr)
                    boundaryDict[xreg][xnbr]=0.0
                    results=False
                else:
                    sumbounds+=boundaryDict[xreg][xnbr]

        if abs(exteriorLength-sumbounds*G.perimeterToMiles)>allowedBoundarySumVariance:
            badreg.append(xreg+'/Sum')
            results=False
            if verbose:
                print('Boundary inconsistency:',xreg)
                print('Perimeter:',round(exteriorLength,4))
                print('   Neighbors: (',len(inList[xreg]['Neighbors']),')')
                for xnbr in inList[xreg]['Neighbors']:
                    print('  ',xnbr,round(boundaryDict[xreg][xnbr]*G.perimeterToMiles,4))
                print('   Sum:',round(sumbounds*G.perimeterToMiles,4))


        Totresults=Totresults and results
    if not Totresults:
        print()
        print(badreg,': Boundary inconsistency')        
        print()
    elif verbose:
        print('OK')

    return Totresults

def boundaryDecomposer(inSet,baseList=None,verbose=True,bndyDetail=False):

    if baseList is None:
        baseList=G.regionMemberList

    # Use 'Group' for regionMemberLists, 'Regions' for districtMemberLists
    if 'Regions' in inSet.keys():
        thisCollection=inSet['Regions']
    elif 'Group' in inSet.keys():
        thisCollection=inSet['Group']
    else:
        print('Error: no Region list or Group list found')
        return

    if verbose:
        thisCollection=set(sorted(thisCollection))
        print('  Data perimeter: ',round(inSet['Data']['Perimeter'],4))

        print('  Collection:',thisCollection)
        print('  Members : length:')
        xPerimeters=0.0
        for xreg in thisCollection:
            xPerimeters+=baseList[xreg]['Data']['Perimeter']
            if bndyDetail: print('    ',xreg,':',round(baseList[xreg]['Data']['Perimeter'],4))
        print('     Total:',round(xPerimeters,4))
        print('  Common boundaries:')
        totCommon=0.0
        for xreg in thisCollection:
            xCommon=sum([G.boundaryLengthDict[xreg][xnbr] for xnbr in baseList[xreg]['Neighbors'] if xnbr in thisCollection])*G.perimeterToMiles
            totCommon+=xCommon
            if bndyDetail: print('    ',xreg,':',round(xCommon,4),'(',end='')
            for xnbr in baseList[xreg]['Neighbors'] & thisCollection:
                if bndyDetail: print(xnbr,round(G.boundaryLengthDict[xreg][xnbr]*G.perimeterToMiles,4),'), ',end='')
            if bndyDetail: print()
        print('     Total',round(totCommon,4))
        print('   Computed perimeter - Members:',round(xPerimeters,4),'; Common:',round(totCommon,4),'; Net:',round(xPerimeters-totCommon,4))

    thisID=inSet['ID']
    sumTotal=0.0
    print('   Computed from Neighbors:')
    for bndyReg in inSet['Neighbors']:
        if verbose and bndyDetail: 
            if bndyReg not in G.boundaryLengthDict[thisID].keys():
                print('    ',bndyReg,':','Missing','  (',end='')
            else:
                print('    ',bndyReg,':',round(G.boundaryLengthDict[thisID][bndyReg]*G.perimeterToMiles,3),'  (',end='')
        if bndyReg!='Outside':
            nbrSet=baseList[bndyReg]['Neighbors'] & thisCollection
            for xreg in nbrSet:
                if verbose and bndyDetail: print(xreg,':',end='')
                if xreg in G.boundaryLengthDict[bndyReg].keys():
                    thisTotal=G.boundaryLengthDict[bndyReg][xreg]*G.perimeterToMiles
                    sumTotal+=thisTotal
                    if verbose and bndyDetail: print(round(thisTotal,3),', ',end='')

                else:
                    if verbose and bndyDetail: print('Missing, ',end='')
            if verbose and bndyDetail: print(')')
        else:
            outsideSet={xreg for xreg in thisCollection if 'Outside' in baseList[xreg]['Neighbors']}
            xCommon=sum([G.boundaryLengthDict[xreg]['Outside'] for xreg in outsideSet])*G.perimeterToMiles
            if bndyDetail: print('    ','Outside',':',round(xCommon,4),'(',end='')
            for xreg in outsideSet:
                thisTotal=G.boundaryLengthDict[xreg]['Outside']*G.perimeterToMiles
                sumTotal+=thisTotal
                if bndyDetail: print(xreg,round(thisTotal,4),'), ',end='')
            if verbose and bndyDetail: print(')')
    if verbose:
        print('     Total',round(sumTotal,4))
    return

def CheckRegionDistrictConsistency(inDist,inReg):
    results=True
    try:
        temp={}
        for x in inDist.keys():
            regSum=[]
            for y in inDist[x]['Regions']:
                regSum.append(inReg[y]['Data']['TAPersons'])
            temp[x]=(inDist[x]['Data']['TAPersons'],sum(regSum))
    except:
        print('  Region/District population calculation failed.')
        return False
    badDistricts={}
    for x in temp.keys():
        if not (temp[x][0]==temp[x][1]):
            badDistricts[x]=temp[x]
            results=False
    if not results:
        print('  Bad (district,region) consistency:')
        districtPrefix=list(badDistricts.keys())[0][0]
        allbadPairs=set()
        for x in badDistricts.keys():
            badPairs={(x,y) for y in inDist[x]['Regions'] if inReg[y]['Districts'][districtPrefix]!=x}
            allbadPairs.update(badPairs)
        print('  Unpaired (district,Region): ',allbadPairs)
        print()

    return results

def checkRegionDistrictAssignmentConsistency(inRegList,inDistList):

    districtPrefix=list(inDistList.keys())[0][0]
    results=True
    for thisDist in inDistList.keys():
        regions=inDistList[thisDist]['Regions']
        for thisReg in regions:
            regDist=inRegList[thisReg]['Districts'][districtPrefix]
            if regDist!=thisDist:
                results=False
                print('  District/Region inconsistency (District,region member,region district):',thisDist,thisReg,regDist)
    return results

def CheckDistrictRegions(inList):
    for thisDist in inList.keys():
        if len(G.districtMemberList[thisDist]['Regions'])==0:
            print('  District has no regional members:',thisDist)
            return False
    return True

def CheckRegionAssignments(inList):
    for thisReg in inList.keys():
        xdist=inList[thisReg]['Districts'][G.districtPrefix]
        if xdist is None:
            print('  Region has no district assignment',thisReg)
            return False
        if len(xdist)<2:
            print('  Regional district assignment malformed:',thisReg,xdist)
            return False
    return True

def CheckForMissing(inVariable,inList,zeroOK=False):

    if 'Outside' in inList.keys():
        del inList['Outside']

    listTotal=list()
    xset=inList[list(inList.keys())[0]]
    if inVariable in G.regionData:
        if not zeroOK:
            listTotal.extend([x for x in inList.keys() if inList[x]['Data'][inVariable]==0])
        listTotal.extend([x for x in inList.keys() if inList[x]['Data'][inVariable] is None])
    elif inVariable in xset.keys():
        listTotal.extend([x for x in inList.keys() if inList[x][inVariable] is None])
    else:
        print('  Missing/zero Check:',inVariable,'not found. Ignoring.')
        return True

    if len(listTotal)>0:
        print('  Missing/zero Check failed for ',inVariable,':',listTotal)
        return False
    
    return True

def CheckRegionNeighbors(inList):
    # All neighbors must be known regions
    thisRegionSet=set(inList.keys())
    badRegions={x:inList[x]['Neighbors']-thisRegionSet for x in inList if len(inList[x]['Neighbors'] - inList[x]['Neighbors'])!=0}
    if len(badRegions)>0:
        print('  Regions:Neighbors not in region list:',badRegions)
        return False

    # Regions may not be their own neighbors
    badRegions={x for x in inList if x in inList[x]['Neighbors']}
    if len(badRegions)>0:
        print('  Regions in their own neighbor set:',badRegions)
        return False

    return True

def CheckPopTotal(inList,popCheck=None,verbose=None):

    if verbose is None: verbose=G.verbose
    if popCheck is None:
        popCheck=G.populationTotal
    listTotal=sum(inList[x]['Data']['TAPersons'] for x in inList.keys() if inList[x]['Data']['TAPersons'] is not None)
    if listTotal!=popCheck:
        if verbose: print('  Population Total Check failed.',listTotal,'(Actual:',popCheck,', Var:',listTotal-popCheck,')')
        return False
    return True

def CheckGeometries(inList,plotBad=False):
    badList=list()
    for x in inList.keys():
        if not checkPolygon(inList[x]['geometry']): badList.append(x)
    if len(badList)==0:
        return True
    print('  Polygon failure(s)')
    for x in badList:
        if len(badList)>1: print('  Bad Geometries:',badList)
        print('    ',x,'-- ',inList[x]['Name'],'-- ',inList[x]['geometry'].geom_type)
        if plotBad:
            QuickPlot({x:inList[x]['geometry']},showPlots=True)
    return False

def CheckConnectedness(distList,regList,verbose=None):

    if verbose is None: verbose=G.verbose
    results=True
    badDict={}
    for xdist in distList.keys():
        distRegs=distList[xdist]['Regions']
        # Make sure all the members are known
        missing=distRegs-set(regList.keys())
        if len(missing)==0:
            pieces=CCFD.getConnected(distRegs,{x:regList[x]['Neighbors'] for x in distRegs})
            if len(pieces)!=1:
                badDict[xdist]=pieces
                results=False
        else:
            if verbose: print('  No connectedness test: Missing regions:',missing)
            return False

    if not results:
        if verbose: print('  Error: some districts not connected:')
        for x in badDict.keys():
            print('    ',x,badDict[x])

    return results

# Check regionMemberList for completeness
# Returns extra & missing regions in subregion list compared with region list
def CheckSubregionCompleteness(xsubregionMemberList,xregionMemberList):
    checkFlag=True
    if isinstance(xregionMemberList,dict):
        baseList=list(xregionMemberList.keys())
    else:
        baseList=list(xregionMemberList)

    if isinstance(xsubregionMemberList,dict):
        checkList=list(xsubregionMemberList.keys())
    else:
        checkList=list(xsubregionMemberList)

    for x in baseList:
        if x not in checkList:
            if x[4]=='G':
                try:
                    for y in xregionMemberList[x]['Group']:
                        if y not in checkList:
                            checkFlag=False
                            print('  Extra G region:',y,'in',x)
                        else:
                            checkList.remove(y)
                except:
                    continue
            elif x[0]=='C':
                for y in {xreg for xreg in checkList if xreg[1:4]==x[1:4]}:
                    if y not in checkList:
                        checkFlag=False
                        print('  Extra subregion:',y,'in',x)
                    else:
                        checkList.remove(y)
            else:
                checkFlag=False
                print('  Unassigned subregion:',x)

        else:
            checkList.remove(x)

    # Remove the outside region references:
    for x in copy(checkList):
        if x[1:4]=='000' or x=='Outside': checkList.remove(x)

    if len(checkList)>0:
        for x in checkList:
            checkFlag=False
            print('  Missing region:',x)
    
    if not checkFlag: print('  Region completeness check failed.')
    return checkFlag

# Check neighbors for consistency
#    (The neighbors of a region should have the region as a neighbor

def CheckNeighborConsistency(inMemberList):
    checkFlag=True
    for thisRegion in inMemberList.keys():
        if len(inMemberList[thisRegion]['Neighbors'])==0:
            checkFlag=False
            print('  No neighbors for region',thisRegion)
        else:
            for nbr in {xnbr for xnbr in inMemberList[thisRegion]['Neighbors'] if xnbr!='Outside'}:
                if nbr not in inMemberList.keys():
                    checkFlag=False
                    print('  Region Not Found:',nbr,'(a neighbor of ',thisRegion,')')
                else:
                    if thisRegion not in inMemberList[nbr]['Neighbors'] and thisRegion!='Outside':
                        checkFlag=False
                        print('  Neighbor inconsistency. ',nbr,'nbr of',thisRegion,' but not reverse.')

    if not checkFlag: print('  Neighbor consistency check failed.')

    return checkFlag

def StageTests(regionMemberList=None,districtMemberList=None,useGeometry=False):

    if regionMemberList is None: regionMemberList=G.regionMemberList
    if districtMemberList is None: districtMemberList=G.districtMemberList

    allResults=True
    results=True
    print('Region checks ...',end=' ')
    results=results and CheckRegionAssignments(regionMemberList)
    results=results and CheckPopTotal(regionMemberList)
    results=results and CheckRegionNeighbors(regionMemberList)
    if useGeometry:
        results=results and CheckGeometries(regionMemberList)
    if results:
        print('OK')
    allResults=results

    results=True
    print('District checks ...',end=' ')
    results=results and CheckDistrictRegions(districtMemberList)
    results=results and CheckPopTotal(districtMemberList)
    results=results and CheckRegionNeighbors(districtMemberList)
    results=results and CheckConnectedness(districtMemberList,regionMemberList)
    if results:
        print('OK')
    allResults=allResults and results


    results=True
    print('Consistency checks ...',end=' ')
    results=results and CheckNeighborConsistency(regionMemberList)
    results=results and checkRegionDistrictAssignmentConsistency(regionMemberList,districtMemberList)
    results=results and CheckRegionDistrictConsistency(districtMemberList,regionMemberList)
    if results:
        print('OK')
    allResults=allResults and results

    return allResults

