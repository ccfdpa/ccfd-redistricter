
import Global_Settings as G

from statistics import mean,pstdev
from copy import deepcopy
from math import sqrt,pi
from shapely.ops import unary_union


import CCFD_Functions as CCFD
import Shape_Handler as CCFD_shapes

# This function computes a data change associated with a swap of regions in a district
# "dataName" is a valid name for a data dictionary member
# "inreg" is a string value or a set of region string values added to the district;
# "outreg" is a string value or a set of regions removed from the district
# "inDistrictSet" is a single district data set;
# "regionMemberList" is a full dictionary of all region data composing districts
def makeAsSet(inreg):
    if isinstance(inreg,str):
        outSet={inreg}
    else:
        outSet=inreg
    return outSet

def chgDistrictRegions(districtSet,inregSet,outregSet):

    xdist=deepcopy(districtSet)
    inSet=makeAsSet(inregSet)
    outSet=makeAsSet(outregSet)

    if len(inSet)>0: xdist['Regions'].update(inSet)
    if len(outSet)>0:
        xdist['Regions']-=outSet
    return xdist

# Tests whether making the specified swap of regions breaks connectedness of the district
def MoveConnectedTest(inReg,outReg,districtSet,baseMemberList=None):
    testSet=districtSet['Regions'].union(makeAsSet(inReg))
    testSet-=makeAsSet(outReg)

    if baseMemberList==None:
        baseMemberList=G.regionMemberList

    # Invalid move if nothing is left in district (should not happen)
    if len(testSet)==0:
        if verbose: print('Program error. After move, no region is member of district',districtSet['ID'])
        return False

    # Test for connectedness if inreg enters district "district" and outreg leaves it
    xtest=CCFD.testConnectedness(testSet,baseMemberList)            
#            print('Test:',districtSet['ID'],xtest,testSet)
#            print('  pt 2',CCFD.getConnected(testSet,{x:G.regionMemberList[x]['Neighbors'] for x in G.regionMemberList.keys()}))
    if not xtest:
        if G.verbose:
            badList=CCFD.getConnected(testSet,{x:baseMemberList[x]['Neighbors'] for x in baseMemberList.keys()})
            badpieces=[len(x) for x in badList]
            print('  Rejected: Moving ',inReg,' in and ',outReg,' out of ',districtSet['ID'],' breaks connectedness (pieces):',badpieces)
        return False

    return True


# General Metric Change from Move
def metricCalc(metric,dataSet={'xxx':None},targetval=1.0,vars=False):

    if metric in dataSet.keys():
        if vars:
            return [metric]
        return dataSet[metric]
    elif metric=='compactTest':
        if vars:
            return list(G.regionData)
        return dataSet['BdgSqrIndex']
    elif metric=='PolsbyPopper':
        if vars:
            return ['Area','Perimeter']
        if abs(dataSet['Perimeter'])<0.0001:
               return 0.0
        return 4.0*pi*dataSet['Area']/dataSet['Perimeter']/dataSet['Perimeter']
    elif metric=='popTest':
        if vars:
            return ['TAPersons']
        return abs(float(dataSet['TAPersons'])/targetval-1.0)
    elif metric=='popVarTest':
        if vars:
            return ['TAPersons']
        return (float(dataSet['TAPersons'])/targetval-1.0)**2.0

    print('Stop at metric calculator:',metric)
    return None
# Creates a (mean,standard deviation) tuple from the data in a memberList
def metricMoments(thisVar,memberList,targetDict):
        if targetDict is None:
            targetDict={x:0.0 for x in memberList.keys()}
        valList=[metricCalc(thisVar,memberList[x]['Data'],targetDict[x]) for x in memberList.keys()]
        return (mean(valList),pstdev(valList))

# Computes the change in the value of a metric before and after a membership change
def metricChange(dataName,inregSet,outregSet,districtSet,regionMemberList=None,targetval=0.0):
    if regionMemberList is None:
        regionMemberList=G.regionMemberList

    current=metricCalc(dataName,districtSet['Data'],targetval)

    xdist=chgDistrictRegions(districtSet,inregSet,outregSet)

    varList=metricCalc(dataName,vars=True)
    rebuilt={var:0.0 for var in varList}
    for xreg in xdist['Regions']:
        rebuilt=CCFD.aggregate_Data(rebuilt,regionMemberList[xreg]['Data'],varList)
#    testOne=CCFD_shapes.checkPolygon(unary_union([regionMemberList[x]['geometry'] for x in xdist['Regions']]))
#    if not testOne:
#        print('Problem with geometry:',xdist)

#    rebuilt=CCFD.build_Data(xdist,regionMemberList)

    indexChange=metricCalc(dataName,rebuilt,targetval)-current 
    return indexChange

def metricMove(metric,inreg,outreg,inDistrict,outDistrict,regionMemberList,targetDict=None):

    if targetDict is None:
        targetDict={inDistrict['ID']:0.0,outDistrict['ID']:0.0}
    moveValue= metricChange(metric,inreg,outreg,inDistrict,regionMemberList,targetDict[inDistrict['ID']])
    moveValue+=metricChange(metric,outreg,inreg,outDistrict,regionMemberList,targetDict[outDistrict['ID']])
    return moveValue

def objectiveTest(goalIDs,inregSet,outregSet,srcDistrict,destDistrict,regionMemberList,districtMemberList,targetDict=None):

    # If a target value for each District is not defined, use 0.0 as the target
    if targetDict is None:
        targetDict={}
        for thisGoal in goalIDs:
            targetDict[thisGoal]={x:0.0 for x in districtMemberList}

    # Calculate the test metric changes
    # Each goal must have a target or a maximum possible value
    # Calculate a mean and variance from each target computed for each district from their deviation from the applicable target.
    # The weights apply to the relative current variances.
    #   Example: if one metric improves by one standard deviation and the other worsens by one standard deviation, the equally weighted goalID is zero.
    #   Proposed changes take place if the goalID increases.

    # Calculate the current distributional values for all metrics
    goalTuples={}
    goalValues={}

    tempSrcDist=deepcopy(districtMemberList[srcDistrict])
    tempDestDist=deepcopy(districtMemberList[destDistrict])

    objectiveValue=0.0
    if G.verbose: print('   Objectives:')
    for goalName in goalIDs.keys():
        goalTuples[goalName]=metricMoments(goalName,districtMemberList,targetDict[goalName])

        # Initial Values for destination
        goalValues[goalName]=[metricCalc(goalName,tempDestDist['Data'],targetDict[goalName][destDistrict])]
        # Initial Values for Source
        goalValues[goalName].append(metricCalc(goalName,tempSrcDist['Data'],targetDict[goalName][srcDistrict]))
        # Total Value change from the move (sum of changes to source and destination)
        goalValues[goalName].append(metricMove(goalName,inregSet,outregSet,tempSrcDist,tempDestDist,regionMemberList,targetDict[goalName]))
        # the mean and standard deviation for current values
        goalValues[goalName].extend(goalTuples[goalName])
        
    # Objective function improvement
    # goalValues list elements [initial Value
        contrib=goalValues[goalName][2]/goalTuples[goalName][1]    # Be careful here. the contribution is hardwired to the third element of goalValues
        goalValues[goalName].append(contrib)
        objectiveValue+=contrib*goalIDs[goalName]

        if G.verbose:
            print('     ',goalName,':',round(objectiveValue,3),
                            ', wt:',round(goalIDs[goalName],3),
                            ', contrib:',round(contrib,3))
            print('      Details:',[round(x,4) for x in goalValues[goalName]])

    if objectiveValue>0.0:
        results=True
    else:
        results=False

    return results,goalValues 

def EvaluateConstraints(constraintsSet,inregSet,outregSet,targDistrict,nbrDistrict,regionMemberList,districtMemberList):
    # Evaluate the constraints -- they apply to the modified districts
    
    if len(constraintsSet)==0:
        return True

    constraintsTest=True
    for thisConstraint in constraintsSet:

        countyCode=list(inregSet)[0][1:4]
        if thisConstraint=='sameCounty':
            for xreg in inregSet|outregSet:
                if xreg[1:4]!=countyCode:
                    # if any one region in the inset or outset is not in the county, the test fails
                    if G.verbose:
                        print('  sameCounty test failed: ',xreg)
                    return False

    # All of the constraints below this point require evaluations
        if thisConstraint in ['isolationTest','populationTest']:
            tempTarget=chgDistrictRegions(districtMemberList[targDistrict],outregSet,inregSet)
            tempNbr=chgDistrictRegions(districtMemberList[nbrDistrict],inregSet,outregSet)

        # Isolation of district (disallows donut-shaped districts)
        if thisConstraint=='isolationTest':
            tempTarget['distNeighbors']=CCFD.build_distNeighbors(tempTarget,regionMemberList)
            tempNbr['distNeighbors']=CCFD.build_distNeighbors(tempNbr,regionMemberList)
            isolationTest=len(tempTarget['distNeighbors'])>1 and len(tempNbr['distNeighbors'])>1
            if G.verbose and len(tempTarget['distNeighbors'])<2:
                print('  Isolation test failed:',tempTarget['ID'])
                return False
            if G.verbose and len(tempNbr['distNeighbors'])<2:
                print('  Isolation test failed:',tempNbr['ID'])
                return False


        # As a constraint, the population test requires that new populations remain either:
        #    1) within the established population variation limits
        #    2) the same as or lower than the input variations from the population target
        # This is weaker than the goal but allows non-compliance solutions to be selected provided they are improvements
        if thisConstraint=='populationTest':
            tempTarget['Data']=CCFD.build_Data(tempTarget,regionMemberList,['TAPersons'])
            tempNbr['Data']=CCFD.build_Data(tempNbr,regionMemberList,['TAPersons'])

            # Compute the 4 variances
            targVar=abs(districtMemberList[targDistrict]['Data']['TAPersons']-G.popTarget)
            tempTargVar=abs(tempTarget['Data']['TAPersons']-G.popTarget)

            nbrVar=abs(districtMemberList[nbrDistrict]['Data']['TAPersons']-G.popTarget)
            tempNbrVar=abs(tempNbr['Data']['TAPersons']-G.popTarget)

            if not (tempTargVar<=targVar or tempTargVar<int(G.popTarget*G.allowedPopVariance)):
                if G.verbose: print('  Population target variance constraint failed: in:',round(targVar/G.popTarget,3),'; out:',round(tempTargVar/G.popTarget,3))
                return False
            if not (tempNbrVar<=nbrVar or tempNbrVar<int(G.popTarget*G.allowedPopVariance)):
                if G.verbose: print('  Population neighbors variance constraint failed: in:',round(nbrVar/G.popTarget,3),'; out:',round(tempNbrVar/G.popTarget,3))
                return False

    return True
