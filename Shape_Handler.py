import Global_Settings as G
import CCFD_Functions as CCFD

import numpy as np
import pandas as pd
import os
import sys
from math import log,sqrt
from copy import copy,deepcopy

# geospatial
import geopandas as gp
from shapely import geometry
from shapely.geometry import Polygon, LineString, MultiPoint, shape,JOIN_STYLE
from shapely.ops import snap,unary_union
from pysal.lib import weights # Needed to constract contiguity matrix
import matplotlib.pyplot as plt

# Projections
import pyproj

global defaultCRS, mercatorCRS, cartesianCRSstrg
defaultCRS=pyproj.CRS.from_epsg(4269)
mercatorCRS=pyproj.CRS.from_epsg(3395)
cartesianCRSstrg=G.projType+G.projUnits

def GroupPlot(inGroupDict,showPlots=None):
    # IngroupDict is a dictionary of {x:geometry} dictionaries
    if showPlots is None:
        showPlots=G.showPlots

    if inGroupDict is None:
        print('No districts to plot.')
        return

    smallFrame=gp.GeoDataFrame()
    smallFrame['subRegion']=None
    smallFrame['color']=None

    colors='rgbcmy'
    modcolors=len(colors)
    seq=0

    for xdist in inGroupDict.keys():
        for xreg in inGroupDict[xdist].keys():

            try:
                smallFrame=smallFrame.append({'subRegion':xdist+'/'+xreg,
                                              'color':colors[seq % modcolors],
                                              'edgecolor':'k',
                                              'geometry':inGroupDict[xdist][xreg]},ignore_index=True)
            except:
                print('Error in GroupPlot')
        seq+=1

    if not smallFrame.empty:
        smallFrame['geomtype']=smallFrame.geom_type
        print(smallFrame[['subRegion','color','geomtype']])
        if showPlots:
            smallFrame.plot(color=smallFrame['color'],edgecolor=smallFrame['edgecolor'], legend=True)
            plt.show()
    else:
        print('GroupPlot: No printable geometries.')

    return

def QuickPlot(inregDict,targetRegions=None,showPlots=None):
    # inregDict is a dictionary of {label:shape}
    # targetRegions is a list subset of the inregDict keys
    if showPlots is None:
        showPlots=G.showPlots
    if inregDict is None:
        print('No regions to plot.')
        return

    # Keep only plotable geometries
    plotList=[x for x in inregDict if not inregDict[x].is_empty]

    smallFrame=gp.GeoDataFrame()
    smallFrame['subRegion']=None
    smallFrame['color']=None

    colors='rgbcmy'
    modcolors=len(colors)
    seq=0

    for subseq in range(len(plotList)):
        try:
            smallFrame=smallFrame.append({'subRegion':plotList[subseq],'color':colors[seq % modcolors],'geometry':inregDict[plotList[subseq]]},ignore_index=True)
            seq+=1
        except:
            print('Plotting skipped: ',plotList[subseq])

    smallFrame.set_index('subRegion',inplace=True)
    if targetRegions is not None:
        for x in targetRegions:
            if x in plotList:
                smallFrame.loc[x,'color']='k'

    if not smallFrame.empty:
        smallFrame['geom_type']=smallFrame.geom_type
        smallFrame['edgecolor']='k'
        print(smallFrame[['color','edgecolor','geom_type']])
        if showPlots:
            smallFrame.plot(color=smallFrame['color'],edgecolor=smallFrame['edgecolor'])            
            plt.show()
    else:
        print('QuickPlot: No printable geometries.')
    
    return

def setEPSG(inEPSG):
    return pyproj.CRS.from_epsg(inEPSG)

def changeCRS(inFrame,crsID,geometry='geometry',zeroLongLatTuple=None):

    thisCRS=defaultCRS
    if inFrame.crs is None:
        inFrame=inFrame.set_crs(defaultCRS)
    if crsID=='default':
        thisCRS=defaultCRS
    elif crsID=='mercator':
        thisCRS=mercatorCRS
    elif crsID=='cartesian':
        if zeroLongLatTuple is None:
            print('Warn: cartesian CRS type requires the zeroLongLatTuple argument. Using mercator CRS')
            thisCRS=mercatorCRS
        else:
            projLat='+lat_0='+str(zeroLongLatTuple[0])+' '
            projLong='+long_0'+str(zeroLongLatTuple[1])+' '
            thisCRS=cartesianCRSstrg+projLat+projLong

    else:
        print('Error: Unknown CRS. Returning with default CRS')
    
    if thisCRS!=inFrame.crs:
        inFrame=inFrame.to_crs(thisCRS)
    return inFrame


def CalcPreciseGeoms(inFrame,geometry='geometry'):
    # Change the crs to constant area located at each VTD centroid
    # This must be done one region at a time because geopandas allows only one crs for the whole matrix
    
    chgFrame=inFrame.copy()
    chgFrame=changeCRS(chgFrame,'default')

    # Compute an approximate centroid for the Frame members.
    chgFrame['centroids']=changeCRS(chgFrame[geometry],'mercator').centroid
    chgFrame['centroids']=changeCRS(chgFrame['centroids'],'default')
    inFrame['Latitude']=chgFrame['centroids'].y
    inFrame['Longitude']=chgFrame['centroids'].x

    if len(inFrame)==0:
        print('Error: no regions to analyze.')
        print('  Check raw subregion selection substring for the CreateLists action.')
        print('Exiting ...')
        sys.exit(0)


    print('  Computing areas and perimeters ',end='')
    areas=[]
    lengths=[]
    xcounter=0.0
    dotcounter=0.0
    print('[                   ]',end='')
    print(chr(8)*20,end='')
    markIterate=1.0/len(inFrame)
    markADot=0.05

    for row in chgFrame.itertuples():

        testGEO=row.Index

        gShape=gp.GeoSeries(row.geometry)

        # the centroid here is used only to locate the zero point for the projection. Doesn't have to be exact.
#        thisCRS=pyproj.crs.CRS(crsString)

#        g2Shape=gShape.to_crs(thisCRS)
        g2Shape=changeCRS(gShape,'cartesian',zeroLongLatTuple=(row.centroids.x,row.centroids.y))
        areas.append(g2Shape.area[0])
        lengths.append(g2Shape.length[0])

        xcounter+=markIterate
        if xcounter>markADot:
            xdots=int(20*(dotcounter+xcounter))-int(20*dotcounter)
            dotcounter+=xcounter
            print(chr(9608)*xdots,end='')
            xcounter=0
    print()

    inFrame['Area']=areas
    inFrame['Perimeter']=lengths

    return inFrame

def IsNeighbor(shape1,shape2):

    # First, check the bounds for adjacency.
    try:
        # In PA, these are 4-tuples of the extrema of the lat-longs (west,south,east,north)
        Onetest=shape1.bounds
        Twotest=shape2.bounds
    except:
        print('IsNeighbor: One or more shapes not accepted.')
        return None
    
    try:
        if CCFD.deltaLONG_to_miles(Onetest[0]-Twotest[2],Onetest[1])>0.01:   # More than 0.01 miles separation
            if CCFD.deltaLONG_to_miles(Twotest[0]-Onetest[2],Onetest[1])>0.01: 
                return False
    except:
        print('IsNeighbor Fail 1: ', Onetest,Twotest)
        return False
    try:
        if CCFD.deltaLAT_to_miles(Onetest[1]-Twotest[3])>0.01:
            if CCFD.deltaLAT_to_miles(Twotest[1]-Onetest[3])>0.01: 
                return False
    except:
        print('IsNeighbor Fail 2: ', Onetest,Twotest)
        return False

    # Indirect test -- based on overlap
    if shape1.overlaps(shape2):
        if shape1.intersection(shape2).area*3632.5>17.0:             # rough size of a city block
            print('IsNeighbor Fail3-overlaps')
            QuickPlot({'Shape 1':shape1,'Shape 2':shape2,'Intersection':shape1.intersection(shape2)})
            return False
        shapex=shape1.difference(shape2)
        sliverLength=shapex.intersection(shape2.buffer(0.000001)).area/0.000001*3.14159/180.0*3960.0*5280.0
    else:
        sliverLength=shape1.intersection(shape2.buffer(0.000001)).area/0.000001*3.14159/180.0*3960.0*5280.0
    if sliverLength>175.0:                       # Roughly 5 times the length in feet of the narrow dimension of a typical urban lot (1/2 city block)
            return True

    return False

# regression function returns a tuple (intercept, slope)
def estimate_coef(x, y):
    # Ref: https://www.geeksforgeeks.org/linear-regression-python-implementation/
    # number of observations/points 
    n = np.size(x) 
  
    # mean of x and y vector 
    m_x, m_y = np.mean(x), np.mean(y) 
  
    # calculating cross-deviation and deviation about x 
    SS_xy = np.sum(y*x) - n*m_y*m_x 
    SS_xx = np.sum(x*x) - n*m_x*m_x 
  
    # calculating regression coefficients 
    b_1 = SS_xy / SS_xx 
    b_0 = m_y - b_1*m_x 
  
    return (b_0, b_1) 

def DensifyPolygon(thisShape,scaleLength=60.27,pieceLength=0.10):

    # scaleLength: delta lat/long to miles conversion (at ESPG:4269 latitude)
    # pieceLength: miles between straightline subsegments
    try:
        if thisShape.geom_type!='Polygon':
            print('Densify: Shape not polygon. Return unchanged.')
            return thisShape
    except:
        print('Densify: Input shape not recognized:',thisShape)
        return thisShape

    # Shape is a polygon. We can proceed.
    theseLongs=thisShape.exterior.coords.xy[0]
    theseLats=thisShape.exterior.coords.xy[1]
    fixedPoints=[]
    for segment in range(1,len(theseLongs)):
        segLength=sqrt((theseLongs[segment]-theseLongs[segment-1])**2 +(theseLats[segment]-theseLats[segment-1])**2)
        numSubSegments=max(1.0,segLength*scaleLength/pieceLength)
        xSubLength=(theseLongs[segment]-theseLongs[segment-1])/numSubSegments
        ySubLength=(theseLats[segment]-theseLats[segment-1])/numSubSegments
        for segSeq in range(int(numSubSegments)+1):

            nextPoint=(theseLongs[segment-1]+segSeq*xSubLength,theseLats[segment-1]+segSeq*ySubLength)
            fixedPoints.append(nextPoint)
 
    fixedShape=Polygon(fixedPoints).buffer(0.0)

    return fixedShape

def rayfromPoints(pt1,pt2):

    if pt1[0] is None or pt1[1] is None or pt2[0] is None or pt2[1] is None:
        return (None,None)

    slope=(pt1[0]-pt2[0])/(pt1[1]-pt2[1])
    intercept=(pt1[0]+pt2[0]-slope*(pt1[1]+pt2[1]))/2.0

    return (intercept,slope)

def AddOutsideShape(inShapes):
    # The input variable is a list of geometries
    # The output is a donut shape with the union of the inDict shapes the cutout.
    
    # Build a big union of the shapes in inDict
    bigshape=unary_union(inShapes)
    if not checkPolygon(bigshape,verbose=False):
        print('Warning: Union of all members not a polygon.')
    outside=bigshape.exterior.buffer(0.1)
    outside=outside.difference(bigshape)
    # We want the outside to be a Polygon with one interior ring
    if not checkPolygon(outside,verbose=False):
        pieces={x: outside.geoms[x] for x in range(len(outside.geoms))}
        bigpiece=max(pieces,key=lambda x:pieces[x].area)
        outside=pieces[bigpiece]
        if not checkPolygon(outside,verbose=False):
            print('Outside region not accepted.')
            return None

    if len(outside.interiors)!=1:
        print('Outside region has too many interior holes.')
        return None
   
    return outside

def makeSplitterRay(inFrame, subRegpopDict):

    subReg_Frame=gp.GeoDataFrame()
    subReg_Frame['subRegion']=None
    subReg_Frame['geometry']=gp.GeoSeries
    #set the geometry
    subReg_Frame=subReg_Frame.set_geometry('geometry')

    # Start with the bounds of the full region. Split by population share
    xbounds=inFrame.bounds

    sPoint=xbounds['miny'][0]
    wPoint=xbounds['minx'][0]
    ePoint=xbounds['maxx'][0]
    nPoint=xbounds['maxy'][0]
    regPop=sum([subRegpopDict[x] for x in subRegpopDict.keys()])
    nsTotal=nPoint-sPoint
    xsPoint=sPoint
    inFrame.drop(['TAPersons'],axis=1,inplace=True)

    #collect the subregions into a geodataframe
    for regionID in subRegpopDict.keys():

        xnPoint=float(subRegpopDict[regionID])/float(regPop)*nsTotal

        p1=(wPoint,xsPoint)
        p2=(wPoint,xsPoint+xnPoint)
        p3=(ePoint,xsPoint+xnPoint)
        p4=(ePoint,xsPoint)
        px=Polygon([p1,p2,p3,p4])
        subReg_Frame=subReg_Frame.append({'subRegion':regionID,'TAPersons':subRegpopDict[regionID],'geometry':px},ignore_index=True)

        # Up to next slice
        xsPoint+=xnPoint

    # The group intersections creates a group of subregion geometries
    new_subRegions=gp.overlay(inFrame,subReg_Frame,how='intersection')

    # Now find the centroids of these regions and regress a straight line through them
    slice_centroids=new_subRegions.centroid
    # We're going to divide in n-s direction, so regress longitude (dependent) agains latitude (independent)
    # In doing so, we avoid a division by zero problem when the centroid regression is straight a north and south line.
    if len(subRegpopDict.keys())>2:
        regr=estimate_coef(slice_centroids.centroid.y,slice_centroids.centroid.x)  # We choose the n-s direction to be the independent variable here.
        intercept=regr[0]
        slope=regr[1]
    else:
        outRay=rayfromPoints((slice_centroids[0].x,slice_centroids[0].y),(slice_centroids[1].x,slice_centroids[1].y))
        if outRay[0] is None:
            print('Problem with centroids.')
        else:
            intercept=outRay[0]
            slope=outRay[1]

    # for debugging purposes, we add the splitter to the geodataframe that is returned
    fitted=LineString([((ePoint+wPoint)/2.0,sPoint),((ePoint+wPoint)/2.0,nPoint)])

#    print(new_subRegions)
    if G.verbose:
        print('Plotting a splitter ray for : ',subRegpopDict)
        new_subRegions=new_subRegions.append({'subRegion':'Fitted','geometry':LineString(fitted)},ignore_index=True)
        QuickPlot({x:new_subRegions.loc[new_subRegions['subRegion']==x]['geometry'].unary_union for x in new_subRegions['subRegion']})
#        new_subRegions.plot(cmap='Set1')
#        plt.show()
        
    return (intercept,slope)

# Returns TRUE is the input geometry is a LineString
def checkLineString(srcGeom):

    if srcGeom is None:
        return False

    elif isinstance(srcGeom.geom_type,str):
        if srcGeom.geom_type!='LineString':
            if G.verbose:
                print('checkLineString False. ',srcGeom.geom_type)
            return False
    else:
        return False

    return True

# Returns TRUE is the input geometry is a polygon
# If srcPolygon is empty, returns True
def checkPolygon(srcPolygon,verbose=True):

    if srcPolygon is None:
        return False

    if not isinstance(srcPolygon.geom_type,str):
        return False

    if srcPolygon.geom_type=='Polygon':
        return True

    if verbose:
        print('checkPolygon False. ',srcPolygon.geom_type)

    return False

# this function analyzes the interior structures of an input polygon geometry
#   eliminating any structure not resolvable to a polygon bigger than about 0.01 acres (about 21 ft square)
def InteriorsRepair(regionGeom):
    if not checkPolygon(regionGeom):
        print('Error. Shape is not a polygon.')
        return None
            
    numInteriors=len(regionGeom.interiors)
    if numInteriors>0:
        keptInteriors=[]
        for x in range(numInteriors):
            thisInt=regionGeom.interiors[x]
            try:
                Int_to_Polygon=Polygon(thisInt)
                xarea=Int_to_Polygon.area*2324800.0
            except:
                xarea=0.0
            if xarea>0.01:
                keptInteriors.append(thisInt)
            else:
                if G.verbose: print('  Removing a small hole of area:',round(xarea,3))
        repairedGeom=Polygon(regionGeom.exterior,[[pt for pt in xhole.coords] for xhole in keptInteriors])

    else:
        repairedGeom=regionGeom

    return repairedGeom

# Remove very small linearstrings from  a subset of a memberList
# Returns a subset of inList with repaired geometries
# Does not recompute the geometry-based values
def ScrubPolygons(inCore,inList=None):
    if inList is None:
        inList=G.regionMemberList
    fixedList={}
    coreLength=len(inCore)
    coreList={x for x in inList.keys() if CCFD.coreID(x)[:coreLength]==inCore}
    fixedCount=0
    fixedDict={}
    for thisReg in coreList:
        thisGeom=inList[thisReg]['geometry']
        inHoles=len(thisGeom.interiors)
        if inHoles>0:
            tryFix=InteriorsRepair(thisGeom)
            if tryFix is None:
                print('      Problem with fixing:',thisReg)
            else:
                if inHoles>len(tryFix.interiors):
                    fixedDict[thisReg]=deepcopy(inList[thisReg])
                    fixedDict[thisReg]['geometry']=copy(tryFix)
                    fixedCount+=1
                    print('      Fixed small holes:',thisReg)
    if fixedCount>0:
        print('    Total fixed:',fixedCount)

    return fixedDict

# Makes neighbors of inreg to touch by adding a bit to inreg then cutting off any overlaps with inList members
# inreg: a region ID
# inList: a dictionary {ID:geometry}
def Fix_Neighbor_Geo(inreg,inList=None):
    if inList is None: 
        Neighbors=G.regionMemberList[inreg]['Neighbors']
        inList={x:G.regionMemberList[x]['geometry'] for x in Neighbors+{inreg}}
    else:
        Neighbors=set(inList.keys())-{inreg}
    # Probably unnecessary action
    Neighbors.discard(inreg)

    # This routine can add a little bit to holes, water boundaries, larger political boundaries, etc.
#    G.verbose=True
    print('    Fix_Neighbor_Geo: ',inreg,Neighbors)

#    if inreg in ['C029.0','C029.1']:
#        print('Interrupt:')

    sigDigits=G.sigDigits

    geo1=inList[inreg]

    # try snapping the region to its neighbors after doing a miter buffer
    snapTolerance=100.0*1.5*10.0**-sigDigits
    touchFlag=True
    for x in Neighbors:
        if not geo1.touches(inList[x]):
            geo1Fixed=geo1.buffer(snapTolerance,1,join_style=JOIN_STYLE.mitre)
            geo1Fixed=snap(geo1,inList[x],snapTolerance)
            geo1Fixed=geo1Fixed.buffer(0.0)
            geo1Fixed=geo1Fixed.intersection(inList[x]).buffer(snapTolerance,1,join_style=JOIN_STYLE.mitre)
            geo1Fixed=geo1Fixed.union(geo1)
        else:
            geo1Fixed=geo1
        if not geo1Fixed.touches(inList[x]): 
            touchFlag=False
            print('    Buffer/snap failed:',inreg,x)
            print('    Is valid?:',geo1Fixed.is_valid)
            print('    Intersects?:',geo1Fixed.intersects(inList[x]))
            print('    Overlaps?:',geo1Fixed.overlaps(inList[x]))
            print('    Disjoint?:',geo1Fixed.disjoint(inList[x]))
            if geo1Fixed.intersects(inList[x]) or geo1Fixed.overlaps(inList[x]):
                print('    Intersection area: ',str(round(geo1Fixed.intersection(inList[x]).area*3632.5*5280.0*5280.0,0)),'sq. ft.')
                xintersection=geo1Fixed.intersection(inList[x]).buffer(snapTolerance,1,join_style=JOIN_STYLE.mitre)
                if G.verbose: QuickPlot({inreg:geo1Fixed,x:inList[x],'Inter':xintersection})
            else:
                if G.verbose: QuickPlot({inreg:geo1Fixed,x:inList[x]})

        # try removing any intersections of inreg with neighbors
        if not geo1Fixed.touches(inList[x]):
            if geo1Fixed.intersects(inList[x]):
                xarea=geo1Fixed.intersection(inList[x]).area*3632.5*5280.0*5280.0
                if xarea<1.0:
                    geo1Fixed=geo1Fixed.difference(inList[x])
                    geo1Fixed.buffer(0.0)
                if geo1Fixed.touches(inList[x]):
                    touchFlag=True and touchFlag
                    print('  Touch success from intersection removal: ',inreg,x)

        if not touchFlag: print('  Touch failure 1: ',inreg,x)
        if checkPolygon(geo1Fixed) and checkPolygon(inList[x]):
            if IsNeighbor(geo1Fixed,inList[x]): print('IsNeighbor success 1:')
        else:
            print('Polygon problem.')

    if G.verbose:
        checkPlot=copy(inList)
        checkPlot[inreg]=geo1Fixed
        QuickPlot(checkPlot,inreg)

    if touchFlag and checkPolygon(geo1Fixed):
        return geo1Fixed

    # If we still have problems try some buffering

    touchFlag=True
    for x in Neighbors:
        touchCheck=round(geo1.distance(inList[x])*60.27*5280.0,0)
        if touchCheck<1.0:
            smallpiece=geo1Fixed.difference(inList[x].buffer(snapTolerance))
            geo1Fixed=snap(geo1Fixed,inList[x],snapTolerance)
            geo1Fixed=geo1Fixed.buffer(0.0)
            if not geo1Fixed.touches(inList[x]):
                print('  Touch failure from infill: ',inreg,x)
                touchFlag=False
        else:
            print('  No action taken. Large gap exists: ',inreg,x,'(',int(touchCheck),'ft. )')
            checkPlot=copy(inList)
            checkPlot[inreg]=geo1
            QuickPlot(checkPlot,inreg)
            return geo1

    if not checkPolygon(geo1Fixed):
        # Try repairing it
        geo1Fixed=Single_from_MultiPolygon(geo1Fixed)
        if not checkPolygon(geo1Fixed):
            print('No action taken. Geometry repair failed for ',inreg)
            checkPlot=copy(inList)
            checkPlot[inreg]=geo1Fixed
            QuickPlot(checkPlot,inreg)
            return geo1
 
    return geo1Fixed

# takes a multipolygon, returns the piece with greatest area
def Single_from_MultiPolygon(newGeoseries,forced=False):

    if newGeoseries is None:
        print('Single_from_MultiPolygon: Geoseries is None.')
        return None

    xlist=[x for x in range(len(newGeoseries)) if isinstance(newGeoseries[x].geom_type=='Polygon',bool)]
    # First, eliminate very small polygons
    keepIndex=[x for x in xlist if newGeoseries[x].area>0.1*10**-G.sigDigits]
    
    if len(keepIndex)==0:
        print('No polygons found.')
        return None

    areadict={x:newGeoseries[x].area for x in keepIndex}
    if forced:
        # Choose the largest and use it
        keepIndex=[max(areadict,key=areadict.get)]

    # Try to reassemble into a single polygon
    # In this fails, return the largest

    # Start with the largest
    keepPolygon=newGeoseries[keepIndex[0]]
    if G.verbose: print('Combining: ',len(keepIndex))
    for x in range(1,len(keepIndex)):
        try:
            keepPolygon=keepPolygon.union(newGeoseries[x][1])
        except:
            keepPolygon=keepPolygon.union(newGeoseries[x])

        if not checkPolygon(keepPolygon):
            print('Reassembly failed.')
            if G.verbose: QuickPlot({'Piece'+str(x):newGeoseries[x] for x in range(len(newGeoseries))},keepIndex[0])

            return None

    # Pick the multipolygon with the largest area. 
    # Make it the actual newRegion.
    if len([x['subRegion'] for x in newGeoseries if x.geom_type!='Polygon'])>0:
        print('problem with explode for ',newRegion)
        print('  ',newGeoseries)

    keepPolygon=newGeoseries[max(areadict,key=areadict.get)]
    if G.verbose: print('Single_from_Multipolygon succeeded: ',keepPolygon.geom_type)
    
    return keepPolygon

# Input is a geodataframe with a single mulitpolygon geometry
# Repair involves dropping very small polygons, combining larger ones into single polygons
def MultiPolygon_Repair(badRegion,srcRegion,inFrame,regionColumnName=None,geometryColumnName=None):
    print('  MultiPolygon_Repair for: ',badRegion,'; srcRegion: ',srcRegion)

    if regionColumnName is None: regionColumnName='subRegion'
    if geometryColumnName is None: geometryColumnName='geometry'

    fixedFrame=copy(inFrame)
    fixedFrame.set_index(regionColumnName)
    # Explode the multipolygons and add them to the frame

    # Begin by merging srcRegion and badRegion. Doing so must create a single polygon geometry
    srcGeo=inFrame.loc[inFrame[regionColumnName]==srcRegion,geometryColumnName].unary_union
    badGeo=inFrame.loc[inFrame[regionColumnName]==badRegion,geometryColumnName]
    holdsrcGeometry=srcGeo.union(badGeo.unary_union)
    if not checkPolygon(holdsrcGeometry): holdsrcGeometry=srcGeo.union(badGeo.unary_union.buffer(0.01))
    if not checkPolygon(holdsrcGeometry):
        print('Union failed.')
        # Try buffering here?
        
    # Fix newRegion and srcRegion if a multipolygon is created
#    newGeoseries=xFrame.loc[inFrame[regionColumnName]==badRegion,geometryColumnName].unary_union
#    if not checkPolygon(newGeoseries):
        # newRegion is a multipolygon. Not acceptable.
    if G.verbose: print('  Fixing multipolygon issue for ',badRegion,'(',badGeo.geom_type,')')
    newGeoseries=badGeo.explode()

    # Simplify by elminating very small slices, if possible
    keepPolygon=Single_from_MultiPolygon(newGeoseries,True)
    # srcRegion becomes the original source region with keepPolygon cut out
    # if it returns "None" then the fix failed.
    if keepPolygon is not None:
        srcPolygon=holdsrcGeometry.difference(keepPolygon)
    else:
        srcPolygon=None

# ****************************** Working here ******************************

    if not checkPolygon(srcPolygon):
        if G.verbose: print('  Fixing problem with srcPolygon:',srcPolygon.geom_type)
        srcPolygon=Single_from_MultiPolygon(srcPolygon)

    if checkPolygon(keepPolygon) and checkPolygon(srcPolygon):
        # this worked. We're done. return.
        fixedFrame.loc[fixedFrame[regionColumnName]==badRegion,geometryColumnName]=Polygon(keepPolygon)
        fixedFrame.loc[fixedFrame[regionColumnName]==srcRegion,geometryColumnName]=Polygon(srcPolygon)
        return fixedFrame

    # Prints below only if failure
    if G.verbose:
        print('Repair failure. No changes.')
        print(fixedFrame)
        QuickPlot({x:fixedFrame.loc[fixedFrame['subRegion']==x,'geometry'] for x in fixedFrame['subRegion']})
 
    return inFrame

# Creates a set of subregions from a dataframe with one region and returns a dataframe with many subregions based on a splitting value (e.g. population share)
#   "splitRay" is a (intercept,slope) tuple defining a straight line. The splits will be perpendicular to this line
def subRegion_Shaper(regionCode,regionList=None,regionFrame=None,splitRay=None):
    if G.verbose: print('    subRegion_Shaper for: ',regionCode)
    if regionList is None: regionList=G.regionMemberList

    codestr=CCFD.coreID(regionCode)
    ctySubregions=list(x for x in regionList.keys() if CCFD.coreID(x)==codestr)
    regPop=sum(regionList[x]['Data']['TAPersons'] for x in ctySubregions)     #Total population for combined region

    sigDigits=G.sigDigits      # required to make sure subregions "touch"

    if regionFrame is None:
        # By default, if the regionFrame is not specified, it is regionCode
        regionFrame=gp.GeoDataFrame()
        regionFrame['subRegion']=None
        regionFrame['TAPersons']=None
        regionFrame['geometry']=gp.GeoSeries
        regionFrame.set_geometry('geometry',inplace=True)

        regionFrame=regionFrame.append({'subRegion':codestr,'TAPersons':regPop,'geometry':regionList[regionCode]['geometry']},ignore_index=True)
        regionFrame.set_index('subRegion',inplace=True)

    if len(ctySubregions)<2:
        print('  subRegion_Shaper: Found only one region. Exiting ...')
        return None

    # Get the extrema for the region
    xbounds=regionFrame.bounds

    sPoint=round(xbounds['miny'][0],sigDigits)
    wPoint=round(xbounds['minx'][0],sigDigits)
    ePoint=round(xbounds['maxx'][0],sigDigits)
    nPoint=round(xbounds['maxy'][0],sigDigits)

    if splitRay is None:
        # if no ray is specified, create subregions using E-W lines, then fit a ray through their centroids
        splitRay=makeSplitterRay(regionFrame,{x:regionList[x]['Data']['TAPersons'] for x in ctySubregions})
    intercept=splitRay[0]
    slope=splitRay[1]

    # To avoid several potential problems, the N-S bounds relevant for this process are those
    #   points where the splitter line intersects with the boundary of the region
    
    # Add a line segment of the splitRay that runs from the north end of the region to the south end
    splitLine=LineString([(intercept+nPoint*slope,nPoint),(intercept+sPoint*slope,sPoint)])
    # Find the intersections of this line with the boundary of the region

#    print(regionFrame)
    regBoundary=regionFrame.loc[codestr,'geometry'].boundary
    xpoints=regBoundary.intersection(splitLine)
    # Limit the N-S range for splits to the max range between the intersections
    nPtUse=0.0
    sPtUse=90.0
    for pts in range(len(xpoints)):
        nPtUse=max(nPtUse,xpoints[pts].y)
        sPtUse=min(sPtUse,xpoints[pts].y)

    # Split along lines perpendicular to "Split"  ray above
    # We follow the n-s direction with population-proportional splits running up the west side of the bounding square
    #
    subReg_Frame=gp.GeoDataFrame()
    subReg_Frame['subRegion']=None
    subReg_Frame['geometry']=gp.GeoSeries
    #set the geometry
    subReg_Frame=subReg_Frame.set_geometry('geometry')
    subReg_Frame.set_index('subRegion',drop=False,inplace=True)

    slice_Frame=gp.GeoDataFrame()
    slice_Frame['subRegion']=None
    slice_Frame['geometry']=gp.GeoSeries
    #set the geometry
    slice_Frame=subReg_Frame.set_geometry('geometry')


    nsplits=len(ctySubregions)

    nsHeight=nPtUse-sPtUse
    popCum=0
    subReg_sPoint=(intercept-wPoint)*slope +(1.0+slope*slope)*sPoint

    if G.verbose: print('    ctySubregions: ',ctySubregions)
    for regionID in ctySubregions:
    
        if G.verbose: print('    This region: ',regionID)
        regionData=regionList[regionID]['Data']
        popSlice=int(regionData['TAPersons'])

        nsSplit=sPtUse+nsHeight*float(popCum+popSlice)/float(regPop)

        # Compute the intercept and slope of the perpendicular line passing through the nsSplit Latitude
        if abs(slope)>1000000.0:
            nsWestSide=nsSplit
            nsEastSide=nsSplit
        else:
            pSlope=-1/slope
            pIntercept=intercept+nsSplit*(slope-pSlope)
            nsWestSide=(pIntercept-wPoint)*slope
            nsEastSide=(pIntercept-ePoint)*slope

        xp1=round(subReg_sPoint,sigDigits)
        xp2=round(nsWestSide,sigDigits)
        xp3=round(nsEastSide,sigDigits)
        xp4=round(nsEastSide-(nsWestSide-subReg_sPoint),sigDigits)

        # We have to make sure no parts of any subregion are cut off completely.
        # So, at the bottom and top, set the extremum points at the appropriate boundaries.
        if popCum<1:
            xp1=min(xp1,sPoint)
            xp4=min(xp4,sPoint)
            slice_Frame=slice_Frame.append({'subRegion':regionID+' Base','geometry':LineString([(wPoint,xp1),(ePoint,xp4)])},ignore_index=True)
        if popCum+popSlice==regPop:
            xp2=max(xp2,nPoint)
            xp3=max(xp3,nPoint)

        # Create the polygon points
        p1=(wPoint,xp1)
        p2=(wPoint,xp2)
        p3=(ePoint,xp3)
        p4=(ePoint,xp4)

        if popCum>0: 
            linexPrior=linex
            pxPrior=px

        px=Polygon([p1,p2,p3,p4])
        if popCum>0:
            tryTouches=px.touches(pxPrior)
            if not tryTouches:
                print('    Touches? False')
                # Find key distance between points
# ************************************************ woking here
                # Compute a scale for the px polygon (tuples of lat,long)
                xdist=CCFD.get_distance(p2,p3)
                # Add points to make resulting polygon more suitable for analysis
                print('xdist for Densify: ',str(xdist))
                px=DensifyPolygon(px,pieceLength=xdist/20.0)
                print('    Touches after densify?',px.touches(pxPrior))

        linex=LineString([p2,p3])

        subReg_Frame=subReg_Frame.append({'subRegion':regionID,'geometry':px},ignore_index=True)
        slice_Frame=slice_Frame.append({'subRegion':regionID+' Slice','geometry':linex},ignore_index=True)

        if G.verbose:
            print('    Plotting slices through ',regionID)
            print('    Pop values: Base:',popCum,' increment: ',int(regionData['TAPersons']))
            print('    Polygon points: ',p1,p2,p3,p4)

        # Set inputs for next loop
        popCum+=popSlice    # add this subregion population to the total
        subReg_sPoint=nsWestSide   # the southern west side point for the next slice is the northern one for this slice

#    if G.verbose:
#        pd.concat([subReg_Frame,regionFrame],ignore_index=True,sort=True).plot(cmap='Set1')
#        plt.show()

    # Create the sliced region from the slices
    subReg_Frame=gp.overlay(subReg_Frame,regionFrame,how='intersection')
    
    # check to make sure that the resulting slices are simple polygons and not multipolygons
    fixList=[]
    for regionID in ctySubregions:
        xGeo=subReg_Frame.loc[subReg_Frame['subRegion']==regionID]['geometry'].unary_union
        if not checkPolygon(xGeo):
            print('    Repairing polygon issue for: ',regionID)
            # We will probably have to merge some regions to fix this.
            # Assume that the subregion with the next higher dot.ID will play this role.
            # If the bad one has the highest ID, then use the one with the next lower ID.
            xreg,xnum=CCFD.region_and_Number(regionID)
            nextRegion=xreg+'.'+str(xnum+1)
            try:
                xtry=regionList[nextRegion]['ID']
            except:
                nextRegion=xreg+'.'+str(xnum-1)
                print('   Set srcPoly as prior:',nextRegion)

            subReg_Frame=MultiPolygon_Repair(regionID,nextRegion,subReg_Frame)
            xGeo=subReg_Frame.loc[subReg_Frame['subRegion']==regionID]['geometry']
            if not checkPolygon(xGeo.unary_union):
                   print('    Repair failed.')
                    
# ******************************************* Working here **************************************
    # Add the slice lines
    subReg_Frame=subReg_Frame.append(slice_Frame,ignore_index=True,sort=False)

    # Add the fitted line to the plot
    fitted=LineString([(intercept+slope*sPoint,sPoint),(intercept+slope*nPoint,nPoint)])
    subReg_Frame=subReg_Frame.append({'subRegion':'Fitted','geometry':LineString(fitted)},ignore_index=True)

    #Plot the results
    if G.verbose:
        print('Finished plot for: ',regionCode)
        print(subReg_Frame)
        subReg_Frame.plot(cmap='Set1')
        plt.show()

    return subReg_Frame

def getNeighbors(inDict,verbose=False,useIsNeighbors=False):

    if inDict is None: return None

    regionList=list(inDict.keys())
    x_Frame=gp.GeoDataFrame({'ID':regionList})
    x_Frame['geometry']=[inDict[val] for val in regionList]
#    for regionID in ctySubregions:
#        subReg_Frame.loc[regionID,'Neighbors']=list(set(rW[regionID].keys()))

# pysal weights approach
    rW = weights.contiguity.Rook.from_dataframe(x_Frame, idVariable = 'ID')
    neighbors=[list(rW[regionID].keys()) for regionID in inDict.keys()]
    x_Frame['Neighbors']=neighbors

# less sophisticated IsNeighbor function
    if useIsNeighbors:
        regCount=len(regionList)
        rW={}
        for x in range(regCount): rW.update({regionList[x]:set()})

        if G.verbose: print('getNeighbors. Checking: ')
        for x in range(regCount):
            regtry=regionList[x]
    #        if G.verbose: print('  ',regtry)
            reg1=inDict[regtry]
            # print(regList[x])
            for y in range(1,regCount):
                if regionList[y]!=regtry:
                    reg2=inDict[regionList[y]]
    #                if regtry=='C001.0' and regionList[y]=='C041.0':
    #                    QuickPlot({regtry:reg1,regionList[y]:reg2})
                    if IsNeighbor(reg1,reg2):
                        if verbose: print('  Adding: ',regionList[y],' <--> ',regtry)
                        rW[regtry].add(regionList[y])
                        rW[regionList[y]].add(regtry)

    neighbors=[list(rW[regionID]) for regionID in inDict.keys()]
    x_Frame['Neighbors']=neighbors
    
    x_Frame.set_index('ID',inplace=True)
    
    return x_Frame

# This is similar to the following but takes a single set element of a Member list rather than a geodataframe.
def calcRegionGeoData(inRegionSet):

    try:
        xbounds=inRegionSet['geometry'].bounds
        xcentroid=inRegionSet['geometry'].centroid
        thisGeo=gp.GeoSeries(inRegionSet['geometry'])
        thisGeo=changeCRS(thisGeo,'cartesian',zeroLongLatTuple=(xcentroid.y,xcentroid.x))
        xArea=thisGeo.area[0]*G.areaToSqMiles
        xlength=thisGeo.boundary.length[0]*G.perimeterToMiles
    except:
        # Reset All to None
        geodata=set(inRegionSet['Data'].keys())-G.regionDemoData
        for xdata in geodata:
            inRegionSet['Data'][xdata]=None
        return inRegionSet


    nLat=xbounds[3]
    sLat=xbounds[1]
    eLong=xbounds[2]
    wLong=xbounds[0]

    inRegionSet['Data']['Area']=xArea
    inRegionSet['Data']['Perimeter']=xlength
    inRegionSet['Data']['North']=nLat
    inRegionSet['Data']['South']=sLat
    inRegionSet['Data']['East']=eLong
    inRegionSet['Data']['West']=wLong
    inRegionSet['Data']['Longitude']=xcentroid.x
    inRegionSet['Data']['Latitude']=xcentroid.y

    inRegionSet['Data']['BdgSqrIndex']=xArea/CCFD.bounding_square(nLat,sLat,eLong,wLong)
    if inRegionSet['Data']['BdgSqrIndex']>1.0:
        print('  BdgSqrIndex problem:', inRegionSet['ID'],xArea,xbounds)

    return inRegionSet            

def recalcRegionData(regSet,inList):
    # The input frame has some members identified with inList, some do not.
    # The adjustments here apply only to those that do.

    for x in regSet:
        thisReg=inList[x]

        thisReg=calcRegionGeoData(thisReg)

        # Pare down neighbors as necessary. (Note: we must always start with a set containing all possible neighbors)
        # Start with all counties and subcounties we know of right now plus their neighbors
        coreSet={y[:4] for y in thisReg['Neighbors']}
        coreSet.add(x[:4])
        for y in inList.keys():
            if y[:4] in coreSet:
                thisReg['Neighbors'].update(inList[y]['Neighbors'])
        # Remove itself
        thisReg['Neighbors'].discard(x)
        for thisNeighbor in copy(thisReg['Neighbors']):
            if not IsNeighbor(thisReg['geometry'],inList[thisNeighbor]['geometry']):
                # Not a neighbor
                thisReg['Neighbors'].discard(thisNeighbor)
                inList[thisNeighbor]['Neighbors'].discard(x)
            else:
                # Is a neighbor. Add this region to its neighbor's neighborlist
                inList[thisNeighbor]['Neighbors'].discard(CCFD.coreID(x)+'.0')   # The neighbors started with with this. If we've split, this may not be correct anymore.
                inList[thisNeighbor]['Neighbors'].add(x)
        if G.verbose:
            plotDict={y:inList[y]['geometry'] for y in thisReg['Neighbors']}
            plotDict[x]=thisReg['geometry']
            QuickPlot(plotDict,[x])

        inList[x]=thisReg
    return inList

def tryUnion(regionList,memberList):

    geoList=[]
    for xreg in regionList:
        if memberList[xreg]['geometry'].exterior is not None:
            geoList.append(memberList[xreg]['geometry'])
            
    try:    
        outGeometry=unary_union(geoList)
        checkFlag=checkPolygon(outGeometry)
    except:
        nameList=list(regionList)
        print('  Error: Union failure:',{x:checkPolygon(x) for x in geoList})
        return None
    # Check to make sure this is a polygon
    if not checkFlag:
        print('  Union problem with regions ',regionList)
        QuickPlot({xreg:memberList[xreg]['geometry'] for xreg in regionList},showPlots=True)
     
        # if not, test connectedness
        xresults=CCFD.testConnectedness(regionList,memberList)
        if not xresults:
            print('  Special union attempts failed: region not connected')
            return Polygon()

        # if connected, look to see if any pieces are not polygons
        xresults=True
        for xreg in regionList:
            checkFlag=checkPolygon(memberList[xreg]['geometry'])
            if not checkFlag:
                print('  Problem with region geometry:',xreg,memberList[xreg]['geometry'].geom_type)
                xresults=False
        if not xresults:
            print ('  Exiting with geometry issues')
            return Polygon()

        # Try assembling one at a time.
        print('  Attempting sequential union ...')
        remainingList=copy(regionList)
        xreg=remainingList.pop()
        outGeometry=memberList[xreg]['geometry']
        while len(regionList)>0:
            # Select a neighbor of xreg
            testSet=memberList[xreg]['Neighbors'] & remainingList
            if len(testSet)==0:
                print('  Special sequential union attempts fails: ',xreg,'Neighbors',memberList[xreg]['Neighbors'],'found no connecting region',remainingList)
                return Polygon()
            xreg=list(testSet)[0]
            remainingList.remove(xreg)
            checkFlag=False
            trial=0
            while not checkFlag and trial<10:
                trial+=1
                outGeometry=unary_union([outGeometry,memberList[xreg]['geometry']])
                checkFlag=checkPolygon(outGeometry)
                if not checkFlag:
                    print('    Attempting fix for ',xreg)
                    geoDict[xreg]=Fix_Neighbor_Geo(xreg)
                else:
                    print('    Added',xreg)
                if trial>9: GroupPlot({{xreg:outGeometry},{x:memberList[x]['geometry'] for x in remainingList}})
                    
        if not checkFlag:
            print('  Region union failed (returning Multipolygon):',regionList)

    return outGeometry

def Region_Merger(srcRegion,badRegion,inList=None,inDist=None):
    # Combines the two regions together, zeroing only selected information about badRegion
    if inList is None: inList=G.regionMemberList
    if inDist is None: inDist=G.districtMemberList

    thisDist=list(inList[srcRegion]['Districts'])[0]
    # Merge the two, making sure beforehand that they are neighbors
    if srcRegion in inList[badRegion]['Neighbors'] and badRegion in inList[srcRegion]['Neighbors']:
        inList[srcRegion]['geometry']=inList[srcRegion]['geometry'].union(inList[badRegion]['geometry'])

        if not checkPolygon(inList[srcRegion]['geometry']):
            print('Merger problem with ',srcRegion,' and ',badRegion)

            if G.verbose: QuickPlot({x:inList[x]['geometry'] for x in [srcRegion,badRegion]})

        # Reset the neighbors
        inList[srcRegion]['Neighbors'].remove(badRegion)           # Done this way as a check: missing elements will cause this to fail
        inList[badRegion]['Neighbors'].remove(srcRegion)           # Done this way as a check: missing elements will cause this to fail
        for x in inList[badRegion]['Neighbors']:
            inList[x]['Neighbors'].remove(badRegion)
        inList[srcRegion]['Neighbors'].update(inList[badRegion]['Neighbors'])

        # Reset the data
        inList[srcRegion]['Data']['TAPersons']+=inList[badRegion]['Data']['TAPersons']

        # Delete badRegion
        del inList[badRegion]

        # Make the appropriate district Changes
        inDist[thisDist]['Regions'].remove(badRegion)
        for x in inDist[thisDist]['distNeighbors']:
            inDist[x]['Neighbors'].discard(badRegion)

        # Reset the pop in the destination region to zero (useful to check)
        inList=recalcRegionData(gp.GeoDataFrame({'subRegion':[srcRegion]},geometry=gp.GeoSeries(inList[srcRegion]['geometry'])),inList)

        print('Merged: ',srcRegion,badRegion,inList[srcRegion]['geometry'].geom_type)
        if G.verbose:
            gp.GeoSeries(inList[srcRegion]['geometry']).plot()
            plt.show()
    else:
        print('Regions disjoint: ',srcRegion,destRegion)


    return inList,inDist


def Region_Splitter(xexcess,srcRegion,destRegion,newRegion,regionMemberList=None,districtMemberList=None):
    if regionMemberList is None: regionMemberList=G.regionMemberList
    if districtMemberList is None: districtMemberList=G.districtMemberList

    srcDistrict=list(regionMemberList[srcRegion]['Districts'])[0]
    destDistrict=list(regionMemberList[destRegion]['Districts'])[0]


    # The splitter ray for this action will be the straight line between the centroids of srcRegion and destRegion
    # The boundary line between srcRegion and newRegion will be perpendicular to this line.

    splitRay=rayfromPoints((regionMemberList[srcRegion]['Data']['Longitude'],regionMemberList[srcRegion]['Data']['Latitude']),
                                        (regionMemberList[destRegion]['Data']['Longitude'],regionMemberList[destRegion]['Data']['Latitude']))

    # One of the two regions will remain a neighbor of the region in the receiving district
    # By convention, we stack the subRegions created by splitting south to north starting with zero.
    # The new region is going to be the northern of the split region.
    # So, if the latitude of destRegion is larger than the latitude of srcRegion, it's the most northerly one (newRegion).
    # Otherwise, it's the southerly one (srcRegion) and we'll have to switch the names after splitting

    regPop=regionMemberList[srcRegion]['Data']['TAPersons']

    flipRegion=False
    # By convention, the srcRegion (.0 region) will be on the south end. This is the region that will stay.
    # If we want the north end to stay, we need to give the south end the move amount of population, then recode the IDs after splitting
    if regionMemberList[destRegion]['Data']['Latitude']>regionMemberList[srcRegion]['Data']['Latitude']:
        print('  Make newRegion the northernmost.')
        regionMemberList[newRegion]['Data']['TAPersons']=xexcess
        regionMemberList[srcRegion]['Data']['TAPersons']=regPop-xexcess
    else:
        flipRegion=True
        print('  Make newRegion the southernmost.')
        regionMemberList[srcRegion]['Data']['TAPersons']=xexcess
        regionMemberList[newRegion]['Data']['TAPersons']=regPop-xexcess

    regionMemberList[newRegion]['Districts']=set([srcDistrict])

    # Remeber the extrema for the unsplit source region for possible plotting purposes
    nPoint=regionMemberList[srcRegion]['Data']['North']
    sPoint=regionMemberList[srcRegion]['Data']['South']
    ePoint=regionMemberList[srcRegion]['Data']['East']
    wPoint=regionMemberList[srcRegion]['Data']['West']
    Pt1=(regionMemberList[srcRegion]['Data']['Longitude'],regionMemberList[srcRegion]['Data']['Latitude'])
    Pt2=(regionMemberList[destRegion]['Data']['Longitude'],regionMemberList[destRegion]['Data']['Latitude'])
    # Also remember the geometry for the unsplit source region in case of failures below
    holdsrcGeometry=regionMemberList[srcRegion]['geometry']

    xFrame=subRegion_Shaper(srcRegion,{x:regionMemberList[x] for x in [srcRegion,newRegion]},splitRay=splitRay)


    if xFrame is None:
        print('Problem in region_shaper. Src: ',srcRegion,' New: ',newRegion)
        print('SplitRay: ',splitRay)

    if flipRegion:
        # This is the case where the most northern of the regions is really the source region.
        # So, we have to flip the subregion names and populations in the frame and the region list before proceeding.
        # (Clumsy: Might be able to modify subRegion_Shaper in the future to remedy this.)
        if G.verbose:
            print('Flipping:')            # when verbose interrupt here.
        else:
            print('Flipping:')
        xFrame.at[xFrame['subRegion']==newRegion,'subRegion']='xxxReg'
        xFrame.at[xFrame['subRegion']==srcRegion,'subRegion']=newRegion
        xFrame.at[xFrame['subRegion']=='xxxReg','subRegion']=srcRegion
        regionMemberList[srcRegion]['Data']['TAPersons']=regPop-xexcess
        regionMemberList[newRegion]['Data']['TAPersons']=xexcess

    # Fix newRegion and srcRegion if a multipolygon is created
    newGeoseries=xFrame.loc[xFrame['subRegion']==newRegion,'geometry'].unary_union
    if not checkPolygon(newGeoseries):
        # newRegion is a multipolygon. Not acceptable.
        print('Fixing multipolygon issue for ',newRegion)
        newGeoseries=xFrame.loc[xFrame['subRegion']==newRegion,'geometry'].explode()

        keepPolygon=Single_from_MultiPolygon(newGeoseries)

        if keepPolygon is None:
            print('  Forcing: ')
            keepPolygon=Single_from_MultiPolygon(newGeoseries,True)

        if not checkPolygon(keepPolygon):
            print('Problem with newRegion')
        else:
            xFrame.loc[xFrame['subRegion']==newRegion,'geometry']=Polygon(keepPolygon)

        # srcRegion becomes the original source region with keepPolygon cut out
        srcPolygon=holdsrcGeometry.difference(keepPolygon)

        if not checkPolygon(srcPolygon):
            print('Fixing problem with srcPolygon:',srcPolygon.geom_type)
            srcPolygon=Single_from_MultiPolygon(srcPolygon)
            #Make sure that srcPolygon

        if not checkPolygon(srcPolygon) or G.verbose:
            print('Rebuild of srcRegion failed.')
            smallDict={'newGeoseries':newGeoseries}
            smallDict.update({srcRegion:srcPolygon})
            smallDict.update({destRegion:regionMemberList[destRegion]['geometry']})
            QuickPlot(smallDict)
  
    regionMemberList=recalcRegionData(xFrame,copy(regionMemberList))


    # Recompute neighbors for all regions in srcDistrict and destDistrict
    # To start, build a small region with the member regions and all of their neighbors.
    shortList=set([srcRegion,destRegion,newRegion])
    shortList.update(districtMemberList[srcDistrict]['Regions'])
    shortList.update(districtMemberList[destDistrict]['Regions'])
    shortList.update(districtMemberList[srcDistrict]['Neighbors'])
    shortList.update(districtMemberList[destDistrict]['Neighbors'])
    shortList.update(regionMemberList[srcRegion]['Neighbors'])
    shortList.update(regionMemberList[destRegion]['Neighbors'])
    shortList.update(regionMemberList[newRegion]['Neighbors'])

    print('Neighbors recalc...')
    neighborFrame=getNeighbors({x:regionMemberList[x]['geometry'] for x in shortList})

    for x in shortList:
        oldNeighbors=regionMemberList[x]['Neighbors']
        regionMemberList[x]['Neighbors']=set(neighborFrame.loc[x,'Neighbors'])
        regionMemberList[x]['Neighbors'].discard(x)
        if len(regionMemberList[x]['Neighbors'])==0:
            print('Error. Region has no neighbors:', x)
        if G.verbose and len(regionMemberList[x]['Neighbors'].difference(oldNeighbors))>0:
            print('Regions added to neighbors of ',x,': ',regionMemberList[x]['Neighbors'].difference(oldNeighbors))
        if G.verbose and len(oldNeighbors.difference(regionMemberList[x]['Neighbors']))>0:
            print('Regions removed as neighbors of ',x,': ',oldNeighbors.difference(regionMemberList[x]['Neighbors']))
        for y in regionMemberList[x]['Neighbors']:
            if regionMemberList[y]['Neighbors'] is None: regionMemberList[y]['Neighbors']=set()
            if x not in regionMemberList[y]['Neighbors'] and G.verbose: print('  ',x,' added as neighbor to ',y)
            regionMemberList[y]['Neighbors'].update([x])
            regionMemberList[y]['Neighbors'].discard(y)


    # Check to make sure that the new region we created lies between srcRegion and destRegion
    if len(regionMemberList[newRegion]['Neighbors'] & set([srcRegion,destRegion]))!=2 and newRegion!=destRegion:
        print('Fixing location issue for ',newRegion)
        print('Neighbors check: ')
        print('src ',srcRegion,regionMemberList[srcRegion]['Neighbors'])
        print('new ',newRegion,regionMemberList[newRegion]['Neighbors'])
        print('dest',destRegion,regionMemberList[destRegion]['Neighbors'])
#        G.verbose=True
        # The new region does not lie between srcRegion and destRegion
        # Reset it
        print('Monitor how this is fixed:')
        if G.verbose: QuickPlot({x:regionMemberList[x]['geometry'] for x in [srcRegion,destRegion,newRegion]})

        for x in [srcRegion,destRegion]:
            if newRegion not in regionMemberList[x]['Neighbors']:
                print('Try fixing neighbor:',newRegion,x)
                regionMemberList[newRegion]['geometry']=Fix_Neighbor_Geo(newRegion,regionMemberList)
        
        nSet=set([newRegion,srcRegion,destRegion])
        for x in [newRegion,srcRegion,destRegion]:
            nSet.update(regionMemberList[x]['Neighbors'])
        neighborFrame=getNeighbors({x:regionMemberList[x]['geometry'] for x in nSet})
        for x in [srcRegion,newRegion,destRegion]:
            newset=set(neighborFrame.loc[x,'Neighbors'])
        if newRegion not in regionMemberList[srcRegion]['Neighbors'] or newRegion not in regionMemberList[destRegion]['Neighbors']:
            print(neighborFrame)
            QuickPlot({x:regionMemberList[x]['geometry'] for x in nSet},newRegion)
            print('Location repairs failed.')


        # Now create a modest intersection between srcRegion and destRegion buffered.
        tryBuffer=0.01
        findGeo=xFrame.loc[xFrame['subRegion']==srcRegion,'geometry']
        if len(findGeo)>1:
            print('Problem (multiple rows in xFrame) with ',srcRegion)
        else:
            srcGeo=findGeo.unary_union
        # Iterate a bit.
        findGeo=xFrame.loc[xFrame['subRegion']==newRegion,'geometry']
        if len(findGeo)>1:
            print('Problem (multiple rows in xFrame) with ',newRegion)
        else:
            newGeo=findGeo.unary_union
        # Iterate a bit.
        for count in range(3):
            sliceGeo=regionMemberList[destRegion]['geometry'].buffer(tryBuffer).intersection(srcGeo)
            restGeo=srcGeo.difference(sliceGeo)
            # Adjust the buffer size based on areas and populations
            if sliceGeo.area>0.00001:
                bufferShare=sliceGeo.area/srcGeo.area
            else:
                bufferShare=0.5
            tryBuffer*=(regionMemberList[newRegion]['Data']['TAPersons']/regionMemberList[srcRegion]['Data']['TAPersons'])/bufferShare
            if restGeo.geom_type!='Polygon':
                print()
        # Make the slice big enough to create a single polygon
        count=0
        tryBuffer=0.01
        while (sliceGeo.geom_type!='Polygon'or restGeo.geom_type!='Polygon') and count<100:
            tryBuffer*=1.15
            sliceGeo=regionMemberList[destRegion]['geometry'].buffer(tryBuffer).intersection(srcGeo)
            count+=1
        if sliceGeo.geom_type!='Polygon' or restGeo.geom_type!='Polygon':
            print('Polygon slice problem found for ',srcRegion)
            print('Src: ',srcGeo.geom_type,'Slice: ',sliceGeo.geom_type,' Rest: ',restGeo.geom_type)
            xFrame=xFrame.append({'subRegion':'restGeo','geometry':restGeo},ignore_index=True)
            xFrame=xFrame.append({'subRegion':'spliceGeo','geometry':sliceGeo},ignore_index=True)
            print(xFrame)
            # Create a focused plot of the first one that is not a polygon
            smallFrame=xFrame.loc[xFrame['geometry'].geom_type!='Polygon']
            smallFrame.plot(cmap='Set1')
            plt.show()


        # Replace the geometries in xFrame:
        xFrame.at[xFrame['subRegion']==srcRegion,'geometry']=Polygon(restGeo)
        xFrame.at[xFrame['subRegion']==newRegion,'geometry']=Polygon(sliceGeo)

        # Recompute the neighbors lists
        neighborFrame=getNeighbors({x:regionMemberList[x]['geometry'] for x in shortList})
        for x in shortList:
            regionMemberList[x]['Neighbors'].update(set(neighborFrame.loc[x,'Neighbors']))
            regionMemberList[x]['Neighbors'].discard(x)
            for y in regionMemberList[x]['Neighbors']:
                regionMemberList[y]['Neighbors'].update([x])
                regionMemberList[y]['Neighbors'].discard(y)

    regionMemberList=recalcRegionData(xFrame,copy(regionMemberList))
    print('District adjustment of: ',xexcess,' from ',srcDistrict,'(split ',srcRegion,':',regionMemberList[srcRegion]['Name'],') to ',destDistrict,'(at ',destRegion,':',regionMemberList[destRegion]['Name'],')')

    if G.verbose:
        # Start with the elements of xFrame we want to keep
        criterion=xFrame['subRegion'].isin([srcRegion,newRegion,srcRegion+' Slice'])
        plot_Frame=xFrame.loc[criterion]

        # Add the relevant surrounding regions
        for x in [x for x in regionMemberList.keys() if x[:4]==srcRegion[:4] and x not in [newRegion,srcRegion]]:
            plot_Frame=plot_Frame.append({'subRegion':x,'geometry':regionMemberList[x]['geometry']},ignore_index=True)
            nPoint=max(nPoint,regionMemberList[x]['Data']['North'])
            sPoint=min(sPoint,regionMemberList[x]['Data']['South'])
            ePoint=max(ePoint,regionMemberList[x]['Data']['East'])
            wPoint=min(wPoint,regionMemberList[x]['Data']['West'])
                          
        for x in [x for x in regionMemberList.keys() if x[:4]==destRegion[:4]]:
            plot_Frame=plot_Frame.append({'subRegion':x,'geometry':regionMemberList[x]['geometry']},ignore_index=True)
            nPoint=max(nPoint,regionMemberList[x]['Data']['North'])
            sPoint=min(sPoint,regionMemberList[x]['Data']['South'])
            ePoint=max(ePoint,regionMemberList[x]['Data']['East'])
            wPoint=min(wPoint,regionMemberList[x]['Data']['West'])

        # Add the splitter ray between the centroids
        plot_Frame=plot_Frame.append({'subRegion':'Splitter','geometry':LineString([Pt1,Pt2])},ignore_index=True)

        # Resize the slice line for plotting purposes
        xpoints=list(xFrame.loc[xFrame['subRegion']==srcRegion+' Slice','geometry'].unary_union.coords)
        intercept,slope=rayfromPoints(xpoints[0],xpoints[1])
 
        boundaryLine=LineString([(intercept+slope*(sPoint-1.0),sPoint-1.0),(intercept+slope*(nPoint+1.0),nPoint+1.0)])
        plotBox=Polygon([(wPoint,sPoint),(wPoint,nPoint),(ePoint,nPoint),(ePoint,sPoint)])
        xplotBox=plotBox.intersection(boundaryLine).bounds
        if slope>0.0:
            sliceLine=LineString([(xplotBox[0],xplotBox[1]),(xplotBox[2],xplotBox[3])])
        else:
            sliceLine=LineString([(xplotBox[0],xplotBox[3]),(xplotBox[2],xplotBox[1])])

        plot_Frame.set_index('subRegion',inplace=True)
        plot_Frame.at[srcRegion+' Slice','geometry']=sliceLine

        print('Used to move: ',str(regionMemberList[newRegion]['ID']),' with ',str(regionMemberList[newRegion]['Data']['TAPersons']))
        print('Remaining in: ',regionMemberList[srcRegion]['ID'],' to ',regionMemberList[srcRegion]['Data']['TAPersons'],' from ',str(regPop))

        plot_Frame.plot(cmap='Set1',legend=True)
        plt.show()


    return regionMemberList

