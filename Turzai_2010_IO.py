import Global_Settings as G
import Global_IO as F
import Shape_Handler as CCFD_shapes

# The part below is saved in individual txt files
# Note: the ID is the first element in the variable list

# Source: "Turzai" files from 2011 PA redistricting effort (Project redmap?)
# Google drive: Shared with me: pa-redistricting/data/Turzai data/shapefiles
# Pulled VTDs on: 8/6/2019
F.dataPath=G.dataPath

# Set the CRS for the input shapefiles
CCFD_shapes.defaultCRS=CCFD_shapes.setEPSG(4269)


# Input Detailed Attributes (at State level)

F.stateinputAreaScaleFactor=1.0/1000000.0            # Converts raw area to sq km equivalent
F.stateShapeFile='tl_2010_42_state10.SHP'
F.stateCoder=('rawRegionCode',1)
F.stateVariableDict={'rawRegionCode':'GEOID10',
                 'Name':'NAME10',
                 'CountyRawCode':'STATEFP10',
                 'ALAND':'ALAND10',
                 'AWATER':'AWATER10',
                 'INTPTLAT':'INTPTLAT10',
                 'INTPTLON':'INTPTLON10',
                 'geometry':'geometry'
                }

# Input Detailed Attributes (at County level)
F.countyinputAreaScaleFactor=1.0/1000000.0            # Converts raw area to sq km equivalent
F.countyShapeFile='Turzai - 01652.SHP'
F.countyCoder=('CountyRawCode',3)
F.countyVariableDict={'rawRegionCode':'GEOID10',
                 'Name':'NAMELSAD10',
                 'CountyRawCode':'COUNTYFP10',
                 'ALAND':'ALAND10',
                 'AWATER':'AWATER10',
                 'INTPTLAT':'INTPTLAT10',
                 'INTPTLON':'INTPTLON10',
                 'Size':'TAPERSONS',
                 'geometry':'geometry'
                }

# Input Detailed Attributes (at Municipal level)
F.municipalinputAreaScaleFactor=1.0/1000000.0            # Converts raw area to sq km equivalent
F.municipalShapeFile='Turzai - 01643.SHP'
F.municipalCoder=('MuniRawCode',5)
F.municipalVariableDict={'rawRegionCode':'GEOID10',
                 'Name':'NAMELSAD10',
                 'CountyRawCode':'COUNTYFP10',
                 'MuniRawCode':'COUSUBFP10',
                 'ALAND':'ALAND10',
                 'AWATER':'AWATER10',
                 'INTPTLAT':'INTPTLAT10',
                 'INTPTLON':'INTPTLON10',
                 'Size':'TAPERSONS',
                 'geometry':'geometry'
                }

# Input Detailed Attributes (at VTD level)
F.VTDinputAreaScaleFactor=9430.              # Converts raw area to sq km equivalent
# Note the raw area in this file is expressed in sq degrees.
F.VTDShapeFile='VTDs 2010 Plus.SHP'
F.VTDCoder=('VTDRawCode',4)
F.VTDVariableDict={'rawRegionCode':'CORRECTSTF',
                 'Name':'NAMELSAD',
                 'CountyRawCode':'COUNTYFP10',
                 'MuniRawCode':'COUSUBFP10',
                 'VTDRawCode':'rawVTDCode',
                 'ALAND':'SHAPE_AREA',
                 'INTPTLAT':'Y_COORD',
                 'INTPTLON':'X_COORD',
                 'Size':'TAPERSONS',
                 'geometry':'geometry'
                }

# Special Adjustments for this data source

# Source regions to exclude
G.regionExclusions=[('M02906544B.0','M02906544A.0'),('C029B.0','C029A.0')]    # (deleted,received demo data) Unpopulated piece of Birmingham Township, Chester county

# Source graph edges to remove
G.neighborExclusions=[('C013','C033'),('C031','C019'),('C043','C067'),('C097','C099'),('C003','C005'),('M04512442','M02906544A'),('M02906544A','000000000'),
                    ('M101WD053','M101WD023'),('M101WD053','M101WD062'),('M00332832','M00527784')]

# Source graph edges to add
G.neighborForced=[('M101WD034.0','M09144976.0')]
G.neighborForced.extend([('V045476161720.0', 'V045476161710.0'),
                         ('V045132120670.0', 'V045033360080.0'),
                         ('V061363680190.0', 'V061363680220.0'),
                         ('V107465840535.0', 'V107700561345.0')])


