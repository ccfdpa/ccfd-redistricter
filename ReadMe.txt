Using the memberList Assemblers

The entry point for this process is the file named MemberList Assemblers.py. Some documentation
for the various packages used in the assembly process appear in the .py code.

The memberList Assemblers code translates raw geographical and demographic information about
political subdivisions for a state or sub-state region into a standard form used by the
rest of the system. Currently, the code translates and formats the data in standard form
then stores them in JSON formatted files.

The code is assembled using Python 3.7. I requires import of the following public python
packages:
	os, sys, time, copy, operator, pathlib, socket
	math, numpy, random, statistics
	matplotlib.pyplot
	fiona, json, pandas, geopandas, pysal.lib
	pyproj, shapely, shapely.ops, shapely.geometry

There are two key files that the user must modify as required to accomplish all the district
creation processes:

	Global_settings.py 	-- sets key parameters and global variables for processing
	[template_IO].py	-- sets key filename and variable maps for IO purposes

Not all settings apply to the Assembler stage. This section descributes only those that will
apply to the raw data input steps.

Global_settings.py

Navigate to the section titled: # Key folder navigation for all files
This section allows for the creation of path information for any computers on the network
that might be used for this process. Doing so requires the articulation of numerous items:

	hostName -- the name of the computer
	thisRoot -- the path to all operating files and subdirectories
	thisPath -- the name of the subdirectory (with trailing \\) holding solutions
	dataPath -- the name of the subdirectory (with trailing \\) holding raw data
	dataInitID -- an ID that uniquely idenfies this dataset

Example:
	# Key folder navigation for all files
	global thisPath,dataPath,thisRoot
	if socket.gethostname()=='Computer 1':
	    thisRoot='N:\\Work\\Districting\\'
	elif socket.gethostname()=='Computer 2':
	    thisRoot='L:\\User1\\Work\\Districting\\'
	else:
 	   thisPath='C:\\'

	thisPath=thisRoot+'MemberList Assemblers Solutions\\'
	dataPath=thisRoot+'MemberList Assemblers Raw Data\\'
	dataInitID='2010'


[template_IO].py

Navigate to the section titled: # IO mapping (small special purpose py file)
This section contains a single line which defines the file name of the .py file that
articulates the specific source files, including geometry, and input variables names for the
region precedence levels (e.g. state, county, municipal, VTD). These global variables carry the F. prefix. 
For each of these precedence levels, the IO mapping file requires the input of:

	[level]inputAreaScaleFactor -- constant factor that converts the input area to sq km equivalent
	[level]ShapeFile            -- the name of the file with geometry containing the information
					This file nmust be located in the "dataPath" directory identified above.
	[level]Coder                -- Not currently used. Reserved for generalizing this process.
	[level]VariableDict         -- A dictionary mapping the raw input variables with the variable nmemonics used internally
					The keys for this dictionary identify the internal variables, and the values identify the raw input variables.

Example:
(Note: a complete template file titled "Turzai_2010_IO.py" is included here.)

	# Input Detailed Attributes (at VTD level)
	F.VTDinputAreaScaleFactor=9430.              # Converts raw area to sq km equivalent
	# Note the raw area in this file is expressed in sq degrees.
	F.VTDShapeFile='VTDs 2010 Plus.SHP'
	F.VTDCoder=('VTDRawCode',4)
	F.VTDVariableDict='rawRegionCode':'CORRECTSTF',
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

