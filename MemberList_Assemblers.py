# Standard Entry
import sys
import json
import importlib
from copy import copy,deepcopy

import Global_Settings as G
import Global_IO as F
import CCFD_Functions as CCFD
import MemberList_Creators as CCFD_inputs


print()
print('Begin MemberList Assemblers ...')

# Sets info about the raw input files.
sys.path.append('L:\\\\HessViao\\Hess Work\\Districting\\MemberList Assemblers Raw Data\\')
importlib.import_module(F.this_IO_File)

# Open boundary dictionary if it exists
try:
    CCFD_inputs.tryOpenBoundaryDict()
except:
    print('No boundary dictionary found. Creating...')
    G.setGlobal('boundaryLengthDict',{})
    G.memberListOpenFlags['boundaryDict']=True

# Input listLevel can be one of: ['state','county','municipal','VTD']
# Note: in this implementation, the countyMemberList must be created before the stateMemberList
# To rerun: select a subset of the list
for x in ['county','municipal']:
    print('Creating',x,'Member List ...')
    thisMemberList=CCFD_inputs.CreateList(x,Exit=False)
    CCFD_inputs.boundaryRepair(x,ctyList=None,inBoundaries=False,inBad=False)
    print(x,'MemberLists and boundaryLengthDict created.')


