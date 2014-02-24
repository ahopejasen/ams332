#!/bin/bash

#insert copyright notice into all files in curdir and below

toolsDir=${0%/*}
#echo ${toolsDir}
copyrightFile="${toolsDir}/copyrightheader.txt"
#copyrightFile="${toolsDir}/copyrightheader.txt"



#testing
#find . -name '*.m' -print0 | xargs -0 sed  "/<copyright>/,/<\/copyright>/ {\%</\?copyright>%! d}
							 #/<copyright>/ r ${copyrightFile}" |less


#inplace:
find . -name '*.m' -print0 | xargs -0 sed  -i.copybak "/<copyright>/,/<\/copyright>/ {\%</\?copyright>%! d}
							 /<copyright>/ r ${copyrightFile}"

