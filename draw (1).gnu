#set cbrange [-10:100]
set pm3d #depthorder
#set isosamples 100
#set samples 100

set view xrot,zrot 
#splot x+2*y with pm3d

splot  "temp.data" using 1:2:3 w pm3d


