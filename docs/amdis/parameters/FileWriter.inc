FileWriter|Parameters for Data output. Typically the label ist \texttt{(name)}:=\texttt{(problem-name)->output}|\renewcommand{\basis}{{(problem-name)->output->}}
filename|string||Used filename prefix.
AMDiS format|bool|0|0: Don't write AMDiS files; 1: Write AMDiS files. 
AMDiS mesh ext|string|.mesh|AMDiS mesh-file extension.
AMDiS data ext|string|.dat|AMDiS solution-file extension.
ParaView format|bool|0|0: Don't write ParaView files; 1: Write ParaView files.
ParaView vector format|bool|0|0: Don't write ParaView vector files; 1: Write ParaView vector files.
ParaView animation|bool|0|0: Don't write ParaView animation file; 1: Write ParaView animation file.
ParaView ext|string|.vtu|VTK file extension.
Periodic format|bool|0|0: Don't write periodic files; 1: Write periodic files.
Periodic ext|string|.per|Periodic file extension.
PNG format|bool|0|0: Don't write png files; 1: Write png image files.
PNG type|int|0|0: Gray color picture; 1: RGB picture.
append index|int|0|0: Don't append time index to filename prefix, 1: Append time index to filename prefix.
index length|int|5|Total length of appended time index.
index decimals|int|3|Number of decimals in time index.
write every i-th timestep|int|1|Timestep modulo: write only every tsModulo-th timestep! 
Povray format|bool|0|0: Don't write Povray scripts; 1: Write Povray scripts
Povray template|string||name of the template file that will be prepended to all created *.pov files 
Povray camera location|string||camera position for povray script files
Povray camera look\_at|string||orientation for camera in povray script files
DOF format|bool|0|0: Don't write DOF files; 1: Write DOF files
ARH format|bool|0|0: Don't write ARH files; 1: Write ARH files
compression|enum \{\texttt{gz, gzip, bz2, bzip2}\}||