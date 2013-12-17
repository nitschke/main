Meshdistributor|The label \texttt{(name)} represents the string that is passed to the constructor of MeshDistributor.
repartitioning|bool|0|En/disables repartitioning of the macro mesh, when derivation of mean number of DOFs exceeds a threshold value.
repartition ith change|int|20|Sets number of mesh changes to wait before threshold check for repartitioning will be performed.
partitioner|enum \{\texttt{parmetis, zoltan, simple}\}|parmetis|Defines the external tool that performs the partitioning of the, e.g.~by graph-partitioning. {\em simple} does not change the initial partitioning, i.e., its a random distribution of the macro elements to the processors.
box partitioning|bool|0|If the macro mesh is globally refined from macro.stand.2d or macro.stand.3d, than the partitioner may compute the partitioning not based on triangled or tetrahedron, but on (composed) rectangles or boxes. Till now this is implemented only for 3D and Zoltan partitioner.
log main rank|bool|0|If set to {\em true}, stdout output will be printed only by the main rank 0. Otherwise, all ranks stdout output will be created.
pre refine|int|-1|If set to -1, the number of pre refinements for the macro mesh will be calculated for the given number of processors. This value can be overwritten by setting a value $>= 0$.
output->serialization filename|string||Name of the parallel serialization file. If at least one stationary problem is serialized, this parameter must be set.
input->serialization filename|string||Name of the parallel deserialization file. If at least one stationary problem is deserialized, this parameter must be set.
debug output dir|string||Path name where debug data should be written to.
write parallel debug file|bool|0|If set to {\em true}, the parallelization will create for each rank a file with the name ``mpi-dbg-\{rank-no\}.dat''. This files contain information about all DOF indices in ranks domain. They may be useful for debugging or some postprocessing steps.