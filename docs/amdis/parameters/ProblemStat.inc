ProblemStat|Standard ProblemStat(\texttt{(name)}). First construktor argument is \texttt{(name)}
components|int|-1|Number of problem components (must be set)
input->read serialization|int|0|
input->serialization with adaptinfo|int|0|
input->serialization filename|string||
mesh|string||Name of the mesh
dim|int|0|problem dimension
refinement set[i]|int|-1|$i=0,\ldots,$\texttt{nComponents},
polynomial degree[i]|int|1|$i=0,\ldots,$\texttt{nComponents},
solver|enum \{\texttt{cg, cgs, bicg, bicgstab, bicgstab2, bicgstab\_ell, qmr, tfqmr, gmres, idr\_s, minres, (umfpack), 0}\}|0|iterative/direct solver for the linear system
estimator[i]|enum \{\texttt{residual, simple-residual, recovery, 0}\}|0|$i=0,\ldots,$\texttt{nComponents}, estimator type for each components
output->write serialization|int|0|write serialization files