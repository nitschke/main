HL\_SignedDistTraverse|Reinitialization class HL\_SignedDistTraverse(\texttt{(name)},\ldots) with \texttt{(name)} as first argument in the constructor.
tolerance|double||Tolerance for Hopf-Lax update iteration loop.
maximal number of iteration steps|int||Maximal number of mesh iterations for Hopf-Lax update.
Gauss-Seidel iteration|bool||Indicates whether Gauss-Seidel or Jacobi iteration is used. 0: Jacobi, 1: Gauss-Seidel
infinity value|double||Initialization value ``inifinity'' for non-boundary vertices. Must be $>$ 1000
boundary initialization|int||Define boundary initialization strategy. 0: BoundaryElementLevelSetDist, 1: BoundaryElementTopDist, 2: BoundaryElementEdgeDist, 3: BoundaryElementNormalDist