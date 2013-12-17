|Parameters for AdaptInfo for each component separately. \texttt{(name):=(name)[$i$]} for $i=0,\ldots,$\texttt{nComponents}$-1$|\renewcommand{\basis}{{(name)[i]->}}
tolerance|double|0.0|Tolerance for the (absolute or relative) error
time tolerance|double|0.0|Time tolerance.
coarsen allowed|int \{0,1\}|0|true if coarsening is allowed, false otherwise.
refinement allowed|int \{0,1\}|1|true if refinement is allowed, false otherwise.
refine bisections|int|1|parameter to tell the marking strategy how many bisections should be performed when an element is marked for refinement; usually the value is 1 or DIM
coarsen bisections|int|1|parameter to tell the marking strategy how many bisections should be undone when an element is marked for coarsening; usually the value is 1 or DIM
sum factor|double|1.0|factors to combine max and integral time estimate
max factor|double|0.0|factors to combine max and integral time estimate