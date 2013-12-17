Mesh|Mesh(\texttt{(name)},$\cdot$)
macro file name|string||Filename for the macrofile
value file name|string||Filename for value file. 
periodic file|string||filename for periodic file
check|int|1|Check the mesh structure
preserve coarse dofs|bool|0|When an element is refined, not all dofs of the coarse element must be part of the new elements. An example are centered dofs when using higher lagrange basis functions. The midpoint dof of the parents element is not a dof of the both children elements. Therefore, the dof can be deleted. In some situation, e.g., when using multigrid techniques, it can be necessary to store this coarse dofs. Then this variable must be set to true. If false, the not required coarse dofs will be deleted.
macro weights|string||file of file that contains weight number for each element. The higher the weigh, the more this element will be refined and thus more processors share this element.
global refinements|int|0|Number of global refinements