AdaptInstationary|\texttt{(name)} ist the first argument of the constructor: AdaptInstationary(\texttt{(name)},\ldots)
strategy|int|0|Strategy for choosing one timestep: strategy 0: Explicit strategy, strategy 1: Implicit strategy.
time delta 1|double|0.7071|Parameter $\delta_1$ used in time step reduction
time delta 2|double|1.4142|Parameter $\delta_2$ used in time step enlargement
info|int|10|Info level (from AdaptBase)
break when stable|int|0|If this parameter is 1 and the instationary problem is stable, hence the number of solver iterations to solve the problem is zero, the adaption loop will stop.
time adaptivity debug mode|bool|0|In debug mode, the adapt loop will print information about timestep decreasing and increasing.
queue->runtime|int|-1|Runtime of the queue (of the servers batch system) in seconds. If the problem runs on a computer/server without a time limited queue, the value is -1.
queue->serialization filename|string|\_\_serialized\_problem.ser|Name of the file used to automatically serialize the problem.