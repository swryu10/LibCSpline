# About
A C++ library to perform cubic spline interpolations up to 3D.
Discrete data points are interpolated with given first derivatives at the boundaries.
Alternatively, natural cubic spline (i.e., vanishing second derivatives at the boundaries) can be carried out.
One can use OpenMP to expedite the initialization process in the case of 2D or 3D interpolations.

# Build
This library can be built with **cmake**. \
In a **Linux/UNIX** system, one can build at a subdirectory with the following commands. \
&ensp;$ mkdir [subdirectory name] \
&ensp;$ cd [subdirectory name] \
&ensp;$ cmake [directory for the LibCSpline local repository] \
&ensp;$ cmake --build .
