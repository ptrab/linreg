project: Linreg
project_github: https://github.com/ptrab/linreg/blob/master/fortran
author: Philipp Traber
author_description: I am a doctoral student from Jena, Germany.
date: April 18, 2019
github: https://github.com/ptrab
graph: true
summary: This program calculates the fitting parameters for a linear least-squares
problem. However, it is not a *simple* regression in a way of `y` over `x`, which
is fitted by `y_i = a + b * x_i`, but a fit of many `y_ijk` over their resp. `x_ik`,
`y_ijk = a + b * c_ij * x_ik`.  
In other words, if you have values to be fitted that can be described by the same
fitting parameters `A` and `B` and the slope only differs by a parameter, you can
fit all together.
