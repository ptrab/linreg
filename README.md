# linreg

This program calculates the fitting parameters for a linear least-squares problem.
However, it is not a "simple" regression in a way of *y* over *x*, which is fitted by `y_i = a + b * x_i`,
but a fit of many *y_ij* over their resp. *x_i*, `y_ij = a + b * c_j * x_i`.

In other words, if you have values to be fitted that can be described by the same fitting parameters A and B
and the slope only differs by a parameter, you can fit all together.

![linear fit of four measure series](/abx1.png?raw=true)
