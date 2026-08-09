# TOV Neutron Star Solver

A Python implementation of the Tolman–Oppenheimer–Volkoff (TOV) equations for modeling the structure of neutron stars.

## Overview

The pressure gradient is given by the Tolman–Oppenheimer–Volkoff equation:

$$
\frac{dP}{dr}
=
-\frac{
G(\epsilon+P)
\left(
m+\frac{4\pi r^3P}{c^2}
\right)
}
{
c^2 r^2
\left(
1-\frac{2Gm}{rc^2}
\right)
}
$$

The equations are solved using numerical integration.

## Features

- Numerical solution of TOV equations
- Polytropic equation of state
- Mass-radius relation calculation
- Maximum neutron star mass estimation
- Visualization using Matplotlib

## Methods

The equations are solved using:

- Python
- SciPy `solve_ivp`
- NumPy
- Matplotlib

## Result

The code produces the mass-radius relation of neutron stars and estimates the maximum stable mass for the chosen EOS.

## Author

Parsa Arabzadeh Asadi
MSc Physics  
University of Tehran
