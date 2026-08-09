# TOV Neutron Star Solver

A Python implementation of the Tolman–Oppenheimer–Volkoff (TOV) equations for modeling the structure of neutron stars.

## Overview

This project solves the relativistic stellar structure equations:

\[
\frac{dm}{dr}=4\pi r^2\frac{\epsilon}{c^2}
\]

and

\[
\frac{dP}{dr}
=
-\frac{G(\epsilon+P)
(m+4\pi r^3P/c^2)}
{c^2r^2(1-2Gm/rc^2)}
\]

using numerical integration.

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
