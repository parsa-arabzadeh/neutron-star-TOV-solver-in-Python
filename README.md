# TOV Neutron Star Solver

A Python implementation of the **Tolman–Oppenheimer–Volkoff (TOV) equations** to model the structure of neutron stars in General Relativity.

This project numerically solves the relativistic hydrostatic equilibrium equations and generates the mass–radius relation of neutron stars using a polytropic equation of state.

---

## Overview

Newtonian hydrostatic equilibrium is not sufficient to describe extremely compact objects such as neutron stars. The TOV equations provide the relativistic extension of stellar structure equations by including General Relativistic effects.

In this project, the TOV equations are solved using numerical integration to obtain:

- Pressure profile inside the star
- Mass distribution as a function of radius
- Neutron star mass–radius relation
- Maximum stable mass for a given equation of state

---

## Physical Model

The stellar structure is determined by solving the TOV equations:

### Mass equation

$$
\frac{dm}{dr}=4\pi r^2\frac{\epsilon}{c^2}
$$


### Pressure equation

The pressure gradient is:

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


where:

- \(m(r)\) is the enclosed gravitational mass
- \(P(r)\) is the pressure
- \(\epsilon\) is the energy density
- \(G\) is the gravitational constant
- \(c\) is the speed of light

---

## Equation of State

A polytropic equation of state is used:

$$
P=K\rho^\gamma
$$

where:

- \(K\) is the polytropic constant
- \(\gamma\) is the adiabatic index
- \(\rho\) is the rest-mass density

The energy density is calculated as:

$$
\epsilon=\rho c^2+\frac{P}{\gamma-1}
$$

---

## Numerical Method

The equations are solved using:

- Python
- NumPy
- SciPy (`solve_ivp`)
- Matplotlib

The integration stops automatically when the pressure reaches zero, corresponding to the surface of the neutron star.

---

## Result

The code produces the mass-radius relation of neutron stars and estimates the maximum stable mass for the chosen EOS.

## Author

Parsa Arabzadeh Asadi

MSc Physics  
University of Tehran
