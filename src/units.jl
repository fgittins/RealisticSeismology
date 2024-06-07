"""
Physical constants and conversion factors.

# Constants
- `c`
- `G`
- `Msol`
- `eV`
- `mₙ`
- `length_geometric_to_CGS`
- `mass_geometric_to_CGS`
- `time_geometric_to_CGS`
- `pressure_geometric_to_CGS`
- `mass_geometric_to_Msol`
- `pressure_geometric_to_natural`

# Notes
Constants are taken from NIST [1] and given in Centimetre-Gram-Second (CGS)
system of units.

Geometric units assume G = c = 1 and measure all variables in powers of km.
Natural units assume hbar = c = kB = 1.

# References
[1] https://physics.nist.gov/cuu/Constants/index.html.
"""

"Speed of light in vacuum [cm s^-1]."
const c = 29979245800.0
"Gravitational constant [cm^3 g^-1 s^-2]."
const G = 6.67430e-8
"Solar mass [g]."
const Msol = 1.98841e33
"Electron volt [erg]."
const eV = 1.602176634e-12
"Neutron mass [g]."
const mₙ = 1.67492749804e-24

"Length conversion [cm km^-1]."
const length_geometric_to_CGS = 1e5
"Mass conversion [g km^-1]."
const mass_geometric_to_CGS = 1e5*c^2/G
"Time conversion [s km-1]."
const time_geometric_to_CGS = 1e5/c

"Pressure conversion [erg cm^-3 km^2]."
const pressure_geometric_to_CGS = (mass_geometric_to_CGS
                                   /length_geometric_to_CGS
                                   /time_geometric_to_CGS^2)

"Mass conversion [Msol km^-1]."
const mass_geometric_to_Msol = mass_geometric_to_CGS/Msol

"Pressure conversion [MeV fm^-3 km^2]."
const pressure_geometric_to_natural = pressure_geometric_to_CGS/(eV*1e45)
