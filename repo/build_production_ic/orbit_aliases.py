"""
Orbit alias configuration for build_production_ic.

Define top-level string variables with exactly six whitespace-separated values.
Tuples/lists with six numbers are also accepted.
When a case's first six orbit values in mergegalaxy.par match one of these
definitions, the output will store the variable name in the orbit_label column.

Example:
RPF1 = "-12.7   220.1   -94.9  -23.7   -67.2   -68.4"
"""

LowOrbitRPF1 = "-12.7   220.1   -94.9  -23.7   -67.2   -68.4"
HighOrbitRPF1 = "-38.7   271.2  -197.9   -8.0  -115.4   -2.0"

MW51 = "median mass MW, low density CGM"
MW23 = "median mass MW, low densityCGM"
MW41 = "low mass MW, high density CGM"
MW43 = "median mass MW, low densityCGM"
MW44 = "median mass MW, high densityCGM"
MW53 = "median mass MW, high densityCGM"

sigma_Fornax = 11
d_Fornax = 139.6
median_fd = 8.3