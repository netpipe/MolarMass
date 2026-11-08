# MolarMass
molar mass chemical formula calculator


Normal formulas: H2O, C6H12O6, NaCl
Nested groups: Ca(OH)2, Al2(SO4)3, K4[Fe(CN)6]
Parentheses, square brackets, braces: (), [], {}
Hydrates/adducts: CuSO4.5H2O, CuSO4·5H2O, CuSO4*5H2O
Common isotopes: D, T
Charge signs are mostly ignored: NH4+, SO4^2-


Element symbols are case-sensitive:
Co = cobalt
CO = carbon + oxygen
Leading coefficients on the first formula are ignored:
2H2O gives the molar mass of H2O, not twice it.
Hydrate coefficients after a dot are counted:
CuSO4.5H2O includes 5 water molecules.
For ions, plain text charge notation can be ambiguous. These work well:
NH4+
SO4^2-
Ca^2+
Plain Ca2+ is ambiguous because 2 could be interpreted as a subscript. If you want calcium ion mass, use Ca or Ca^2+.