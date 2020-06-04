# A package that parse mathematica into sympy
# Worked with sympy community. Python 3.6 or higher is required
from sympy.parsing import mathematica

# Usage 
# put your exprseeion into mathematica.parse('<expression>')
# I am struggling with how index notation will be understood as index 
# not a power or subscript

res1=mathematica.parse('Sin[x]^2 + \[Delta] + x^3')
print(res1)
