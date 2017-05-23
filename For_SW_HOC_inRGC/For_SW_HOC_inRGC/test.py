import sympy
sympy.init_printing()
x = sympy.Symbol('x')
sympy.series((1-sympy.exp(-x))/(1+sympy.exp(-x)),x)



