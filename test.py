from escvv.equation.helper import parse_equation, from_gibbs

if __name__ == "__main__":    
    eq = parse_equation('AgNO3(s) === Ag(s) + NO2(g) + 0.5O2(g)')
    eq.equilibrium_K = from_gibbs(156.9, 0.2447, 621)
    print(eq)
