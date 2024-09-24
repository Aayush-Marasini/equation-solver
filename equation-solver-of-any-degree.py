import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import math  

def solve_linear(y):
    if y[1] == 0:
        return None  
    soln = (-y[0]) / y[1]
    return soln

def solve_quadratic(y):
    a, b, c = y
    discriminant = b**2 - 4*a*c
    
    if discriminant > 0:
        x1 = (-b + math.sqrt(discriminant)) / (2*a)
        x2 = (-b - math.sqrt(discriminant)) / (2*a)
        return f"Two real solutions: x1 = {x1}, x2 = {x2}"
    elif discriminant == 0:
        x = -b / (2*a)
        return f"One real solution: x = {x}"
    else:
        real_part = -b / (2*a)
        imaginary_part = (abs(discriminant)**0.5) / (2*a)
        return f"Complex solutions: x1 = {real_part} + {imaginary_part}i, x2 = {real_part} - {imaginary_part}i"

def solve_higher_degree_newton_raphson(coefficients, initial_guess, tolerance=1e-6, max_iterations=100):
    x = sp.symbols('x')
    equation = sum(sp.S(coeff) * x**i for i, coeff in enumerate(coefficients[::-1]))
    derivative = sp.diff(equation, x)
    
    x_0 = initial_guess
    for _ in range(max_iterations):
        f_x_0 = equation.subs(x, x_0)
        f_prime_x_0 = derivative.subs(x, x_0)
        
        x_1 = x_0 - f_x_0 / f_prime_x_0
        
        if abs(x_1 - x_0) < tolerance:
            return x_1
        
        x_0 = x_1
    
    return None

def solve_higher_degree_simpsons_rule(coefficients, lower_limit, upper_limit, num_intervals):
    x = sp.symbols('x')
    equation = sum(sp.S(coeff) * x**i for i, coeff in enumerate(coefficients[::-1]))
    result = sp.integrate(equation, (x, lower_limit, upper_limit))
    interval_width = (upper_limit - lower_limit) / num_intervals
    
    integral_estimate = result
    for i in range(num_intervals + 1):  # Fix the loop range
        x_i = lower_limit + i * interval_width
        if i == 0 or i == num_intervals:
            integral_estimate += equation.subs(x, x_i)
        elif i % 2 == 1:
            integral_estimate += 4 * equation.subs(x, x_i)
        else:
            integral_estimate += 2 * equation.subs(x, x_i)
    
    integral_estimate *= interval_width / 3
    return integral_estimate

def plot_curve(coefficients):
    x = sp.symbols('x')
    equation = sum(sp.S(coeff) * x**i for i, coeff in enumerate(coefficients[::-1]))

    # Find the real roots of the equation
    real_roots = sp.solve(equation, x, domain=sp.S.Reals)

    if not real_roots:
        print("No real roots found. Unable to plot the curve.")
        lower_limit = -10  
        upper_limit = 10   
    else:
        # Convert SymPy expressions to numerical values
        lower_limit = float(min(real_roots))
        upper_limit = float(max(real_roots))

    x_values = np.linspace(lower_limit - 1, upper_limit + 1, 400)
    y_values = [float(equation.subs(x, val)) for val in x_values]  # Convert the result to a numerical float

    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, label=f'Equation: {equation}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Curve of the Equation')
    plt.grid(True)
    plt.legend()
    plt.show()


degree = int(input("What is the degree of the equation? "))

equation = []
for i in range(degree + 1):
    b = int(input(f"What is the coefficient of x^{i}? "))
    equation.append(b)

if degree == 1:
    solution = solve_linear(equation)
    print(f"The solution of the given equation is {solution}")
    plot_curve(equation)
elif degree == 2:
    solution = solve_quadratic(equation)
    print(f"The solution of the given equation is {solution}")
    plot_curve(equation)
else:
    method = input("Choose a method to solve the equation (Newton-Raphson or Simpson's rule): ").strip().lower()

    if method == "newton-raphson":
        initial_guess = float(input("Enter an initial guess for the root: "))
        solution = solve_higher_degree_newton_raphson(equation, initial_guess)
        if solution is not None:
            print(f"The solution of the given equation is approximately x = {solution}")
            plot_curve(equation)
        else:
            print("No solution found within the specified tolerance.")
    elif method == "simpson's rule":
        lower_limit = float(input("Enter the lower limit of integration: "))
        upper_limit = float(input("Enter the upper limit of integration: "))
        num_intervals = int(input("Enter the number of intervals: "))
        
        integral_estimate = solve_higher_degree_simpsons_rule(equation, lower_limit, upper_limit, num_intervals)
        print(f"The approximate integral value is {integral_estimate}")
        plot_curve(equation, lower_limit, upper_limit)
    else:
        print("This program only solves 1st-degree (linear) and 2nd-degree (quadratic) equations or uses the Newton-Raphson method or Simpson's rule for higher-degree equations.")
