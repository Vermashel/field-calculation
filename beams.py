#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 15:27:30 2020
Предполагается, что это что-то по научьке
@author: verma
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import cmath


def integrate_2D(target_func, limits_x, limits_y, Nx=100, Ny=100):    
    x_min, x_max = limits_x
    y_min, y_max = limits_y
    delta_x = (x_max - x_min) / (Nx - 1)
    delta_y = (y_max - y_min) / (Ny - 1)
    
    probe_value = target_func(x_min, y_min)
    
    dtype = None
    if isinstance(probe_value, complex):
        dtype = np.complex                  # check the type of function
    else:
        dtype = np.float
    
    grid = np.ndarray(shape=(Nx, Ny), dtype=dtype)
    
    x = x_min
    for i in range(0, Nx):
        y = y_min
        for j in range(0, Ny):              # fill the grid with function values
            grid[i,j] = target_func(x,y)
            y += delta_y
        x += delta_x
        
    summ = 0
    for i in range(0, Nx-1):
       for j in range(0, Ny-1):
           summ += 0.25 * (grid[i,j] + grid[i+1,j] + grid[i,j+1] + grid[i+1,j+1]) * delta_x * delta_y
    
    #print(grid)
    #print(summ)
    return summ


def get_field(ampl_func, phase_func, limits_X, limits_Y, point, N_points, vawe_lambda=600e-9):
    x, y, z = point
    k = 2 * math.pi / vawe_lambda
    
    def podyntegralnaya_funkcia(X, Y):
        l = math.sqrt((X-x)**2 + (Y-y)**2 + z**2)
        result = 1/(4*math.pi) * ampl_func(X,Y) * cmath.exp(1j * phase_func(X,Y) + 1j * k * l) / l * \
                 (1j * k * z / l)
        
#        result = 1/(4*math.pi) * ampl_func(X,Y) * cmath.exp(1j * phase_func(X,Y) + 1j * k * l) / l * \
#                 ( - z / l**2 + 1j * k * z / l - 1j * k)
        return result
    
    E = integrate_2D(podyntegralnaya_funkcia, limits_X, limits_Y, Nx=N_points, Ny=N_points)
    
    return abs(E), cmath.phase(E)


def test_integrate_2D():
    r = integrate_2D( lambda x, y : x + 2*y, (-5, 5), (-1000, 1000) )
    if math.fabs(r) > 1e-5:
        print("test failed on f = x + 2y")
    
    r = integrate_2D( lambda x, y : 1.0, (-5, 5), (-1000, 1000) )
    if math.fabs(r - 2000 * 10) > 1e-5:
        print("test failed on f = 1.0")
        
    r = integrate_2D( lambda x, y : math.sin(x)+math.cos(x), (-10*math.pi, 10*math.pi), (-5*math.pi, 5*math.pi) )
    if math.fabs(r) > 1e-5:
        print("test failed on f = sin x + cos y, result is %f" % r)
    
    r = integrate_2D( lambda x, y : math.exp(x), (0, 5), (10, 11), Nx=500 )
    r0 = (11-10) * (math.exp(5) - math.exp(0))
    if math.fabs((r - r0) / r0) > 1e-5:
        print("test failed on f = exp(x), result is %f, ground truth is %f" % (r, r0))
    
    r = integrate_2D( lambda x, y : complex(1.0, 0.0), (-5.0, 5.0), (-1000, 1000))
    if not isinstance(r, complex):
        print("test failed: result is not complex")
    if abs(r - 2000 * 10) > 1e-5:
        print("test failed on f = 1.0")


def circle(X,Y):
    R = 0.05
    if X**2 + Y**2 < R**2:
        return 1
    else:
        return 0

N = 0    # number of integration points
D = 0.05    # length of integration segment
            # integration area is a square

L = 1
vawe_lambda=1000e-9
N_points = D**2 / L / vawe_lambda

x = np.linspace(-0.5, 0.5, 300) 
Amp = np.ndarray(x.shape) 
ph = np.ndarray(x.shape) 

for i in range(x.shape[0]):
    print(i)
    Amp[i], ph[i] = get_field(circle, lambda x,y: 0, (-D/2,D/2), (-D/2,D/2), (x[i],0,L), N_points=4*int(N_points), vawe_lambda=vawe_lambda)

print(x)
print(Amp)
print(ph)

plt.plot(x, Amp)

#test_integrate_2D()








