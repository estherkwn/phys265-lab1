#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 15:26:27 2022
Modified on Fri Mar 11 15:26:27 2022
@author: Your name

Description
------------
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from lab1a_utilities import calculate_force
from lab1a_utilities import calculate_potential

# Create the source charges
#q1
Q1, Q1x, Q1y = -2.5, -75., 50.
Q2, Q2x, Q2y = 2.5, -55., 75.
Q3, Q3x, Q3y = 2.8, 120., -75.
Q4, Q4x, Q4y = -2, 20., -50.
Q5, Q5x, Q5y = -2.5, 90., 100.
Q = np.array([Q1, Q1x, Q1y,\
              Q2, Q2x, Q2y,\
              Q3, Q3x, Q3y,\
              Q4, Q4x, Q4y,\
              Q5, Q5x, Q5y]).reshape(5,3)

# Set the default initial conditions for v0, angle, and y0
v0 = 2 #m/s
angle = 60 #degrees
y0 = 50 #ypos


# Keep x0 fixed at -100
x0 = -100.

def clear():
    
    # NO NEED TO EDIT THIS FUNCTION

    fig = plt.figure('Game Window')
    ax = fig.axes[0]
    ax.cla()
    ax.axis('square')
    ax.set_xlim(-200,200)
    ax.set_ylim(-200,200)
    ax.set_title('Electrostatic Projectile Game',fontsize=16)
    ax.set_xlabel('x position (meters)',fontsize=16)
    ax.set_ylabel('y position (meters)',fontsize=16)
    ax.grid(visible=True)
    fig.tight_layout()    
    fig.show()
    return

def create_game_window():
    fig = plt.figure('Game Window')
    fig.clf()
    ax = fig.add_subplot()
    clear()
    return


############################################################

def play():

    # NO NEED TO EDIT THIS FUNCTION

    global v0, angle, y0
    
    print("Starting x location is -100")
    v0 = float(input("Enter the initial speed between zero and 100.\n"))
    assert(v0 >= 0 and v0 <= 100), "Initial velocity should > 0 and < 100"
    angle = float(input("Enter the initial angle in degrees.\n"))
    assert(angle >= -180.00 and angle <= 180.00), \
        "Angle should be between -180 and +180"
    y0 = float(input("Enter the initial y position.\n"))
    assert(y0 >= -200 and y0 <= 200), "y0 should be between -200 and +200"
    
    plot_trajectory()
    
    return
 


############################################################

def reveal_potential():

    fig = plt.figure('Game Window')
    # for a 3D wireframe or surface plot, comment out this line:
    ax = fig.axes[0]
# Uncomment these lines to create a 3D wireframe or surface plot
#    fig.clf()
#    ax = fig.add_subplot(projection='3d')

    # your code goes here
    #Contour Plot
    nx = 101
    ny = 101
    x = np.linspace(-200.,200.,nx)
    y = np.linspace(-200.,200.,ny)
    zq = np.zeros(nx*ny).reshape(ny,nx)
    for row in range(ny):
        for col in range(nx):
            zq[row][col] = calculate_potential(x[col],y[row],Q)
    X, Y = np.meshgrid(x,y)
    import matplotlib.cm as cm
    levels = np.linspace(-1800., 25000., 21) 
    #should i use isclose to prevent crowding?
    ax.contour(X,Y,zq,levels, cmap=cm.cool)
    #print('zq: min,max',zq.min(),zq.max())
    
    #Surface Plot
#    import matplotlib.cm as cm
#    ax.plot_surface(X,Y,zq, cmap=cm.cool)

    return

############################################################

def reveal_forcefield():
    fig = plt.figure('Game Window')
    ax = fig.axes[0]
    
    # your code goes here
    

    return

############################################################


############################################################

def plot_trajectory():

    fig = plt.figure('Game Window')
    ax = fig.axes[0]
    
    # Your code goes here
    vx0 = v0* np.cos(angle*(np.pi/180))
    vy0 = v0* np.sin(angle*(np.pi/180))
    #print(vx0,vy0)
    t0, t1 = 0.0, 10.0
    tt = np.linspace(t0,t1,81)
    def derivatives(t,s):
        #s[0] = x
        #s[1] = vx
        #s[2] = y
        #s[3] = vy
        D = np.zeros(4)
        Fx, Fy = calculate_force(s[0],s[2],Q)
        D[0] = s[1]
        D[2] = s[3]
        D[1] = Fx
        D[3] = Fy
        return D
    soltraj = solve_ivp(derivatives, (t0,t1), (x0,vx0,y0,vy0), t_eval=tt)
    xtraj = soltraj.y[0]
    ytraj = soltraj.y[2]
    ax.plot(xtraj,ytraj, color='gray')
    ax.plot(xtraj[::2],ytraj[::2], marker='o', color='k', markersize=3, linestyle='')
    
    fig.show()
    
    return

############################################################

def solve_it():
    quad1 = int(input('Does quadrant 1 have:\n(1) positive charge; (2) negative charge; (3) neither; or (4) both?\n'))
    if quad1 == 2:
        print('Correct')
    else:
        print('Incorrect')
    quad2 = int(input('How about quadrant 2?\n'))
    if quad2 == 4:
        print('Correct')
    else:
        print('Incorrect')
    quad3 = int(input('How about quadrant 3?\n'))
    if quad3 == 3:
        print('Correct')
    else:
        print('Incorrect')
    quad4 = int(input('And finally, how about quadrant 4?\n'))
    if quad4 == 4:
        print('Correct')
    else:
        print('Incorrect')
    if quad1 == 2 and quad2 == 4 and quad3 ==3 and quad4 == 4:
        print('Congratulations! You have won!')
        reveal_potential()
    else:
        print('Sorry, try again.')
    return









