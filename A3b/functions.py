import numpy as np
import math
from matplotlib import pyplot as plt

def sinewave(t, f:float=1.0, d:float=0.0):
    """
    function to create a sinewave

    params:
    t: float - the time value
    f: float - the frequency of the wave
    d: float - the delay in time

    return:
    float - the function value for the specified time
    """

    if(type(t) == list):
        t = np.array(t)
    phi = 2 * np.pi * d * f
    return np.sin(2 * np.pi * f * t + phi)

def cosinewave(t, f:float=1, d:float=0):
    t = np.array(t) if type(t) is list else t
    
    return np.cos((2 * np.pi * f * t) + (2 * np.pi * d * f))

def gabor(t, a=1.0, sigma=1, f=1.0, phi=0.0):
    return a * np.exp((-t**2)/(2.0 * (sigma**2))) * np.cos(2*np.pi*f*t+phi)


def gabore(t, a=1.0, sigma=1, f=1.0):
    return gabor(t, a, sigma, f=f, phi=0.0)


def gaboro(t, a=1.0, sigma=1, f=1.0):
    return gabor(t, a, sigma, f=f, phi=(np.pi/2))


def gabor_norm(fs, sigma=1, f=1.0, phi=0.0):
    vanish_point = math.sqrt(-math.log(0.01)*2*(sigma)**2)
    time_val = -vanish_point
    t = []
    while time_val < vanish_point:
        t.append(time_val)
        time_val += 1/fs

    gabor_vals = []
    for t_val in t:
        gabor_vals.append(gabor(t_val, a=1.0, sigma=sigma, f=f, phi=phi))
    return np.linalg.norm(gabor_vals, ord=2)


def gabore_norm(fs, sigma=1, f=1.0):
    vanish_point = math.sqrt(-math.log(0.01)*2*(sigma)**2)
    time_val = -vanish_point
    t = []
    while time_val < vanish_point:
        t.append(time_val)
        time_val += 1/fs

    gabor_vals = []
    for t_val in t:
        gabor_vals.append(gabore(t_val, a=1.0, sigma=sigma, f=f))
    return np.linalg.norm(gabor_vals, ord=2)


def gaboro_norm(fs, sigma=1, f=1.0):
    vanish_point = math.sqrt(-math.log(0.01)*2*(sigma)**2)
    time_val = -vanish_point
    t = []
    while time_val < vanish_point:
        t.append(time_val)
        time_val += 1/fs

    gabor_vals = []
    for t_val in t:
        gabor_vals.append(gaboro(t_val, a=1.0, sigma=sigma, f=f))
    return np.linalg.norm(gabor_vals, ord=2)


def plot_gabor(t, sigma=4.0, f=1.0, a=1.0):
    gabore_vals = []
    for t_val in t:
        gabore_vals.append(gabore(t_val, sigma=sigma, f=f, a=a))

    plt.figure(figsize=(8, 6), dpi=80)
    plt.plot(t, gabore_vals)
    plt.ylabel("Amplitude")
    plt.xlabel("Time (s)")
    plt.show()

# 1c

def gammatone(t, a:float=1.0, f:float=1.0, n:float=4.0, phi:float=0.0, norm:bool=True):
    """
    function to implement the gammatone
    
    params:
    t: list or np.array or float - the input time value
    f: float - frequency
    n: float - shape parameter
    phi: float - phase
    """

    # check that t is the right form
    t = np.array(t) if type(t) is list else t

    # calculate other params, math
    b = 1.019*(24.7*(((4.37*f)/1000) + 1))
    G = a * (t**(n-1))*np.exp(-2*np.pi*b*t)*np.cos(2*np.pi*f*t + phi)

    # normalize if necessary and return
    G = ( G / np.sum(np.absolute(G)) ) if norm is True else G
    return G

def u(t=0, f=1):
    # function of tlim
    c = 0
    T = np.array(t)

    # unit step
    return np.array([1 if t_ >= c else 0 for t_ in T]) if T.size > 1 else 1 if t >= c else 0

def delta(t=0, fs:float=1):
    # function of tlim
    c = 0
    T = np.array(t)

    # delta and return
    if T.size > 1:
        return np.array([1 if abs(t_-c) <= abs(c-1/(fs*2)) else 0 for t_ in T]) 
    else:
        return 1 if abs(t) < abs(c-1/(fs*2)) else 0

def power(x):
    X = np.array(x)
    return energy(x)/X.size

def energy(x):
    # raise all elements to 2nd power and sum all elements
    X = np.array(x)
    return np.sum(np.power(X, 2))
