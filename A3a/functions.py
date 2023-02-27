import numpy as np
import math

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

    # check that t is of type np.array
    t = np.array(t) if type(t) is list else t

    # math and return
    return np.sin((2 * np.pi * f * t) + (2 * np.pi * d * f))

def gabor(t, a:float=1.0, sigma:float=1, f:float=1.0, phi:float=0.0):
    """
    method to define the gabor function
    
    params:
    t: float or np.array or list - the input x values
    a: float - the amplitude of the wave
    sigma: float - the gaussian width (standard distrubution)
    f: float - frequency
    phi: float - the phase to give gabore or gaboro functions
    """

    # math and return
    return a * np.exp((-t**2)/(2.0 * (sigma**2))) * np.cos(2*np.pi*f*t+phi)

def gabore(t, a:float=1.0, sigma:float=1, f:float=1.0):
    """
    a helper function to perform the even gabor function

    params:
    t: float or np.array or list - the input x values
    a: float - the amplitude of the wave
    sigma: float - the gaussian width (standard distrubution)
    f: float - frequency
    phi: float - the phase to give gabore or gaboro functions

    return:
    the even function of a gabor
    """

    # call gabor with correct phase and return
    return gabor(t, a, sigma, f, phi=0.0)

def gaboro(t, a:float=1.0, sigma:float=1, f:float=1.0):
    """
    a helper function to perform the odd gabor function

    params:
    t: float or np.array or list - the input x values
    a: float - the amplitude of the wave
    sigma: float - the gaussian width (standard distrubution)
    f: float - frequency
    phi: float - the phase to give gaboro functions

    return:
    the even function of a gabor
    """

    # call gabor with correct phase and return
    return gabor(t, a, sigma, f, phi=(np.pi/2))

def gabor_norm(fs, sigma=1, f=1.0, phi=0.0):
    """
    function to define the normalized gabor
    
    params:
    fs: float - the space between x-values
    sigma: float - the gaussian width
    f: float - the frequency
    phi: the gabor phase (0 or np.pi/2)
    
    return:
    np.array - the normalized gabor values
    """

    # define time values
    T = np.arange(-math.sqrt(-math.log(0.01)*2*(sigma)**2), math.sqrt(-math.log(0.01)*2*(sigma)**2), fs)

    # gabor and return
    G = gabor(T, a=1.0, sigma=sigma, f=f, phi=phi)
    return np.linalg.norm(G, ord=2)

def gabore_norm(fs, sigma=1, f=1.0):
    """
    function to define the normalized gabor
    
    params:
    fs: float - the space between x-values
    sigma: float - the gaussian width
    f: float - the frequency
    
    return:
    np.array - the normalized gabor values
    """

    # define time values
    T = np.arange(-math.sqrt(-math.log(0.01)*2*(sigma)**2), math.sqrt(-math.log(0.01)*2*(sigma)**2), fs)

    # gabor and return
    G = gabore(T, a=1, sigma=sigma, f=f)
    return np.linalg.norm(G, ord=2)

def gaboro_norm(fs, sigma=1, f=1.0):
    """
    function to define the normalized gabor
    
    params:
    fs: float - the space between x-values
    sigma: float - the gaussian width
    f: float - the frequency
    phi: the gabor phase (0 or np.pi/2)
    
    return:
    np.array - the normalized gabor values
    """

    # define time values
    T = np.arange(-math.sqrt(-math.log(0.01)*2*(sigma)**2), math.sqrt(-math.log(0.01)*2*(sigma)**2), fs)

    # gabor and return
    G = gaboro(T, a=1.0, sigma=sigma, f=f)
    return np.linalg.norm(G, ord=2)

def gammatone(t, a:float=1, f:float=200, n:float=4.0, phi:float=0.0, norm:bool=True):
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
    G = gaboro_norm(0.01, 100, 10000) * (t**(n-1))*np.exp(-2*np.pi*b*t)*np.cos(2*np.pi*f*t + phi)

    # normalize if necessary and return
    G = ( G / max(G) ) if norm is True else G
    return G

def u(t=0):
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
