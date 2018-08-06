# -*- coding: utf-8 -*-
"""
Created on Sun May  7 12:18:52 2017

@author: Theodor
"""
import numpy as np
__all__=["bregProj"]

def J(x,p,norm=np.linalg.norm):    
    """
    Duality function 
    
    
    Parameters
    ----------
    x : array_like - The point which dual is searched
    p : value - The .....
    """
    if norm(x,ord=2)>0:
        return x*norm(x,ord=2)**(p-2)
    else:
        return x

def quasi_newton(args):
    """
    Quasi Newtonian method for finding roots.
    starts at x=1 and assumes as gradient the secant going througn 1,f(1) 
    and 0,f(0). Afterwards it assumes as gradient the secant going though the
    points defined by the last and current iteration.
    
    Parameters
    ----------
    
    
    args : dict    - The dictionary containing all of the below
    f : function   - The function for which we want to find the root
    epsilon : float - The wished precision
    maxiter : int   - The maximum allowed iterationcount
    
    Returns
    -------
    float     - The root
    
    """
    std_args={"epsilon":10**(-9),"maxiter":10**3}
    std_args.update(args)
    f=args["f"]
    epsilon=args["epsilon"]
    maxiter=args["maxiter"]
    last_x=0
    f_x=f(last_x)
    f_prime_guess=f(1)-f_x
    last_x=1
    next_x=last_x-f_x/f_prime_guess
    f_x=f(last_x)
    iteration=0
    while iteration<maxiter and abs(next_x-last_x)>epsilon:     
        f_nx=f(next_x)
        f_prime_guess=(f_x-f_nx)/(last_x-next_x)
        tmp=next_x
        next_x=last_x-f_x/f_prime_guess
        last_x=tmp
        f_x=f_nx
        iteration+=1
        print(abs(last_x-next_x))
    return next_x


def bregProj(args):
    """
    Bregman projection operator onto a hyperplane given defined by 
    all x such that u*x=a.
    
    The bregman projection is to the standard bregman function ||.||**p
    
    
    Parameters
    ----------
    args : dict    - The dictionary containing all of the below
    u : array_like - The orthogonal vector defining the hyperplane
    a : value -      The offset needed to define the hyperplane
    p : value -      The parameter needed to define the bregman function  
    
    x : array_like - The point which we want to project onto the hyperplane
    
    method : function- The rootfinding algorithm
    norm : function - The underlying norm
    Returns
    -------
    array_like     - The projection
    """ 
    #parse the arguments
    std_args={"x":np.zeros(1),
              "u":np.zeros(1),
              "a":0,
              "p":2,
              "method":quasi_newton,
              "norm":np.linalg.norm}
    std_args.update(args)
    x=std_args["x"]
    u=std_args["u"]
    a=std_args["a"]
    p=std_args["p"]
    method=std_args["method"]
    norm=std_args["norm"]
    
    
    hp= lambda t: -np.dot(u,J(J(x,p)-t*u,p/(p-1)))+a
    if not np.any(u):
        t0=0
    elif a==0:
        t0=np.dot(J(x,p),u)/norm(u,ord=2)**2    
    else:
        arguments=dict(f=hp)
        t0=method(arguments)
    return J(J(x,p)-t0*u,p/(p-1))


