#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 11:00:02 2021

@author: yavar001
"""
import numpy as np
import matplotlib.pyplot as plt

def enumarate11(sentence):

    '''
    Parameters
    ----------
    sentence : string
        input sent of user.

    Returns
    -------
    n : int
        number of words in the sentence.

    '''
    
    # tmp is a list of strings
    tmp=sentence.split()
    
    n=len(tmp)
    
    return n



def mypower(m, n):
    '''
    Parameters
    ----------
    m : int
        DESCRIPTION.
    n : int
        DESCRIPTION.

    Returns
    -------
    c : int 
        DESCRIPTION.

    '''
    
    
    
    
    c=m
    for i in range (1,n):
        c=c*m
        
    
    return c



def plotter1(fname):
    '''
    

    Parameters
    ----------
    fname : string
        the file name.

    Returns
    -------
    None.

    '''
    
    data = np.loadtxt(fname)
    plt.xlabel("Energy (eV)")
    plt.plot(data[:,0], np.sin(data[:,1]), lw=3, label='test')
    plt.legend()
    plt.savefig("test.png")
    
    return 










    
    
    
    
    
    
    
    