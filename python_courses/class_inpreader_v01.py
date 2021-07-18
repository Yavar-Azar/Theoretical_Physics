#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:04:13 2021

@author: yavar001


Here we can explain difference between property and method
"""

import numpy as np
import matplotlib.pyplot as plt

def dataextract(filename):
    """
    
        Parameters
        ----------
        filename : string
            filename.

        Returns
        -------
        dataarray : np.array
            data in array format.

     """
    dataarray=np.loadtxt(filename)
    return dataarray


class analyser1:
    
    def __init__(self, log):
        self.log=log
        
#  @staticmethod    

    
    @property
    def arrayform(self):
        temp=dataextract(self.log)
        return temp
     
    def plotter(self):
        
        data=self.arrayform
        plt.plot(data[:,0], data[:,1])
        
        
        plt.savefig("test1.png")
        
        
        return
        
