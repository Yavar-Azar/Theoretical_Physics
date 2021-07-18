#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:04:13 2021

@author: yavar001
"""

import numpy as np
import matplotlib.pyplot as plt




class analyser1:
    
    def __init__(self, log):
        self.log=log
        
        
        
        
    
        
    def plotter(self):
        
        
        
        data=np.loadtxt(self.log)
        plt.plot(data[:,0], data[:,1])
        
        
        plt.savefig("test.png")
        
        
        return
        
        