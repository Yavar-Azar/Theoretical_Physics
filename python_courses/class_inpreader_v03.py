#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:04:13 2021

@author: yavar001


Here we can explain difference between property and method
"""

import numpy as np
import matplotlib.pyplot as plt


class analyser1:
    
    def __init__(self, log):
        self.log=log
         

    @staticmethod   
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

    
    @property
    def arrayform(self):
        temp=self.dataextract(self.log)
        return temp

        

# I want to linewidth to plotter argument
# IN HERFEYI NIST KE DIALOGUE INJAD KONIM
    def plotter(self, linewidth=1):
        
        
        colorid=int(input("""please enter number of color for your plot\n
                    red   1
                    blue  2
                    green 3
                    magenta  4
                    yellow 5
                    """
                    ))
        coldict={1:"red", 2:"blue", 3:"green", 4:"magenta", 5:"yellow"}
        
        mycolor=coldict[colorid]
        
        data=self.arrayform
        plt.plot(data[:,0], data[:,1], lw=linewidth, c=mycolor)
        
        
        plt.savefig("test1.png")
        
        
        return
        
