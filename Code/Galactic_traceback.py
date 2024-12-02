#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:04:32 2024

@author: karan
"""

import numpy as np
import pandas as pd
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import matplotlib.pyplot as plt
from astropy.table import Column, join, Table, vstack, hstack
from astropy.io import ascii
from datetime import datetime

import os

class GalacticTraceback:
    def __init__(self,table):
        
        # self.l = long
        # self.b = lat
        # self.mu_l = mu_l_cosb
        # self.mu_b = mu_b
        self.table = table
        self.color_map = {
        "xkcd:blue": "O I-III",
        "xkcd:red": "B I-III",
        "xkcd:bright Blue": "Oe I-III",
        "xkcd:dark blue": "Oe IV-V",
        "xkcd:green": "B0e I-III",
        "xkcd:bright green": "B1e I-III",
        "xkcd:grass green": "B2e I-III",
        "xkcd:black": "OB IV-V",
        "xkcd:grey": "M,A,None",
        "xkcd:purple": "B0e IV-V",
        "xkcd:light purple": "B1e IV-V",
        "xkcd:dark purple": "B2e IV-V",
    }
    def traceback_time(self):
        '''
        Calculate the traceback time for a star
        The time is takes to return to the galactic midplane in years
        
        input:
        self
        b - galactic longitude degrees
        mu_b- proper motion in b mas/yr

        Returns
        -------
        None.

        '''

        lat = self.table['b']
        mu_b = self.table['pm_b_poleski']
        
        #convert mu_b to degree/yr 1deg = 3.6 million milliarcseconds
        mu_b_deg = mu_b /(3.6e6) 
        trace_time = lat/mu_b_deg
        
        self.table['Trace Time'] = np.array(trace_time)/float(1e6)
        self.table['Trace Time'].unit = 'Million years'
        
        
        return self.table
        
    def trace_linear_path(self, source_id, time_step=1000):
        """
        Trace the path of the star in Galactic coordinates until b = 0, using the Euler method.

        Parameters:
        - time_step (float): Step size in years for tracing the path.
        - max_steps (int): Maximum number of steps for tracing.

        Returns:
        - path (list of tuples): List of (l, b) pairs representing the path.
        """

 
        row = self.table[self.table['source_id'] == source_id]
        #Convert to float because row[x] is a numpy array of length 1
        #makes issues for plotting
        l = float(row['l'][0])
        b = float(row['b'][0])
        mu_l = float(row['pm_l_poleski'][0])
        mu_b = float(row['pm_b_poleski'][0])
        dist = float(row['distance'][0])
        # l = 121.22
        # b = -1.46
        # mu_l = -1.8178140519374264
        # mu_b = -0.4422414945716403
        # dist = 3.677948996847052
        mu_l_deg_per_year = mu_l/ 3.6e6
        mu_b_deg_per_year = mu_b / 3.6e6
        #initalize path
        long_path, lat_path, z_path = [l],[b],[dist*np.sin(np.radians(b))]
        
        max_steps = int(1e5)
    # Trace back in time using Euler's method
        current_l = l
        current_b = b
        for _ in range(max_steps-1):
            # Update coordinates using Euler method
            #fixed timestep
            current_l -= mu_l_deg_per_year * time_step  # Adjust longitude
            current_b -= mu_b_deg_per_year * time_step #adjust latitiude
            # Adjust latituded

            #Wrap longitude to [0, 360)
            current_l %= 360
            #calculate the height
            current_z = dist*np.sin(np.radians(current_b))*1000
            
            #append to path
            long_path.append(current_l)
            lat_path.append(current_b)
            z_path.append(current_z)
            
        # Stop if b crosses zero 
            if current_b * lat_path[0] < 0:
                break
    #add to table
    #self.table['Galactic Path'] = traced_paths
        
        return long_path, lat_path, z_path

    def plot_trace(self,savefig=False):
        #parllax mask
        #only plot data with confident parallax
        prlx_mask = self.table['parallax']/self.table['parallax_error'] >=5.0
        table = self.table[prlx_mask]
        plt.figure(figsize=(10,5))
        #plot zero line
        plt.axhline(y=0, color= 'xkcd:black',linestyle='--')
        source_ids = table['source_id']
        for ids in source_ids:
            try:
                data = self.table[self.table['source_id']==ids]
                z_naught = data['distance']*np.sin(np.radians(data['b']))*1000
                l_naught = data['l']
                b_naught = data['b']
                arrow_pml, arrow_pmb = data['pm_l_poleski'], data['pm_b_poleski']
                long_path, lat_path, z_path = self.trace_linear_path(ids, time_step=1000)
                N = len(long_path)
                #plot path and color plot by specrtral type
                sp_type = data['Mod_SpType'][0] #dumb
                path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'
                plt.scatter(long_path[1:N], lat_path[1:N],s=3,color=path_color)
                
                #plot star current position
                plt.scatter(l_naught,b_naught,color='xkcd:black')
                #plot pm vectors
                plt.quiver(l_naught,b_naught, arrow_pml,arrow_pmb,color='xkcd:grey',angles='uv',width=0.002)
            except Exception as e:
                print(f'{e}')
        
        plt.xlabel('Galactic Longitude (deg)')
        plt.ylabel('Galactic Latitude (deg)')
        plt.title('Traceback Paths in Galactic Coordinates')
       # plt.legend()
        plt.grid(True)
        
        mydir = os.path.dirname(os.path.realpath(__file__))
        today = datetime.now().strftime("%Y%m%d")
        if savefig == True:
            plt.savefig(mydir+'/Figures/'+f"Tracepath_{today}.png")
        plt.show()
        return None

test_table = ascii.read('/home/karan/Documents/UvA/Thesis/HMXB_practice_analysis.csv',format='csv')
mydir = os.path.dirname(os.path.realpath(__file__))
if __name__ == "__main__":
    GalacticTraceback(test_table).plot_trace(savefig=True)

