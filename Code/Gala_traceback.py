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
import gala.potential as gp
import gala.dynamics as gd
import astropy.units as u
import os
import astropy.coordinates as coord

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
        self.table['Trace Time'].unit = u.Myr
        
        
        return self.table
        
    def trace_galactic_path(self, source_id, time_step=100):
        """
        Steps:
            define the position in galactic coodinate system
            convert to cartesian coordinates
            integrate orbit
            
        """

 
        row = self.table[self.table['source_id'] == source_id]
        #Convert to float because row[x] is a numpy array of length 1
        #makes issues for plotting
        l = float(row['l'][0]) *u.deg
        b = float(row['b'][0]) * u.deg
        mu_l = float(row['pm_l_poleski'][0]) * u.mas/u.yr
        mu_b = float(row['pm_b_poleski'][0]) *u.mas/u.yr
        dist = float(row['distance'][0]) *u.kpc
        k = 4.74 * (u.km/u.s)/(u.mas *u.kpc/u.yr) #km/s per mas/yr 
        
        #transform to galactic frame
        with coord.galactocentric_frame_defaults.set('v4.0'):
            galcen_frame = coord.Galactocentric()
        galactic_rep = coord.SkyCoord(l=l,b=b,pm_l_cosb=mu_l,pm_b=mu_b,distance=dist,frame='galactic')
        #transform frame
        #star_galacto = galactic_rep.transform_to(galcen_frame)
        cartesian_rep  = galactic_rep.cartesian
        
        # x, y, z = cartesian_rep.x.value, cartesian_rep.y.value, cartesian_rep.z.value
        # vx, vy, = k*mu_l*dist,  k*mu_b*dist
        # if 'radial_velocity' in row.colnames and not np.isnan(row['radial_velocity'][0]):
        #     vz = row['radial_velocity'][0]  # Use value from table
        # else:
        #     vz =-60 *u.km/u.s  # Default to 0 km/s if not available
        # print(vz)
        initial_pos = gd.PhaseSpacePosition(galactic_rep.data)
        #gd.PhaseSpacePosition(pos=(x,y,z)*u.kpc,vel=(vx,vy,vz)*(u.km/u.s))
        total_time = 3 *1e6 *u.yr
        print(total_time)
        potential = gp.MilkyWayPotential2022()  # Or another potential model  
        hamiltonian = gp.Hamiltonian(potential)  
        orbit = potential.integrate_orbit(initial_pos, dt=-time_step, t1=0, t2=-total_time)
        
        #tick_times = np.arange(-1,-total_time.value -1, -1) * u.Myr
        #ticks = orbit.evaluate_at(tick_times)
        
        
        #tick_positions = [(t.represent_as('spherical').lon.value, t.z.value) for t in ticks]            
        #limit time stepsto 3 million years

        return orbit #tick_positions

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
            #try:
                data = self.table[self.table['source_id']==ids]
                star_name = str(data['Name'].value[0])
                z_naught = data['distance']*np.sin(np.radians(data['b']))*1000
                l_naught = data['l']
                b_naught = data['b']
                arrow_pml, arrow_pmb = data['pm_l_poleski'], data['pm_b_poleski']
                
                orbit  = self.trace_galactic_path(ids)
                print(orbit.t.shape)
                x_vals, y_vals, z_vals = orbit.pos.xyz  

                # Convert to Galactic coordinates
                galactic_coords = coord.SkyCoord(
                    x=x_vals, y=y_vals, z=z_vals, 
                    frame='galactocentric', representation_type='cartesian'
                ).transform_to('galactic')

                # Extract Galactic longitude (l) and latitude (b)
                l_vals = galactic_coords.l.deg  # Galactic longitude in degrees
                b_vals = galactic_coords.b.deg
                print(f"Initial position: l={l_naught}, b={b_naught}")
                print(f"Orbit start: l={l_vals[0]}, b={b_vals[0]}")
# Galactic latitude in degrees
                #plot path and color plot by specrtral type
                #print(len(l_vals),len(b_vals))
                sp_type = data['Mod_SpType'][0] #dumb
                path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'
                
                plt.plot(l_vals,b_vals, marker='*', color='xkcd:vermillion',label='Gala path')
                #plt.scatter(long_path[1:N], lat_path[1:N],s=2,color=path_color,alpha=0.7)
                #plot star current position
                
                #plt.scatter(l_naught,b_naught,color='xkcd:black',s=50)
                plt.scatter(l_naught,b_naught,color=path_color,label=star_name)
                #plot pm vectors
                arrow_pmz = data['distance']*np.sin(np.radians(arrow_pmb))*1000
                #two arrows for b or z
                #plt.quiver(l_naught,b_naught, arrow_pml,arrow_pmb,color='xkcd:grey',angles='xy',width=0.002)
                plt.quiver(l_naught,b_naught, arrow_pml,arrow_pmb,color='xkcd:shit brown',angles='xy',width=0.002)
               
                # for tick_l, tick_b in ticks:
                #     tick_z = data['distance']*np.sin(np.radians(tick_b))*1000
                    
                #     plt.scatter(tick_l, tick_z, color='xkcd:electric purple', s=20)
            # except Exception as e:
            #      print(type(e))
        
        plt.xlabel('Galactic Longitude (deg)')
        plt.gca().invert_xaxis() # resverse the x-axis, standrd to show longitude
        plt.ylabel('Galactic Latitiude (deg)')
        plt.title('Traceback Paths in Galactic Coordinates')
        plt.legend()
        plt.grid(True)
        
        mydir = os.path.dirname(os.path.realpath(__file__))
        today = datetime.now().strftime("%Y%m%d")
        if savefig == True:
            plt.savefig(mydir+'/Figures/'+f"Gala_tracepath_{today}.png")
        plt.show()
        return None

test_table = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/1700_37.csv',format='csv')
mydir = os.path.dirname(os.path.realpath(__file__))
source_id = test_table['source_id']
if __name__ == "__main__":
    GalacticTraceback(test_table).plot_trace(savefig=False)
#     x = GalacticTraceback(test_table).trace_galactic_path(source_id=source_id)
    
# # y = x.spherical.long[0]
# # z = x.z[0]  
# # Extract Cartesian coordinates from the orbit
# x_vals, y_vals, z_vals = x.pos.xyz  

# # Convert to Galactic coordinates
# galactic_coords = SkyCoord(
#     x=x_vals, y=y_vals, z=z_vals, 
#     frame='galactocentric', representation_type='cartesian'
# ).transform_to('galactic')

# # Extract Galactic longitude (l) and latitude (b)
# l_vals = galactic_coords.l.deg  # Galactic longitude in degrees
# b_vals = galactic_coords.b.deg  # Galactic latitude in degrees

# Plot the orbit in Galactic coordinates
#plt.plot(l_vals, b_vals, '-', color='k')
