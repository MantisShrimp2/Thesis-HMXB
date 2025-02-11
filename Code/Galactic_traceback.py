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
        self.table['Trace Time'].unit = 'Million years'
        
        
        return self.table
        
    def trace_linear_path(self, source_id, time_step=1000):
        """
        Trace the path of the star in Galactic coordinates until b = 0, using the Euler method.

        Parameters:
        - time_step (float): Step size in years for tracing the path.
        - max_steps (int): Maximum number of steps for tracing.

        Returns:
        - path  3 lists, for longitiude, latitiude and height. 
        Each path is for one star
        """

 
        row = self.table[self.table['source_id'] == source_id]
        #Convert to float because row[x] is a numpy array of length 1
        #makes issues for plotting
        l = float(row['l'][0])
        b = float(row['b'][0])
        mu_l = float(row['pm_l_poleski'][0])
        mu_b = float(row['pm_b_poleski'][0])
        dist = float(row['distance'][0])

        mu_l_deg_per_year = mu_l/ 3.6e6
        mu_b_deg_per_year = mu_b / 3.6e6
        #initalize path
        long_path, lat_path, z_path = [l],[b],[dist*np.sin(np.radians(b))]
        
        
        #limit time stepsto 3 million years
        max_steps = int(3e6/1000)
        
        #add ticks to the lines every 1 million years
        ticks = []
        
    # Trace back in time using Euler's method
        current_l = l
        current_b = b
        current_time= 0
        for _ in range(max_steps):
            # Update coordinates using Euler method
            #fixed timestep
            current_time += time_step
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
            # Add tick position at each 1 Myr
            if current_time % 1e6 == 0:
                ticks.append((current_l, current_b))
        # Stop if b crosses zero -optional
            # if current_b * lat_path[0] < 0:
            #     break
    #add to table
    #self.table['Galactic Path'] = traced_paths
        
        return long_path, lat_path, z_path, ticks
    def trace_galactic_path(self, source_id, time_step=1000):
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
        radial_velocity = float(row['RV'][0])*u.km/u.s
        k = 4.74 * (u.km/u.s)/(u.mas *u.kpc/u.yr) #km/s per mas/yr 
        print(radial_velocity)
        #transform to galactic frame
        with coord.galactocentric_frame_defaults.set('v4.0'):
            galcen_frame = coord.Galactocentric()
        galactic_rep = coord.SkyCoord(l=l,b=b,pm_l_cosb=mu_l,pm_b=mu_b,distance=dist,
                                      radial_velocity =radial_velocity, frame='galactic')
        #transform frame
        star_galacto = galactic_rep.transform_to(galcen_frame)


        initial_pos = gd.PhaseSpacePosition(star_galacto.data)
        total_time = -3.0 *u.Myr
        dt = -0.2 *u.Myr
        n_steps = int(np.abs(total_time/dt))
        potential = gp.MilkyWayPotential2022()  
       # hamiltonian = gp.Hamiltonian(potential)
        #orbit = potential.integrate_orbit(initial_pos, dt=-time_step, t1=0, t2=-total_time)
        orbit = potential.integrate_orbit(initial_pos, dt=dt,n_steps=n_steps)
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
            try:
                data = self.table[self.table['source_id']==ids]
                z_naught = data['distance']*np.sin(np.radians(data['b']))*1000
                l_naught = data['l']
                b_naught = data['b']
                arrow_pml, arrow_pmb = data['pm_l_poleski'], data['pm_b_poleski']
                long_path, lat_path, z_path,ticks = self.trace_linear_path(ids, time_step=1000)
                N = len(long_path)
                #plot path and color plot by specrtral type
                sp_type = data['Mod_SpType'][0] #dumb
                path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'
                plt.scatter(long_path[1:N], z_path[1:N],s=3,color='xkcd:black',alpha=0.7)
                #plt.scatter(long_path[1:N], lat_path[1:N],s=2,color=path_color,alpha=0.7)
                #plot star current position
                
                #plt.scatter(l_naught,b_naught,color='xkcd:black',s=50)
                plt.scatter(l_naught,z_naught,color=path_color)
                #plot pm vectors
                arrow_pmz = data['distance']*np.sin(np.radians(arrow_pmb))*1000
                #two arrows for b or z
                #plt.quiver(l_naught,b_naught, arrow_pml,arrow_pmb,color='xkcd:grey',angles='xy',width=0.002)
                plt.quiver(l_naught,z_naught, arrow_pml,arrow_pmz,color='xkcd:shit brown',angles='xy',width=0.002)
               
                for tick_l, tick_b in ticks:
                    tick_z = data['distance']*np.sin(np.radians(tick_b))*1000
                    
                    plt.scatter(tick_l, tick_z, color='xkcd:canary', s=20)
            except Exception as e:
                print(f'{e}')
        
        plt.xlabel('Galactic Longitude (deg)')
        plt.gca().invert_xaxis() # resverse the x-axis, standrd to show longitude
        plt.ylabel('height (pc)')
        plt.title('Traceback Paths in Galactic Coordinates')
       # plt.legend()
        plt.grid(True)
        
        mydir = os.path.dirname(os.path.realpath(__file__))
        today = datetime.now().strftime("%Y%m%d")
        if savefig == True:
            plt.savefig(parentdir+'/Figures/'+f"Tracepath_{today}.png")
        plt.show()
        return None
    def plot_with_cluster(self,source_id, clustertable=None,savefig=False):
        
        single_star = self.table[self.table['source_id']==source_id]
        long_path, lat_path, _, ticks = self.trace_linear_path(source_id=source_id)
        #convert to numpy arrays
        long_path, lat_path, ticks = np.array(long_path), np.array(lat_path),np.array(ticks)
        
        #inital positions
        
        l_naught = single_star['l']
        b_naught = single_star['b']
        
        #proper motion vectors
        arrow_pml, arrow_pmb = single_star['pm_l_poleski'], single_star['pm_b_poleski']
        #member clusters
        if clustertable is not None:
            member_long = clustertable['l']
            member_lat = clustertable['b']
        
        #include gala
        orbit  = self.trace_galactic_path(source_id)
        #fig = orbit.plot()
        # Convert to Galactic coordinates
        galcen_frame = coord.Galactocentric()
        galactic_coords = orbit.to_coord_frame(galcen_frame).transform_to(coord.Galactic())
        # galactic_coords = coord.SkyCoord(
        #     x=x_vals, y=y_vals, z=z_vals, 
        #     frame='galactocentric', representation_type='cartesian'
        # ).transform_to('galactic')

        # Extract Galactic longitude (l) and latitude (b)
        gala_l_vals = galactic_coords.l.deg % 360 # Galactic longitude in degrees
        gala_b_vals = galactic_coords.b.deg
        
        # print(f"Initial position: l={l_naught}, b={b_naught}")
        # print(f"Orbit start: l={gala_l_vals[0]}, b={gala_b_vals[0]}")
        plt.figure(figsize=(10,5))
        if clustertable is not None:
            plt.scatter(member_long,member_lat, s=50,marker='*',label='Members',color='xkcd:grey')
        
        N = len(long_path)
        star_name = str(single_star['Name'].value[0])
        #cluster_name = clustertable['Name']
        #plot path and color plot by specrtral type
        sp_type = single_star['Mod_SpType'][0] #dumb
        path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'
        #linear path
        plt.scatter(long_path[1:N], lat_path[1:N],s=10,color='xkcd:black',alpha=0.7)
        plt.plot([],[],color='xkcd:black',label='Linear Path')
        #inistal position
        plt.scatter(l_naught,b_naught,color=path_color,label=star_name)
        #gala
        plt.plot(gala_l_vals,gala_b_vals, lw=1.1, marker='*', color='xkcd:vermillion',label='Gala path')
        #pm arrows
        plt.quiver(l_naught,b_naught, arrow_pml,arrow_pmb,color='xkcd:shit brown',angles='xy',width=0.002)
        
        #tick markers
        for tick_l, tick_b in ticks:
            plt.scatter(tick_l,tick_b,color='xkcd:canary',s=20,marker='x')
        plt.xlabel('Galactic Longitude(deg)')
        plt.ylabel('Galactic Latitude (deg)')
        plt.gca().invert_xaxis()
        plt.title(f"Retraced path of {star_name}")
        plt.grid(True)
        plt.legend(loc='best')
        plt.show()
        return None

test_table = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/HMXB_20250210_.ecsv',format='ecsv')
test_170037 = test_table[test_table['Name']=='4U 1700-377']
scoob1 = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/SCO OB1-result.ecsv')
mydir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(mydir)
if __name__ == "__main__":
    GalacticTraceback(test_table).plot_with_cluster(test_170037['source_id'],clustertable=scoob1,savefig=False)
#    GalacticTraceback(test_table).plot_trace(savefig=False)
    

