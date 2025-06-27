#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:04:32 2024

@author: Karan Kumar
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
from datetime import datetime
import gala.potential as gp
import gala.dynamics as gd
import gala.integrate as gi
from gala.coordinates import reflex_correct
import astropy.units as u
import os
import astropy.coordinates as coord
from astropy.coordinates import CartesianRepresentation, CartesianDifferential
class GalacticTraceback:
    def __init__(self,table,int_time):
        
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
        self.int_time = int_time
        self.k = 4.74047 #yr/kms convert mas/yr to km/yr
        self.sun_params = Table()
        self.sun_params['l'] = [0.0*u.deg]               # Galactic longitude in degrees
        self.sun_params['b'] = [0.0*u.deg]               # Galactic latitude in degrees
        self.sun_params['distance_bj'] = [8.15*u.kpc]        # Distance in kpc (typical Galactocentric distance of the Sun)
        self.sun_params['RV'] = [0.0*u.km/u.s]              # Radial velocity in km/s (assumed zero for the Sun in this frame)
        self.sun_params['pm_l_poleski'] = [0.0*u.mas/u.yr]    # Proper motion in l (mas/yr)
        self.sun_params['pm_b_poleski'] = [0.0*u.mas/u.yr]    # Proper motion in b (mas/yr)
    def traceback_time(self):
        '''
        Calculate the traceback time for a star
        The time is takes to return to the galactic midplane in years
        
        input:
        self.table - astopy table to put result back into
        b - galactic longitude degrees
        mu_b- proper motion in b mas/yr

        Returns
        -------
        traceback time, into table.

        '''

        lat = self.table['b']
        mu_b = self.table['pm_b_poleski']
        
        #convert mu_b to degree/yr 1deg = 3.6 million milliarcseconds
        mu_b_deg = mu_b /(3.6e6) 
        trace_time = lat/mu_b_deg
        
        self.table['Trace Time'] = np.array(trace_time)/float(1e6)
        self.table['Trace Time'].unit = 'Million years'
        
        
        return self.table
    def cluster_linear_path(self, cluster_params, int_time, time_step=1000):
        '''
        Same as trace_linear_path only for cluster, lazy way to do orbit integration
        

        Parameters
        ----------
        cluster_params : TYPE
            DESCRIPTION.
        int_time : TYPE
            DESCRIPTION.
        time_step : TYPE, optional
            DESCRIPTION. The default is 1000.

        Returns
        -------
        None.

        '''
        row = cluster_params
        l = float(row['l'][0])
        b = float(row['b'][0])
        mu_l = float(row['pm_l_poleski'][0])
        mu_b = float(row['pm_b_poleski'][0])
        dist = float(row['distance_bj'][0])
        
        mu_l_deg_per_year = (mu_l)/ 3.6e6
        mu_b_deg_per_year =(mu_b )/ 3.6e6
        #initalize path
        long_path, lat_path, z_path = [l],[b],[dist*np.sin(np.radians(b))]
        
        
        #limit time steps to 3 million years
        total_int_time = int_time* 1e6
        max_steps = abs(int(total_int_time/time_step))
        
        total_int_time = int_time* 1e6
        max_steps = abs(int(total_int_time/time_step))
        
        #add ticks to the lines every 1 million years
        ticks = []
        current_l = l
        current_b = b
        current_time= 0
        for _ in range(max_steps):
            # Update coordinates using Euler method
            #fixed timestep
            current_time += time_step
            current_l += mu_l_deg_per_year * time_step * np.sign(int_time)  # Adjust longitude
            current_b += mu_b_deg_per_year * time_step * np.sign(int_time) #adjust latitiude
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
        return long_path, lat_path, z_path, ticks
        
        
    def trace_linear_path(self, source_id, cluster_params, int_time, time_step=1000):
        """
        Trace the linear path of a source in galactic coordinates for {int_time} years.

        Parameters:
        - time_step (float): Step size in years for tracing the path.
        - max_steps (int): Maximum number of steps for tracing.
        - cluster_params (table): l,b,pm_l_cosb, pm_b, distance, rv of cluster
        -int_time (float): how many Myrs of integration ex if int_time = -3.0,
        integration will be for -3e6 years

        Returns:
        - path  3 lists, for longitiude, latitiude and height
        - ticks: dots to represent every 1 million years of integration. 
        Each path is for one star
        """

 
        row = self.table[self.table['source_id'] == source_id]
        #Convert to float because row[x] is a numpy array of length 1
        #makes issues for plotting
        l = float(row['l'][0])
        b = float(row['b'][0])
        mu_l = float(row['pm_l_poleski'][0])
        mu_b = float(row['pm_b_poleski'][0])
        dist = float(row['distance_bj'][0])
        
        #plot path wrt cluster
        if cluster_params is not None:
            mu_cl_l = cluster_params['pm_l_poleski'][0]
            mu_cl_b = cluster_params['pm_b_poleski'][0]
            dist_cl  = cluster_params['distance_bj'][0]
            
            mu_l_deg_per_year = (mu_l-mu_cl_l)/ 3.6e6
            mu_b_deg_per_year =(mu_b-mu_cl_b) / 3.6e6
        else:
            mu_l_deg_per_year = (mu_l)/ 3.6e6
            mu_b_deg_per_year =(mu_b )/ 3.6e6
        #initalize path
        long_path, lat_path, z_path = [l],[b],[dist*np.sin(np.radians(b))]
        
        
        #limit time steps to 3 million years
        total_int_time = int_time* 1e6
        max_steps = abs(int(total_int_time/time_step))
        
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
            current_l += mu_l_deg_per_year * time_step * np.sign(int_time)  # Adjust longitude
            current_b += mu_b_deg_per_year * time_step * np.sign(int_time) #adjust latitiude
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
        return long_path, lat_path, z_path, ticks

        
    def trace_galactic_path(self,source_info,int_time):
        """
        Steps:
            source_info- table of object's diatance, l, b, proper motion, radial velocity
            think a star or cluster'
            int_time float of integreation time
            define the position in galactic coodinate system
            convert to cartesian coordinates
            integrate orbit
            
        """
 
 
        row = source_info #self.table[self.table['source_id'] == source_id
        # l = row['l']
        # b = row['b'] 
        # mu_l = row['pm_l_poleski']
        # mu_b = row['pm_b_poleski']

        # dist = row['distance_bj'] #kpc
        # radial_velocity = row['RV'] 

   
        # mu_l_total  = mu_l
        # mu_b_total = mu_b
        # dist_total  = (dist.value) * u.kpc
        # radial_velocity_total = radial_velocity
        
        ra = row['ra']
        dec = row['dec']
        pmra = row['pmra']
        pmdec = row['pmdec']
        radial_velocity = row['RV'] 
        # if row['parallax'].value/row['parallax_error'].value >5.0:
            
        #     distance_icrs = row['distance_bj']
        # else:
        if row['parallax'][0] <= 0.0:
            distance_icrs = row['distance_bj']
        #distance_icrs = row['distance_bj']
        else:
            distance_icrs = (1/row['parallax']).value *u.kpc
        
        
            
        with coord.galactocentric_frame_defaults.set('v4.0'):
            galcen_frame = coord.Galactocentric()
        # galactic_rep = coord.SkyCoord(l=l,
        #                               b=b,
        #                               pm_l_cosb=mu_l_total,
        #                               pm_b=mu_b_total,
        #                               distance=dist_total,
        #                               radial_velocity=radial_velocity_total,
        #                               frame="galactic")
        #galactic_rep = coord.SkyCoord(ra=ra,dec=dec,pm_ra_cosdec=pmra, pm_dec=pmdec,radial_velocity=radial_velocity)
        
        galactic_rep = coord.SkyCoord(ra=ra, dec=dec, distance=distance_icrs,
              pm_ra_cosdec=pmra, pm_dec=pmdec,
              radial_velocity=radial_velocity)
        #transform frame
        star_galacto = galactic_rep.transform_to(galcen_frame)
        # star_galacto = reflex_correct(star_galacto)
        #correct for solar motion
       # print(star_galacto)

        initial_pos = gd.PhaseSpacePosition(star_galacto.data)
        total_time = int_time *u.Myr
        dt = np.sign(int_time)* 0.1 *u.Myr
        n_steps = int(abs(total_time/dt))

        integrator = gi.LeapfrogIntegrator
        potential = gp.MilkyWayPotential2022()  
        orbit = potential.integrate_orbit(initial_pos, dt=dt, t1=0*u.Myr, t2=total_time)
        #orbit = potential.integrate_orbit(initial_pos, 
                                        # dt=dt,n_steps=n_steps,
                                        #  Integrator=integrator)
                

        return orbit 

    def plot_trace(self,savefig=False, cluster_params=None):
     
        table = self.table
        plt.figure(figsize=(10,5))
        #plot zero line
        plt.axhline(y=0, color= 'xkcd:black',linestyle='--')
        source_ids = table['source_id']
        for ids in source_ids:
            try:
                data = self.table[self.table['source_id']==ids]
                z_naught = data['distance_bj']*np.sin(np.radians(data['b']))*1000
                l_naught = data['l']
                b_naught = data['b']
                arrow_pml, arrow_pmb, arrow_pmz, _ = self.trace_linear_path(ids, None, 0.5,time_step=1000)
                long_path, lat_path, z_path,ticks = self.trace_linear_path(ids, None,int_time = self.int_time, time_step=1000)
                N = len(long_path)
                #plot path and color plot by specrtral type
                sp_type = data['Mod_SpType'][0] #dumb
                path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'
                #plt.scatter(long_path[1:N], z_path[1:N],s=3,color='xkcd:black',alpha=0.7)
                plt.scatter(long_path[1:N], lat_path[1:N],s=2,color='xkcd:black',alpha=0.5)
                #plot star current position
                
                plt.scatter(l_naught,b_naught,color=path_color,s=50,label='')
                #plt.scatter(l_naught,z_naught,color=path_color)
                #plot pm vectors
               # arrow_pmz = data['distance_para']*np.sin(np.radians(arrow_pmb))*1000
                #two arrows for b or z
                delta_l = arrow_pml[-1] - l_naught
                delta_b = arrow_pmb[-1] - b_naught
                plt.quiver(l_naught,b_naught, delta_l, delta_b,color='xkcd:brown',angles='xy',width=0.002)
                #causing issues if i dont convert to np array, something about masking in astropy
                # l_naught, z_naught = np.array(l_naught), np.array(z_naught)
                # plt.quiver(l_naught,z_naught, arrow_pml[:-1],arrow_pmz,color='xkcd:brown',angles='xy',width=0.002)
               
                for tick_l, tick_b in ticks:
                    #tick_z = data['distance_bj']*np.sin(np.radians(tick_b))*1000
                    
                    plt.scatter(tick_l, tick_b, color='xkcd:orange', s=20,alpha=0.5)
            except Exception as e:
                print(f'{e}')
        if cluster_params is not None:
            cl_long, cl_lat, cl_z, cl_ticks = self.cluster_linear_path(cluster_params, int_time=-3.0, time_step=1000)
            cl_l_arrow, cl_b_arrow,_,_ = self.cluster_linear_path(cluster_params, int_time=0.5,time_step=1000)
            cl_l, cl_b = cluster_params['l'], cluster_params['b']
            plt.scatter(cl_long, cl_lat,s=2,color='xkcd:grey',alpha=0.7)
            plt.scatter(cl_l, cl_b, color='xkcd:grey', label='Cluster Position')
            delta_cl_l = cl_l_arrow[-1] - cl_l
            delta_cl_b = cl_b_arrow[-1] - cl_b
            
            #proper motion arrows
            plt.quiver(cl_l,cl_b, delta_cl_l, delta_cl_b,color='xkcd:brown', angles='xy', width=0.002)
            for cl_tick_l, cl_tick_b in cl_ticks:
                plt.scatter(cl_tick_l, cl_tick_b, color='xkcd:orange', s=20)
            
        plt.xlabel('Galactic Longitude (deg)')
        #plt.gca().invert_xaxis() # resverse the x-axis, standrd to show longitude
        plt.ylabel('Galactic Latitude (deg)')
        plt.title('Traceback Paths in Galactic Coordinates')
        #plt.legend()
        plt.grid(True)
        
        today = datetime.now().strftime("%Y%m%d")
        if savefig == True:
            plt.savefig(parentdir+'/Figures/'+f"Tracepath_{today}.png")
        plt.show()
        return None
    
    def plot_trace_aitoff(self,savefig=False, cluster_params=None):
        table = self.table
        plt.figure(figsize=(10, 6))
        ax = plt.subplot(111, projection="aitoff")
        
        # Convert and center longitudes around 0 radians (Aitoff needs -π to π)
        def convert_coords(l, b):
            l_rad = np.radians(l)
            l_rad = ((l_rad + np.pi) % (2 * np.pi)) - np.pi  # Wrap to [-π, π]
            b_rad = np.radians(b)
            return l_rad, b_rad
        
        source_ids = table['source_id']
        for ids in source_ids:
            try:
                data = table[table['source_id'] == ids]
                l_naught, b_naught = data['l'][0], data['b'][0]
                sp_type = data['Mod_SpType'][0]
                path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'
        
                long_path, lat_path, _, ticks = self.trace_linear_path(ids, None, -10, time_step=1000)
                long_path_rad, lat_path_rad = convert_coords(long_path, lat_path)
                ax.plot(long_path_rad[1:], lat_path_rad[1:], lw=1.0, color='xkcd:black', alpha=1.0)
                
                # Star current position
                l0_rad, b0_rad = convert_coords(l_naught, b_naught)
                ax.scatter(l0_rad, b0_rad, color=path_color, s=50)
        
                # Proper motion arrow
                arrow_pml, arrow_pmb, _, _ = self.trace_linear_path(ids, None, 0.5, time_step=1000)
                delta_l = arrow_pml[-1] - l_naught
                delta_b = arrow_pmb[-1] - b_naught
                delta_l_rad, delta_b_rad = np.radians(delta_l), np.radians(delta_b)
                ax.quiver(l0_rad, b0_rad, delta_l_rad, delta_b_rad, color='xkcd:brown', angles='xy', width=0.002)
        
                # for tick_l, tick_b in ticks:
                #     tick_l_rad, tick_b_rad = convert_coords(tick_l, tick_b)
                #     ax.scatter(tick_l_rad, tick_b_rad, color='xkcd:orange', s=20, alpha=0.5)
            except Exception as e:
                print(f'{e}')
        
        if cluster_params is not None:
            cl_long, cl_lat, _, cl_ticks = self.cluster_linear_path(cluster_params, int_time=-3.0, time_step=1000)
            cl_l_arrow, cl_b_arrow, _, _ = self.cluster_linear_path(cluster_params, int_time=0.5, time_step=1000)
            cl_l, cl_b = cluster_params['l'], cluster_params['b']
            cl_long_rad, cl_lat_rad = convert_coords(cl_long, cl_lat)
            cl_l_rad, cl_b_rad = convert_coords(cl_l, cl_b)
            ax.scatter(cl_long_rad, cl_lat_rad, s=2, color='xkcd:grey', alpha=0.7)
            ax.scatter(cl_l_rad, cl_b_rad, color='xkcd:grey', label='Cluster Position')
        
            delta_cl_l = cl_l_arrow[-1] - cl_l
            delta_cl_b = cl_b_arrow[-1] - cl_b
            delta_cl_l_rad, delta_cl_b_rad = np.radians(delta_cl_l), np.radians(delta_cl_b)
            ax.quiver(cl_l_rad, cl_b_rad, delta_cl_l_rad, delta_cl_b_rad, color='xkcd:brown', angles='xy', width=0.002)
        
            # for cl_tick_l, cl_tick_b in cl_ticks:
            #     cl_tick_l_rad, cl_tick_b_rad = convert_coords(cl_tick_l, cl_tick_b)
            #     ax.scatter(cl_tick_l_rad, cl_tick_b_rad, color='xkcd:orange', s=20)
        
        ax.set_title("Traceback Paths in Galactic Coordinates (Aitoff)", pad=30)
        ax.grid(True)
        plt.tight_layout()
        
        if savefig:
            today = datetime.now().strftime("%Y%m%d")
            plt.savefig(f"traceback_aitoff_{today}.png", dpi=300)
        else:
            plt.show()
            
        return None
    def plot_separation(self,star,cluster,savefig=False):
        star_orbit = self.trace_galactic_path(star, self.int_time)
        cluster_orbit = self.trace_galactic_path(cluster, self.int_time)
        
        star_vel_x, star_vel_y, star_vel_z = self.velocity_check(star_orbit,False)
        cluster_vel_x, cluster_vel_y, cluster_vel_z = self.velocity_check(cluster_orbit,False)
        
        rel_velocity = np.sqrt((star_vel_x - cluster_vel_x)**2 + (star_vel_y - cluster_vel_y)**2 + (star_vel_z - cluster_vel_z)**2)
        
        std_rel_vel = np.std(rel_velocity)
        print(f'relative velocity {np.mean(rel_velocity)}, std {std_rel_vel}')
        
        separation = np.linalg.norm(star_orbit.xyz - cluster_orbit.xyz, axis=0)
        
        separation = np.array(separation)
        fig, ax = plt.subplots(1, figsize=(10,10))
        t_array = np.array(star_orbit.t)
        print(np.amin(separation)*1e3, 'pc')
        print(t_array[np.argmin(separation)]/1e6, 'Myr ago')    
        star_name = star['Name'][0]
        cluster_name = cluster['Name'][0]
        
        sep_where_min = np.amin(separation)*1e3
        time_where_min = t_array[np.argmin(separation)]
        #for labeling
        sep_where_min_round = round(np.amin(separation)*1e3,2)
        time_where_min_round = round(t_array[np.argmin(separation)],2)
        
        ax.set_title(f'Separation Between {cluster_name} with {star_name}')
        ax.set_ylabel('Separation pc')
        ax.set_xlabel('Time (Myr)')
        ax.scatter(time_where_min, sep_where_min, label=f'{time_where_min_round} Myr at {sep_where_min_round} pc',color='xkcd:purple', s=100)
        ax.plot(t_array[:-1], separation[:-1]*1e3, color='dodgerblue', ls='-', marker='', lw=2, zorder=-1)
        ax.legend(title='Kinematical Age and Minimum Separation')
        
        if savefig:
            mydir = os.path.dirname(os.path.realpath(__file__))
            parentdir = os.path.dirname(mydir)
            today = datetime.now().strftime("%Y%m%d")
            save_path = os.path.join(parentdir, 'Figures', 'Traceback/Separations', f"separation_{star_name}_with_{cluster_name}_{today}.png")
            plt.savefig(save_path)
    
        plt.show()
        return None
    def plot_separation_with_uncertainty(self, star, cluster, N=100, savefig=False):
        """
        Plot separation between star and cluster using N samples from star parameter uncertainties.
        Includes confidence intervals and global minimum marker.
    
        Parameters
        ----------
        star : astropy Table (length=1)
            The star data row with columns: ra, dec, distance_bj, pmra, pmdec, RV and their *_err.
        cluster : astropy Table (length=1)
            The cluster data row, same format as star.
        N : int
            Number of Monte Carlo samples for the star orbit.
        savefig : bool
            Whether to save the figure.
        """
    
        int_time = self.int_time
    
        # Get cluster orbit (only once)
        cluster_orbit = self.trace_galactic_path(cluster, int_time)
        cluster_xyz = cluster_orbit.xyz.to_value(u.pc)
        time_array = cluster_orbit.t.to_value(u.Myr)
    
        dist = np.array(star['distance_bj'])[0]
        dist_up = np.array(star['distance_bj_high'])[0] - dist
        dist_down = dist - np.array(star['distance_bj_low'][0])
        dist_err = (dist_up + dist_down)/2.0
        # Sample from star uncertainties
        d     = np.random.normal(dist,  dist_err,  N)
        parallax = np.random.normal(star['parallax'][0], star['parallax_error'][0], N)
        pmra  = np.random.normal(star['pmra'][0],         star['pmra_error'][0],         N)
        pmdec = np.random.normal(star['pmdec'][0],        star['pmdec_error'][0],        N)
        rv    = np.random.normal(star['RV'][0],           star['RV_err'][0],           N)
        ra    = star['ra'][0]
        dec   = star['dec'][0]
    
        sep_curves = []
    
        for i in range(N):
            star_sample = Table(
                [[d[i]] * u.kpc, [ra]*u.deg, [dec]*u.deg, [pmra[i]] * u.mas/u.yr, [pmdec[i]] * u.mas/u.yr, [rv[i]] * u.km/u.s, [parallax[i]]*u.mas],
                names=["distance_bj", "ra", "dec", "pmra", "pmdec", "RV",'parallax']
            )
    
            star_orbit = self.trace_galactic_path(star_sample, int_time)
            star_xyz = star_orbit.xyz.to_value(u.pc)
    
            sep = np.linalg.norm(star_xyz - cluster_xyz, axis=0)
            sep_curves.append(sep)
    
        sep_curves = np.array(sep_curves)
    
        # Percentiles
        p16 = np.percentile(sep_curves, 16, axis=0)
        p50 = np.percentile(sep_curves, 50, axis=0)
        p84 = np.percentile(sep_curves, 84, axis=0)
    
        # Global min
      
        min_idx = np.argmin(p50)
        min_time = time_array[min_idx]
        min_sep = p50[min_idx]

    
        # Plotting
        fig, ax = plt.subplots(1, figsize=(10, 6))
    
        for sep in sep_curves:
            ax.plot(time_array, sep, color='gray', alpha=0.1)
    
        ax.plot(time_array, p50, color='xkcd:navy', lw=2, label='Median separation')
        ax.fill_between(time_array, p16, p84, color='xkcd:emerald', alpha=0.3, label='68% CI')
    
        ax.scatter(min_time, min_sep, color='xkcd:red', s=80,
                   label=f'Min: {min_time:.2f} Myr, {min_sep:.1f} pc')
        ax.axvline(min_time, ls=':', color='xkcd:red')
        ax.axhline(min_sep, ls=':', color='xkcd:red')
    
        ax.set_title(f'Separation: {star["Name"][0]} vs {cluster["Name"][0]}')
        ax.set_xlabel('Time (Myr)')
        ax.set_ylabel('Separation (pc)')
        ax.legend()
        ax.grid(True)
    
        if savefig:
            mydir = os.path.dirname(os.path.realpath(__file__))
            parentdir = os.path.dirname(mydir)
            today = datetime.now().strftime("%Y%m%d")
            save_path = os.path.join(parentdir, 'Figures', 'Traceback/Separations',
                                     f"sep_{star['Name'][0]}_{cluster['Name'][0]}_{today}.png")
            plt.savefig(save_path, dpi=150)
    
        plt.tight_layout()
        plt.show()
        return None

        
    def plot_orbit_planes(self, star_x, star_y, star_z, cluster_x, cluster_y, cluster_z, star_label="Star", cluster_label="Cluster"):
        """
        Plot orbit projections of the star and cluster in (Y, X), (X, Z), and (Y, Z) planes, with motion direction arrows.
    
        Parameters:
        - star_x, star_y, star_z: Star's Cartesian coordinates
        - cluster_x, cluster_y, cluster_z: Cluster's Cartesian coordinates
        - star_label: Label for the star (default: "Star")
        - cluster_label: Label for the cluster (default: "Cluster")
        """
    
        fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    
        # Compute velocity vectors (differences between consecutive positions)
        star_dx = np.diff(star_x)
        star_dy = np.diff(star_y)
        star_dz = np.diff(star_z)
        
        cluster_dx = np.diff(cluster_x)
        cluster_dy = np.diff(cluster_y)
        cluster_dz = np.diff(cluster_z)
    
        # Ensure same length for plotting (remove last element from position arrays)
        star_x_q, star_y_q, star_z_q = star_x[:-1], star_y[:-1], star_z[:-1]
        cluster_x_q, cluster_y_q, cluster_z_q = cluster_x[:-1], cluster_y[:-1], cluster_z[:-1]
    
        # (Y, X) Plane
        axs[0].scatter(star_y, star_x, s=50, label=star_label, color='blue', marker="*")
        axs[0].scatter(star_y[0], star_x[0], s=100, color='green', edgecolor='black', label="Start")  # Initial position
        
        axs[0].plot(star_y, star_x, color='blue', alpha=0.7)
        axs[0].quiver(star_y_q, star_x_q, star_dy, star_dx, color="blue", angles="xy", scale_units="xy", scale=1)
    
        axs[0].scatter(cluster_y, cluster_x, s=10, label=cluster_label, color='red', marker="o")
        axs[0].scatter(cluster_y[0], cluster_x[0], s=50, color='purple', edgecolor='black', label="Cluster Start")  # Initial position
        axs[0].plot(cluster_y, cluster_x, color='red', alpha=0.5)
        axs[0].quiver(cluster_y_q, cluster_x_q, cluster_dy, cluster_dx, color="red", angles="xy", scale_units="xy", scale=1)
    
        axs[0].set_xlabel('Y')
        axs[0].set_ylabel('X')
        axs[0].set_title('(Y, X) Plane')
        axs[0].legend()
    
        # (X, Z) Plane
        axs[1].scatter(star_x, star_z, s=50, label=star_label, color='blue', marker="*")
        axs[1].scatter(star_x[0], star_z[0], s=100, color='green', edgecolor='black', label="Start")  # Initial position
        
        axs[1].plot(star_x, star_z, color='blue', alpha=0.7)
        axs[1].quiver(star_x_q, star_z_q, star_dx, star_dz, color="blue", angles="xy", scale_units="xy", scale=1)
    
        axs[1].scatter(cluster_x, cluster_z, s=10, label=cluster_label, color='red', marker="o")
        axs[1].scatter(cluster_x[0], cluster_z[0], s=50, color='purple', edgecolor='black', label="Cluster Start")  # Initial position
        axs[1].plot(cluster_x, cluster_z, color='red', alpha=0.5)
        axs[1].quiver(cluster_x_q, cluster_z_q, cluster_dx, cluster_dz, color="red", angles="xy", scale_units="xy", scale=1)
    
        axs[1].set_xlabel('X')
        axs[1].set_ylabel('Z')
        axs[1].set_title('(X, Z) Plane')
        axs[1].legend()
    
        # (Y, Z) Plane
        axs[2].scatter(star_y, star_z, s=50, label=star_label, color='blue', marker="*")
        axs[2].scatter(star_y[0], star_z[0], s=100, color='green', edgecolor='black', label="Start")  # Initial position
        
        axs[2].plot(star_y, star_z, color='blue', alpha=0.7)
        axs[2].quiver(star_y_q, star_z_q, star_dy, star_dz, color="blue", angles="xy", scale_units="xy", scale=1)
    
        axs[2].scatter(cluster_y, cluster_z, s=10, label=cluster_label, color='red', marker="o")
        axs[2].scatter(cluster_y[0], cluster_z[0], s=50, color='purple', edgecolor='black', label="Cluster Start")  # Initial position
        axs[2].plot(cluster_y, cluster_z, color='red', alpha=0.5)
        axs[2].quiver(cluster_y_q, cluster_z_q, cluster_dy, cluster_dz, color="red", angles="xy", scale_units="xy", scale=1)
    
        axs[2].set_xlabel('Y')
        axs[2].set_ylabel('Z')
        axs[2].set_title('(Y, Z) Plane')
        axs[2].legend()
    
        plt.tight_layout()
        plt.show()

    def plot_comoving_cluster(self,star_params,cluster_params,plotting=False):
        
        
        star_orbit = self.trace_galactic_path(star_params, int_time=-3.0)
        cluster_orbit = self.trace_galactic_path(cluster_params,int_time=-3.0)
        
        relative_orbit =  (star_orbit.xyz - cluster_orbit.xyz).to(u.kpc)

        time = star_orbit.t
    
        separation = np.linalg.norm(star_orbit.xyz - cluster_orbit.xyz, axis=0)

        min_sep = np.argmin(separation)
        time_min_sep = time[min_sep]
        
        
   
        rel_x,rel_y,rel_z = np.array(relative_orbit[0]), np.array(relative_orbit[1]), np.array(relative_orbit[2])
        
        #star_xyz = star_orbit.xyz
        star_x, star_y,star_z = np.array(star_orbit.x),np.array(star_orbit.y),np.array(star_orbit.z)

        #cluster_xyz = cluster_orbit.xyz
        cluster_x,cluster_y,cluster_z = np.array(cluster_orbit.x),np.array(cluster_orbit.y),np.array(cluster_orbit.z)
        

        #offset the star 
        star_x_shifted = star_x - cluster_x
        star_y_shifted = star_y - cluster_y
        star_z_shifted = star_z - cluster_z
        
        rel_dx = np.diff(star_x_shifted)
        rel_dy = np.diff(star_y_shifted)
        rel_dz = np.diff(star_z_shifted)
        rel_x_q, rel_y_q, rel_z_q = star_x_shifted[:-1], star_y_shifted[:-1], star_z_shifted[:-1]

        if plotting==True:
            # plt.figure(figsize=(10,5))
            # plt.scatter(star_y, star_z,label=f"{star_params['Name'][0]}",s=100,alpha=1.0)
            # plt.plot(star_y,star_z,alpha=1.0)
            
            # plt.scatter(cluster_y,cluster_z,s=50,label='Cluster')
            # plt.ylabel('Z')
            # plt.xlabel('Y')
            # plt.title('No Offset')
            # plt.gca().invert_yaxis()
            # plt.legend()
            # plt.show()
            
            plt.figure(figsize=(10,5))
          
            plt.scatter(rel_y,rel_z,label=f"{star_params['Name'][0]}",s=100,alpha=1.0)
            plt.quiver(rel_y_q, rel_z_q, rel_dy, rel_dz, color="xkcd:black", angles="xy", scale_units="xy", width=0.005)
            plt.scatter(rel_y[0], rel_z[0], color="red", label="Initial Position", s=100)
            crossing_indices = np.argwhere(star_z_shifted[1:] * star_z_shifted[:-1] <  0)
            #distance = np.sqrt(star_y_shifted**2 + star_z_shifted**2)
            #crossing_indices = np.argmin(distance)
            plt.plot(star_y_shifted,star_z_shifted,alpha=1.0)
            plt.scatter(star_y_shifted[crossing_indices], star_z_shifted[crossing_indices],color='xkcd:purple',s=50,label='Crossing')
            for idx in np.atleast_1d(crossing_indices):
                plt.axvline(star_y_shifted[idx], color='gray', linestyle='dashed', alpha=0.7)            
                #vertical line
                plt.axhline(star_z_shifted[idx], color='gray', linestyle='dashed', alpha=0.7) 
            #plt.scatter(cluster_y,cluster_z,s=50,label='Cluster')
            plt.ylabel('Z kpc')
            plt.xlabel('Y kpc')
            plt.title(f"{star_params['Name'][0]} WRT cluster")
            plt.gca().invert_yaxis()
            plt.legend()
            plt.show()
            #self.plot_orbit_planes(star_x_shifted,star_y_shifted,star_z_shifted,cluster_x,cluster_y,cluster_z)
            self.plot_orbit_planes(star_x,star_y,star_z,cluster_x,cluster_y,cluster_z)
            fig_3d = plt.figure()
            ax = fig_3d.add_subplot(projection='3d')
            #star
            ax.scatter(star_x, star_y, star_z, marker='*',label=f"{star_params['Name'][0]}")
            ax.scatter(cluster_x, cluster_y, cluster_z, marker='o',label='cluster')
            ax.set_xlabel('X kpc')
            ax.set_ylabel('Y kpc')
            ax.set_zlabel('Z kpc ')
            ax.legend()
            
            
            fig_3d = plt.figure(figsize=(8, 6))
            ax = fig_3d.add_subplot(projection='3d')
            ax.scatter(rel_x[0], rel_y[0], rel_z[0], color="red", s=100, label="Initial Position")  # Highlight start
            ax.plot(rel_x, rel_y, rel_z, color="blue", label="Relative Trajectory")
    
            ax.set_xlabel("X (kpc)")
            ax.set_ylabel("Y (kpc)")
            ax.set_zlabel("Z (kpc)")
            ax.set_title(f"3D Motion of {star_params['Name'][0]} WRT Cluster")
            ax.legend()
            plt.show()
            

        return rel_x,rel_y,rel_z, time_min_sep #star_x_shifted ,star_y_shifted, star_z_shifted
    def calc_cluster_radius(self,members,cluster_params):
        ra = members['ra']
        dec = members['dec']
        
        #compute center as mean
        ra_cen, dec_cen = np.mean(ra), np.mean(dec)
        
        positions = coord.SkyCoord(ra=ra,dec=dec,frame='icrs')
        center = coord.SkyCoord(ra=ra_cen*u.deg,dec=dec_cen*u.deg,frame='icrs')
        
        theta = positions.separation(center).to(u.rad)
        
        #distance
        if 'parallax' in cluster_params.colnames:
            distance = 1/cluster_params['parallax'].value * u.kpc
        else:
            distance = cluster_params['distance_bj']
        projected_radius = distance * theta.value
        return projected_radius, np.max(projected_radius), np.std(projected_radius)
        
    def velocity_check(self, orbit,values=False):
        # Extract Galactocentric Cartesian velocity components
        vx1 = np.array(orbit.v_x.to_value()) #kpc/myr
        vy1 = np.array(orbit.v_y.to_value())  #kpc/myr
        vz1 = np.array(orbit.v_z.to_value()) # kpc/myr
    
        to_km_per_s = 977.8
        
        vx1, vy1, vz1 = vx1 * to_km_per_s * u.km/u.s, vy1 * to_km_per_s * u.km/u.s, vz1 * to_km_per_s * u.km/u.s
        
        # Define Galactocentric position
        x1 = np.array(orbit.x) * u.kpc
        y1 = np.array(orbit.y) * u.kpc
        z1 = np.array(orbit.z) * u.kpc
    
        # Create the combined coordinate and velocity representation
        galacto_rep = coord.SkyCoord(
            x=x1,
            y=y1,
            z=z1,
            v_x=vx1,
            v_y=vy1,
            v_z=vz1,
            representation_type="cartesian",
            frame="galactocentric"
        )
    
        # Transform to Galactic coordinates
        galactic_vel = galacto_rep.transform_to(coord.Galactic())
        pml = galactic_vel.pm_l_cosb
        pmb = galactic_vel.pm_b
        vrad = galactic_vel.radial_velocity
        if values:
            print(pml, pmb)
            # Compute proper motions in mas/yr
    
            # Print results
            print('means')
            print(f"Vx: {np.mean(vx1)}")
            print(f"Vy: {np.mean(vy1)}")
            print(f"Vz: {np.mean(vz1)}")
            print(f"pm_l: {np.mean(pml)}")
            print(f"pm_b: {np.mean(pmb)}")
            print(f"Vrad: {np.mean(vrad)}")
        return vx1, vy1 , vz1
    def plot_with_cluster(self,clustername, cluster_params=None, clustertable=None, savefig=False):
        '''
        plot the integreated motion of a star wrt to a host cluster

        Parameters
        ----------
        clustername : str - name of cluster to plot against
        cluster_params : astropy table. parameter of cluster (distance, l,b, radial velocity, proper motions)
        clustertable : astropy table, optional
         position of cluster mebers in galactic coorindates
        savefig : boolean optional
            DESCRIPTION. save the figure for plotting

        Note on the proper motion arrows
        i need the angle between the proper motions but later  i invert the x-axis so i should flip the 
        flip the pm_l arrow
        Returns
        -------
        None

        '''
        
        single_star = self.table
        source_id = single_star['source_id']
        long_path, lat_path, _, ticks = self.trace_linear_path(source_id,None,int_time=self.int_time)
        #convert to numpy arrays
        long_path, lat_path, ticks = np.array(long_path), np.array(lat_path),np.array(ticks)
        
        
        #inital positions
        
       # l_naught = single_star['l']
       # b_naught = single_star['b']
        
        #proper motion vectors
   
        #member clusters
        if clustertable is not None:
            member_long = clustertable['l']
            member_lat = clustertable['b']
        
        #include gala
        
        orbit  = self.trace_galactic_path(single_star,int_time=self.int_time)
      #  self.velocity_check(orbit)
        
       # sun_orbit = self.trace_galactic_path(self.sun_params, self.int_time)
        #print(orbit)
        x1 = np.array(orbit.x)
        y1 = np.array(orbit.y)
        z1 = np.array(orbit.z)
        
        rel_x1 = x1# - sun_x
        rel_y1 = y1# - sun_y
        rel_z1 = z1 #- sun_z

        coord_rep = coord.SkyCoord(x=rel_x1*u.kpc, y=rel_y1*u.kpc, z=rel_z1*u.kpc, frame='galactocentric')
        xyz_galactic = coord_rep.transform_to(coord.GalacticLSR())
        gala_l_vals = xyz_galactic.l.deg # Galactic longitude in degrees
        gala_b_vals = xyz_galactic.b.deg

        
        plt.figure(figsize=(10,5))
        if clustertable is not None:
            plt.scatter(member_long,member_lat, s=50,marker='*',label=f'{clustername}',color='xkcd:grey')
        
        #plot the orbit integrated path of the cluster
        if cluster_params is not None:
            cluster_orbit = self.trace_galactic_path(cluster_params,int_time=self.int_time)
            #print('cluster velocity check')
           # self.velocity_check(cluster_orbit)
            cluster_x1  = np.array(cluster_orbit.x)
            cluster_y1 = np.array(cluster_orbit.y)
            cluster_z1 = np.array(cluster_orbit.z)
            
            cl_rel_x = cluster_x1# - sun_x
            cl_rel_y = cluster_y1# - sun_y
            cl_rel_z = cluster_z1# - sun_z
            cluster_galactocentric  = coord.SkyCoord(x=cl_rel_x*u.kpc, y=cl_rel_y*u.kpc,z=cl_rel_z*u.kpc, frame = coord.Galactocentric())
            
            cluster_galactic = cluster_galactocentric.transform_to(coord.GalacticLSR())
            
            cluster_l_vals = cluster_galactic.l
            cluster_b_vals = cluster_galactic.b
            #huuh?
            # rel_x, rel_y,rel_z, time_min_sep  = self.plot_comoving_cluster(single_star, cluster_params, plotting=False)
            # print(f'Minimum separation time {time_min_sep}')
            # relative_sky = coord.SkyCoord(x=rel_x, y=rel_y, z=rel_z, unit=u.kpc,
            #                 representation_type='cartesian', frame='galactocentric').transform_to('galactic')
            
            plt.scatter(cluster_l_vals,cluster_b_vals,label='Cluster Path',color='xkcd:dark grey')
            
            plt.plot(cluster_l_vals,cluster_b_vals,color='xkcd:grey')
            
            # cl_l_arrow, cl_b_arrow,_,cl_ticks = self.cluster_linear_path(cluster_params, int_time=0.5,time_step=1000)
            # cl_l, cl_b = cluster_params['l'], cluster_params['b']
            cl_arrow_orbit = self.trace_galactic_path(cluster_params, 0.5)
            cl_arrow_x, cl_arrow_y, cl_arrow_z = np.array(cl_arrow_orbit.x), np.array(cl_arrow_orbit.y), np.array(cl_arrow_orbit.z)
            #delta_cl_l = cl_l_arrow[-1] - cl_l
            #delta_cl_b = cl_b_arrow[-1] - cl_b
            cl_arrow_galactocentric = coord.SkyCoord(x=cl_arrow_x*u.kpc, y=cl_arrow_y*u.kpc,
                                                     z=cl_arrow_z*u.kpc,frame=coord.Galactocentric())
            cl_arrow_galactic = cl_arrow_galactocentric.transform_to(coord.GalacticLSR())
            cl_arrow_l = cl_arrow_galactic.l
            cl_arrow_b  = cl_arrow_galactic.b
            
            cl_arrow_plot_l, cl_arrow_plot_b = cl_arrow_l[-1].value, cl_arrow_b[-1].value
            #proper motion arrows

            plt.plot(cl_arrow_l, cl_arrow_b, color='xkcd:orange')

        N = len(long_path)
        #print(abs(ttotal))
        star_name = str(single_star['Name'].value[0])
        #plot path and color plot by specrtral type
        sp_type = single_star['Mod_SpType'][0] #dumb
        path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'

        plt.scatter(gala_l_vals[0],gala_b_vals[0],color=path_color,label=star_name,marker='*')
        #gala
        plt.plot(gala_l_vals,gala_b_vals, lw=1.1, marker='o', color=path_color,label=f'Retraced Path {self.int_time} Myrs ago')
        #pm arrows
        # i need the angle between the proper motions but later  i invert the x-axis so i should 
        star_arrow_orbit = self.trace_galactic_path(single_star, 0.5)
        
        star_arrow_x, star_arrow_y, star_arrow_z = np.array(star_arrow_orbit.x), np.array(star_arrow_orbit.y), np.array(star_arrow_orbit.z)

        star_arrow_galactocentric = coord.SkyCoord(x=star_arrow_x*u.kpc, y=star_arrow_y*u.kpc,
                                                 z=star_arrow_z*u.kpc,frame=coord.Galactocentric())
        star_arrow_galactic = star_arrow_galactocentric.transform_to(coord.GalacticLSR())
        star_arrow_l = star_arrow_galactic.l
        star_arrow_b  = star_arrow_galactic.b
        plt.plot(star_arrow_l, star_arrow_b,color='xkcd:light blue')
        
        #plot cluster integrated path

        #tick markers
        # for tick_l, tick_b in ticks:
        #     plt.scatter(tick_l,tick_b,color='xkcd:canary',s=20,marker='x')
        plt.xlabel('Galactic Longitude(deg)')
        plt.ylabel('Galactic Latitude (deg)')
        plt.gca().invert_xaxis()
        plt.title(f"Retraced path of {star_name}")
        plt.grid(True)
        plt.legend(loc='best')
        today = datetime.now().strftime("%Y%m%d")
        if savefig == True:
            if clustertable == None:
                fig_clus_name = ' '
            else:
                fig_clus_name = clustername
                
            plt.savefig(parentdir+'/Figures/Traceback/'+f"{star_name}_with_{fig_clus_name}_{today}.png")
        plt.show()
        return orbit
    def plot_in_galactocentric(self, stars,figname, cluster=None, cluster_members=None, savefig=False):
        """
        Plot the integrated orbits of the star and cluster in the Galactocentric frame
        along with their 2D projections (XY, YZ, ZX).
    
        Parameters
        ----------
        star : astropy Table or dict
            The star's parameters (including 'l', 'b', 'distance_bj', 'RV', etc.).
        cluster : astropy Table or dict
            The cluster's parameters (including 'l', 'b', 'distance_bj', 'RV', etc.).
        savefig : bool, optional
            If True, save the figure to disk.
    
        Returns
        -------
        star_orbit, cluster_orbit : Orbit objects
            The integrated orbits in the Galactocentric frame.
        """
        #initialize figure
        
        fig = plt.figure(figsize=(16, 12))
        ax3d = fig.add_subplot(221, projection='3d')
        
        ax_xy = fig.add_subplot(222)
        ax_yz = fig.add_subplot(223)
        ax_zx = fig.add_subplot(224)
        
        cluster_color = 'xkcd:steel blue'
        arrow_color= 'xkcd:orange'
        
        import matplotlib as mpl
        mpl.rcParams.update({
        'font.size': 20,        # Adjust as needed
        'axes.titlesize': 22,
        'axes.labelsize': 20,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        'legend.fontsize': 10,
        'lines.linewidth': 2,   # Thicker lines for better visibility
    })
        # Integrate the orbits for the star and the cluster.
        for star in stars:
            star = Table(star)
            star_orbit = self.trace_galactic_path(star, int_time=self.int_time)
            ts = star_orbit.t  
            star_short = self.trace_galactic_path(star, int_time=0.5)
            plot_time = str(abs(round(ts[-1].value,2)))
            # Extract Cartesian positions (assumed in kpc) from the orbits.
            # We assume star_orbit.x, .y, .z are Quantities.
            star_x = np.array(star_orbit.x.to(u.kpc).value)
            star_y = np.array(star_orbit.y.to(u.kpc).value)
            star_z = np.array(star_orbit.z.to(u.kpc).value)
            
            x1s, y1s, z1s = star_short.x[0].to(u.kpc).value, star_short.y[0].to(u.kpc).value, star_short.z[0].to(u.kpc).value
            
            x2s, y2s, z2s = star_short.x[-1].to(u.kpc).value, star_short.y[-1].to(u.kpc).value, star_short.z[-1].to(u.kpc).value
            star_color = star['Mod_SpType'][0]
            star_label= f"{star['Name'][0]}"
            plot_time = str(abs(round(ts[-1].value,2)))
            
            font_dict = {'fontsize':12,'color':star_color}
            
            ax3d.plot(star_x, star_y, star_z, c=star_color, label="Path of "+star_label)
            ax3d.scatter(star_x, star_y, star_z, c=star_color)
            ax3d.plot([x1s, x2s], [y1s, y2s], [z1s, z2s], color=arrow_color, linestyle='--')
            ax3d.scatter(star_x[0], star_y[0], star_z[0], c=star_color, marker='*', s=200, label="Current position of "+star_label)
            ax3d.scatter(star_x[-1], star_y[-1], star_z[-1], c=star_color, marker='x', s=150, label="Position of "+star_label+f' {plot_time} Myr ago')
            
            #arrow tip
            ax3d.plot(x2s,y2s, z2s, color=arrow_color, marker='^')
            ax3d.tick_params(labelsize=11)
            
            ax_xy.plot(star_x, star_y, color=star_color, label='path of '+star_label)

            ax_xy.text(star_x[0],star_y[0],s=star_label, **font_dict)
            ax_xy.scatter(star_x, star_y, color=star_color)
            ax_xy.scatter(star_x[0], star_y[0], color=star_color, marker='*', s=200,label='Current position of '+star_label)
            ax_xy.scatter(star_x[-1], star_y[-1], color=star_color, marker='x', s=150, label="Position of "+star_label+f' {plot_time} Myr ago')
            ax_xy.plot([x1s, x2s], [y1s, y2s], color=arrow_color, linestyle='--')
            
            
            # YZ projection
            ax_yz.plot(star_y, star_z, color=star_color, label='path of '+star_label)

            
           
            # Mark the first and last positions for the cluster
            ax_yz.text(star_y[0],star_z[0],s=star_label, **font_dict)
            ax_yz.scatter(star_y[0], star_z[0], color=star_color, marker='*', s=200,label='Current position of '+star_label)
            ax_yz.scatter(star_y[-1], star_z[-1], color=star_color, marker='x', s=150, label="Position of "+star_label+f' {plot_time} Myr ago')

            ax_yz.scatter(star_y, star_z, color=star_color)
        
            ax_yz.plot([y1s, y2s], [z1s, z2s], color=arrow_color, linestyle='--')
            
            #ZX Projection
            ax_zx.text(star_x[0],star_z[0],s=star_label, **font_dict)
            ax_zx.scatter(star_x[0], star_z[0], color=star_color, marker='*', s=200,label='Current position of '+star_label)
            ax_zx.scatter(star_x[-1], star_z[-1], color=star_color, marker='x', s=150, label="Position of "+star_label+f' {plot_time} Myr ago')
            ax_zx.scatter(star_x, star_z, color=star_color)
            ax_zx.plot(star_x, star_z, color=star_color)
   
            #arrows
            ax_zx.plot([x1s, x2s], [z1s, z2s], color=arrow_color, linestyle='--')
     
        if cluster is not None:
            plt.suptitle(f'Galactocentric Orbit of {stars['Name'][0]} with {cluster['Name'][0]}', fontsize=25)
            cluster_orbit = self.trace_galactic_path(cluster, int_time=self.int_time)
            cluster_label = f"{cluster['Name'][0]}"
            
            cluster_x = np.array(cluster_orbit.x.to(u.kpc).value)
            cluster_y = np.array(cluster_orbit.y.to(u.kpc).value)
            cluster_z = np.array(cluster_orbit.z.to(u.kpc).value)
            #arrow of motion 
            # Integrate both orbits over a short 0.5 Myr interval for motion vector plotting
            #star_short = self.trace_galactic_path(star, int_time=0.5)
            cluster_short = self.trace_galactic_path(cluster, int_time=0.5)

            x1c, y1c, z1c = cluster_short.x[0].to(u.kpc).value, cluster_short.y[0].to(u.kpc).value, cluster_short.z[0].to(u.kpc).value
            x2c, y2c, z2c = cluster_short.x[-1].to(u.kpc).value, cluster_short.y[-1].to(u.kpc).value, cluster_short.z[-1].to(u.kpc).value

            ax3d.plot(cluster_x, cluster_y, cluster_z, color=cluster_color, label="path of "+cluster_label)
       
            ax3d.plot([x1s, x2s], [y1s, y2s], [z1s, z2s], color=arrow_color, linestyle='--')
       #arrow tip
            ax3d.plot(x2s,y2s, z2s, color=arrow_color, marker='^')
            ax3d.plot([x1c, x2c], [y1c, y2c], [z1c, z2c], color=arrow_color, linestyle='--')
            ax3d.plot(x2c,y2c, z2c, color=arrow_color, marker='^')
            
            ax3d.plot(cluster_x, cluster_y, cluster_z, color=cluster_color, label="path of "+cluster_label)
            
            # Mark the first and last position for the cluster
            ax3d.scatter(cluster_x[0], cluster_y[0], cluster_z[0], c=cluster_color, marker='o', s=100, label=f"Current Position of {cluster_label}")
            ax3d.scatter(cluster_x[-1], cluster_y[-1], cluster_z[-1], c=cluster_color, marker='x', s=100, label=f"Position of {cluster_label} {plot_time} Myr ago")
            
            #xy projection 
            ax_xy.plot(cluster_x, cluster_y, color=cluster_color, label='path of '+cluster_label)
            
           
            ax_xy.scatter(cluster_x, cluster_y, color=cluster_color)

            
            cluster_font_dict= {'fontsize':12,'color':cluster_color}
            ax_xy.text(cluster_x[0],cluster_y[0],s=cluster_label, **cluster_font_dict)
            ax_xy.scatter(cluster_x[0], cluster_y[0], color=cluster_color, marker='o', s=100)
            ax_xy.scatter(cluster_x[-1], cluster_y[-1], color=cluster_color, marker='x', s=100)
            
            #arrows 
            ax_xy.plot([x1c, x2c], [y1c, y2c], color=arrow_color, linestyle='--')

            #yz projection 
            ax_yz.plot(cluster_y, cluster_z, color=cluster_color, label='path of '+cluster_label)
            ax_yz.text(cluster_y[0],cluster_z[0],s=cluster_label, **cluster_font_dict)
           
       
            ax_yz.scatter(cluster_y[0], cluster_z[0], color=cluster_color, marker='o', s=100)
            ax_yz.scatter(cluster_y[-1], cluster_z[-1], color=cluster_color, marker='x', s=100,zorder=1)
            
            # ax_yz.scatter(star_y, star_z, color=star_color)
            ax_yz.scatter(cluster_y, cluster_z, color=cluster_color)
        
        
            ax_yz.plot([y1c, y2c], [z1c, z2c], color=arrow_color, linestyle='--')
            #xz projection
            ax_zx.plot(cluster_x, cluster_z, color=cluster_color, label='path of '+cluster_label)
            ax_zx.text(cluster_x[0],cluster_z[0],s=cluster_label, **cluster_font_dict)
            

            ax_zx.scatter(cluster_x, cluster_z, color=cluster_color)
            
           

            ax_zx.scatter(cluster_x[0], cluster_z[0], color=cluster_color, marker='o', s=100)
            ax_zx.scatter(cluster_x[-1], cluster_z[-1], color=cluster_color, marker='x', s=100)
            #arrows

            ax_zx.plot([x1c, x2c], [z1c, z2c], color=arrow_color, linestyle='--')
            


        
       
        #for labeling

        
        # star_label= f"{star['Name'][0]}"

        plot_time = str(abs(round(ts[-1].value,2)))
        # Extract l and b from cluster_members and directly apply units
        
        if cluster_members is not None:
            ra_mem = cluster_members['ra'] 
            dec_mem = cluster_members['dec'] 
            cl_dist = cluster['distance_bj']
            cl_rv = cluster['RV'] 
            cl_pmra = cluster['pmra'] 
            cl_pmdec = cluster['pmdec']
            cl_parallax = cluster['parallax']
            
            
            distance_mem = np.ones(len(ra_mem)) * cl_dist * u.kpc
            pmra_mem = np.ones(len(ra_mem)) * cl_pmra * u.mas/u.yr
            pmdec_mem = np.ones(len(ra_mem)) * cl_pmdec * u.mas/u.yr
            rv_mem = np.ones(len(ra_mem)) * cl_rv * u.km/u.s
            parallax_mem = np.ones(len(ra_mem)) * cl_parallax * u.mas
    
            # Create the member_table with properly formatted columns
            member_table = Table(
        [ra_mem, dec_mem, distance_mem, pmra_mem, pmdec_mem, rv_mem,parallax_mem],
        names=('ra', 'dec', 'distance_bj', 'pmra', 'pmdec', 'RV','parallax')
    )
 
            
             # Extract the integrated Cartesian positions (assumed to be stored as Quantities in kpc)
            member_orbit  = self.trace_galactic_path(member_table, self.int_time)
            x_mem = member_orbit.x.to(u.kpc)
            y_mem = member_orbit.y.to(u.kpc)
            z_mem = member_orbit.z.to(u.kpc)
            
            # First and last positions
            mem_curr_color = 'xkcd:dark'
            mem_prev_color = 'xkcd:light grey'
            mem_alpha = 0.8
            
            zorder = -1
            ax_yz.scatter(y_mem[0], z_mem[0], color=mem_curr_color, label=f"Current position of {cluster_label} members",alpha=mem_alpha,zorder=zorder)
            ax_yz.scatter(y_mem[-1], z_mem[-1],color=mem_prev_color, label=f'Position of {cluster_label} members {plot_time} Myrs ago',zorder=zorder)
            ax_xy.scatter(x_mem[0], y_mem[0], color=mem_curr_color,label=f"Current position of {cluster_label} members",alpha=mem_alpha,zorder=zorder)
            ax_xy.scatter(x_mem[-1], y_mem[-1], color=mem_prev_color,label=f'Position of {cluster_label} members {plot_time} Myrs ago',alpha=mem_alpha,zorder=zorder)
 
            
            ax_zx.scatter(x_mem[-1], z_mem[-1],color=mem_prev_color, label=f'Position of {cluster_label} members {plot_time} Myrs ago',zorder=zorder)
            ax_zx.scatter(x_mem[0], z_mem[0], color=mem_curr_color,label=f"Current position of {cluster_label} members",alpha=mem_alpha,zorder=zorder)
            
    
 
        ax3d.set_xlabel("X (kpc)")
        ax3d.set_ylabel("Y (kpc)")
        ax3d.set_zlabel("Z (kpc)")
        ax3d.set_title("3D Galactocentric Orbit")
        #ax3d.legend(loc = 'best')
    


    
        ax_xy.set_xlabel("X (kpc)")
        ax_xy.set_ylabel("Y (kpc)")
        ax_xy.set_title("XY Projection")
        ax_xy.legend()
    

        ax_yz.set_xlabel("Y (kpc)")
        ax_yz.set_ylabel("Z (kpc)")
        ax_yz.set_title("yz Projection")
        ax_yz.legend()
    
    

    
        ax_zx.set_xlabel("x (kpc)")
        ax_zx.set_ylabel("z (kpc)")
        ax_zx.set_title("zx Projection")
        ax_zx.legend()
        
        # xlow, xhigh = -6.63,-6.35
        # ylow, yhigh = -1.5, 0
        # zlow, zhigh  = 0.020, 0.085
        
        # ax3d.view_init(elev=27, azim=40)
        # ax3d.set_xlim(xlow, xhigh)
        # ax3d.set_ylim(ylow, yhigh)
        # ax3d.set_zlim(zlow, zhigh)
        
        # # XY projection
        # ax_xy.set_xlim(xlow, xhigh)
        # ax_xy.set_ylim(ylow, yhigh)
        
        # # # YZ projection
        # ax_yz.set_xlim(ylow, yhigh)
        # ax_yz.set_ylim(zlow, zhigh)
        
        # # # ZX projection
        # ax_zx.set_xlim(xlow, xhigh)
        # ax_zx.set_ylim(zlow, zhigh)
    
        # Show the plot
        plt.tight_layout()

        if savefig:
            import os
            from datetime import datetime
            mydir = os.path.dirname(os.path.realpath(__file__))
            parentdir = os.path.dirname(mydir)
            today = datetime.now().strftime("%Y%m%d")
            save_path = os.path.join(parentdir, 'Figures', 'Traceback/Galactocentric', f"orbits_{figname}_{today}.png")
            plt.savefig(save_path)
    
        plt.show()
    

        if savefig:
            mydir = os.path.dirname(os.path.realpath(__file__))
            parentdir = os.path.dirname(mydir)
            today = datetime.now().strftime("%Y%m%d")
            save_path = os.path.join(parentdir, 'Figures', 'Traceback/Galactocentric', f"Galactocentric_cylindrical_{figname}_{today}.png")
            plt.savefig(save_path)
        plt.show()
        
    def plot_multiple_stars_in_galactocentric(self, stars_table):
        """
        Plot 3D Galactocentric orbits and XY projection for multiple stars from an Astropy Table.
    
        Parameters
        ----------
        stars_table : astropy Table
            Table containing multiple stars with required kinematic columns.
        """
        fig = plt.figure(figsize=(14, 6))
    
        # Subplots
        ax3d = fig.add_subplot(121, projection='3d')
        ax_xy = fig.add_subplot(122)
    
        prlx_mask = stars_table['parallax']/stars_table['parallax_error'] >=5.0
        stars_table = stars_table[prlx_mask]
        for star in stars_table:
            #star loses units, very DUMB
            #add units back with Table(star), then it rememebers the units
            star_temp = Table(star)
            star_orbit = self.trace_galactic_path(star_temp, int_time=self.int_time)
    
            # Positions in kpc
            x = np.array(star_orbit.x.to(u.kpc).value)
            y = np.array(star_orbit.y.to(u.kpc).value)
            z = np.array(star_orbit.z.to(u.kpc).value)
    
            # Color from spectral type or fallback
            color = star.get('Mod_SpType', 'black')
    
            # 3D plot
            ax3d.plot(x, y, z, color=color)
            ax3d.scatter(x, y, z, color=color, s=1, alpha=0.7)
    
            # XY projection
            ax_xy.plot(x, y, color=color)
            ax_xy.scatter(x[0], y[0], color=color, s=20, alpha=0.7)
    
        # Correct Galactocentric Sun position
        sun_x = -8.15  # kpc
        sun_y = 0.0
        sun_z = 0.0208
    
        # Mark the Sun in both projections
        ax3d.scatter(sun_x, sun_y, sun_z, color='yellow', s=100, marker='*')
        ax_xy.scatter(sun_x, sun_y, color='yellow', s=300,edgecolors='xkcd:black', marker='o',label='Sun',zorder=5)
    
        # Labels and layout
        ax3d.set_xlabel("X (kpc)")
        ax3d.set_ylabel("Y (kpc)")
        ax3d.set_zlabel("Z (kpc)")
        ax3d.set_title("3D Galactocentric Orbits")
    
        ax_xy.set_xlabel("X (kpc)")
        ax_xy.set_ylabel("Y (kpc)")
        ax_xy.set_title("XY Projection")
        ax_xy.legend()
        plt.tight_layout()
        plt.grid(True)
        plt.show()
        
        # fig_cyl = plt.figure(figsize=(8, 6))
        # ax_rz = fig_cyl.add_subplot(111)
        
        # for star in stars_table:
        #     # Restore units to the star data
        #     star_temp = Table(star)
        #     star_orbit = self.trace_galactic_path(star_temp, int_time=self.int_time)
        
        #     # Convert orbit to cylindrical representation
        #     cyl_orbit = star_orbit.represent_as('cylindrical')
        
        #     # Extract cylindrical components
        #     rho = cyl_orbit.rho.to(u.kpc).value
        #     z = cyl_orbit.z.to(u.kpc).value
        
        #     # Determine color based on spectral type or default to black
        #     color = star.get('Mod_SpType', 'black')
        
        #     # 2D plot in cylindrical coordinates (ρ vs. z)
        #     ax_rz.plot(rho, z, color=color)
        #     ax_rz.scatter(rho[0], z[0], color=color, s=20, alpha=0.7)
        
        # # Mark the Sun's position in the plot
        # sun_x = -8.15  # kpc
        # sun_y = 0.0
        # sun_z = 0.0208  # kpc
        # sun_rho = np.sqrt(sun_x**2 + sun_y**2)
        
        # ax_rz.scatter(sun_rho, sun_z, color='yellow', s=100, marker='o', label='Sun')
        
        # # Set labels and title
        # ax_rz.set_xlabel("ρ (kpc)")
        # ax_rz.set_ylabel("Z (kpc)")
        # ax_rz.set_title("Meridional Plane (ρ vs. Z)")
        # ax_rz.legend()
        
        # plt.tight_layout()
        # plt.grid(True)
        # plt.show()

        return None


test_table = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/HMXB_20250527_.ecsv',format='ecsv')
test_170037 = test_table[test_table['Name']=='4U 1700-377']
ngc6231_members = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/Ngc6231_members_galacitc.ecsv',format='ecsv')
ngc6231_params = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/NGC2631_param.ecsv')


plot_galactocentric_table = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/HMXB_20250415_.ecsv',format='ecsv')
mydir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(mydir)
if __name__ == "__main__":
   # GalacticTraceback(test_170037,-3.0).plot_trace_aitoff(savefig=False, cluster_params=None)
   #GalacticTraceback(test_170037,-3.0).plot_with_cluster(clustername='NGC 6231',cluster_params=ngc6231_params, clustertable=ngc6231_members,savefig=False)
  #GalacticTraceback(test_170037,-3.0).plot_with_cluster(clustername='NGC 6231',cluster_params=None, clustertable=None,savefig=False)
   
    x= GalacticTraceback(test_170037, -3.0).plot_in_galactocentric(test_170037, figname='170037_NGC6231',cluster=ngc6231_params,cluster_members=ngc6231_members,savefig=True)
    #GalacticTraceback(test_170037, -5.0).plot_separation(test_170037, ngc6231_params,savefig=True)
    #GalacticTraceback(None, 5.0).plot_multiple_stars_in_galactocentric(plot_galactocentric_table)
   #x,y,z,t_min_sep = GalacticTraceback(test_170037,int_time=-3.0).plot_comoving_cluster(test_170037, ngc6231_params,plotting=True)
#    GalacticTraceback(test_table).plot_trace(savefig=False
