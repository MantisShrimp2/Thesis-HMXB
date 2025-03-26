#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:04:32 2024

@author: karan
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
    def trace_linear_path_vdm(self,table,time_step=1000):
        '''
        from van der meji 2021
        convert  position to cartisian coordinates
        
        integrate with distance d = vr*dt 
        
        convert back to galactic longitude and galactic latitude
        
        
        

        Parameters
        ----------
        table : pd dataframe or astropy table.

        Returns
        -------
        None.

        '''
        ra = np.radians(table['ra'])[0]
        dec = np.radians(table['dec'])[0]
        prlx = table['parallax'][0]
        pmra= table['pmra'][0]
        pmdec = table['pmdec'][0]
        dist = table['distance_bj'][0]
        RV = table['RV'][0]
        
        long = np.radians(table['l'][0])
        lat = np.radians(table['b'][0])

        
        #inital conditions 
        xo = dist*np.cos(lat)*np.cos(long)
        yo = dist*np.cos(lat)*np.sin(long)
        zo  = dist* np.sin(lat)

        p_hat =np.array([-np.sin(ra), np.cos(ra), 0.0])
        q_hat = np.array([-np.sin(dec)*np.cos(ra), - np.sin(dec)*np.sin(ra), np.cos(dec)])
        r_hat= np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra),np.sin(dec)])
        
        #transformation matrix
        
        D_A  = dist*self.k # im gonna use this alot 
        matrix_t  = [
            [p_hat[0]*D_A, q_hat[0]*D_A, r_hat[0]],
                     [p_hat[1]*D_A,q_hat[1]*D_A,r_hat[1]],
                     [p_hat[2]*D_A,q_hat[2]*D_A,r_hat[2]]
                     ]
        galactic_vel = np.array([pmra, pmdec, RV])
        # vra, vdec = D_A*pmra, D_A*pmdec
        # v_vec = vra * p_hat + vdec * q_hat + RV * r_hat
        cartesian_vel = matrix_t @ galactic_vel
        #calculate velocities in cartesian coordinates    # 1 km/s ≈ 1.0227e-6 kpc/yr.
        conv = 1.02201e-9 #km/s to kpc/yr
       # v_vec[0], v_vec[1], v_vec[2]
  
        vx, vy, vz = cartesian_vel[0], cartesian_vel[1], cartesian_vel[2] 
        vx = vx*conv
        vy = vy*conv
        vz = vz*conv
        total_int_time =np.abs(self.int_time* 1e6)
        max_steps =int(total_int_time/time_step)
    
        current_time = 0
        #initalize postion
        xi, yi, zi = xo, yo, zo
        
        x_path, y_path, z_path = [xo], [yo], [zo]   
        for step  in range(max_steps):
                # Update coordinates using Euler method
                #fixed timestep
        
                current_time += time_step
                
                xi -=   vx*time_step
                yi -=  vy* time_step
                zi -=   vz* time_step
                #Wrap longitude to [0, 360]
                x_path.append(xi)
                y_path.append(yi)
                z_path.append(zi)
                step +=1
                if step >1e6:
                    print('broke loop')
                    break
 
        #convert to galactic coordinates
        
        x_path = np.array(x_path)
        y_path = np.array(y_path)
        z_path = np.array(z_path)
        #atan2 is verified by GAIA
        l_path = (np.degrees(np.arctan2(y_path,x_path))) % 360
        b_path = np.degrees(np.arctan2(z_path, np.sqrt(x_path**2 + y_path**2)))
        dist_path = np.sqrt(x_path **2 + y_path **2 + z_path **2)
       # b_path = np.degrees(np.arcsin(z_path/dist_path))
        return l_path, b_path, dist_path
        
        
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
        l = row['l']
        b = row['b'] 
        mu_l = row['pm_l_poleski']
        mu_b = row['pm_b_poleski']
        dist = row['distance_bj'] #kpc
        radial_velocity = row['RV'] 

       # k = 4.74 * (u.km/u.s)/(u.mas *u.kpc/u.yr) #km/s per mas/yr 
        #print(radial_velocity)
        #transform to galactic frame
        #from carreto-castrillo 2023
        mu_l_total  = mu_l
        mu_b_total = mu_b
        dist_total  = dist 
        radial_velocity_total = radial_velocity
            
        with coord.galactocentric_frame_defaults.set('v4.0'):
            galcen_frame = coord.Galactocentric()
        galactic_rep = coord.SkyCoord(l=l,b=b,pm_l_cosb=mu_l_total,pm_b=mu_b_total,distance=dist_total,
                                      radial_velocity =radial_velocity_total,frame=coord.Galactic())
        #transform frame
        star_galacto = galactic_rep.transform_to(galcen_frame)

        #correct for solar motion
        #star_galacto = reflex_correct(star_galacto)
       # print(star_galacto)

        initial_pos = gd.PhaseSpacePosition(star_galacto.data)
        total_time = int_time *u.Myr
        dt = np.sign(int_time)* 0.1 *u.Myr
        n_steps = int(total_time.to_value(u.Myr) / dt.to_value(u.Myr))

        integrator = gi.LeapfrogIntegrator
        potential = gp.MilkyWayPotential2022()  
        #orbit = potential.integrate_orbit(initial_pos, dt=dt, t1=0, t2=total_time)
        orbit = potential.integrate_orbit(initial_pos, 
                                          dt=dt,n_steps=n_steps,
                                          Integrator=integrator)
                

        return orbit 

    def plot_trace(self,savefig=False, cluster_params=None):
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
                z_naught = data['distance_bj']*np.sin(np.radians(data['b']))*1000
                l_naught = data['l']
                b_naught = data['b']
                arrow_pml, arrow_pmb, arrow_pmz, _ = self.trace_linear_path(ids, None, 0.5,time_step=1000)
                long_path, lat_path, z_path,ticks = self.trace_linear_path(ids, None,int_time = -3.0, time_step=1000)
                N = len(long_path)
                #plot path and color plot by specrtral type
                sp_type = data['Mod_SpType'][0] #dumb
                path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'
                #plt.scatter(long_path[1:N], z_path[1:N],s=3,color='xkcd:black',alpha=0.7)
                plt.scatter(long_path[1:N], lat_path[1:N],s=2,color=path_color,alpha=0.7)
                #plot star current position
                
                plt.scatter(l_naught,b_naught,color='xkcd:black',s=50,label='')
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
                    
                    plt.scatter(tick_l, tick_b, color='xkcd:orange', s=20)
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
        plt.gca().invert_xaxis() # resverse the x-axis, standrd to show longitude
        plt.ylabel('Galactic Latitude (deg)')
        plt.title('Traceback Paths in Galactic Coordinates')
        #plt.legend()
        plt.grid(True)
        
        today = datetime.now().strftime("%Y%m%d")
        if savefig == True:
            plt.savefig(parentdir+'/Figures/'+f"Tracepath_{today}.png")
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

    def velocity_check(self, orbit):
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
        long_path, lat_path, _, ticks = self.trace_linear_path(source_id,None,int_time=-3.0)
        #convert to numpy arrays
        long_path, lat_path, ticks = np.array(long_path), np.array(lat_path),np.array(ticks)
        
       # long_vdm_path, lat_vdm_path, _  = self.trace_linear_path_vdm(single_star)
        
        
        
        #inital positions
        
        l_naught = single_star['l']
        b_naught = single_star['b']
        
        #proper motion vectors
   
        #member clusters
        if clustertable is not None:
            member_long = clustertable['l']
            member_lat = clustertable['b']
        
        #include gala
        
        orbit  = self.trace_galactic_path(single_star,int_time=-3.0)
        self.velocity_check(orbit)
        #print(orbit)
        
        x1 = np.array(orbit.x)
        y1 = np.array(orbit.y)
        z1 = np.array(orbit.z)
        
        coord_rep = coord.SkyCoord(x=x1*u.kpc, y=y1*u.kpc, z=z1*u.kpc, representation_type="cartesian", frame='galactocentric')
        xyz_galactic = coord_rep.transform_to(coord.Galactic())
        
        gala_l_vals = xyz_galactic.l.deg # Galactic longitude in degrees
        gala_b_vals = xyz_galactic.b.deg

        
        plt.figure(figsize=(10,5))
        if clustertable is not None:
            plt.scatter(member_long,member_lat, s=50,marker='*',label=f'{clustername}',color='xkcd:grey')
        
        #plot the orbit integrated path of the cluster
        if cluster_params is not None:
            cluster_orbit = self.trace_galactic_path(cluster_params,int_time=-3.0)
            print('cluster velocity check')
            self.velocity_check(cluster_orbit)
            cluster_x1  = np.array(cluster_orbit.x)
            cluster_y1 = np.array(cluster_orbit.y)
            cluster_z1 = np.array(cluster_orbit.z)
            
            
            cluster_galactocentric  = coord.SkyCoord(x=cluster_x1*u.kpc, y=cluster_y1*u.kpc,z=cluster_z1*u.kpc, frame = coord.Galactocentric())
            
            cluster_galactic = cluster_galactocentric.transform_to(coord.Galactic())
            
            cluster_l_vals = cluster_galactic.l
            cluster_b_vals = cluster_galactic.b
            #huuh?
            rel_x, rel_y,rel_z, time_min_sep  = self.plot_comoving_cluster(single_star, cluster_params, plotting=False)
            print(f'Minimum Seperation time {time_min_sep}')
            relative_sky = coord.SkyCoord(x=rel_x, y=rel_y, z=rel_z, unit=u.kpc,
                            representation_type='cartesian', frame='galactocentric').transform_to('galactic')
            
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
            cl_arrow_galactic = cl_arrow_galactocentric.transform_to(coord.Galactic())
            cl_arrow_l = cl_arrow_galactic.l
            cl_arrow_b  = cl_arrow_galactic.b
            
            cl_arrow_plot_l, cl_arrow_plot_b = cl_arrow_l[-1].value, cl_arrow_b[-1].value
            #proper motion arrows

            plt.plot(cl_arrow_l, cl_arrow_b, color='xkcd:orange')
      
            # plt.quiver(cluster_params['l'], cluster_params['b'], -cluster_params['pm_l_poleski'], cluster_params['pm_b_poleski'], angles='uv', width=0.002)
            
        #     arrow_pml= single_star['pm_l_poleski']        
        #     arrow_pmb  = single_star['pm_b_poleski']   
        #     comove_l = relative_sky.l.deg % 360
        #     comove_b =relative_sky.b.deg
        # else:
        #     arrow_pml= single_star['pm_l_poleski'] 
        #     arrow_pmb  = single_star['pm_b_poleski']
            

      
        N = len(long_path)
        #print(abs(ttotal))
        star_name = str(single_star['Name'].value[0])
        #plot path and color plot by specrtral type
        sp_type = single_star['Mod_SpType'][0] #dumb
        path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'
        #linear path
       # plt.scatter(long_path[0:N], lat_path[0:N],s=10,color='xkcd:black',alpha=0.7)
        #integration from van der meji 2021

        plt.scatter(l_naught,b_naught,color=path_color,label=star_name)
        #gala
        plt.plot(gala_l_vals,gala_b_vals, lw=1.1, marker='*', color=path_color,label='Gala path')
        #pm arrows
        # i need the angle between the proper motions but later  i invert the x-axis so i should 
        star_arrow_orbit = self.trace_galactic_path(single_star, 0.5)
        
        star_arrow_x, star_arrow_y, star_arrow_z = np.array(star_arrow_orbit.x), np.array(star_arrow_orbit.y), np.array(star_arrow_orbit.z)

        star_arrow_galactocentric = coord.SkyCoord(x=star_arrow_x*u.kpc, y=star_arrow_y*u.kpc,
                                                 z=star_arrow_z*u.kpc,frame=coord.Galactocentric())
        star_arrow_galactic = star_arrow_galactocentric.transform_to(coord.Galactic())
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

test_table = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/HMXB_20250301_.ecsv',format='ecsv')
test_170037 = test_table[test_table['Name']=='4U 1700-377']
scoob1 = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/SCO OB1-result.ecsv')
ngc6231_params = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/NGC2631_param.ecsv')


mydir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(mydir)
if __name__ == "__main__":
   GalacticTraceback(test_170037,-3.0).plot_trace(savefig=False, cluster_params=ngc6231_params)
   GalacticTraceback(test_170037,-3.0).plot_with_cluster(clustername='NGC 6231',cluster_params=ngc6231_params, clustertable=scoob1,savefig=False)

  # x,y,z,t_min_sep = GalacticTraceback(test_170037,int_time=-3.0).plot_comoving_cluster(test_170037, ngc6231_params,plotting=True)
#  GalacticTraceback(test_170037, -3.0).trace_linear_path_vdm(test_170037)
#    GalacticTraceback(test_table).plot_trace(savefig=False)
