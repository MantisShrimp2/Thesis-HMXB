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
import gala.integrate as gi
from gala.coordinates import reflex_correct
import astropy.units as u
import os
import astropy.coordinates as coord
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
        
    def trace_linear_path(self, source_id, cluster_params, time_step=1000):
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
        
        
        #limit time stepsto 3 million years
        total_int_time = abs(self.int_time)* 1e6
        max_steps = int(total_int_time/time_step)
        
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
        
        return long_path, lat_path, z_path, ticks
    def trace_linear_path_vdm(self,table,time_step=1000):
        '''
        from van der meji 2021
        convert  position to cartisian coordinates
        
        integrate with distance d = vr*dt 
        
        convert back to galactic longitude and galactic latitude
        
        
        

        Parameters
        ----------
        table : TYPE
            DESCRIPTION.

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
        # if cluster_params is not None:
        #    cluster_coords = coord.SkyCoord(l=cluster_params['l'],b = cluster_params['b'],     
        #                                 pm_l_cosb=cluster_params['pm_l_poleski'],
        #                                 pm_b = cluster_params['pm_b_poleski'],
        #                                 distance = cluster_params['distance_bj'],
        #                                 radial_velocity = cluster_params['RV'], frame='galactic')
        #    cluster_coords = cluster_coords.transform_to(coord.Galactocentric)
        #    cl_vx, cl_vy, cl_vz = cluster_coords.v_x.value[0], cluster_coords.v_y.value[0], cluster_coords.v_z.value[0]
        #    vx, vy, vz = cartesian_vel[0] - cl_vx, cartesian_vel[1] - cl_vy, cartesian_vel[2] - cl_vz # km/s
        # else:    
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
        l_path = (np.degrees(np.arctan2(y_path,x_path))) %360
        b_path = np.degrees(np.arctan2(z_path, np.sqrt(x_path**2 + y_path**2)))
        dist_path = np.sqrt(x_path **2 + y_path **2 + z_path **2)
       # b_path = np.degrees(np.arcsin(z_path/dist_path))
        return l_path, b_path, dist_path
        
        
    def trace_galactic_path(self,source_info):
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
        # if cluster_params is not None:
        #     mu_cl_l = cluster_params['pm_l_poleski']
        #     mu_cl_b = cluster_params['pm_b_poleski']
        #     dist_cl  = cluster_params['distance_bj']
        #     radial_velocity_cl = cluster_params['RV']
        #     #shit name
        #     mu_l_total  = (mu_l - mu_cl_l) *u.mas/u.yr
        #     mu_b_total = (mu_b - mu_cl_b) * u.mas/u.yr
        #     dist_total  = dist 
        #     radial_velocity_total =( radial_velocity- radial_velocity_cl)*u.km/u.s
        #     print('wrt cluster')
        # else:
        mu_l_total  = mu_l 
        mu_b_total = mu_b 
        dist_total  = dist 
        radial_velocity_total = radial_velocity
            
        with coord.galactocentric_frame_defaults.set('v4.0'):
            galcen_frame = coord.Galactocentric()
        galactic_rep = coord.SkyCoord(l=l,b=b,pm_l_cosb=mu_l_total,pm_b=mu_b_total,distance=dist_total,
                                      radial_velocity =radial_velocity_total, frame='galactic')
        #transform frame
        star_galacto = galactic_rep.transform_to(galcen_frame)

        #correct for solar motion
        star_galacto = reflex_correct(star_galacto)
        

        initial_pos = gd.PhaseSpacePosition(star_galacto.data)
        total_time = self.int_time *u.Myr
        dt = -0.1 *u.Myr
        n_steps = int(abs(total_time.to_value(u.Myr) / dt.to_value(u.Myr)))
        integrator = gi.LeapfrogIntegrator
        potential = gp.MilkyWayPotential2022()  
        #orbit = potential.integrate_orbit(initial_pos, dt=dt, t1=0, t2=total_time)
        orbit = potential.integrate_orbit(initial_pos, 
                                          dt=dt,n_steps=n_steps,
                                          Integrator=integrator)
        #orbit = potential.integrate_orbit(initial_pos,  dt=-0.5 * u.Myr, t1=0, t2=-4 * u.Myr)


        return orbit 

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
                z_naught = data['distance_bj']*np.sin(np.radians(data['b']))*1000
                l_naught = data['l']
                b_naught = data['b']
                arrow_pml, arrow_pmb = data['pm_l_poleski'], data['pm_b_poleski']
                long_path, lat_path, z_path,ticks = self.trace_linear_path(ids, None,time_step=1000)
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
                arrow_pmz = data['distance_para']*np.sin(np.radians(arrow_pmb))*1000
                #two arrows for b or z
                #plt.quiver(l_naught,b_naught, arrow_pml,arrow_pmb,color='xkcd:grey',angles='xy',width=0.002)
                #causing issues if i dont convert to np array, something about masking in astropy
                l_naught, z_naught = np.array(l_naught), np.array(z_naught)
                plt.quiver(l_naught,z_naught, arrow_pml,arrow_pmz,color='xkcd:shit brown',angles='xy',width=0.002)
               
                for tick_l, tick_b in ticks:
                    tick_z = data['distance_bj']*np.sin(np.radians(tick_b))*1000
                    
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
    

    def plot_orbit_planes(self,star_x, star_y, star_z, cluster_x, cluster_y, cluster_z, star_label="Star", cluster_label="Cluster"):
            """
            Plot orbit projections of the star and cluster in (Y, X), (X, Z), and (Y, Z) planes.
        
            Parameters:
            - star_x, star_y, star_z: Star's Cartesian coordinates
            - cluster_x, cluster_y, cluster_z: Cluster's Cartesian coordinates
            - star_label: Label for the star (default: "Star")
            - cluster_label: Label for the cluster (default: "Cluster")
            """
        
            fig, axs = plt.subplots(1, 3, figsize=(18, 5))
        
            # (Y, X) Plane
            star_marker = "*"
            star_color= 'blue'
            star_size = 50
            cluster_marker='o'
            cluster_color = 'red'
            cluster_size = 10
            axs[0].scatter(star_y, star_x, s=star_size, label=star_label, alpha=1.0, color=star_color,marker=star_marker)
            axs[0].plot(star_y, star_x, alpha=1.0, color=star_color)
            axs[0].scatter(cluster_y, cluster_x, s=10, label=cluster_label, color=cluster_color, marker=cluster_marker)
            axs[0].plot(cluster_y, cluster_x, alpha=0.5, color=cluster_color)
            axs[0].set_xlabel('Y')
            axs[0].set_ylabel('X')
            axs[0].set_title('(Y, X) Plane')
            axs[0].legend()
        
            # (X, Z) Plane
            axs[1].scatter(star_x, star_z, s=star_size, label=star_label, alpha=0.8, color=star_color)
            axs[1].plot(star_x, star_z, alpha=0.5, color=star_color)
            axs[1].scatter(cluster_x, cluster_z, s=cluster_size, label=cluster_label, color=cluster_color, marker=cluster_marker)
            axs[1].plot(cluster_x, cluster_z, alpha=0.5, color=cluster_color)
            axs[1].set_xlabel('X')
            axs[1].set_ylabel('Z')
            axs[1].set_title('(X, Z) Plane')
            axs[1].legend()
        
            # (Y, Z) Plane
            axs[2].scatter(star_y, star_z, s=star_size, label=star_label, alpha=0.8, color=star_color)
            axs[2].plot(star_y, star_z, alpha=0.5, color=star_color)
            axs[2].scatter(cluster_y, cluster_z, s=cluster_size, label=cluster_label, color=cluster_color, marker=cluster_marker)
            axs[2].plot(cluster_y, cluster_z, alpha=0.5, color=cluster_color)
            axs[2].set_xlabel('Y')
            axs[2].set_ylabel('Z')
            axs[2].set_title('(Y, Z) Plane')
            axs[2].legend()
        
            plt.tight_layout()
            plt.show()
            return None
    def plot_comoving_cluster(self,star_params,cluster_params,plotting=False):
        
        
        star_orbit = self.trace_galactic_path(star_params)
        
        cluster_orbit = self.trace_galactic_path(cluster_params)
        
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
        
        #offset the cluster
        # cluster_x_shifted = cluster_x - cluster_x[0]
        # cluster_y_shifted = cluster_y - cluster_y[0]
        # cluster_z_shifted = cluster_z - cluster_z[0]
        
        #offset the star 
        star_x_shifted = star_x - cluster_x
        star_y_shifted = star_y - cluster_y
        star_z_shifted = star_z - cluster_z
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
            plt.scatter(star_y_shifted,star_z_shifted,label=f"{star_params['Name'][0]}",s=100,alpha=1.0)
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
            plt.ylabel('Z')
            plt.xlabel('Y')
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
            ax.set_title('no offset')
            ax.set_xlabel('X ')
            ax.set_ylabel('Y ')
            ax.set_zlabel('Z')
            ax.legend()
            
            fig_3d = plt.figure()
            # ax = fig_3d.add_subplot(projection='3d')
            # #star
            # ax.set_title('Cluster shifted')
            # ax.scatter(star_x_shifted, star_y_shifted, star_z_shifted, marker='*',label=f"{star_params['Name'][0]}")
            
            # ax.scatter(cluster_x, cluster_y, cluster_z, marker='o',label='cluster')
            # #ax.scatter(rel_x,rel_y,rel_z,label='relative_motion')
            # ax.set_xlabel('X')
            # ax.set_ylabel('Y')
            # ax.set_zlabel('Z')
            # ax.legend()
            
        return rel_x,rel_y,rel_z, time_min_sep #star_x_shifted ,star_y_shifted, star_z_shifted
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
        long_path, lat_path, _, ticks = self.trace_linear_path(source_id,None)
        #convert to numpy arrays
        long_path, lat_path, ticks = np.array(long_path), np.array(lat_path),np.array(ticks)
        
        long_vdm_path, lat_vdm_path, _  = self.trace_linear_path_vdm(single_star)
        
        
        
        #inital positions
        
        l_naught = single_star['l']
        b_naught = single_star['b']
        
        #proper motion vectors
        # arrow_pml= single_star['pm_l_poleski'] - cluster_params['pm_l_poleski']        
        # arrow_pmb  = single_star['pm_b_poleski']   - cluster_params['pm_b_poleski']
        #member clusters
        if clustertable is not None:
            member_long = clustertable['l']
            member_lat = clustertable['b']
        
        #include gala
        
        orbit  = self.trace_galactic_path(single_star)
       #print(orbit)
        
        x1 = np.array(orbit.x)
        y1 = np.array(orbit.y)
        z1 = np.array(orbit.z)
        
        xyz_galctocentric = coord.SkyCoord(x=x1*u.kpc,y=y1*u.kpc,z=z1*u.kpc,frame=coord.Galactocentric())
        xyz_galactic = xyz_galctocentric.transform_to(coord.Galactic())
        
        plt.figure(figsize=(10,5))
        if clustertable is not None:
            plt.scatter(member_long,member_lat, s=50,marker='*',label=f'{clustername}',color='xkcd:grey')
        
        #plot the orbit integrated path of the cluster
        if cluster_params is not None:
            cluster_orbit = self.trace_galactic_path(cluster_params)
            cluster_x1  = np.array(cluster_orbit.x)
            cluster_y1 = np.array(cluster_orbit.y)
            cluster_z1 = np.array(cluster_orbit.z)
            
            cluster_galactocentric  = coord.SkyCoord(x=cluster_x1*u.kpc, y=cluster_y1*u.kpc,z=cluster_z1*u.kpc, frame = coord.Galactocentric())
            
            cluster_galactic = cluster_galactocentric.transform_to(coord.Galactic())
            
            cluster_l_vals = cluster_galactic.l.deg
            cluster_b_vals = cluster_galactic.b.deg
            #huuh?
            rel_x, rel_y,rel_z, time_min_sep  = self.plot_comoving_cluster(single_star, cluster_params, plotting=False)
            print(f'Minimum Seperation time {time_min_sep}')
            relative_sky = coord.SkyCoord(x=rel_x, y=rel_y, z=rel_z, unit=u.kpc,
                            representation_type='cartesian', frame='galactocentric').transform_to('galactic')
            plt.scatter(cluster_l_vals,cluster_b_vals,label='Cluster Path',color='xkcd:dark grey')
            
            plt.plot(cluster_l_vals,cluster_b_vals,color='xkcd:grey')
            arrow_pml= single_star['pm_l_poleski'] - cluster_params['pm_l_poleski']        
            arrow_pmb  = single_star['pm_b_poleski']   - cluster_params['pm_b_poleski']
            comove_l = relative_sky.l.deg % 360
            comove_b =relative_sky.b.deg
        else:
            arrow_pml= single_star['pm_l_poleski'] 
            arrow_pmb  = single_star['pm_b_poleski']
            


        gala_l_vals = xyz_galactic.l.deg # Galactic longitude in degrees
        gala_b_vals = xyz_galactic.b.deg

        N = len(long_path)
        #print(abs(ttotal))
        star_name = str(single_star['Name'].value[0])
        #plot path and color plot by specrtral type
        sp_type = single_star['Mod_SpType'][0] #dumb
        path_color = sp_type if sp_type in self.color_map else 'xkcd:grey'
        #linear path
        plt.scatter(long_path[0:N], lat_path[0:N],s=10,color='xkcd:black',alpha=0.7)
       # plt.scatter(long_vdm_path, lat_vdm_path,s=50,color='xkcd:pink',alpha=0.7,label='vdm')
        #plt.plot(long_vdm_path, lat_vdm_path,color='xkcd:pink')

        plt.plot([],[],color='xkcd:black',label='Linear Path')
        #inistal position
        plt.scatter(l_naught,b_naught,color=path_color,label=star_name)
        #gala
        plt.plot(gala_l_vals,gala_b_vals, lw=1.1, marker='*', color=path_color,label='Gala path')
        #pm arrows
        # i need the angle between the proper motions but later  i invert the x-axis so i should flip the 
        #flip the pm_l arrow
        plt.quiver(l_naught,b_naught, -arrow_pml,arrow_pmb,color='xkcd:shit brown',angles='uv',width=0.002)
        
        #plot cluster integrated path

        #tick markers
        for tick_l, tick_b in ticks:
            plt.scatter(tick_l,tick_b,color='xkcd:canary',s=20,marker='x')
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
        return None

test_table = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/HMXB_20250301_.ecsv',format='ecsv')
test_170037 = test_table[test_table['Name']=='4U 1700-377']
scoob1 = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/SCO OB1-result.ecsv')
ngc6231_params = ascii.read('/home/karan/Documents/UvA/Thesis/DATA/NGC2631_param.ecsv')

mydir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(mydir)
if __name__ == "__main__":
   # GalacticTraceback(test_170037,-3.0).plot_trace(savefig=False)
   GalacticTraceback(test_170037,-3.0).plot_with_cluster(clustername='NGC 6231',cluster_params=ngc6231_params, clustertable=scoob1,savefig=False)
   x,y,z,t_min_sep = GalacticTraceback(test_170037,int_time=-3.0).plot_comoving_cluster(test_170037, ngc6231_params,plotting=True)
#  GalacticTraceback(test_170037, -3.0).trace_linear_path_vdm(test_170037)
#    GalacticTraceback(test_table).plot_trace(savefig=False)
    
#ngc6231_params
#test_170037['pmra']
