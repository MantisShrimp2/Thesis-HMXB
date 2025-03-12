#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 19:51:02 2025

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
class galactic_path_MCMC:
    def __init__(self,int_time):
        self.int_time = int_time
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
        
        
        orbit_integrator = gi.LeapfrogIntegrator
        potential = gp.MilkyWayPotential2022()  
        orbit_params = {"dt": dt, "n_steps": n_steps, "Integrator": orbit_integrator}
        #orbit = potential.integrate_orbit(initial_pos, dt=dt, t1=0, t2=total_time)
        orbit = potential.integrate_orbit(initial_pos, **orbit_params)


        return orbit 