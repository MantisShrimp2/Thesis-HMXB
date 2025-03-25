#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 19:51:02 2025

@author: karan
"""
import numpy as np
import pandas as pd
from astroquery.gaia import Gaia
from astropy.io import ascii
import emcee
import gala.potential as gp
import gala.dynamics as gd
import gala.integrate as gi
from gala.coordinates import reflex_correct
import astropy.coordinates as coord
import os
import astropy.units as u
from multiprocessing import Pool 
import h5py
import matplotlib.pyplot as plt
from astropy.table import Table
class galactic_path_MCMC:
    def __init__(self,int_time):
        self.int_time = int_time
    #     self.frac_of_light = 0.1
    #     #10% of light speed  in mas/yr at 1kpc 
    #     self.speed_of_light  = self.frac_of_light*(3e6/(4.74 * 1.0))
    #     self.star_params = star_params
    #     self.cluster_params = cluster_params
    #     self.cluster_radius = cluster_radius
    #     self.star_name= star_name
    #     self.cluster_name = cluster_name
    # def make_theta(self):
    #     '''
    #     Make the theta object for mcme as the inital conditions of the star
    #     in galactic coordinates
    #     Returns
    #     -------
    #     theta : List, for emcee.

    #     '''
    #     star_l   = np.array(self.star_params['l'])[0]
    #     star_b   = np.array(self.star_params['b'])[0]
    #     star_d   = np.array(self.star_params['distance_bj'])[0]
    #     star_pml = np.array(self.star_params['pm_l_poleski'])[0]
    #     star_pmb = np.array(self.star_params['pm_b_poleski'])[0]
    #     star_rv  = np.array(self.star_params["RV"])[0]

    #     theta = [star_l, star_b, star_d, star_pml, star_pmb, star_rv]
    #     return theta
    # def make_star_stds(self):
    #     '''
    #     Compute 1-simga error on each theta paramter in galactic units

    #     Returns
    #     -------
    #     star_stds : list.

    #     '''
    #     dist_err = np.array(self.star_params['distance_bj_high'])[0] - np.array(self.star_params['distance_bj'])[0]
    #     std_l   = np.array(self.star_params['l_err'])[0]
    #     std_b   = np.array(self.star_params['b_err'])[0]
    #     std_d   = dist_err
    #     std_pml = np.array(self.star_params['pm_l_err'])[0]
    #     std_pmb = np.array(self.star_params['pm_b_err'])[0]
    #     std_rv  = np.array(self.star_params['RV_err'])[0]
    
    #     star_stds = [std_l, std_b, std_d, std_pml, std_pmb, std_rv]
    #     return star_stds
    # def calc_cluster_radius(self,cluster_distance,angular_diameter):
    #     '''
    #     radius of a cluster from angular diameter

    #     Parameters
    #     ----------
    #     cluster_distance : TYPE
    #         DESCRIPTION.
    #     angular_diameter : TYPE
    #         DESCRIPTION.

    #     Returns
    #     -------
    #     TYPE
    #         DESCRIPTION.

    #     '''
    #     rad_diameter = angular_diameter* (np.pi/ (60*180)) # diamter of cluster in radians
    #     return cluster_distance.value*np.tan(rad_diameter/2)
    
    # def make_log_gauss(self,x,mu,sigma):
    #     '''
    #     gaussian in log space for prior

    #     Parameters
    #     ----------
    #     x : walker guess for parameter.
    #     mu : expected value from GAIA/ HMXB table.
    #     sigma : 1-sigma deviation for paramter.

    #     Returns
    #     -------
    #     TYPE
    #         DESCRIPTION.

    #     '''
    #     return -0.5 * ((x - mu)/sigma)**2
    # def make_log_uniform(self,low, high):
    #     '''
    #     logarithim uniform distribution

    #     Parameters
    #     ----------
    #     low : paramter lower bound.
    #     high : paramter upper .

    #     Returns
    #     -------
    #     log of a uniform distribution.

    #     '''
    #     return np.log(np.random.uniform(low=low, high = high))
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
        #star_galacto = reflex_correct(star_galacto)
        

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
    # def log_likelihood(self,theta, cluster_params, cluster_radius):
    #     '''theta  =['l','b','dist','pml','pmb','rv']
    #     log likelihood is the time seperation of the star and cluster after orbit integration
        
    #     calculate the orbit of the star and cluster in galactocentric
    #     Very stupid, MCMC shouldn't accept theta with units so i have to remove the units and then put them back into the orbit integration
    
    #     Also very stupid/ annoying adding units back into star_table for integration some reason doing the integration with Gala in the function
    #     caused troubles. So i used my GalacticTrackback.trace_galactic_path function and that works
        
    #     Returns a kinematic age
    #     '''
    #     star_l, star_b,star_d = theta[:3]
    #     star_pml, star_pmb, star_rv = theta[3:6]
    
        
    #     star_table = Table([[star_d]*u.kpc, [star_l]*u.deg, [star_b]*u.deg, [star_pml]*u.mas/u.yr, [star_pmb]*u.mas/u.yr, [star_rv]*u.km/u.s],
    #                 names=["distance_bj", "l", "b", "pm_l_poleski", "pm_b_poleski", "RV"])
    
    #     gp_mcmc = self.trace_galactic_path(star_table)
    
    #     star_orbit = gp_mcmc.trace_galactic_path(star_table)
    #     cluster_orbit = gp_mcmc.trace_galactic_path(cluster_params)
    
    #     time = star_orbit.t
    
    
    #     seperation = np.linalg.norm(star_orbit.xyz - cluster_orbit.xyz, axis=0)
    #     min_sep = np.argmin(seperation)
    
    #     #gaussian penalty
    #     sigma = cluster_radius/2.0
    #     # if seperation[min_sep] > cluster_radius:
    #     #     return -np.inf
    #     if seperation[min_sep].value > cluster_radius:
    #         penalty = -0.5 * ((seperation[min_sep].value - cluster_radius) / sigma)**2
    #     else:
    #         penalty = 0.0
    #     #to maximize liklehood
    #     return -float(time[min_sep].value)
    # def log_normal_prior(self,theta,star_params,star_stds):
    #     '''Assume each parameter of the star comes from a normal distriubtion
    #     Calcualte the prior for each parameter in logspace
        
    #     cannot do self.star_params or self. theta because log probability needs those are arguements'''
    #     #these are sampled from walker
    #     star_l, star_b, star_d = theta[:3]
    #     star_pml, star_pmb, star_rv = theta[3:6]
    
    #     #true value of stars
    #     true_l, true_b, true_d= star_params['l'], star_params['b'], star_params['distance_bj']
    #     true_pml, true_pmb, true_rv = star_params['pm_l_poleski'], star_params['pm_b_poleski'], star_params["RV"]
    
    #     #true standard deviations 
    #     #star_stds = [std_l, std_b, std_d, std_pml, std_pmb, std_rv]
    #     star_stds = self.make_star_stds()
    #     dist_err = star_stds[2] 
    #     std_l, std_b, std_d = star_stds[0], star_stds[1], dist_err
    #     std_pml, std_pmb, std_rv = star_stds[3], star_stds[4], star_stds[5]
    
    #     if star_b > 90.0 or star_b < -90:
    #         return -np.inf
    #     if star_d < 0.0 or np.sqrt(star_pml**2 + star_pmb**2) >self.speed_of_light: 
    #         return -np.inf
    #     #make_log_uniform(true_l - std_l, true_l + std_l
    #     #make_log_uniform(true_b - std_b, true_b + std_b)
    #     log_l = self.make_log_gauss(star_l,true_l,std_l).value
    #     log_b = self.make_log_gauss(star_b,true_b,std_b).value
    #     #print(log_l, log_b)
    #     log_d = self.make_log_gauss(star_d, true_d, std_d).value
        
    #     log_pml = self.make_log_gauss(star_pml,true_pml, std_pml).value
    #     log_pmb = self.make_log_gauss(star_pmb, true_pmb, std_pmb).value
    #     log_rv = self.make_log_gauss(star_rv, true_rv, std_rv).value
    #     return log_l + log_b + log_d + log_pml + log_pmb + log_rv
    
    # def log_probability(self,theta,star_params, star_stds, cluster_params,cluster_radius):
    #     '''Calculate the posterior and the result of the log likliehood'''
    #     lp = self.log_normal_prior(theta,star_params,star_stds)
    #     if not np.all(np.isfinite(lp)):
    #         return -np.inf, np.nan
    #     kinematic_age = self.log_likelihood(theta, cluster_params, cluster_radius)
    #     return lp + kinematic_age, kinematic_age
    
    # def start_mcmc(self,backend_name, cluster_radius):
    #     '''
    #     run the MCMC 
    #     cluster radius will either be from angular diameter or a given

    #     Parameters
    #     ----------
    #     backend_name : TYPE
    #         DESCRIPTION.
    #     cluster_radius : TYPE
    #         DESCRIPTION.

    #     Returns
    #     -------
    #     None.

    #     '''
    #     star = self.star_params
    #     star_stds = self.make_star_stds()
    #     cluster = self.cluster_params
    #     theta = self.make_theta()
        
    #     args = (star,star_stds,cluster,cluster_radius)
    #     backend = emcee.backends.HDFBackend(backend_name+'.h5')
    #     nwalkers = 24
    #     ndim = 6
    #     nsample = 3000
    #     initial_pos = np.hstack([
    # np.random.normal(theta[0:3], star_stds[0:3], size=(nwalkers, 3)),
    # np.random.normal(theta[3:6], star_stds[3:6], size=(nwalkers, 3)),])  
    #     print(f'Compling MCMC for {self.star_name} and {self.cluster_name} now')
    #     with Pool() as pool:
    #             sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_probability,args = args,
    #                                         pool=pool, backend=backend)
    #             state = sampler.run_mcmc(initial_pos,nsample,progress=True)

    #     return sampler, state
       
           
              
           