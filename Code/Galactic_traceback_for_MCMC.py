#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 19:51:02 2025

@author: Karan Kumar
"""
import numpy as np
from datetime import datetime
import emcee
import gala.potential as gp
import gala.dynamics as gd
import gala.integrate as gi
import astropy.coordinates as coord
import os
import astropy.units as u
from multiprocessing import Pool 
import h5py
import matplotlib.pyplot as plt
from astropy.table import Table
class galactic_path_MCMC:
    def __init__(self):
        self.frac_of_light = 0.1
        #10% of light speed  in mas/yr at 1kpc 
        self.speed_of_light  = self.frac_of_light*(3e6/(4.74 * 1.0))

    def make_theta(self, star_params):
        '''
        Make the theta object for mcme as the inital conditions of the star
        in galactic coordinates
        Returns
        -------
        theta : List, for emcee.

        '''
        #star_l   = np.array(star_params['l'])[0]
        #star_b   = np.array(star_params['b'])[0]
        # star_d   = np.array(star_params['distance_bj'])[0]
        # star_pml = np.array(star_params['pm_l_poleski'])[0]
        # star_pmb = np.array(star_params['pm_b_poleski'])[0]
        # star_rv  = np.array(star_params["RV"])[0]
        
        star_d   = np.array(star_params['distance_bj'])[0]
        # star_d   = (1/star_params['parallax'])
        star_pmra = np.array(star_params['pmra'])[0]
        star_pmdec = np.array(star_params['pmdec'])[0]
        star_rv  = np.array(star_params["RV"])[0]

        #theta = [star_l, star_b, star_d, star_pml, star_pmb, star_rv]
       # theta = [star_d, star_pml, star_pmb, star_rv]
        theta = [star_d, star_pmra, star_pmdec, star_rv]
        return theta
    def make_star_stds(self, star_params):
        '''
        Compute 1-simga error on each theta paramter in galactic units

        Returns
        -------
        star_stds : list.

        '''

        dist = np.array(star_params['distance_bj'])[0]
        dist_up = np.array(star_params['distance_bj_high'])[0] - dist
        dist_down = dist - np.array(star_params['distance_bj_low'][0])
        dist_err = (dist_up + dist_down)/2.0
   
        
        
        std_d   = dist_err
        std_pmra = np.array(star_params['pmra_error'])[0]
        std_pmdec = np.array(star_params['pmdec_error'])[0]
        std_rv  = np.array(star_params['RV_err'])[0]
    
   
        star_stds  = [std_d, std_pmra, std_pmdec, std_rv]
        return star_stds
    def calc_cluster_radius(self,cluster_distance,angular_diameter):
        '''
        radius of a cluster from angular diameter

        Parameters
        ----------
        cluster_distance : float in kpc or pc.
        angular_diameter : float radians.

        Returns
        -------
        radius of cluster in pc or kpc (distance dependent).

        '''
        rad_diameter = angular_diameter* (np.pi/ (60*180)) # diameter of cluster in radians
        return cluster_distance.value*np.tan(rad_diameter/2)
    
    def make_log_gauss(self,x,mu,sigma):
        '''
        gaussian in log space for prior

        Parameters
        ----------
        x : walker guess for parameter.
        mu : expected value from GAIA/ HMXB table.
        sigma : 1-sigma deviation for paramter.

        Returns
        -------
        log probability of prior.

        '''
        return -0.5 * ((x - mu)/sigma)**2
    def make_log_uniform(self,low, high):
        '''
        logarithim uniform distribution

        Parameters
        ----------
        low : paramter lower bound.
        high : paramter upper .

        Returns
        -------
        log of a uniform distribution.

        '''
        return np.log(np.random.uniform(low=low, high = high))
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
  

        ra = row['ra']
        dec = row['dec']
        pmra = row['pmra']
        pmdec = row['pmdec']
        radial_velocity = row['RV'] 
        distance_icrs = row['distance_bj'] #(1/row['parallax']).value *u.kpc
       # k = 4.74 * (u.km/u.s)/(u.mas *u.kpc/u.yr) #km/s per mas/yr 

            
        with coord.galactocentric_frame_defaults.set('v4.0'):
            galcen_frame = coord.Galactocentric()

        galactic_rep = coord.SkyCoord(ra=ra, dec=dec, distance=distance_icrs,
              pm_ra_cosdec=pmra, pm_dec=pmdec,
              radial_velocity=radial_velocity)
        #transform frame
        star_galacto = galactic_rep.transform_to(galcen_frame)

        initial_pos = gd.PhaseSpacePosition(star_galacto.data)
        total_time = int_time *u.Myr
        dt = -0.1 *u.Myr
        n_steps = int(abs(total_time.to_value(u.Myr) / dt.to_value(u.Myr)))
        
        
        orbit_integrator = gi.LeapfrogIntegrator
        potential = gp.MilkyWayPotential2022()  
        orbit_params = {"dt": dt, "n_steps": n_steps, "Integrator": orbit_integrator}
        orbit = potential.integrate_orbit(initial_pos, **orbit_params)


        return orbit 
    def log_likelihood(self,theta, cluster_params, cluster_radius,int_time,star_ra,star_dec):
        '''theta  =['l','b','dist','pml','pmb','rv']
        log likelihood is the time separation of the star and cluster after orbit integration
        
        calculate the orbit of the star and cluster in galactocentric
        Very stupid, MCMC shouldn't accept theta with units so i have to remove the units and then put them back into the orbit integration
    
        Also very stupid/ annoying adding units back into star_table for integration some reason doing the integration with Gala in the function
        caused troubles. So i used my GalacticTrackback.trace_galactic_path function and that works
        
        Returns a kinematic age
        '''
        
        #star_d = theta[0]
        star_d, star_pmra, star_pmdec, star_rv = theta[0:4]
    
        
        # star_table = Table([[star_d]*u.kpc, star_l, star_b, [star_pml]*u.mas/u.yr, [star_pmb]*u.mas/u.yr, [star_rv]*u.km/u.s],
        #             names=["distance_bj", "l", "b", "pm_l_poleski", "pm_b_poleski", "RV"])
        
        star_table = Table([[star_d]*u.kpc, star_ra, star_dec, [star_pmra]*u.mas/u.yr, [star_pmdec]*u.mas/u.yr, [star_rv]*u.km/u.s],
                    names=["distance_bj", "ra", "dec", "pmra", "pmdec", "RV"])
    
    

    
        star_orbit = self.trace_galactic_path(star_table, int_time)
        cluster_orbit = self.trace_galactic_path(cluster_params,int_time)
    
        time = star_orbit.t
    
    
        separation = np.linalg.norm(star_orbit.xyz - cluster_orbit.xyz, axis=0)
        
        min_sep = np.argmin(separation)
        
        min_sep_radius = float(separation[min_sep].to(u.kpc).value)
        
        #to maximize liklehood
        return -float(time[min_sep].value), min_sep_radius
    def log_normal_prior(self,theta,star_params):
        '''Assume each parameter of the star comes from a normal distriubtion
        Calcualte the prior for each parameter in logspace
        
        cannot do self.star_params or self. theta because log probability needs those are arguements'''
        #these are sampled from walker
        star_d,star_pmra, star_pmdec, star_rv = theta[0:4]

        true_d= star_params['distance_bj']
        true_pmra, true_pmdec, true_rv = star_params['pmra'], star_params['pmdec'], star_params["RV"]
    
        #true standard deviations 
        
        star_stds = self.make_star_stds(star_params) 
    
        std_d, std_pmra,std_pmdec,std_rv = star_stds[0:4]
        
    

        if star_d < 0.0 or np.sqrt(star_pmra**2 + star_pmdec**2) >self.speed_of_light: 
            return -np.inf


        log_d = self.make_log_gauss(star_d, true_d, std_d).value

        log_pmra = self.make_log_gauss(star_pmra,true_pmra, std_pmra).value
        log_pmdec = self.make_log_gauss(star_pmdec, true_pmdec, std_pmdec).value
        log_rv = self.make_log_gauss(star_rv, true_rv, std_rv).value
        #log_d + log_pml + log_pmb + log_rv
        return log_d + log_pmra + log_pmdec + log_rv
    
    def log_probability(self,theta,star_params, cluster_params,cluster_radius,int_time):
        '''Calculate the posterior and the result of the log likliehood'''
        lp = self.log_normal_prior(theta,star_params)
        if not np.all(np.isfinite(lp)):
            return -np.inf, np.nan, np.nan

        kinematic_age, min_sep = self.log_likelihood(theta, cluster_params, cluster_radius,int_time,star_params['ra'], star_params['dec'])
        return lp + kinematic_age, np.array([kinematic_age, min_sep], dtype=np.float64)
    
    def start_mcmc(self,backend_name,star_params,cluster_params, cluster_radius,int_time,names):
        '''
        

        Parameters
        ----------
        backend_name : to store MCMC simulation without rerunning it.
        star_params : l,b,distance, pml, pmb rv in that order for the star.
        cluster_params : same as star_params only for cluster.
        cluster_radius : float in pc.
        int_time : integration time for orbit (-) for backwards in time. (+) for forwards.
        names : list, name of star and name of cluster.

        Returns
        -------
        sampler, state.

        '''

        star_name, cluster_name = names[0],names[1]
        star_stds = self.make_star_stds(star_params)
        theta = self.make_theta(star_params)
        
        args = (star_params,cluster_params,cluster_radius,int_time)
        backend = emcee.backends.HDFBackend(backend_name+'.h5')
        nwalkers = 24
        ndim = len(theta) #6 if l and b change
        nsample = 1000
        initial_pos = np.hstack([
        np.random.normal(theta[0:ndim], star_stds[0:ndim], size=(nwalkers, ndim)),  # Sample the remaining 4 parameters
    ])

        backend.reset(nwalkers, ndim)
        print(f'Compling MCMC for {star_name} and {cluster_name}')
        with Pool() as pool:
                sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_probability,args = args,
                                            pool=pool, backend=backend)
                state = sampler.run_mcmc(initial_pos,nsample,progress=True)

        return sampler, state
       
        
    def compute_full_separations_and_plot(self,flat_chain, cluster_params, int_time, star_ra, star_dec,cluster_radius, title, cluster_age=None, n_samples=100,savefig=False,ylim=1000):
        '''
        

        Parameters
        ----------
       flat_chain : ndarray
        MCMC samples, shape (n_samples_total, 4) [distance, pmra, pmdec, rv].
    cluster_params : astropy Table
        Cluster parameters for orbit integration.
    int_time : float
        Integration time (Myr).
    star_ra, star_dec : float
        Fixed ra and dec of the star (deg).
    n_samples : int, optional
        Number of random MCMC samples to use for separation computation.
    cluster_radius : float or None, optional
        Radius of the cluster in parsecs to plot on the graph.
        -------
        Separtion over time plot
        with age of cluster
        radius of cluster
        1 sigma, 2 sigma confidence intervals.

        '''
        
        total_samples = len(flat_chain)
        n_samples = min(n_samples, total_samples)
        
        indices = np.random.choice(total_samples, size=n_samples, replace=False)
        chosen_samples = flat_chain[indices]
        
        
        separation_curves = []
        for theta in chosen_samples:
            star_d, pmra, pmdec, rv = theta
            #star_d = 1/star_parallax.value * u.kpc
            star_table = Table(
    [
        [star_d] * u.kpc,
        star_ra,
        star_dec,
        [pmra] * u.mas / u.yr,
        [pmdec] * u.mas / u.yr,
        [rv] * u.km / u.s
    ],
    names=["distance_bj", "ra", "dec", "pmra", "pmdec", "RV"]
)

            star_orbit = self.trace_galactic_path(star_table, int_time)
            cluster_orbit = self.trace_galactic_path(cluster_params, int_time)
        
            time_array = star_orbit.t.to_value(u.Myr)
            separation = np.linalg.norm(star_orbit.xyz - cluster_orbit.xyz, axis=0)
            separation_curves.append(separation*1e3)
        
        separation_curves = np.array(separation_curves)
        
        # Compute percentiles for confidence interval
        p16 = np.percentile(separation_curves, 16, axis=0)
        p50 = np.percentile(separation_curves, 50, axis=0)
        p84 = np.percentile(separation_curves, 84, axis=0)
        
        p2 = np.percentile(separation_curves,2.5,axis=0)
        p98 = np.percentile(separation_curves, 97.5, axis=0)
        
        where_min_index = np.argmin(p50) -1 
        time_where_min = time_array[where_min_index]
        sep_where_min = np.min(p50[where_min_index])
        
        # Plot all individual samples (faint)
        plt.figure(figsize=(10,6))
        if cluster_age[0]:
            age_median = cluster_age[0] 
            age_std =cluster_age[1]
            plt.vlines(x=-age_median,ymin=0,ymax=5000,color='xkcd:grey', 
                       linestyle='--',linewidth=3.0,label=f'Cluster age {age_median:.2f} $\pm$ {age_std:.2f} Myr')
            if age_std:
                age_min = -age_median - age_std
                age_max = -age_median + age_std
                print(age_min, age_max)
                plt.axvspan(age_min, age_max,  color='xkcd:grey',alpha=0.4)
        # for sep in separation_curves[0:100]:
        #     plt.plot(time_array, sep, color='gray', alpha=0.1)
            
        # Plot median and confidence interval
        plt.plot(time_array, p50, color='xkcd:purple',lw=2,label='Median')
        #'xkcd:neon green' xkcd:orangee
        plt.fill_between(time_array, p16, p84, color='xkcd:neon green', alpha=0.5, label='68% confidence interval')
        plt.fill_between(time_array, p2, p98, color='xkcd:orange', alpha=0.5, label='95% confidence interval',zorder=-1)
        
        where_color = 'xkcd:black'
        plt.scatter(time_where_min,sep_where_min,color=where_color, 
                    label=f'Min {time_where_min:.2f} Myr, {sep_where_min:.2f} pc')
        plt.hlines(sep_where_min, xmin=time_array[0], xmax=time_array[-1], ls=':', color=where_color,zorder=2)
        plt.axvline(time_where_min, ymax=max(separation_curves[where_min_index]),  ls=':', color=where_color)
        plt.xlim(time_array[0],time_array[-1])
        plt.ylim(0,ylim) #pc
        
        #plot radius interval



        # Cluster radius reference line
        if cluster_radius is not None:
            cl_radi = np.max(cluster_radius)*1000
            cl_radi_std = np.std(cluster_radius)*1000
            radi_upper = cl_radi+ cl_radi_std
            radi_lower = cl_radi - cl_radi_std
            plt.axhline(cl_radi, color='xkcd:hot pink', ls='--', label=f'Cluster radius = {cl_radi:.2f} pc')
            plt.fill_between(time_array, radi_lower,radi_upper,color='xkcd:pink',alpha=0.5)

        if title:
            plt.title(title)
        else:
            plt.title('Separation of Star and Cluster Over Time')
        
        plt.xlabel('Time (Myr)')
        plt.ylabel('Separation (pc)')
  
        plt.gca().invert_xaxis()
        plt.legend()
        plt.tight_layout()
        if savefig ==True:
            mydir = os.path.dirname(os.path.realpath(__file__))
            parentdir = os.path.dirname(mydir)
            today = datetime.now().strftime("%Y%m%d")
            save_path = os.path.join(parentdir, 'Figures', 'MCMC',
                                     f"MCMCsep_{title}_{today}.png")
            plt.savefig(save_path,dpi=300)
        plt.show()
        return None



  


           
              
           