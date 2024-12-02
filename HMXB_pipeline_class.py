#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 14:39:14 2024

@author: Karan Kumar 
14906619
"""

import numpy as np
import pandas as pd
from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
from astropy.table import Column, Table
from astropy.io import ascii

import astropy.coordinates as coords
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import astropy_mpl_style
import os 
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
class pipeline:
    def __init__(self, table_name,fmt):
     # define constant
     self.RUWE = 1.4
    # self.list_ID = list_ID
     self.R0= 8.5 # solar galactic distance to center
     self.sun_curve = 220 #km/s galactic velocity of the sun
     self.ra_np = np.radians(192.25) # deg to radians right acsension of north pole
     self.dec_np = np.radians(27.4) # deg to radians decliatoin of north pole
     self.theta_o = np.radians(123) # deg to radiansradians
     self.k  = 4.74 #km/s per mas/yr 
      
     #solar values
      # solar motion km/s 
    # M. Carretero-Castrillo 2023 and Ried 2019
     self.U_sun = 10.8
     self.V_sun = 13.6
     self.W_sun = 7.6
     self.table_name = table_name #string
     self.fmt = fmt
    def make_table(self):
        table = ascii.read(self.table_name,format=self.fmt)
        return table
    def modify_parllax_add_distance(self,table):
        '''Calculate the distance as 1/parallax, rememeber its an estimator
        add parallax error to negative parallax to get lower limit on distance
        some parallax still may be negative'''
        #offset Parallax
        p_offset = 0.00 # from GAIA EDR3
        table['parallax'] = (table['parallax'] - p_offset)
        table['parallax'].unit = u.mas# from GAIA EDR3
        #add parallax units
        #adjust the negative values
        negative_prlx_mask = table['parallax'] <=0 
        table['parallax'][negative_prlx_mask] += table['parallax_error'][negative_prlx_mask]
        #add the distance to the table
        table['distance']= (1/table['parallax'])*u.kpc # kpc
        
        # distance from object to sun in kpc
        table['distance'].unit = u.kpc
        return table
    def Galaxy_dist(self,table):
        '''find the distance an object is to the center of the galaxy based on law of cosine
        input:
        long (l)- galactic longitutde in degree
        lat(b)- galactic latitute in degrees 
        object_dist (d) - distance from the sun to the object in kpc
        return 
        R - distance from the object to the galactic centre in kpc
        '''
        long = table['l']
        lat = table['b']
        long_rad = np.radians(long) #convert to radians
        lat_rad = np.radians(lat)
        obj_dist = table['distance'] # kpc
        R_sqrd = self.R0**2 + (obj_dist**2 * np.cos(lat_rad)**2) - 2*(self.R0*obj_dist*np.cos(long_rad)*np.cos(lat_rad))
        galactic_dist = np.sqrt(R_sqrd)
       # table.add_column(galactic_dist, name = 'galactic distance')
        table['galactic distance'] = galactic_dist*u.kpc
        table['galactic distance'].unit = u.kpc
        return table
    def rotation_curve(self,table):
        '''From brand 1993 and fich 1988
        Calculate the circular velocity of a star based on its galactocentric distance
        return in km/s'''
        #best fit constants for rotational curve  fit
        a1 = 1.00767
        a2 = 0.0394
        a3 = 0.00712
        # fit from Brand 1988
        gal_dist = table['galactic distance']
        theta  = a1*((gal_dist/self.R0)**a2) + a3
        theta = self.sun_curve*theta #km/s
        #table.add_column(theta, name='circular velocity')
        table['circular velocity']= theta* u.km/u.s
        table['circular velocity'].unit = u.km/u.s
        return table
    def vlsr_model(self,table):
        '''Calculate the Velocity as Local Standard of Rest based on 
        1) Galactocentric distance
        2) Circular Velocity due to galactic rotation (see rotation_curve function)
        3) Galactic Latitude and Longitude
        
        Compare with a model from Brand 1993 which reduces solar motion
        input:
        table- astropy table
        plot- boolean'''
        long = table['l']
        lat = table['b']
        gal_dist = table['galactic distance']
        theta = table['circular velocity']
        
        long_rad = np.radians(long)
        lat_rad = np.radians(lat)
        V_lsr = (theta * (self.R0/gal_dist) - self.theta_o)*np.sin(long_rad)*np.cos(lat_rad)
        #table.add_column(V_lsr, name='LSR velocity')
        table['LSR velocity']=  V_lsr
        table['LSR velocity'].unit = u.km/u.s
        return table
    def make_sky(self,table):
        '''Create a galactic coordinate table from astropy table proper motion, ra,dec and distance to object'''
        ra = table['ra'] # deg
        dec= table['dec']#deg
        pmra= table['pmra'] #mas/yr
        pmdec = table['pmdec']
       # distance = table['distance']
        #skycoord doesn't like ra = ra*u.deg, predefine the units somewhere else
        sky = SkyCoord(ra = ra, dec =dec, pm_ra_cosdec= pmra, pm_dec=pmdec, frame='icrs')
        galactic= sky.transform_to('galactic')
        pm_l = galactic.pm_l_cosb
        pm_b = galactic.pm_b
        table['pm_l'] = pm_l
        table['pm_b'] = pm_b
        return table
    def transform_pm(self,table):
        '''Poleski'''
        dec = table['dec']
        ra = table['ra']
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)
        pmra = table['pmra']
        pmdec = table['pmdec']

        C1 = np.sin(self.dec_np)*np.cos(dec_rad) - np.cos(self.dec_np)*np.sin(dec_rad)*np.cos(ra_rad - self.ra_np)
        C2 = np.cos(self.dec_np)*np.sin(ra_rad - self.ra_np)
        cosb = np.sqrt(C1**2 + C2**2)
        #C_mtrx = 1/(cosb)*np.array([C1, C2],[-C2, C1])
        
        mu_l_cosb = (1/cosb)* (C1 * pmra + C2 *pmdec)
        mu_b = (1/cosb) * (-C2 * pmra + C1*pmdec)
        
        table['pm_l_poleski'] = mu_l_cosb 
        table['pm_b_poleski'] = mu_b 
        #set units
        table['pm_l_poleski'].unit = u.mas/u.yr
        table['pm_b_poleski'].unit = u.mas/u.yr
        
        return table  

    def transform_space_velocity(self,ra, dec):
        '''Johnson 1986 calculate space velocities
        Use a transformation Matrix to convert ra and dec corrdinates into space velocity components
                                                    UVW
    
        Notes:
        Transform is correct for Johnson 1986 constants
        input:
        ra- right ascension of star in radians
        dec - declinatoin of star in radians
    
        Transform - see johnson 1986
        A- coordinate matrix
        Return the dot product of Transform and A
    
        
        '''
        
        T1= np.array([[np.cos(self.theta_o), np.sin(self.theta_o), 0.0],
                            [np.sin(self.theta_o), -np.cos(self.theta_o), 0.0],
                            [0.0,0.0,1.0]])
        T2 =  np.array([[-np.sin(self.dec_np), 0.0, np.cos(self.dec_np)],
                            [0.0,-1.0,0.0],
                            [np.cos(self.dec_np), 0.0, np.sin(self.dec_np)]])
        T3 = np.array([[np.cos(self.ra_np), np.sin(self.ra_np), 0.0],
                            [np.sin(self.ra_np), -np.cos(self.ra_np), 0.0],
                            [0.0,0.0,1.0]])
        # @ is a matrix operator
        #this is correct with Johnson 1986 values 
        Transform = T1 @ T2 @ T3 # transformation matrix
        A = np.array([[np.dot(np.cos(ra),np.cos(dec)), -np.sin(ra), np.dot(np.cos(ra),np.sin(dec))],
                     [np.dot(np.sin(ra), np.cos(dec)), np.cos(ra), -np.dot(np.sin(ra), np.sin(dec))],
                     [np.sin(dec), 0.0, np.cos(dec)]])
        #each star will have a unique matrix
        B = Transform @ A
        return B
    
    def calculate_space_velocity(self,table):
        '''Caulate the space velocity of a star WRT the local standard of rest, subtracting solar motion Johnson 1986'''
    
        UVW = []
        for row in table:
            ra = np.radians(row['ra'])
            dec = np.radians(row['dec'])
            
            pmra = row['pmra'] * 1e-3  # mas/yr to as/yr
            pmdec = row['pmdec'] * 1e-3
            
            prlx = row['parallax'] * 1e-3  # mas to as
            radial = row['radial_velocity']  # km/s
            B= self.transform_space_velocity(ra, dec) # calculate the transform tensor
        
            comp_array = np.array([radial, self.k*pmra/prlx, self.k*pmdec/prlx])
            UVW_val = np.dot(B,comp_array)
            #subtract solar motions
            UVW_val[0] = UVW_val[0] - self.U_sun
            UVW_val[1] = UVW_val[1] - self.V_sun
            UVW_val[2] = UVW_val[2] - self.W_sun
            UVW.append(UVW_val)
        UVW = np.array([UVW])
        table['U'] = UVW[:,0]*u.km/u.s
        table['V'] = UVW[:,1]*u.km/u.s
        table['W'] = UVW[:,2]*u.km/u.s
       # table.add_columns([space_U,space_V,space_W],names=['U','V','W'])
        
        return table
    def solar_proper_motion(self,table):
        '''Calculate the proper motion compoent due to solar motion in the galactic plane
        moffat 1998'''
        lat = table['b']
        long = table['l']
        dist = table['distance']
        #convert to radians
        long_rad = np.radians(long)
        lat_rad = np.radians(lat)
        
        Kr_mul_sol = self.U_sun*np.sin(long_rad) - self.V_sun*np.cos(long_rad)
        table['mu_l_sol'] = (Kr_mul_sol/(self.k*dist))
        table['mu_l_sol'].unit = u.mas/u.yr
        
        Kr_mub_sol = self.U_sun*np.cos(long_rad)*np.sin(lat_rad) + self.V_sun*np.sin(long_rad)*np.sin(lat_rad) - self.W_sun*np.cos(lat_rad)
        table['mu_b_sol'] = (Kr_mub_sol/(self.k*dist))
        table['mu_b_sol'].unit = u.mas/u.yr
        #table.add_columns([mul_sol,mub_sol],names=['mu_l_sol','mu_b_sol'])
        return table
    def flat_rotation_curve(self,table):
        '''Based off moffat 1998
        Model the proper motion in l and b due to the galactic rotation curve
        Uses flat rotation curve model good for 3 < R < 18 Kpc ( 2*R0)'''
        R = table['galactic distance']
        lat = table['b']
        long = table['l']
        dist = table['distance']
        omega_0 = self.sun_curve/self.R0 # km/s per kpc
        omega  =self.sun_curve/R
        long_rad = np.radians(long)
        lat_rad = np.radians(lat)
        #calculate proper motion due to galactic rotation
        K_mul_rot = ((self.R0/(dist*np.cos(lat_rad)))*(omega-omega_0)*np.cos(long_rad))-omega
        table['mu_l_rot'] = K_mul_rot/(self.k)
        table['mu_l_rot'].unit = u.mas/u.yr
        
        #do the same for latitudal proper motion
        K_mub_rot = -(self.R0/dist)*(omega-omega_0)*np.sin(lat_rad)*np.sin(long_rad)
        table['mu_b_rot'] = K_mub_rot/self.k
        table['mu_b_rot'].unit = u.mas/u.yr

        return table
    def peculiar_velocity(self, table):
        '''Calculate peculiar velocity of a star based on its proer motion 
        Follows Moffat 1998 
        account for peculiar motion as observed motion minus expeected motion
        account for solar bias
        See solar_proper_motion
        See flat_rotation_curve
        returns peculiar velocity in km/s
        peculiar proper motion components in mas/yr
        works now'''
        long = table['l'] # deg
        lat = table['b'] #deg
        long_rad =  np.radians(long)
        lat_rad = np.radians(lat)
        dist = table['distance'] # kpc
        # observed proper Motions
        obs_mu_l = table['pm_l_poleski']
        obs_mu_b = table['pm_b_poleski']
        #motion due to galactic rotation 
        
        rotation_mu_l = table['mu_l_rot']
        rotation_mu_b = table['mu_b_rot']
        
        #calculate solar motion for each object 
        mul_sol = table['mu_l_sol']
        mub_sol = table['mu_b_sol']
    
        pec_mu_l = obs_mu_l - mul_sol - rotation_mu_l
        pec_mu_b = obs_mu_b - mub_sol - rotation_mu_b
        #convert to km/s
        V_pec = self.k*dist*np.sqrt(pec_mu_l**2 + pec_mu_b**2)
        table['Peculiar Velocity'] = V_pec
        table['Peculiar Velocity'].unit = u.km/u.s
        #include the proper motion in mas/yr
        table['peculiar_mu_l']= pec_mu_l
        table['peculiar_mu_l'].unit = u.mas/u.yr
        table['peculiar_mu_b'] = pec_mu_b
        table['peculiar_mu_b'].unit = u.mas/u.yr
        return table
    def lay_pipe(self,filename,filetype):
        '''Combine everything'''
        table = self.make_table()
        table = self.modify_parllax_add_distance(table)
        table =self.transform_pm(table)
        #table = self.make_sky(table)
        table = self.Galaxy_dist(table)
        table = self.rotation_curve(table)
        table = self.vlsr_model(table)
        #table = self.calculate_space_velocity(table)
        table = self.solar_proper_motion(table)
        table = self.flat_rotation_curve(table)
        table = self.peculiar_velocity(table)
        table.write(filename, format=filetype,overwrite=True)
        return table
test_table = "GAIA_HMXB_DNE.ecsv"
if __name__ == "__main__":
    test_massive = pipeline(test_table,fmt='ecsv')
    my_table = test_massive.make_table()
    x = test_massive.lay_pipe(filename='test_lay_pipe.ecsv',filetype='ascii.ecsv')