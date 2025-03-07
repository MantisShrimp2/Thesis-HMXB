#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 11:21:00 2025

@author: karan
"""

import numpy as np
import os 
from astropy.io import ascii
import astropy.units as u
from datetime import datetime
from astropy.coordinates import SkyCoord
rotation_matrix = np.array([[-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
                            [0.4941094278755837, -0.4448296299600112, 0.7469822444972189],
                            [-0.8676661490190047, -0.1980763734312015, 0.4559837761750669]])


class calc_errors:
    def __init__(self):
        #this is already transposed
        self.rotation_matrix = np.array([[-0.0548755604162154, -0.8734370902348850,-0.4838350155487132],
                                [0.4941094278755837, -0.4448296299600112, 0.7469822444972189],
                                [-0.8676661490190047, -0.1980763734312015, 0.4559837761750669]])
        self.k = 4.74
        
    def proper_motion_errors(self,table):
        #define cordinates
        pm_l_err_list = []
        pm_b_err_list = []
        for row in table:
            ra = np.radians(row['ra'])
            dec = np.radians(row['dec'])
            l = np.radians(row['l'])
            b = np.radians(row['b'])
            # erros of coordinates
            sigma_pmra  = np.radians(row['pmra_error'])
            sigma_pmdec = np.radians(row['pmdec_error'])
            # correlation between pmra and pmdec
            correl_coef = row['pmra_pmdec_corr'] #unitless
            
            diag_term = correl_coef*sigma_pmra*sigma_pmdec # the (0,1) and (1,0) term of C_naught
            
            #define transfomation matricies
            p_icrs = np.array([-np.sin(ra), np.cos(ra),0])
            q_icrs = np.array([np.cos(ra)*np.sin(dec), -np.sin(ra)*np.sin(dec), np.cos(dec)])
            icrs_matrix = np.array([p_icrs, q_icrs])
            #same for galactic
            p_gal = np.array([-np.sin(l), np.cos(l),0])
            q_gal = np.array([-np.cos(l)*np.sin(b), -np.sin(l)*np.sin(b), np.cos(b)])
            gal_matrix = np.array([p_gal, q_gal])
   
            G_matrix = gal_matrix @ self.rotation_matrix @ icrs_matrix.T
            
            C_naught = np.array([[sigma_pmra**2 , diag_term],
                                [diag_term, sigma_pmdec**2]])
            Covariance_matrix = G_matrix @ C_naught @ G_matrix.T
            
            #report the standard deviations
            pm_l_err  = np.sqrt(Covariance_matrix[0,0])
            pm_b_err = np.sqrt(Covariance_matrix[1,1])
            pm_l_err_list.append(pm_l_err)
            pm_b_err_list.append(pm_b_err)            
            

        table['pm_l_err'] = pm_l_err_list
        table['pm_l_err'].unit = u.mas/u.yr
        table['pm_b_err'] = pm_b_err_list
        table['pm_b_err'].unit = u.mas/u.yr
    
        
        return table
    def galactic_coord_errs(self,table):
        ra_err = table['ra_error'].to(u.deg) #already as ra*cos(dec)
        
        dec_err = table['dec_error'].to(u.deg)
        err_coord = SkyCoord(ra=ra_err, dec=dec_err, frame='icrs')
        gal_coord = err_coord.galactic
        
        l_err = gal_coord.l
        b_err = gal_coord.b
        
        table['l_err'] = l_err
        table['b_err'] = b_err
        return table
    def calc_all_errors(self,table):
        table = self.galactic_coord_errs(self.proper_motion_errors(table))
        return table
        
cwd = os.getcwd()
home_files = os.path.dirname(cwd)
csv_files  = cwd + '/Documents/UvA/Thesis/DATA/'
today = datetime.now().strftime("%Y%m%d")
test_table = ascii.read(csv_files+'HMXB_pm_errs-result.ecsv',format='ecsv')

if __name__ == "__main__":
    test_table = calc_errors().calc_all_errors(test_table)
    test_table.write(csv_files+f'HMXB_all_errors_{today}.ecsv',format='ascii.ecsv',overwrite=True)
test_table['b_err']