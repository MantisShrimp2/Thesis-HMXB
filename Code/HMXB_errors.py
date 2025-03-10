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
    def calculate_gaia_cov(self,sigma_vector, corr_vector):
        '''Calculate the correlaction matrix for error propogation in galctic coordinates
        does this for one star
        needs:
            sigma_vector- errors in ra, dec, parallax, pmra, pmdec in that order
            corr_vector - correlations between sigmas in the order defined by gaia archive
            '''
        n = int(np.sqrt(len(corr_vector) * 2 + 0.25) + 0.5)  # Calculate the size of the matrix
        rho = np.eye(n)  # Initialize with 1s on the diagonal

        # Fill the upper triangle
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                rho[i, j] = corr_vector[k]
                rho[j, i] = corr_vector[k]  # Enforce symmetry
                k += 1
                
        n = len(sigma_vector)
        C = np.zeros((n,n))
        for i in range(n):
            for j in range(i,n):
                if i ==j:
                    C[i,j] = sigma_vector[i]**2
                else:
                    C[i,j] = sigma_vector[i]* sigma_vector[j] * rho[i,j]
                    C[j,i] = C[i,j]
                    
        return C
            
    def gaia_jacobian(self,table):
        pm_l_err_list = []
        pm_b_err_list = []
        l_err_list = []
        b_err_list = []
        for row in table:
            ra = row['ra'] * u.mas
            dec = row['dec'] * u.mas
            
            ra_rad, dec_rad = ra.to(u.rad), dec.to(u.rad)
            
            l = np.radians(row['l'])
            b = np.radians(row['b'])
            # erros of coordinates
            sigma_ra = row['ra_error']
            sigma_dec = row['dec_error']
            sigma_parallax = row['parallax_error']
            sigma_pmra  = row['pmra_error'] #mas/yr
            sigma_pmdec = row['pmdec_error'] #mas.yr
            
            sigma_vector = np.array([sigma_ra, sigma_dec, sigma_parallax, sigma_pmra, 
                                   sigma_pmdec])
            
            
            # correlation coefficents 
            ra_dec_corr = row['ra_dec_corr']
            ra_parallax_corr = row['ra_parallax_corr']
            ra_pmra_corr = row['ra_pmra_corr']
            ra_pmdec_corr = row['ra_pmdec_corr']
            
            dec_parallax_corr = row['dec_parallax_corr']
            dec_pmra_corr = row['dec_pmra_corr']
            dec_pmdec_corr = row['dec_pmdec_corr']
            
            parallax_pmra_corr = row['parallax_pmra_corr']
            parallax_pmdec_corr = row['parallax_pmdec_corr']
            
            
            pmra_pmdec_corr= row['pmra_pmdec_corr'] #unitless
            
            corr_vector = np.array([ra_dec_corr,ra_parallax_corr,ra_pmra_corr,
                                    ra_pmdec_corr, dec_parallax_corr, dec_pmra_corr,
                                    dec_pmdec_corr, parallax_pmra_corr, parallax_pmdec_corr,
                                    pmra_pmdec_corr])
            Covariance_matrix = self.calculate_gaia_cov(sigma_vector, corr_vector)
            #define e_matrix
            
            e_matrix = np.array([sigma_ra, sigma_dec, sigma_parallax, sigma_pmra, sigma_pmdec])
            
            #define transfomation matricies
            p_icrs = np.array([-np.sin(ra_rad), np.cos(ra_rad),0])
            q_icrs = np.array([np.cos(ra_rad)*np.sin(dec_rad), -np.sin(ra_rad)*np.sin(dec_rad), np.cos(dec_rad)])
            icrs_matrix = np.array([p_icrs, q_icrs])
            #same for galactic
            p_gal = np.array([-np.sin(l), np.cos(l),0])
            q_gal = np.array([-np.cos(l)*np.sin(b), -np.sin(l)*np.sin(b), np.cos(b)])
            gal_matrix = np.array([p_gal, q_gal])
   
            G_matrix = gal_matrix @ self.rotation_matrix @ icrs_matrix.T
            J = np.block([
                [G_matrix, np.zeros((2,1)), np.zeros((2,2))],
                [np.zeros((1,2)), np.array([[1]]), np.zeros((1,2))],
                [np.zeros((2,2)), np.zeros((2,1)), G_matrix]
                ])
            C_gal = np.dot(np.dot(J,Covariance_matrix),J.T)
            l_err_list.append(np.sqrt(C_gal[0,0]))
            b_err_list.append(np.sqrt(C_gal[1,1]))
            #skip parllax
            pm_l_err_list.append(np.sqrt(C_gal[3,3]))
            pm_b_err_list.append(np.sqrt(C_gal[4,4]))
            
        table['l_err'] = l_err_list
        table['b_err'] = b_err_list
        table['pm_l_err'] = pm_l_err_list
        table['pm_b_err'] = pm_b_err_list
    
        
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
    test_table = calc_errors().gaia_jacobian(test_table)
    test_table.write(csv_files+f'HMXB_all_errors_{today}.ecsv',format='ascii.ecsv',overwrite=True)
#test_table['b_err']