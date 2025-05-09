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
from astropy.table import Table

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
    def calc_V_pec_errors(self,table):
        '''
        Error propogation on moffat 1998 with 
        galactic rotation curve brand and blitz 1993

        mu_pec = mu_obs - mu_gal_rot - mu_solar
        calculate errors on each parameter of mu_pec in l and b
        
        calaculate total mu_pec error
        convert to 1 -sigma V_pec_tan error
        
        Also calculate for 3d pec velocity
        
        
        
        propogate with error in cartesian vleocity
        Parameters
        ----------
        table : with errors and paraemters of HMXB.

        Returns
        -------
        results table : with v_pec_tan_errors.

        '''
        #make a result table
        results = Table()
        results['source_id'] = table['source_id'] # for cross matching
        #constants
        a1 = 1.00767
        a2 = 0.0394
        a3 = 0.00712
        #km/s solar values
        U_sun = 10.8
        V_sun = 13.6
        W_sun = 7.6
        # galactic constants
        R0 =8.15
        omega0 =236.0/R0 #km/s
        #variables
        dist = table['distance_bj']
        dist_high = table['distance_bj_high']
        dist_low = table['distance_bj_low']
        
        
        long = table['l']
        lat = table['b']
        long_rad = np.radians(long)
        lat_rad = np.radians(lat)
        sigma_l = table['l_err'] # one sigma
        sigma_b = table['b_err']
        #symetric sigma_d
        #have 1 sigma confidence interval for distances, convert to  symmetric error
        sigma_d= 0.5*((dist_high - dist) + (dist - dist_low)) # one sigma
        
        #errors in omega
        
        theta = table['circular velocity']
        gal_dist = table['galactic distance']
        omega = theta/gal_dist
        sigma_theta = (a1*a2*(gal_dist/R0)**(a2-1) * (1/R0))
        
        dRdd = ((dist- R0*np.cos(long_rad))/gal_dist)**2 
        dRdl = ((R0*dist*np.sin(lat_rad))/gal_dist)**2
        
        sigma_R_squrd = ((dRdd**2 * sigma_d**2 ) + (dRdl**2 * sigma_l**2))
        sigma_R = np.sqrt(sigma_R_squrd)
        
        sigma_omega = omega* np.sqrt((sigma_theta/theta)**2 + (sigma_R/gal_dist)**2)
     
        #this is where the fun begins
        #solar errors 
        dmu_l_sol_dd = -(U_sun*np.sin(long_rad) - V_sun*np.cos(long_rad)/dist)
        dmu_l_sol_dl = (U_sun*np.cos(long_rad) + V_sun*np.sin(long_rad))
        
        sigma_mu_sol_l = np.sqrt((1/dist**2) *(dmu_l_sol_dd**2 * sigma_d**2) + (dmu_l_sol_dl**2 * sigma_l**2))
        #sigma_mu_sol_b
        dmu_b_sol_dd = (U_sun*np.cos(long_rad)*np.sin(lat_rad) + V_sun*np.sin(long_rad)*np.sin(lat_rad) + W_sun*np.cos(lat_rad))/dist
        
        dmu_b_sol_db = (U_sun*np.cos(long_rad)*np.cos(lat_rad) + V_sun*np.sin(long_rad)*np.cos(lat_rad) - W_sun*np.sin(lat_rad))
        
        dmu_b_sol_dl = (-U_sun*np.sin(long_rad)*np.sin(lat_rad) + V_sun*np.sin(long_rad)*np.sin(lat_rad))
        
        sigma_mu_sol_b = np.sqrt((1/dist**2)* ((dmu_b_sol_dd**2 *sigma_d**2) +
                                               (dmu_b_sol_db**2 * sigma_b**2) +
                                               (dmu_b_sol_dl**2 * sigma_l**2)))
        
        
        #sigma mu_rot_b
        dmu_rot_b_dd = +(omega - omega0)*np.sin(lat_rad)*np.sin(long_rad)
        
        dmu_rot_b_domega = -np.sin(lat_rad)*np.sin(long_rad)
        
        dmu_rot_b_db = (omega - omega0)*np.sin(long_rad)*np.cos(lat_rad)
        
        dmu_rot_b_dl = (omega - omega0)*np.sin(lat_rad)*np.cos(long_rad)
        
        sigma_mu_rot_b = np.sqrt(((R0/dist)**2)*((dmu_rot_b_dd**2 * sigma_d**2) + 
                                                 (dmu_rot_b_db**2 * sigma_b**2) + 
                                                 (dmu_rot_b_dl**2  * sigma_l**2)+
                                                 (dmu_rot_b_domega**2 * sigma_omega**2)))
        #sgima mu_rot_l
        dmu_rot_l_dd = (omega - omega0)*np.cos(long_rad)/dist 
        
        dmu_rot_l_dl = -(omega-omega0)*np.sin(long_rad)
        
        dmu_rot_l_domega = (np.cos(long_rad) - (dist*np.cos(lat_rad)/R0))
        
        dmu_rot_l_db = ((omega - omega0)*np.cos(long_rad)*np.sin(lat_rad))/np.cos(lat_rad)
        
        sigma_mu_rot_l =np.sqrt((R0 / (dist * np.cos(lat_rad)))**2 * ((dmu_rot_l_dd**2 *sigma_d**2) + 
                                                                (dmu_rot_l_db**2 *sigma_b**2) +
                                                                (dmu_rot_l_dl**2 * sigma_l**2) +
                                                                (dmu_rot_l_domega**2 * sigma_omega**2)))
        
    
        #now calaculate mu_pec_errors
        sigma_mu_obs_l = table['pm_l_err']
        sigma_mu_obs_b = table['pm_b_err']
        
        mu_pec_l = table['peculiar_mu_l']
        mu_pec_b = table['peculiar_mu_b']
        
        mu_pec = np.sqrt(mu_pec_l**2 + mu_pec_b**2)
        
        #errors in peculair proper motion
        sigma_mu_l_pec = np.sqrt(sigma_mu_obs_l**2 + sigma_mu_rot_l**2 + sigma_mu_sol_l**2)
        sigma_mu_b_pec = np.sqrt(sigma_mu_obs_b**2 + sigma_mu_rot_b**2 + sigma_mu_sol_b**2)
        
        sigma_mu_pec = np.sqrt((mu_pec_l/mu_pec)**2 * sigma_mu_l_pec**2 + (mu_pec_b/mu_pec)**2 * sigma_mu_b_pec**2)
        
        #sigma_tan_pec 
        V_pec_tan = table['Peculiar Velocity']
        V_pec_rad = table['Peculiar Radial Velocity']
        sigma_v_pec_tan = V_pec_tan*np.sqrt((sigma_d/dist)**2 + (sigma_mu_pec/mu_pec)**2)
        results['V_pec_tan_err'] = sigma_v_pec_tan
        results['V_pec_tan_err'].unit = u.km/u.s
        
        #again for radial peculiar
        rv = table['RV']
        rv_err = table['RV_err']
        v_pec_3d = table['Peculiar Velocity 3D']
        
        dv_sol_dl = U_sun*np.sin(long_rad)*np.cos(lat_rad) - V_sun*np.cos(long_rad)*np.cos(lat_rad)
        
        dv_sol_db = U_sun*np.cos(long_rad)*np.sin(lat_rad) + V_sun*np.sin(long_rad)*np.sin(lat_rad) -  W_sun*np.cos(lat_rad)

        sigma_dvr_sol = np.sqrt((dv_sol_db**2 * sigma_b**2) + (dv_sol_dl**2 *sigma_l**2))
        
        #for vr rotational galaxy
        dv_rot_domega = R0*np.sin(long_rad)*np.cos(lat_rad)
        
        dv_rot_dl = R0*(omega - omega0)*np.cos(lat_rad)*np.cos(long_rad)
        
        dv_rot_db = -R0*(omega-omega0)*np.sin(long_rad)*np.sin(lat_rad)
    
        
        sigma_dvr_rot = np.sqrt((dv_rot_db**2 *sigma_b**2) + (dv_rot_dl**2 *sigma_l**2) + (dv_rot_domega**2 *sigma_omega**2))
        
        sigma_rv_pec = np.sqrt(rv_err**2 + sigma_dvr_rot**2 + sigma_dvr_sol**2)
        
        
        #total peculair error
        sigma_vpec_3d = v_pec_3d*np.sqrt((sigma_v_pec_tan/V_pec_tan)**2 + (sigma_rv_pec/V_pec_rad)**2)
        results['V_pec_3d_err'] = sigma_vpec_3d
        results['V_pec_3d_err'].unit = u.km/u.s        
        
        return results
        
        
        
    def calc_all_errors(self,table):
        table = self.galactic_coord_errs(self.proper_motion_errors(table))
        return table
        
cwd = os.getcwd()
home_files = os.path.dirname(cwd)
csv_files  = cwd+ '/Documents/UvA/Thesis/DATA/'
today = datetime.now().strftime("%Y%m%d")

test_table = ascii.read(csv_files+'HMXB_pm_errs-result.ecsv',format='ecsv')
hd7636 = ascii.read(csv_files+'hd7636_complete_analysis.ecsv',format='ecsv')
HMXB_table = ascii.read(csv_files+'HMXB_20250326_.ecsv',format='ecsv')
orion_table = ascii.read(csv_files+'ori_table_20250417.ecsv', format='ecsv')

iota_table = ascii.read(csv_files+'iota_table_20250417.ecsv', format='ecsv')
if __name__ == "__main__":
    #test_table = calc_errors().gaia_jacobian(test_table)
    #test_table.write(csv_files+f'HMXB_all_errors_{today}.ecsv',format='ascii.ecsv',overwrite=True)
    # v_pec_errors = calc_errors().calc_V_pec_errors(HMXB_table)
    # v_pec_errors.write(csv_files+f"HMXB_vpec_errs_{today}.ecsv",format='ascii.ecsv',overwrite=True)
    
    hd7636_jacobian = calc_errors().gaia_jacobian(hd7636)
    hd7636_jacobian.write(csv_files+'hd7636_jacobian.ecsv',format='ascii.ecsv',overwrite=True)
    hd7636_vpec_errors = calc_errors().calc_V_pec_errors(hd7636)
    hd7636_vpec_errors.write(csv_files+'hd7636_vpec_errs.ecsv',format='ascii.ecsv',overwrite=True)
    # orion_jacobian = calc_errors().gaia_jacobian(orion_table)
    # orion_jacobian.write(csv_files+'orion_jacobian.ecsv',format='ascii.ecsv',overwrite=True)
    # orion_v_pec_errors = calc_errors().calc_V_pec_errors(orion_jacobian)
    # orion_v_pec_errors.write(csv_files+'orion_v_pec_errors.ecsv',format='ascii.ecsv',overwrite=True)
    
    # iota_jacobian = calc_errors().gaia_jacobian(iota_table)
    # iota_jacobian.write(csv_files+'iota_jacobian.ecsv',format='ascii.ecsv',overwrite=True)
    # iota_v_pec_errors = calc_errors().calc_V_pec_errors(iota_jacobian)
    # iota_v_pec_errors.write(csv_files+'iota_v_pec_errors.ecsv',format='ascii.ecsv',overwrite=True)
    
    