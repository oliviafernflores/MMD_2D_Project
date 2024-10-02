import dadi
import matplotlib.pyplot as plt
import numpy as np
import csv
# import inbreeding_models as mods
import demes
import demesdraw
import pandas as pd

def get_popt(fname):
    '''
    get best fit parameters from a file
    
    returns full list, including theta0 and log likelihood
    '''
    popt = []
    fits_file = open(fname, 'r')
    fits = fits_file.readlines()
    popt = fits[0]
    popt = popt.split('\t')
    for f in fits:
        if f.split('\t')[0] < popt[0]:
            popt = f.split('\t')
    for i in range(len(popt)):
        popt[i] = float(popt[i])
    return popt
def convert(model, popt, factor):
    if 'im_pre' in model:
        # get raw param values
        nupre = popt[1]
        tpre = popt[2]
        s = popt[3]
        one_minus_s = 1 - s
        nu1 = popt[4]
        nu2 = popt[5]
        t = popt[6]
        m12 = popt[7]
        m21 = popt[8]
        theta = popt[-1]
        # convert param values
        nref = theta / (4 * factor)
        npre = nupre * nref
        time_pre_in_gen = 2 * tpre * nref
        split_1 = s * nref
        split_2 = one_minus_s * nref
        final_1 = nu1 * nref
        final_2 = nu2 * nref
        time_in_gen = 2 * t * nref
        mig_12 = m12 / (2 * nref)
        mig_21 = m21 / (2 * nref)
        # add converted values to dictionary
        params = {'Log likelihood: ':popt[0], 'Nref: ':nref, 'Npre: ':npre, 'Tpre (generations): ':time_pre_in_gen, 
                  'S1: ':split_1, 'S2: ':split_2, 'Nu1: ':final_1, 'Nu2: ':final_2, 
                  'T (generations): ':time_in_gen, 'M12: ':mig_12, 'M21: ':mig_21}
        # acccount for inbreeding or not, add misid and inbreding (if applicable) to dictionary
        if len(popt[9:-1]) > 1:
            params['F1: '] = popt[9]
            params['F2: '] = popt[10]
            params['Misid: '] = popt[11]
        else:
            params['Misid: '] = popt[9]
        return params
    elif 'im' in model:
        # get raw param values
        s = popt[1]
        one_minus_s = 1 - s
        nu1 = popt[2]
        nu2 = popt[3]
        t = popt[4]
        m12 = popt[5]
        m21 = popt[6]
        theta = popt[-1]
        # convert param values
        nref = theta / (4 * factor)
        split_1 = s * nref
        split_2 = one_minus_s * nref
        final_1 = nu1 * nref
        final_2 = nu2 * nref
        time_in_gen = 2 * t * nref
        mig_12 = m12 / (2 * nref)
        mig_21 = m21 / (2 * nref)
        # add converted values to dictionary
        params = {'Log likelihood: ':popt[0], 'Nref: ':nref, 'S1: ':split_1, 'S2: ':split_2, 'Nu1: ':final_1, 'Nu2: ':final_2, 
                  'T (generations): ':time_in_gen, 'M12: ':mig_12, 'M21: ':mig_21}
        # acccount for inbreeding or not, add misid and inbreding (if applicable) to dictionary
        if len(popt[7:-1]) > 1:
            params['F1: '] = popt[7]
            params['F2: '] = popt[8]
            params['Misid: '] = popt[9]
        else:
            params['Misid: '] = popt[7]
        return params
    elif 'split_mig' in model:
        # get raw param values
        nu1 = popt[1]
        nu2 = popt[2]
        t = popt[3]
        m = popt[4]
        theta = popt[-1]
        # convert param values
        nref = theta / (4 * factor)
        final_1 = nu1 * nref
        final_2 = nu2 * nref
        time_in_gen = 2 * t * nref
        mig = m / (2 * nref)
        # add converted values to dictionary
        params = {'Log likelihood: ':popt[0], 'Nref: ':nref, 'Nu1: ':final_1, 'Nu2: ':final_2, 
                  'T (generations): ':time_in_gen, 'Migration: ':mig}
        # acccount for inbreeding or not, add misid and inbreding (if applicable) to dictionary
        if len(popt[5:-1]) > 1:
            params['F1: '] = popt[5]
            params['F2: '] = popt[6]
            params['Misid: '] = popt[7]
        else:
            params['Misid: '] = popt[5]
        return params
    elif 'split_asym_mig' in model:
        # get raw param values
        nu1 = popt[1]
        nu2 = popt[2]
        t = popt[3]
        m12 = popt[4]
        m21 = popt[5]
        theta = popt[-1]
        # convert param values
        nref = theta / (4 * factor)
        final_1 = nu1 * nref
        final_2 = nu2 * nref
        time_in_gen = 2 * t * nref
        mig_12 = m12 / (2 * nref)
        mig_21 = m21 / (2 * nref)
        # add converted values to dictionary
        params = {'Log likelihood: ':popt[0], 'Nref: ':nref, 'Nu1: ':final_1, 'Nu2: ':final_2, 
                  'T (generations): ':time_in_gen, 'M12: ':mig_12, 'M21: ':mig_21}
        # acccount for inbreeding or not, add misid and inbreding (if applicable) to dictionary
        if len(popt[6:-1]) > 1:
            params['F1: '] = popt[6]
            params['F2: '] = popt[7]
            params['Misid: '] = popt[8]
        else:
            params['Misid: '] = popt[6]
        return params
    elif 'split_delay_mig' in model:
        # get raw param values
        nu1 = popt[1]
        nu2 = popt[2]
        tpre = popt[3]
        t = popt[4]
        m12 = popt[5]
        m21 = popt[6]
        theta = popt[-1]
        # convert param values
        nref = theta / (4 * factor)
        final_1 = nu1 * nref
        final_2 = nu2 * nref
        time_pre_in_gen = 2 * tpre * nref
        time_in_gen = 2 * t * nref
        mig_12 = m12 / (2 * nref)
        mig_21 = m21 / (2 * nref)
        # add converted values to dictionary
        params = {'Log likelihood: ':popt[0], 'Nref: ':nref, 'Nu1: ':final_1, 'Nu2: ':final_2, 'Tpre (generations): ':time_pre_in_gen,
                  'Tmig (generations): ':time_in_gen, 'M12: ':mig_12, 'M21: ':mig_21}
        # acccount for inbreeding or not, add misid and inbreding (if applicable) to dictionary
        if len(popt[7:-1]) > 1:
            params['F1: '] = popt[7]
            params['F2: '] = popt[8]
            params['Misid: '] = popt[9]
        else:
            params['Misid: '] = popt[7]
        return params
    else:
        return None
def make_table(dict, fname):
    '''
    {'split mig' : {'log likelihood:' -900, 'Nref': 100000, etc}, 'split mig + inbreeding': ....}
    '''
    cols = ['split mig', 'split mig + inbreeding', 
            'split asym mig', 'split asym mig + inbreeding', 
            'split delay mig', 'split delay mig + inbreeding', 
            'im', 'im + inbreeding', 
            'im pre', 'im pre + inbreeding']
    indx = ['Log likelihood', 'Npre', 'Tpre (generations)', 
            'S1', 'S2', 'Nu1', 'Nu2', 
            'T (generations)', 'Tmig (generations)', 
            'Migration', 'M12', 'M21', 
            'F1', 'F2', 'Misid']
    df = pd.DataFrame(np.nan, columns = cols, index = indx)
    for key in dict.keys():
        for param in dict[key].keys():
            p = param.split(':')[0]
            df.loc[p, key] = (dict[key][param])
    df.to_csv(fname, sep='\t', index=True)
    return df

def main():
    # IRA vs. FRA
    split_mig = convert('split_mig', get_popt('IRA_FRA/demo_results/IRA_FRA_split_mig_demo_fits_combined.txt'), 0.0756)
    split_mig_inbreeding = convert('split_mig_inbreeding', get_popt('IRA_FRA/demo_results/IRA_FRA_split_mig_inbreeding_demo_fits_combined.txt'), 0.0756)
    split_asym_mig = convert('split_asym_mig', get_popt('IRA_FRA/demo_results/IRA_FRA_split_asym_mig_demo_fits_combined.txt'), 0.0756)
    split_asym_mig_inbreeding = convert('split_asym_mig_inbreeding', get_popt('IRA_FRA/demo_results/IRA_FRA_split_asym_mig_inbreeding_demo_fits_combined.txt'), 0.0756)
    split_delay_mig = convert('split_delay_mig', get_popt('IRA_FRA/demo_results/IRA_FRA_split_delay_mig_demo_fits_combined.txt'), 0.0756)
    split_delay_mig_inbreeding = convert('split_delay_mig_inbreeding', get_popt('IRA_FRA/demo_results/IRA_FRA_split_delay_mig_inbreeding_demo_fits_combined.txt'), 0.0756)
    im = convert('im', get_popt('IRA_FRA/demo_results/IRA_FRA_im_demo_fits_combined.txt'), 0.0756)
    im_inbreeding = convert('im_inbreeding', get_popt('IRA_FRA/demo_results/IRA_FRA_im_inbreeding_demo_fits_combined.txt'), 0.0756)
    im_pre = convert('im_pre', get_popt('IRA_FRA/demo_results/IRA_FRA_im_pre_demo_fits_combined.txt'), 0.0756)
    im_pre_inbreeding = convert('im_pre_inbreeding', get_popt('IRA_FRA/demo_results/IRA_FRA_im_pre_inbreeding_demo_fits_combined.txt'), 0.0756)
    IRA_FRA_dict = {'split mig': split_mig, 'split mig + inbreeding': split_mig_inbreeding, 
                    'split asym mig': split_asym_mig, 'split asym mig + inbreeding': split_asym_mig_inbreeding,
                    'split delay mig': split_delay_mig, 'split delay mig + inbreeding': split_delay_mig_inbreeding,
                    'im': im ,'im + inbreeding': im_inbreeding,
                    'im pre': im_pre,'im pre + inbreeding': im_pre_inbreeding}
    table = make_table(IRA_FRA_dict, 'IRA_FRA_demography_converted_parameter_results.csv')
    print(table)
    
    # FRA vs. GER
    split_mig = convert('split_mig', get_popt('FRA_GER/demo_results/FRA_GER_split_mig_demo_fits_combined.txt'), 0.0756)
    split_mig_inbreeding = convert('split_mig_inbreeding', get_popt('FRA_GER/demo_results/FRA_GER_split_mig_inbreeding_demo_fits_combined.txt'), 0.0756)
    split_asym_mig = convert('split_asym_mig', get_popt('FRA_GER/demo_results/FRA_GER_split_asym_mig_demo_fits_combined.txt'), 0.0756)
    split_asym_mig_inbreeding = convert('split_asym_mig_inbreeding', get_popt('FRA_GER/demo_results/FRA_GER_split_asym_mig_inbreeding_demo_fits_combined.txt'), 0.0756)
    split_delay_mig = convert('split_delay_mig', get_popt('FRA_GER/demo_results/FRA_GER_split_delay_mig_demo_fits_combined.txt'), 0.0756)
    split_delay_mig_inbreeding = convert('split_delay_mig_inbreeding', get_popt('FRA_GER/demo_results/FRA_GER_split_delay_mig_inbreeding_demo_fits_combined.txt'), 0.0756)
    im = convert('im', get_popt('FRA_GER/demo_results/FRA_GER_im_demo_fits_combined.txt'), 0.0756)
    im_inbreeding = convert('im_inbreeding', get_popt('FRA_GER/demo_results/FRA_GER_im_inbreeding_demo_fits_combined.txt'), 0.0756)
    im_pre = convert('im_pre', get_popt('FRA_GER/demo_results/FRA_GER_im_pre_demo_fits_combined.txt'), 0.0756)
    im_pre_inbreeding = convert('im_pre_inbreeding', get_popt('FRA_GER/demo_results/FRA_GER_im_pre_inbreeding_demo_fits_combined.txt'), 0.0756)
    FRA_GER_dict = {'split mig': split_mig, 'split mig + inbreeding': split_mig_inbreeding, 
                    'split asym mig': split_asym_mig, 'split asym mig + inbreeding': split_asym_mig_inbreeding,
                    'split delay mig': split_delay_mig, 'split delay mig + inbreeding': split_delay_mig_inbreeding,
                    'im': im ,'im + inbreeding': im_inbreeding,
                    'im pre': im_pre,'im pre + inbreeding': im_pre_inbreeding}
    table = make_table(FRA_GER_dict, 'FRA_GER_demography_converted_parameter_results.csv')
    print(table)
    
    # GER vs. HEL
    split_mig = convert('split_mig', get_popt('GER_HEL/demo_results/GER_HEL_split_mig_demo_fits_combined.txt'), 0.0756)
    split_mig_inbreeding = convert('split_mig_inbreeding', get_popt('GER_HEL/demo_results/GER_HEL_split_mig_inbreeding_demo_fits_combined.txt'), 0.0756)
    split_asym_mig = convert('split_asym_mig', get_popt('GER_HEL/demo_results/GER_HEL_split_asym_mig_demo_fits_combined.txt'), 0.0756)
    split_asym_mig_inbreeding = convert('split_asym_mig_inbreeding', get_popt('GER_HEL/demo_results/GER_HEL_split_asym_mig_inbreeding_demo_fits_combined.txt'), 0.0756)
    split_delay_mig = convert('split_delay_mig', get_popt('GER_HEL/demo_results/GER_HEL_split_delay_mig_demo_fits_combined.txt'), 0.0756)
    split_delay_mig_inbreeding = convert('split_delay_mig_inbreeding', get_popt('GER_HEL/demo_results/GER_HEL_split_delay_mig_inbreeding_demo_fits_combined.txt'), 0.0756)
    im = convert('im', get_popt('GER_HEL/demo_results/GER_HEL_im_demo_fits_combined.txt'), 0.0756)
    im_inbreeding = convert('im_inbreeding', get_popt('GER_HEL/demo_results/GER_HEL_im_inbreeding_demo_fits_combined.txt'), 0.0756)
    im_pre = convert('im_pre', get_popt('GER_HEL/demo_results/GER_HEL_im_pre_demo_fits_combined.txt'), 0.0756)
    im_pre_inbreeding = convert('im_pre_inbreeding', get_popt('GER_HEL/demo_results/GER_HEL_im_pre_inbreeding_demo_fits_combined.txt'), 0.0756)
    GER_HEL_dict = {'split mig': split_mig, 'split mig + inbreeding': split_mig_inbreeding, 
                    'split asym mig': split_asym_mig, 'split asym mig + inbreeding': split_asym_mig_inbreeding,
                    'split delay mig': split_delay_mig, 'split delay mig + inbreeding': split_delay_mig_inbreeding,
                    'im': im ,'im + inbreeding': im_inbreeding,
                    'im pre': im_pre,'im pre + inbreeding': im_pre_inbreeding}
    table = make_table(GER_HEL_dict, 'GER_HEL_demography_converted_parameter_results.csv')
    print(table)

if __name__ == '__main__':
    main()