import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import camb as cb
import healpy as hp
import warnings
import os as os
import shutil as sh
import sys
import time
sys.path.insert(0, '/global/cscratch1/sd/nishant/lenspix')
warnings.filterwarnings('ignore')

#camb runs a cosmological simulation based off start parameters
#returns a bunch of stuff in results, but we choose to return powers
#powers contains lensed unlensed power spectra etc
def generateCAMB():

    #Set up a new set of parameters for CAMB
    pars = cb.CAMBparams()
    #This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
    pars.InitPower.set_params(ns=0.965, r=0)
    pars.set_for_lmax(2500, lens_potential_accuracy=0);

    #calculate results for these parameters
    results = cb.get_results(pars)

    #get dictionary of CAMB power spectra (total, lens_potential, lensed_scalar, unlensed_scalar, unlensed_total, tensor)
    powers = results.get_cmb_power_spectra(pars, CMB_unit='muK', raw_cl=True)

    return powers

#takes a power spectrum and creates a gaussian random field based off that spectrum
def generateGRF(cl):
    return hp.sphtfunc.synfast(cl, 1024)

#generates a power spectrum given a map
def generatePowerSpectrum(InputMap):
    cl = hp.anafast(InputMap, iter=10)
    return cl

#generates a cross power spectrum given two maps
def generateCrossPowerSpectrum(Input1, Input2):
    cl = hp.anafast(Input1, Input2)
    return cl

#converts fits file to map
def fitsTOmap(fits):
    fits_Read = hp.read_map(fits, verbose=False)
    mask = hp.read_map(fits, verbose=False).astype(np.bool)
    fits_Read_masked = hp.ma(fits_Read)

    return fits_Read_masked.filled()

#converts map to fits file
def mapTOfits(fits_name, map):
    hp.write_map(fits_name, map, overwrite=True)

#plots a spectrum
def plotSpectrum(cl, label):
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r"$\ell$", fontsize=15); plt.ylabel(r'$\mathcal{C}_l$', fontsize=15); plt.grid()
    ell = np.arange(len(cl))
    plt.plot(ell, cl, label=label)
    plt.legend()

#takes the ratio between two power spectra
def plotRatioSpectrum(cl1, cl2, n, label):
    cl_ratio = cl1[1:n]/cl2[1:n]
    plt.yscale('linear')
    plt.xscale('log')
    plt.xlabel(r"$\ell$", fontsize=15); plt.ylabel('ratio', fontsize=15); plt.grid()
    ell = np.arange(len(cl_ratio))
    plt.plot(ell, cl_ratio, label=label)
    plt.legend()

#takes multiple power spectra and averages for reduced error at a rate of 1/sqrt(n)
def averagedSpectrum(n, Start_cl):
    total_cl = np.zeros(3072)
    for num in range(0,n):
        map = generateGRF(Start_cl)
        cl = generatePowerSpectrum(map)
        total_cl = total_cl + cl

    final_cl = total_cl/n

    return final_cl

#moves into the scratch directory

def reMake():
    os.system("make clean")
    os.system("make")



#change parameters inside the params.ini, make sure the nside matches your other maps
#def SimLens():
#    reMake()
    os.system('srun -n 1 ./simlens params.ini')






##########^^^^^^^^^BASIC FUNCTIONS^^^^^^^^^^^^^^#################
##########_________FUNCTIONS THAT ARE DIRECTLY RUN___________#############



def main():
    localtime = time.asctime( time.localtime(time.time()) )
    print("Local current time :", localtime)
    num = sys.argv[1]
    print("working")
    powers = generateCAMB()
    print("working2")
    unlensed_cl_CAMB = powers['unlensed_scalar'][:,0]
    lensing_potential_cl_CAMB = powers['lens_potential'][:,0]
#    lensed_cl_CAMB = powers['lensed_scalar'][:,0]
    print("working3")

    unlensed_cl_CAMB_map = generateGRF(unlensed_cl_CAMB)
    lensing_potential_cl_CAMB_map = generateGRF(lensing_potential_cl_CAMB)
#    lensed_cl_CAMB_map = generateGRF(lensed_cl_CAMB)
    print("working4")


    str4 = 'unlensed_cl_CAMB_map' + str(num) +'.fits'
    str5 = 'lensing_potential_cl_CAMB_map' + str(num) +'.fits'
    mapTOfits(str4 , unlensed_cl_CAMB_map)
    mapTOfits(str5 , lensing_potential_cl_CAMB_map)
#    mapTOfits('lensed_cl_CAMB_map.fits', lensed_cl_CAMB_map)
    print("working5")

    move_params = 'mv /global/cscratch1/sd/nishant/params_numbered/params' + str(num) + ".ini  /global/cscratch1/sd/nishant/lenspix/"
    back_params = 'mv /global/cscratch1/sd/nishant/lenspix/params' + str(num) + ".ini  /global/cscratch1/sd/nishant/params_numbered/"

    os.system(move_params)
    command_str = 'srun -n 1 ./simlens params' + str(num) + '.ini'
    os.system(command_str)
    os.system(back_params)


    print("working6")

    str6 = 'lensed_map_LENSPIX' + str(num) + '.fits'

    lensed_cl_LENSPIX_map = fitsTOmap(str6)

    print("working7")

    unlensed_cl_CAMB_map_cl = generatePowerSpectrum(unlensed_cl_CAMB_map)
    lensing_potential_cl_CAMB_map_cl = generatePowerSpectrum(lensing_potential_cl_CAMB_map)
    lensed_cl_LENSPIX_map_cl = generatePowerSpectrum(lensed_cl_LENSPIX_map)
    print("working8")

    str1 = 'unlensed_cl_CAMB_map_cl' + str(num) + '.txt'
    str2 = 'lensing_potential_cl_CAMB_map_cl' + str(num) + '.txt'
    str3 = 'lensed_map_LENSPIX_cl' +str(num) + '.txt'
    print("working9")



    np.savetxt(str1, unlensed_cl_CAMB_map_cl)
    np.savetxt(str2, lensing_potential_cl_CAMB_map_cl)
    np.savetxt(str3, lensed_cl_LENSPIX_map_cl)
    print("working10")

    str7 = 'mv /global/cscratch1/sd/nishant/lenspix/' + str1 + "  /global/cscratch1/sd/nishant/unlensed_cl_CAMB_map_cl"
    str8 = 'mv /global/cscratch1/sd/nishant/lenspix/' + str2 + "  /global/cscratch1/sd/nishant/lensing_potential_cl_CAMB_map_cl"
    str9 = 'mv /global/cscratch1/sd/nishant/lenspix/' + str3 + "  /global/cscratch1/sd/nishant/lensed_map_LENSPIX_cl"

    str10 = 'mv /global/cscratch1/sd/nishant/lenspix/' + str4 + "  /global/cscratch1/sd/nishant/unlensed_cl_CAMB_map"
    str11 = 'mv /global/cscratch1/sd/nishant/lenspix/' + str5 + "  /global/cscratch1/sd/nishant/lensing_potential_cl_CAMB_map"
    str12 = 'mv /global/cscratch1/sd/nishant/lenspix/' + str6 + "  /global/cscratch1/sd/nishant/lensed_map_LENSPIX"
    print("working11")

    os.system(str7)
    os.system(str8)
    os.system(str9)
    os.system(str10)
    os.system(str11)
    os.system(str12)
    print("working12")

    localtime = time.asctime( time.localtime(time.time()) )
    print("Local current time :", localtime)

main()

