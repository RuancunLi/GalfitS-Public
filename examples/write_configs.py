import sys
import os
from astropy.table import Table, Column
import matplotlib.pyplot as plt
import numpy as np
import time
import shutil
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM,z_at_value
from astropy.io import fits, ascii
from astropy.wcs import WCS

sample =Table.read('./gs_test.txt', format='ascii')
fucata = ascii.read('./phot_sub.mrt')
tmassBand = ['j','h','ks']
Band2mass = ['J','H','K']
mag2AB = {'nuv': 0, 'fuv': 0, 'u': -0.04,'g': 0., 'r': 0., 'i': 0., 'z': 0.02,
                  'J': 0.941, 'H': 1.38, 'K': 1.86, 'w1': 2.699, 'w2': 3.339, 'w3': 5.174, 'w4': 6.620 }
Bands = [ 'sloan_u','sloan_g', 'sloan_r', 'sloan_i', 'sloan_z', 'j','h','ks', 'galex_nuv', 'galex_fuv','wise_ch1', 'wise_ch2', 'wise_ch3', 'wise_ch4','herschel_psw','herschel_pmw','herschel_plw']
cosmo = FlatLambdaCDM(H0=67.8 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.308)
loop =0
while loop < (len(sample)):
    name = sample['ID'][loop]
    ebv = sample['e_B_V'][loop]
    ra,dec = sample['RAdeg'][loop],sample['DEdeg'][loop]
    z = sample['z'][loop]
    age = cosmo.age(z).value
    cutr = float (3 * sample['R90'][loop])
    outpath = './result/{0}/'.format(name)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    outname = outpath + 'SEDfit_2exp.lyric'
    standard = open('./phot.lyric', 'r')
    fobj=open(outname,'w')
    imgpath = '/path/to/your/photometric_data/'

    for liindex,line in enumerate(standard):
        if liindex == 4:
            fobj.write("R1) {0} \n".format(name))
        elif liindex == 5:
            fobj.write("R2) [{0},{1}] \n".format(ra,dec))
        elif liindex == 6:
            fobj.write("R3) {0} \n".format(z))
        elif liindex == 10:
            img = imgpath+ '{0}'.format(name)+Bands[0]+'.fits'
            fobj.write("Ia1) [{0},0] \n".format(img))
        elif liindex == 12:
            img = imgpath+ '{0}'.format(name)+Bands[0]+'.fits'
            fobj.write("Ia3) [{0},2] \n".format(img))
        elif liindex == 13:
            img = imgpath+ '{0}'.format(name)+Bands[0]+'.fits'
            fobj.write("Ia4) [{0},3] \n".format(img))
        # same for g band
        elif liindex == 27:
            img = imgpath+ '{0}'.format(name)+Bands[1]+'.fits'
            fobj.write("Ib1) [{0},0] \n".format(img))
        elif liindex == 29:
            img = imgpath+ '{0}'.format(name)+Bands[1]+'.fits'
            fobj.write("Ib3) [{0},2] \n".format(img))
        elif liindex == 30:
            img = imgpath+ '{0}'.format(name)+Bands[1]+'.fits'
            fobj.write("Ib4) [{0},3] \n".format(img))
        # same for r band
        elif liindex == 44:
            img = imgpath+ '{0}'.format(name)+Bands[2]+'.fits'
            fobj.write("Ic1) [{0},0] \n".format(img))
        elif liindex == 46:
            img = imgpath+ '{0}'.format(name)+Bands[2]+'.fits'
            fobj.write("Ic3) [{0},2] \n".format(img))
        elif liindex == 47:
            img =imgpath+ '{0}'.format(name)+Bands[2]+'.fits'
            fobj.write("Ic4) [{0},3] \n".format(img))
        elif liindex == 61:
            img = imgpath+ '{0}'.format(name)+Bands[3]+'.fits'
            fobj.write("Id1) [{0},0] \n".format(img))
        elif liindex == 63:
            img = imgpath+ '{0}'.format(name)+Bands[3]+'.fits'
            fobj.write("Id3) [{0},2] \n".format(img))
        elif liindex == 64:
            img = imgpath+ '{0}'.format(name)+Bands[3]+'.fits'
            fobj.write("Id4) [{0},3] \n".format(img))
        # same for z band
        elif liindex == 78:
            img = imgpath+ '{0}'.format(name)+Bands[4]+'.fits'
            fobj.write("Ie1) [{0},0] \n".format(img))
        elif liindex == 80:
            img = imgpath+ '{0}'.format(name)+Bands[4]+'.fits'
            fobj.write("Ie3) [{0},2] \n".format(img))
        elif liindex == 81:
            img = imgpath+ '{0}'.format(name)+Bands[4]+'.fits'
            fobj.write("Ie4) [{0},3] \n".format(img))
        # same for J band
        elif liindex == 95:
            img = imgpath+ '{0}'.format(name)+Bands[5]+'.fits'
            fobj.write("If1) [{0},0] \n".format(img))
        elif liindex == 97:
            img = imgpath+ '{0}'.format(name)+Bands[5]+'.fits'
            fobj.write("If3) [{0},2] \n".format(img))
        elif liindex == 98:
            img =imgpath+ '{0}'.format(name)+Bands[5]+'.fits'
            fobj.write("If4) [{0},3] \n".format(img))
        # same for H band
        elif liindex == 112:
            img = imgpath+ '{0}'.format(name)+Bands[6]+'.fits'
            fobj.write("Ig1) [{0},0] \n".format(img))
        elif liindex == 114:
            img = imgpath+ '{0}'.format(name)+Bands[6]+'.fits'
            fobj.write("Ig3) [{0},2] \n".format(img))
        elif liindex == 115:
            img = imgpath+ '{0}'.format(name)+Bands[6]+'.fits'
            fobj.write("Ig4) [{0},3] \n".format(img))
        # same for K band
        elif liindex == 129:
            img = imgpath+ '{0}'.format(name)+Bands[7]+'.fits'
            fobj.write("Ih1) [{0},0] \n".format(img))
        elif liindex == 131:
            img =  imgpath+ '{0}'.format(name)+Bands[7]+'.fits'
            fobj.write("Ih3) [{0},2] \n".format(img))
        elif liindex == 132:
            img =  imgpath+ '{0}'.format(name)+Bands[7]+'.fits'
            fobj.write("Ih4) [{0},3] \n".format(img))
        # same for fuv band
        elif liindex == 146:
            img =  imgpath+ '{0}'.format(name)+Bands[8]+'.fits'
            fobj.write("Ii1) [{0},0] \n".format(img))
        elif liindex == 148:
            img = imgpath+ '{0}'.format(name)+Bands[8]+'.fits'
            fobj.write("Ii3) [{0},2] \n".format(img))
        elif liindex == 149:
            img = imgpath+ '{0}'.format(name)+Bands[8]+'.fits'
            fobj.write("Ii4) [{0},3] \n".format(img))
        # same for nuv band
        elif liindex == 163:
            img = imgpath+ '{0}'.format(name)+Bands[9]+'.fits'
            fobj.write("Ij1) [{0},0] \n".format(img))
        elif liindex == 165:
            img = imgpath+ '{0}'.format(name)+Bands[9]+'.fits'
            fobj.write("Ij3) [{0},2] \n".format(img))
        elif liindex == 166:
            img = imgpath+ '{0}'.format(name)+Bands[9]+'.fits'
            fobj.write("Ij4) [{0},3] \n".format(img))
        # same for w1 band
        elif liindex == 180:
            img = imgpath+ '{0}'.format(name)+Bands[10]+'.fits'
            fobj.write("Ik1) [{0},0] \n".format(img))
        elif liindex == 182:
            img = imgpath+ '{0}'.format(name)+Bands[10]+'.fits'
            fobj.write("Ik3) [{0},2] \n".format(img))
        elif liindex == 183:
            img = imgpath+ '{0}'.format(name)+Bands[10]+'.fits'
            fobj.write("Ik4) [{0},3] \n".format(img))
        # same for w2 band
        elif liindex == 197:
            img = imgpath+ '{0}'.format(name)+Bands[11]+'.fits'
            fobj.write("Il1) [{0},0] \n".format(img))
        elif liindex == 199:
            img = imgpath+ '{0}'.format(name)+Bands[11]+'.fits'
            fobj.write("Il3) [{0},2] \n".format(img))
        elif liindex == 200:
            img = imgpath+ '{0}'.format(name)+Bands[11]+'.fits'
            fobj.write("Il4) [{0},3] \n".format(img))
        # same for w3 band
        elif liindex == 214:
            img = imgpath+ '{0}'.format(name)+Bands[12]+'.fits'
            fobj.write("Im1) [{0},0] \n".format(img))
        elif liindex == 216:
            img = imgpath+ '{0}'.format(name)+Bands[12]+'.fits'
            fobj.write("Im3) [{0},2] \n".format(img))
        elif liindex == 217:
            img = imgpath+ '{0}'.format(name)+Bands[12]+'.fits'
            fobj.write("Im4) [{0},3] \n".format(img))
        # same for w4 band
        elif liindex == 231:
            img = imgpath+ '{0}'.format(name)+Bands[13]+'.fits'
            fobj.write("In1) [{0},0] \n".format(img))
        elif liindex == 233:
            img = imgpath+ '{0}'.format(name)+Bands[13]+'.fits'
            fobj.write("In3) [{0},2] \n".format(img))
        elif liindex == 234:
            img = imgpath+ '{0}'.format(name)+Bands[13]+'.fits'
            fobj.write("In4) [{0},3] \n".format(img))
        # same for 250um band
        elif liindex == 248:
            img = imgpath+ '{0}'.format(name)+Bands[14]+'.fits'
            fobj.write("Io1) [{0},0] \n".format(img))
        elif liindex == 250:
            img = imgpath+ '{0}'.format(name)+Bands[14]+'.fits'
            fobj.write("Io3) [{0},2] \n".format(img))
        elif liindex == 251:
            img = imgpath+ '{0}'.format(name)+Bands[14]+'.fits'
            fobj.write("Io4) [{0},3] \n".format(img))
        # same for 350um band
        elif liindex == 265:
            img = imgpath+ '{0}'.format(name)+Bands[15]+'.fits'
            fobj.write("Ip1) [{0},0] \n".format(img))
        elif liindex == 267:
            img = imgpath+ '{0}'.format(name)+Bands[15]+'.fits'
            fobj.write("Ip3) [{0},2] \n".format(img))
        elif liindex == 268:
            img = imgpath+ '{0}'.format(name)+Bands[15]+'.fits'
            fobj.write("Ip4) [{0},3] \n".format(img))
        # same for 500um band
        elif liindex == 282:
            img = imgpath+ '{0}'.format(name)+Bands[16]+'.fits'
            fobj.write("Iq1) [{0},0] \n".format(img))
        elif liindex == 284:
            img = imgpath+ '{0}'.format(name)+Bands[16]+'.fits'
            fobj.write("Iq3) [{0},2] \n".format(img))
        elif liindex == 285:
            img = imgpath+ '{0}'.format(name)+Bands[16]+'.fits'
            fobj.write("Iq4) [{0},3] \n".format(img))


        elif liindex == 319:
            fobj.write("Pb10) [[{0},0.01,{1},0.1,1]] \n".format(0.5*age, age-0.2))
        elif liindex == 338:
            fobj.write("Ga3) [{0},{1},{2},0.01,0] \n".format(z,z-0.05,z+0.05))
        elif liindex == 339:
            fobj.write("Ga4) {0} \n".format(ebv.value[0]))
        else:
            fobj.write(line)
    fobj.close()
    #break
    loop += 1
