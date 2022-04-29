#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 17:44:36 2022

@author: aaa
"""
import os
import shutil
import numpy as np
from PIL import Image as im
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
from scipy.stats import chisquare as cs
import math
import pandas as pd

sorted_fold_path="/home/aaa/Desktop/Thesis/Script/Trial/Sorted data/" #insert folder of sorted meausements files
allmeasurements = sorted_fold_path+"All measurements/"
allrenamed = allmeasurements +"All renamed/"
allmatrixes = allmeasurements + "All matrixes/"
allpictures = allmeasurements + "All pictures/"
allrawpictures = allmeasurements + "All raw pictures/"
alldata_analysis = allmeasurements + "All Data Analysis/"
allcropped_pictures = alldata_analysis + "All Cropped Pictures/"
allcontrolplots = alldata_analysis + "All Control plots/"
allcontrolfits = alldata_analysis + "All Control Fits/"
tiltangles=[0,40,48,61,69,71,79,80,81]
foldername=[]
for i in range(len(tiltangles)):
    foldername.append(str(tiltangles[i])+"deg")
foldername.append("79-4U_77c88deg")
foldername.append("79-8U_76c76deg")
foldername.append("79-12U_75c64deg")
foldername.append("79-16U_74c52deg")
n_theta=[26,46,28,17,16,20,21,20,19,48,43,59,24]  #number of measurements files for each folder (no flat, no compromised data)
n_pixel = 16384 #number of pixels in one measurement

for k in range(len(foldername)):
    bckg=0.
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
    stack = np.stack(matrixes,axis=2)
    xyzabsmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
    yabsmax = xyzabsmax[0][0]
    xabsmax = xyzabsmax[1][0]
    zabsmax = xyzabsmax[2][0]
    roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI+Peaks.mpa',skiprows=1).astype(int)
    print(foldername[k])
    data_and_fit = np.zeros((len(stack[0,0,:])*len(roi[:,0]), 21))
    xmin= int(np.amin(roi[:,1]))
    xmax= int(np.amax(roi[:,2]))
    noise=np.zeros(len(stack[0,0,:]))
    if(len(stack[0,0,:])>25):
        for z in range(len(stack[0,0,:])):
            value=0
            for y in range(len(roi[:][0])):
                if(k==6 and (z==3 or z==4)):
                    continue
                value+=sum(stack[roi[y][0],xmin:xmax,z])
            noise[z] = sum(stack[yabsmax,xmin:xmax,z])
            print(noise) 
        if(k==6):
            noise=np.delete(noise,np.where(noise==0)[0])
            plt.hist(noise[0:-2], bins=4,rwidth=0.95)
        else:
            plt.hist(noise, bins=5,rwidth=0.95)
            plt.title(foldername[k])
            plt.show()