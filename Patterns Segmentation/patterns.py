# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
def clean_img_exp(original):
    img = ndimage.binary_fill_holes(original, structure=np.ones((60, 60)))
    img = ndimage.binary_fill_holes(img, structure=np.ones((30, 30)))
    
    return img

def clean_img_sim(original):
    img = ndimage.binary_opening(original)
    img = img.astype(int)
    dilate = ndimage.grey_dilation(img, size = (2,2), structure=np.ones((2, 2)))
    
    return (2 - dilate)

def count_regions(img):
    distance = ndimage.distance_transform_edt(img)
    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((100, 100)), labels=img)
    markers = ndimage.label(local_maxi)[0]
    labels = watershed(-distance, markers, mask = img)
    
    #plt.figure()
    #plt.imshow(labels, cmap='gist_ncar')
    
    #plt.figure()
    #plt.imshow(distance, cmap='gist_ncar')
    return labels.max()

def count_regions_alt(img, ind):
    distance = ndimage.distance_transform_edt(img)
    arr = np.asarray(distance).reshape(-1)
    (y,x,patches) = plt.hist(arr)
    mask = (distance > x[ind])
    label_im, nb_labels = ndimage.label(mask)
    
    #plt.figure()
    #plt.imshow(label_im, cmap='gist_ncar')
    
    #plt.figure()
    #plt.imshow(distance, cmap='gist_ncar')
    return nb_labels

import numpy as np
import csv
from skimage.morphology import watershed
from skimage.feature import peak_local_max
import matplotlib.pyplot as plt
import matplotlib.image as mpimg 
#from skimage.measure import label, regionprops
from scipy import ndimage
import os

path = os.getcwd()
simulation = []
sim_type = []
sim_holes = []
sim_area = []

experiment = []
exp_type = []
exp_holes = []
exp_area = []

n = 1000

for root, dirs, files in os.walk(path):
    for name in files:
        #print(name)
        if name.endswith('png'):
            img = mpimg.imread(os.path.join(root, name)) 
            w,h = img.shape
            #cimg = clean_img_sim(1 - img)
            simulation.append(img)
            sim_holes.append(count_regions_alt(img, 1))
            mnozh = 100.0/w**2
            sim_area.append(mnozh*np.sum(img))
        
        
        if name.endswith('tif'):
            img = mpimg.imread(os.path.join(root, name))
            w,h = img.shape
            #print(w,h)
            for i in range(0,w - n,200):
                for j in range(0,h - n,100):
                    subimg = img[i:i+n, j:j+n]
                    subimg = np.where(subimg > 254, 1, 0)
                    cimg = clean_img_exp(subimg)
                    experiment.append(cimg)
                    exp_holes.append(count_regions_alt(cimg, 2))
                    mnozh = 100.0/n**2
                    exp_area.append(mnozh*np.sum(cimg))
                    
                    if 'H1152' in name:
                        if 'WT' in name:
                            exp_type.append('0R')
                        if 'CCM1' in name:
                            exp_type.append('1R')
                        if 'CCM2' in name:
                            exp_type.append('2R')
                        if 'CCM3' in name:
                            exp_type.append('3R')
                    else:
                        if 'WT' in name:
                            exp_type.append('0')
                        if 'CCM1' in name:
                            exp_type.append('1')
                        if 'CCM2' in name:
                            exp_type.append('2')
                        if 'CCM3' in name:
                            exp_type.append('3')
                  
