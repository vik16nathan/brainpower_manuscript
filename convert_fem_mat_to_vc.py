#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 23:13:24 2024

@author: vikramn
"""

labels=[('brain', "GM + WM + CSF"), ('skull', 'skull'), ('scalp', 'scalp')]
import numpy as np
import pandas as pd


subject_list = ['sub-0002','sub-0004', 'sub-0006', 'sub-0007']
#folder="C:/Users/vik16/OneDrive/Documents/Baillet Lab/"
folder='/export02/data/vikramn/brainstorm3/'
output_folder='./'
for subject in subject_list:
    # Convert MATLAB data to NumPy arrays
    fem_numpy = {}
    fem_numpy['nodes'] = np.array(pd.read_csv(folder+subject+'_nodes.csv', header=None))
    fem_numpy['elements'] = np.array(pd.read_csv(folder+subject+'_elements.csv', header=None)) - 1
    fem_numpy['labels'] = np.array(pd.read_csv(folder+subject+'_labels.csv', header=None)).flatten() - 1; #switch from matlab 1-based to python 0-based indexing
    fem_numpy['conductivities'] = [0.33, 0.0042, 0.33]
    fem_numpy['label_info'] = [(0, 'brain'), (1, 'skull'), (2, 'scalp')]
    np.savez(output_folder+subject+'_vc_3layer.npz', nodes=fem_numpy['nodes'], 
             elements=fem_numpy['elements'], labels=fem_numpy['labels'], label_info=fem_numpy['label_info'], 
             conductivities=fem_numpy['conductivities'])

