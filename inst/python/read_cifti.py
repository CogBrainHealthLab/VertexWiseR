"""Python function to read CIFTI files without needing to install connectome workbench
   Adapted and modified from https://nbviewer.org/github/neurohackademy/nh2020-curriculum/blob/master/we-nibabel-markiewicz/NiBabel.ipynb
"""

import nibabel as nib
import numpy as np

def read_cifti(file):
    '''
    file: filename of dtseries or dscalar file
    '''
    
    cifti = nib.load(file)
    cifti_hdr = cifti.header
    cifti_data = cifti.get_fdata(dtype=np.float32)
    axes = [cifti_hdr.get_axis(i) for i in range(cifti.ndim)]
    
    for name, data_indices, model in axes[1].iter_structures():  # Iterates over volumetric and surface structures
        if name == "CIFTI_STRUCTURE_CORTEX_LEFT":                                 # Just looking for a surface
            data = cifti_data.T[data_indices]                       # Assume brainmodels axis is last, move it to front
            vtx_indices = model.vertex                        # Generally 1-N, except medial wall vertices
            LH_surf_data = np.zeros((vtx_indices.max() + 1,) + data.shape[1:], dtype=data.dtype)
            LH_surf_data[vtx_indices] = data


    for name, data_indices, model in axes[1].iter_structures():  # Iterates over volumetric and surface structures
        if name == "CIFTI_STRUCTURE_CORTEX_RIGHT":                                 # Just looking for a surface
            data = cifti_data.T[data_indices]                       # Assume brainmodels axis is last, move it to front
            vtx_indices = model.vertex                        # Generally 1-N, except medial wall vertices
            RH_surf_data = np.zeros((vtx_indices.max() + 1,) + data.shape[1:], dtype=data.dtype)
            RH_surf_data[vtx_indices] = data

    return np.concatenate((LH_surf_data, RH_surf_data))
