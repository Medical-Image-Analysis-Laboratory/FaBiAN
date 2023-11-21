# -----------------------------------------------------------
# Python script to crop the white matter volume from the 
# high resolution reference model and apply on it a 2 classes 
# unsupervised segmentation through Gaussian Mixture Models.
#
# 2022-03-23 Andrés le Boeuf
# andres.le.boeuf@estudiantat.upc.edu
# -----------------------------------------------------------

import os
from os import listdir
import numpy as np
from sklearn.mixture import GaussianMixture
import nibabel as nib

def _gmm(seg_path, path, K, covariance='full', warm_start=True):
    """
    This function apply GMM algorithm to the image forwarded as 
    its path with K classes. Saves the clustered nifti image
    and the posteriori probability maps of each class.

            inputs: - seg_path: Path to the white matter volume that we want to segment
                    - path: Folders main path to save the segmentation and posteriori probability maps
                    - K: Number of classes
                    - covariance: - full: each component has its own general covariance matrix.
                                  - Other options see sklearn.mixture.GaussianMixture documentation
                    - warm_start: If is True, the solution of the last fitting is used as 
                                  initialization for the next call of fit(). Can speed up the convergence

    """
    Image_nii = nib.load(seg_path)
    Image = Image_nii.get_fdata()
    
    #Converts input data into 1D array (eack row = data point, each column (1) = feature (intensity))
    if(len(Image.shape)<3):
        Z = Image.reshape((-1,1))
    elif (len(Image.shape) == 3 and Image.shape[2]==3):
        Z = Image.reshape((-1,3)) 
    elif (len(Image.shape) == 3 and Image.shape[2]>3):
        Z = Image.reshape((-1,1))
    
    #Creates GMM model
    gmm = GaussianMixture(n_components=K, covariance_type=covariance, init_params = 'kmeans',  n_init = 1, warm_start=warm_start) 
    #Trains the GMM model
    gmm.fit(Z)
    #Sort means of each class to match labels between segmentations and with FAST method
    gmm.means_ = _sort(gmm.means_)
    #Applies the predictions from the model to the image
    labels = gmm.predict(Z)
    prob_maps = gmm.predict_proba(Z)
    cluster = labels.reshape(Image.shape)
    print(gmm.means_)
    clustered_image = nib.Nifti1Image(dataobj=cluster, affine=Image_nii.affine)
    clustered_image._header = Image_nii.header

    gmm_path = path + "\\" + path.split("\\")[-1] + "_gmm.nii"
    print("Saving images into " + path)
    nib.save(img=clustered_image, filename=gmm_path)

    for i in range(prob_maps.shape[-1]):
        map_path = path + "\\" + path.split("\\")[-1] + "_map" + str(i) + ".nii"
        prob_map = nib.Nifti1Image(dataobj=prob_maps[:,i].reshape(Image.shape), affine=Image_nii.affine)
        prob_map._header = Image_nii.header
        nib.save(img=prob_map, filename=map_path)

#Extracts the WM from all the nifti images allocated into Gholipour or FeTA Dataset.
def extract_WM(folder_path):
    """
    Extracts the white matter from all the nifti 
    images allocated into the dataset.

            inputs: - folder_path: Path were reference 
                                   volumes' folders are
    """
    segmentation_name = "tissue"                                                                    
    wm =  np.array([91, 110, 111, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 125]) #This WM label selection was performed through the annotations.
                                                                                          #without cerebellum and brainstem

    for folders in listdir(folder_path):
        path= folder_path + "\\" + folders 
        for files in listdir(path):
            if segmentation_name in files:
                
                tissue_path = path + "\\" + files
                replacement = "_" + segmentation_name
                brain_path = tissue_path.replace(replacement, "")
                print("Extracting WM from " + brain_path.split("\\")[-1] )

                brain = nib.load(brain_path).get_fdata()
                tissues_nii = nib.load(tissue_path)
                tissues = tissues_nii.get_fdata()
                
                tissues[~(np.isin(tissues,wm))] = 0 #Every voxel non classified as white matter will not be considered
                tissues[tissues!=0] = 1
                wm_tissue = brain*tissues #WM extraction

                wm_image = nib.Nifti1Image(dataobj=wm_tissue, affine=tissues_nii.affine)
                wm_image._header = tissues_nii.header
                save_path = path + "\\" + path.split("\\")[-1] + "_WM.nii"
                if os.path.exists(save_path):
                    os.remove(save_path)
                nib.save(img=wm_image,filename=save_path)  

def segment(folder_path, classes):
    """
    Applies GMM 2 classes segmentation
    to each extracted WM nifti image from a dataset

            inputs: - folder_path: Path were reference 
                                   volumes' folders are
    """
    for folders in listdir(folder_path):
        path= folder_path + "\\" + folders 
        for files in listdir(path):
            if "WM" in files:
                seg_path = path + "\\" + files 

        print("Performing GMM clustering on " + seg_path.split("\\")[-1])
        _gmm(seg_path, path, classes, covariance='full', warm_start=True)
        
def _sort(array):
    """
    Sorts an array to match the order of labels between segmentations and with FAST method

            inputs: - array: Array of mean values of each class.
            
            outputs: - rev_sorted: Array of mean values sorted in the desired order.
    """
    rev_sorted = np.sort(array, axis=0)[::-1]
    rev_sorted=rev_sorted[rev_sorted!=0]
    rev_sorted = np.insert(rev_sorted, [0],0).reshape(-1,1)
    
    return rev_sorted
def main():
    #Define paths to dataset
    folder_path = r'..\data\atlas_gmm_clustering'
    
    #Extract the WM from each subject 
    extract_WM(folder_path)

    #Segments the WM of each subject 
    segment(folder_path, 3) #2 classes for WM and one for background

if __name__ == "__main__":
    main()