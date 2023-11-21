# -----------------------------------------------------------
# Python script which retrieves the clipping values to use 
# on fetal brain subjects from a dataset based on the contrast
# of their white matter volume
#
# 2022-04-22 Andrés le Boeuf
# andres.le.boeuf@estudiantat.upc.edu
# -----------------------------------------------------------

from os import listdir
import nibabel as nib
import numpy as np
from sklearn.mixture import GaussianMixture


def variances(folder_path):
    """
    Gathers the variance from WM volumes

            inputs: - folder_path: Path were the cropped WM 
                                   reference volumes' folders are
            outputs:- variance: Python list with the standard 
                                deviation of WM intensity per subject
    """
    gmm = GaussianMixture(n_components=2) # GMM of two classes: (1) Black background; (2) WM Volume
    variance = []
    for folders in listdir(folder_path):
        path= folder_path + "\\" + folders 
        for files in listdir(path):
            if "WM" in files:
                WM_path = path + "\\" + files
                WM_data = nib.load(WM_path).get_fdata()
                gmm.fit(WM_data.reshape(-1,1))
                variance.append(gmm.covariances_[-1])
                print("Appended WM variance from ", files, " to the list. That is: ", gmm.covariances_[-1])
                
    return variance

def clipping_values(variance):
    """
    Returns the clipping values to use on the current dataset

            inputs: - variance: Python list with the standard 
                                deviation of WM intensity per subject
            outputs:- clip_values: Python list with the clipping values to use for each subject
    """
    var_aux = np.array(variance)
    init_clip_value = 0.25 # Define an initial clipping vaue for the first subject.
    var_perc = []
    clip_values = []
    clip_values.append(init_clip_value)

    for i in range(len(var_aux)-1):
        j = i+1
        perc = ((var_aux[j] - var_aux[i])/var_aux[i])*100
        var_perc.append(perc)

    var_perc = np.array(var_perc)
    var_perc[var_perc>40]=40 # To don't let the clipping values get stuck neaar to 0.
    var_perc = (-1)*var_perc

    clip_change = var_perc/100 + 1

    for i in range(len(clip_change)):
        clip_values.append(clip_values[-1]*clip_change[i])

    return clip_values

def main():
    folder_path = r'..\data\atlas_fast_clustering'   
    var = variances(folder_path)
    values = clipping_values(var)
    print("The clipping values for the spatiotemporal normative atlas are: ")
    for i in range(len(values)):
        print(i+21, ":\t", values[i])

if __name__ == "__main__":
    main()