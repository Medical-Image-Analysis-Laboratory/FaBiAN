White matter maturation processes implementation in FaBiAN simulator
===========================================================================
----

__Copyright (c) 2019-2023__  

__Medical Image Analysis Laboratory (MIAL), Lausanne, Switzerland__  
__Department of Radiology, Lausanne University Hospital (CHUV) and University of Lausanne (UNIL), Lausanne, Switzerland__  
__CIBM Center for Biomedical Imaging, Switzerland__

----
Magnetic resonance imaging is commonly used as a complement to the ultrasound gold-standard to investigate equivocal patterns in the developing fetal brain. However, stochastic motion of the fetus may alter the image quality and subsequent diagnostic. Several post-processing techniques have been developed to mitigate the resultant motion artefacts, but their evaluation requires a ground truth. FaBiAN is an open source Fetal Brain magnetic resonance Acquisition Numerical phantom that simulates clinical T2-weighted magnetic resonance images of the in utero fetal brain throughout maturation. However, this controlled environment was originally built on a three-class model of the fetal brain (white matter, gray matter and cerebrospinal fluid) which does not take into account key maturation processes that occur in the white matter throughout gestation. 
To overcome this limitation, we present in this work a new numerical fetal brain model that accounts for magnetic resonance signal changes in the white matter over development. We show that unsupervised segmentation methods such as Gaussian Mixture Models or FMRIB’s Automated Segmentation Tool are effective in identifying local variations of water content within white matter, thus making it possible to fine-tune the relaxometric properties of the tissues accordingly.

The code from this folder is adapted to work on the public Gholipour's normative spatiotemporal fetal brain [atlas](http://crl.med.harvard.edu/research/fetal_brain_atlas/). If you desire to change to another dataset you should check the white matter labels in *fsl_clustering.py* or *gmm_clustering.py* and make sure that all the paths which refers to the high-resolution reference model are adapted to your dataset.


>Below is described how the white matter implementation should be done:
>
> - First apply the the *White_Matter_improvement/fsl_clustering.py* or *White_Matter_improvement/gmm_clustering.py* script on the desired dataset. Each high resolution brain volume should be in a different folder along with its segmentation map. If you are going to  apply this script to a different dataset, make sure you have changed the white matter labels to recover a good cropping of the white matter volume.
> - Now that FAST has been applied, in each folder you will have available the high resolution brain volume, the segmentation map, the white matter volume and the outputs from FAST tool. But we are only interested in the partial volume maps which will be take as inputs for the FaBiAN simulator in *Utilities/WM_maturation.m* file.
> - If adaptive clipping values is desired for the current dataset, apply the *White_Matter_improvement/clipping_value_building.py* to (ONLY) the the white matter volumes of each subject. The result can be then pasted in *Utilities/clip_function.m*. An example of the clipping values adquisition is displayed through the jupyter notebook *histogram_distr_for_Gholipour.ipynb*.
>
Note: Make sure that the partial volume map paths are correct in *Utilities/WM_maturation.m*

## Contact
For any questions, please contact:\
Andrés le Boeuf\
andres.le.boeuf@estudiantat.upc.edu

OR

Hélène Lajous  
helene.lajous@unil.ch 
