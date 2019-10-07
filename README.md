# pulsed_qMT
This site is currently under construction. This software can be used for analysing qMT data using the single point qMT method originally outlined in [Yarnykh's Paper](https://www.ncbi.nlm.nih.gov/pubmed/22190042) and modified for use in the [spinal cord](https://www.ncbi.nlm.nih.gov/pubmed/24632465) and [Parkinson's disease](https://www.ncbi.nlm.nih.gov/pubmed/28986653).

## Requirements:

* [NIFTI Tools](https://uk.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)

## Running in the Spinal Cord:

I've found the following scripts work well when processing data in the spinal cord. To run them, you'll need the [Spinal Cord Toolbox](https://github.com/neuropoly/spinalcordtoolbox):

### Register MT-weighted Images to MEDIC/mcFFE:
  sct_register_multimodal -i {MT_Image} -d {MEDIC} -param step=1,type=im,algo=slicereg,metric=MeanSquares:step=2,type=im,algo=bsplinesyn,metric=MeanSquares,iter=5,shrink=2
### B0 Image to MEDIC
  sct_register_multimodal -i {B0_Image} -d {MEDIC} -param step=1,type=im,algo=slicereg,metric=MeanSquares
### B1 Image to MEDIC
  sct_register_multimodal -i {B1ref_Image} -d {MEDIC} -param step=1,type=im,algo=slicereg,metric=MeanSquares
### Create mask
  sct_create_mask -i {MEDIC} -p center
### Apply transforms
  sct_apply_transfo -i {B1_Map} -d {MEDIC} -w warp_MTw2MEDIC.nii.gz 