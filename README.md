# AASP-Graph

Source code for the paper "An Adaptive Affinity Graph with Subspace Pursuit for Natural Image Segmentation". The AASP-Graph paper can be found [here](https://ieeexplore.ieee.org/document/8784904).

AASP-Graph is modified from [the offcial GL-Graph implementation](https://github.com/xiaofanglegoc/global-local-affinity-graph).


### Requirements
The code requires the version of Matlab2018a, Ubuntu 16.04.


### Data
The BSDS300 dataset and predata (extracted features) can be downloaded from [GL-Graph](https://github.com/xiaofanglegoc/global-local-affinity-graph) and place them in `BSD/` and `bsd_300_feat/` folder, respectively.


### Demo
Run the demo `demo_AASP_BSDS300.m`.


### Results
The detailed results can be found in `BSDS300_evaluation.txt`. More relevant results of natural image segmentation can be found in our [projects](https://github.com/Yangzhangcst/Natural-color-image-segmentation).


### Citing
If you find this repository useful in your research, please consider citing:
```
@INPROCEEDINGS{AASP-Graph,  
  author={Y. {Zhang} and H. {Zhang} and Y. {Guo} and K. {Lin} and J. {He}},  
  booktitle={IEEE International Conference on Multimedia and Expo (ICME)},   
  title={An Adaptive Affinity Graph with Subspace Pursuit for Natural Image Segmentation},   
  year={2019}，  
  pages={802-807},}
```
