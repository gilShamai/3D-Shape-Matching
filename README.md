# 3D-Shape-Matching
A Matlab implementation of the paper "geodesic distance descriptors" for 3D non-rigid shape correspondence. 

[[Paper 1]](https://ieeexplore.ieee.org/abstract/document/8509134) [[Paper 2]](https://docs.wixstatic.com/ugd/28cd82_91f41197b793480ab25b1f97f10f818a.pdf) [[Main paper]](https://docs.wixstatic.com/ugd/28cd82_bb48e8cf06984a18b6016997beda5e4f.pdf)

The code works on 3D triangle meshes. It can be modified to arbitrary graphs, e.g. creating graphs from 3D cloud points but connecting near neighbors.

<img align="right" img src="Images/Viz_pic.png" width="400px">

## Assumptions
(1) Deformation between shapes doesn't have to be rigid.
(2) Deformation between shapes should be nearly isometry (i.e. one shape is some nearly elastic bending of the other one).
(3) Shapes can have different number of points.
(4) Shapes must have the same topology (cannot deal with partial correspondence, or a hole in one shape that is not in the other one).

## Setup & Usage
The project was tested on OSX with Matlab R2019a, and should work on windows 64 as well (maybe need a few modifications of paths), to run it: 
1) Download files
2) unzip ann_mwrapper.zip
3) run build_mex.m if needed inside laplace_beltrami folder
4) run ann_compile_mex.m if needed inside ann_mwrapper folder
5) Run DEMO.m


## Citation
If you use these ideas, please cite the paper <a href="https://docs.wixstatic.com/ugd/28cd82_bb48e8cf06984a18b6016997beda5e4f.pdf"> Geodesic Distance Descriptors </a>. This paper is using the pairwise geodesic distances. For this code Please cite also and (2) <a href="https://ieeexplore.ieee.org/abstract/document/8509134"> Efficient Inter-Geodesic Distance Computation and Fast Classical Scaling</a> and (3) <a href="https://docs.wixstatic.com/ugd/28cd82_91f41197b793480ab25b1f97f10f818a.pdf"> Accelerating the computation of canonical forms for 3D nonrigid objects using Multidimensional Scaling </a>:

```
@inproceedings{shamai2017geodesic,
  title={Geodesic distance descriptors},
  author={Shamai, Gil and Kimmel, Ron},
  booktitle={Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition},
  pages={6410--6418},
  year={2017}
}
```

```
@article{shamai2018efficient,
  title={Efficient Inter-Geodesic Distance Computation and Fast Classical Scaling},
  author={Shamai, Gil and Zibulevsky, Michael and Kimmel, Ron},
  journal={IEEE transactions on pattern analysis and machine intelligence},
  year={2018},
  publisher={IEEE}
}
```

```
@inproceedings{shamai2015accelerating,
  title={Accelerating the computation of canonical forms for 3D nonrigid objects using multidimensional scaling},
  author={Shamai, Gil and Zibulevsky, Michael and Kimmel, Ron},
  booktitle={Proceedings of the 2015 Eurographics Workshop on 3D Object Retrieval},
  pages={71--78},
  year={2015},
  organization={Eurographics Association}
}
```
