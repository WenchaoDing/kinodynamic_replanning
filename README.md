# kinodynamic_replanning
Materials for the paper "Trajectory Replanning for Quadrotors Using Kinodynamic Search and Elastic Optimization", which is submitted to ICRA2018.
## File Discription
* contain_linear_segments.nb

The file is to solve the max-min problem in Appendix C in closed form, which is written in [Mathematica](https://www.wolfram.com/mathematica/) and has been tested in Mathematica 10.0.

* rbk_inflation.m

The file is to numerically check the inflation needed to guarantee the RBK trajectory to be collision-free, which is written in Matlab and has been tested in Matlab R2016b.

* rbk_inflation_3d.m

The file is to show the 3d case, which provides the complementary results. 

* permn.m

The file is to enumerate all the patterns, which is from [Jos@Mathworks](https://www.mathworks.com/matlabcentral/fileexchange/7147-permn-v--n--k-).

* grid2coord.m

Some basic coordinate transformation. 

* print.pdf

For the users who do not have [Mathematica](https://www.wolfram.com/mathematica/) installed, we print the script and results from [Mathematica](https://www.wolfram.com/mathematica/) for them.


## Authors
* Wenchao DING
* William Wu

HKUST Robotics Institute
