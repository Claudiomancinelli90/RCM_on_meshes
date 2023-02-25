# Computing the Riemannian Center of Mass on Meshes

This repository contains the implementation of the algorithm described in "Computing the Riemannian Center of Mass on Meshes", C. Mancinelli and E. Puppo. The code consists of one GUI (RCM_GUI) supporting all the algorithms for curve tracing described in the paper.

## Compilation
You can compile the project using the standard CMake routine

`mkdir build`<br/>
`cd build`<br/>
`cmake ../`<br/>
`make`<br/>

For mac users, there is also the possibility of building the Xcode project by compiling a python script. Namely

`python3 scripts/build.py xcode`

This will generate the Xcode project "RCM_on_meshes" in the repository build/xcode.

## Dependencies
Mac users only need to have C++ 17 or later installed. Linux users may have to install some additional packages, as

`libxrandr-dev`<br/>
`libxinerama-dev`<br/>
`libxcursor-dev`<br/>
`libxi-dev.`<br/>

Similarly, Windows users may need to install some OS-related packages, 
but there are no dependencies that prevent the use of the code on such operating system. 

## Run
Once the projects are built, you can run the app from within the project repository by issuing

`./bin/RCM_GUI <path_to_a_model>`


## Using the GUI

                                 ---Basics---
Rotate: SHIFT + dragging.
Panning: ALT + dragging.                                 
                                 ---Control Points Selection/Editing---
Control points can be picked with left click and moved on the mesh by dragging them with right click

                                ---Method for geodesic distances---
The method used to compute geodesic distances can be selected with the combobox "Method". The 
available choices are the one analyzed in the paper, i.e VTP and a graph based geodesic solver.

                                ---Curve Tracing---
Once the control points are picked, you can choose which curve you want to trace: rational Bézier, rational B-spline or 
Bézier Interpolants (see paper). By default, all the weights are set to 1, but you can edit them in the section "Parameters": choose the weight you want to edit, and then use the slider to set its value. The curve will be udpated in real time.

If you choose just 1 control point, you can trace a circle made of Rational Bézier segment.






                         






