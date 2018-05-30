evosoro: soft robot simulator
=======================================

<div class="row">
<a href=https://youtu.be/EXuR_soDnFo>
<img src="https://github.com/skriegman/img/blob/master/nick.png" height="135" width="135">
</a>

<a href=https://youtu.be/HgWQ-gPIvt4>
<img src="https://github.com/skriegman/img/blob/master/electro.png" height="135" width="135">
</a>

<a href=https://youtu.be/4ZqdvYrZ3ro>
<img src="https://github.com/skriegman/img/blob/master/swimming.png" height="135" width="135">
</a>

<a href=https://youtu.be/Cw2SwPNwcfM>
<img src="https://github.com/skriegman/img/blob/master/plant1.png" height="135" width="135">
</a>

<a href=https://youtu.be/XqIUJcuOgmw>
<img src="https://github.com/skriegman/img/blob/master/teeth1.png" height="135" width="135">
</a>

<a href=https://youtu.be/r_SL8VUt-wA>
<img src="https://github.com/skriegman/img/blob/master/cage.png" height="135" width="135">
</a>

</div>

Evosoro is a Python soft robot simulation library based on the Voxelyze physics engine. It provides a high-level interface for the dynamic simulation and automated design of soft multimaterial robots.
<!-- evolutionary design of soft multimaterial robots. -->

Evosoro was designed and developed by the [Morphology, Evolution & Cognition Laboratory](http://www.meclab.org), University of Vermont. 
The library is built on top of the open source [VoxCAD](https://github.com/jonhiller/VoxCAD
) and the underlying voxel physics engine ([Voxelyze](https://github.com/jonhiller/Voxelyze)) which were both developed by the [Creative Machines Lab](http://www.creativemachineslab.com/), Columbia University.



(1) Citing
------

If using this code for academic purposes please consider citing the following two papers.

Voxelyze physics engine and VoxCAD GUI:

>Hiller, J., & Lipson, H. (2014). 
>*Dynamic simulation of soft multimaterial 3d-printed objects.*
>Soft Robotics, 1(1), 88-101.

Encoding and optimization (genetic algorithm):

>Cheney, N., MacCurdy, R., Clune, J., & Lipson, H. (2013). 
>*Unshackling evolution: evolving soft robots with multiple materials and a powerful generative encoding.* 
>In Proceedings of the 15th annual conference on Genetic and evolutionary computation (pp. 167-174). ACM.



(2) Installation
------------

It is recommended that you install [Anaconda](https://docs.continuum.io/anaconda/install#) as your Python (2.7) distribution. Anaconda is a free package manager and Python distribution that includes all of the dependencies required for evosoro. However if you instead choose to manually install Python 2.7, the following packages are required: scipy, numpy, networkx, decorator.

Also networkx must be <2.0. When networkx updated 1.0-->2.0 some function changed and I haven't updated the python code to reflect this change.

    pip install networkx==1.11


Install Qt and QMake if you have not already done so, specifically these packages: "libqt4-dev", "qt4-qmake", "libqwt-dev", "freeglut3-dev" and "zlib1g-dev".

    sudo apt-get install libqt4-dev qt4-qmake libqwt-dev freeglut3-dev zlib1g-dev


Install git if you have not already done so.

    sudo apt-get install git

Navigate to your working directory (e.g. your home).

    cd ~

Clone the repo.

    git clone https://github.com/skriegman/evosoro.git

There are different well documented examples (evosoro/examples) and custom versions of VoxCad/Voxelyze included in this repository (evosoro/_voxcad* folders).
Let's try running an example in which soft robots are optimized to locomote in a terrestrial environment, using an evolutionary algorithm and a basic version (_voxcad) of the physics engine (the procedure is the same for all the examples). 

Navigate to the _voxcad directory:

    cd evosoro/evosoro/_voxcad/

The following command compiles both VoxCad and Voxelyze, installing the library at the same time:

    ./rebuild_everything.sh

If you happen to modify VoxCad or Voxelyze in the future, you can call the same script to be sure to clean and recompile everything. 

    make

Install the voxelyze library.

    cd Voxelyze
    make
    make installusr
    cd ../voxelyzeMain/
    make

Navigate back out to the examples folder and run basic.py
    
    cd ../examples
    python basic.py

You should start seeing some output being produced in your console, and a new directory being created (evosoro/evosoro/basic_data), which contains the results of the simulation.



(3) Examples
--------

After allowing basic.py to run for a few generations, you can view the evolved morphologies and behaviors by opening up the generated .vxa files within the VoxCAD GUI. A .vxa file is just an XML file representing a robot that can be simulated by VoxCad/Voxelyze. Different versions of the physics engine can play slightly different .vxa files.
Navigate to evosoro/evosoro/_voxcad/release:
    
    cd ../_voxcad/release
    
Open VoxCad:

    ./VoxCad

Then select the desired .vxa file from 

    "File -> Import -> Simulation"

The .vxa files for the best performing individuals will be saved in 

    evosoro/evosoro/examples/basic_data/bestSoFar/fitOnly.

Once the design is loaded, you can start the physics simulation by clicking the <img src="https://github.com/skriegman/evosoro/blob/master/evosoro/_voxcad/VoxCad/Icons/Sandbox.png" height="25" width="25"> icon in the top bar ("Physics Sandbox").  The robot should start moving: if it doesn't, please check the following section (Known issues).


(4) Known issues
--------

Some versions of VoxCad utilize the qhull library (http://www.qhull.org/). In these versions of the GUI, the 'qhull' executable must be present in your binary path.

    sudo apt-get install qhull-bin

If the robot does not move, disappears, or seems to behave in an unexpected manner when running a .vxa file in VoxCad (GUI), you may be affected by a known problem observed on some non-US machines.
The problem is due to an unexpected behavior of the <a href="http://www.cplusplus.com/reference/cstdlib/atof/">atof</a> function when the system's numeric <a href="https://en.wikipedia.org/wiki/Locale_(computer_software)">locale</a> differs from en_US.UTF-8, which entails loading wrong parameters from the .vxa file (in some cases it was observed how the atof function was approximating all double and floating point values to their integer part, which was the cause of the unexpected behavior).

While we work on a better solution, you can fix this problem by making sure that your machine is configured according to a US numeric locale.
Open the following file:

    sudo gedit /etc/default/locale

Make sure that LC_NUMERIC is set as follows:

    LC_NUMERIC="en_US.UTF-8"

Save, close the file, and reboot.


(5) Documentation
-------------

Although the code included in this repository diverged from the main VoxCad/Voxelyze development branch some time ago, useful indications could be find in the online Voxelyze documentation, available [here](http://jonhiller.github.io/Voxelyze/annotated.html).


(6) License
-------

Released under a MIT License (MIT)


