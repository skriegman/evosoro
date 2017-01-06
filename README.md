evosoro: soft robot simulator
=======================================

<div class="row">
<a href=https://youtu.be/EXuR_soDnFo>
<img src="https://github.com/skriegman/img/blob/master/nick.png" height="135" width="135">
</a>

<a href=https://www.youtube.com/watch?v=HgWQ-gPIvt4>
<img src="https://github.com/skriegman/img/blob/master/electro.png" height="135" width="135">
</a>

<a href=https://youtu.be/4ZqdvYrZ3ro>
<img src="https://github.com/skriegman/img/blob/master/swimming.png" height="135" width="135">
</a>

<a href=https://youtu.be/Cw2SwPNwcfM>
<img src="https://github.com/skriegman/img/blob/master/plant1.png" height="135" width="135">
</a>

<a href=https://www.youtube.com/watch?v=XqIUJcuOgmwl>
<img src="https://github.com/skriegman/img/blob/master/teeth1.png" height="135" width="135">
</a>

<a href=https://www.youtube.com/watch?v=r_SL8VUt-wA>
<img src="https://github.com/skriegman/img/blob/master/cage.png" height="135" width="135">
</a>

</div>

Evosoro is a Python soft robot simulation library based on the Voxelyze physics engine. It provides a high-level interface for the dynamic simulation and automated design of soft multimaterial robots.
<!-- evolutionary design of soft multimaterial robots. -->

Evosoro was designed and developed by the [Morphology, Evolution & Cognition Laboratory](http://www.meclab.org), University of Vermont. 
The library is built on top of the open source [VoxCAD](https://github.com/jonhiller/VoxCAD
) and the underlying voxel physics engine ([Voxelyze](https://github.com/jonhiller/Voxelyze)) which were both developed by the [Creative Machines Lab](http://www.creativemachineslab.com/), Columbia University.



<!--
<a href="https://skriegman.github.io/" target="_blank">Sam Kriegman</a>, <a href="http://sssa.bioroboticsinstitute.it/user/1507" target="_blank">Francesco Corucci</a> and <a href="http://www.ncheney.com/" target="_blank">Nick Cheney</a> at 
the <a href="http://www.meclab.org" target="_blank">Morphology, Evolution & Cognition Laboratory</a>, 
<a href="http://www.uvm.edu/~cmplxsys/" target="_blank">Vermont Complex Systems Center</a>,
University of Vermont, USA.
-->

1. Citing
------

A full list of papers which contributed to the development of this code base may be found in Section 6. 
If using this code for academic purposes please cite the following two papers and consider citing any relevant publications from Section 6. 

Physical simulation:

>Hiller, J., & Lipson, H. (2014). 
>*Dynamic simulation of soft multimaterial 3d-printed objects.*
>Soft Robotics, 1(1), 88-101.

Encoding and optimization:

>Cheney, N., MacCurdy, R., Clune, J., & Lipson, H. (2013). 
>*Unshackling evolution: evolving soft robots with multiple materials and a powerful generative encoding.* 
>In Proceedings of the 15th annual conference on Genetic and evolutionary computation (pp. 167-174). ACM.

When using specific versions of the physics engine (as indicated by VOXELYZE_VERSION in the Python examples), please check the accompanying ReadMe.txt (within the relevant "_voxcad_*" folder) for relevant additional references to be cited.


<!--
Dependencies
------------

- Python 2.7

### Mandatory

- [numpy](http://www.numpy.org/)

- [scipy](http://www.scipy.org/)

- [networkx](http://networkx.github.io/)

### Recommended

- [pandas](http://pandas.pydata.org/)

- [matplotlib](http://matplotlib.org/)

- [seaborn](http://seaborn.pydata.org/)
-->


2. Installation
------------

<!--To install the released version, just do-->
    
<!--    pip install evosoro-->

<!--You may instead want to use the development version from Github, by running-->

<!--    pip install git+git://github.com/skriegman/evosoro.git#egg=evosoro-->

It is recommended that you install [Anaconda](https://docs.continuum.io/anaconda/install#) as your Python (2.7) distribution. Anaconda is a free package manager and Python distribution that includes all of the dependencies required for evosoro. However if you instead choose to manually install Python 2.7,

    sudo apt-get install python-dev python-pip
    sudo pip install scipy numpy networkx decorator  

If you experience an error installing scipy it might be due to the incompatibility of different fortran compilers (see [scipy installation](https://docs.scipy.org/doc/numpy-1.10.1/user/install.html)). In most cases, you must build numpy/scipy with the same fortran compiler used to build blas/lapack/atlas. Try the following:

	sudo apt-get install libatlas-base-dev gfortran
	sudo pip install scipy numpy networkx decorator


Install Qt and QMake if you have not already done so, specifically these packages: "libqt4-dev", "qt4-qmake", "libqwt-dev", "freeglut3-dev" and "zlib1g-dev".

    sudo apt-get install libqt4-dev qt4-qmake libqwt-dev freeglut3-dev zlib1g-dev


Install git if you have not already done so.

    sudo apt-get install git

Navigate to your working directory (e.g. your home).

    cd ~

Clone the repo.

    git clone https://github.com/skriegman/evosoro.git

There are different well documented examples (evosoro/examples) and custom versions of VoxCad/Voxelyze included in this repository (evosoro/_voxcad folders, see the corresponding readme file).
Let's try running an example in which soft robots are optimized to locomote in a terrestrial environment, using an evolutionary algorithm and a basic version (_voxcad) of the physics engine (the procedure is the same for all the examples). 

Navigate to the _voxcad directory:

    cd evosoro/evosoro/_voxcad/

The following command compiles both VoxCad and Voxelyze, installing the library at the same time:

    ./rebuild_everything.sh

If you happen to modify VoxCad or Voxelyze in the future, you can call the same script to be sure to clean and recompile everything. 

<!--
    make

Install the voxelyze library.

    cd Voxelyze
    make
    make installusr
    cd ../voxelyzeMain/
    make -->

Navigate back out to the examples folder and run basic.py
    
    cd ../../examples
    python basic.py

You should start seeing some output being produced in your console, and a new directory being created (evosoro/evosoro/basic_data), which contains the results of the simulation.
    
<!--
------------------------------------
**Installing from scratch on a virtual machine (graphics may not work properly)**
- Install the latest version of VirtualBox
- Download Ubuntu 14.04 x64 (ISO) and install it
- Install VirtualBox Guest Additions (openGL support and other useful features)
- Then run the virtual machine, and follow the instructions below
..*a ready to use VirtualBox image is available

Increasing the video memory could be useful: go in your VirtualBox installation folder and run:

    vboxmanage modifyvm "VIRTUAL_MACHINE_NAME" --vram 256

If you are having difficulty building scipy try

    sudo apt-get install libatlas-base-dev gfortran

------------------------------------
-->


3. Examples
--------

After running basic.py for some time, you can start having a look at some of the evolved morphologies and behaviors by opening up some of the generated .vxa files within the VoxCAD GUI. A .vxa file is just an XML file representing a robot that can be simulated by VoxCad/Voxelyze. Different custom versions of the physics engine can play slightly different .vxa files.
Navigate to evosoro/evosoro/_voxcad/release:
    
    cd ../_voxcad/release
    
Open VoxCad:

    ./VoxCad

Then select the desired .vxa file from 

    "File -> Import -> Simulation"

The .vxa files for the best performing individuals will be saved in 

    evosoro/evosoro/basic_data/bestSoFar/fitOnly.

Once the design is loaded, you can start the physics simulation by clicking the <img src="https://github.com/skriegman/evosoro/blob/master/evosoro/_voxcad/VoxCad/Icons/Sandbox.png" height="25" width="25"> icon in the top bar ("Physics Sandbox").  The robot should start moving: if it doesn't, please check the following section (Known issues).

4. Known issues
--------
If the robot does not move, disappears, or seems to behave in an unexpected manner when running a .vxa file in VoxCad (GUI), you may be affected by a known problem observed on some non-US machines.
The problem is due to an unexpected behavior of the <a href="http://www.cplusplus.com/reference/cstdlib/atof/">atof</a> function when the system's numeric <a href="https://en.wikipedia.org/wiki/Locale_(computer_software)">locale</a> differs from en_US.UTF-8, which entails loading wrong parameters from the .vxa file (in some cases it was observed how the atof function was approximating all double and floating point values to their integer part, which was the cause of the unexpected behavior).

While we work on a better solution, you can fix this problem by making sure that your machine is configured according to a US numeric locale.
Open the following file:

    sudo gedit /etc/default/locale

Make sure that LC_NUMERIC is set as follows:

    LC_NUMERIC="en_US.UTF-8"

Save, close the file, and reboot.

<!--
---------------------------------------------

The examples:

1. **basic.py** evolve running soft robots in a terrestrial environment. 

2. **swimming_basic.py** a simple mesh-based fluid model, that allows you to observe how soft morphologies evolve in a fluid environment.

3. **swimming_complex.py** evolve soft robots in a simple fluid environment along with additional attribtues that allow for more complexe behaviors.

4. **growth_basic.py** plants that grow towards a light source.

--------------------------------------------
 
-->

5. Documentation
-------------

Although the code included in this repository diverged from the main VoxCad/Voxelyze development branch some time ago, useful indications could be find in the online Voxelyze documentation, available [here](http://jonhiller.github.io/Voxelyze/annotated.html).


6. License
-------

Released under a MIT License (MIT)


7. References
--------------------

The difficult of co-optimizing brain and body:

>Cheney, N., Bongard, J., SunSpiral, V., & Lipson, H. (2016). *On the Difficulty of Co-Optimizing Morphology and Control in Evolved Virtual Creatures.* In Proceedings of The Fifteenth International Conference on the Synthesis and Simulation of Living Systems, ALIFE XV

Evolution of growing soft robots:

>Corucci, F., Cheney, N., Lipson, H., Laschi, C., & Bongard, J. (2016). *Material properties affect evolution’s ability to exploit morphological computation in growing soft-bodied creatures.* In Proceedings of The Fifteenth International Conference on the Synthesis and Simulation of Living Systems, ALIFE XV (pp. 234-241).


Evolution of swimming soft robots:

>Corucci, F., Cheney, N., Lipson, H., Laschi, C., & Bongard, J. (2016). *Evolving swimming soft-bodied creatures.* In Late Breaking Proceedings of The Fifteenth International Conference on the Synthesis and Simulation of Living Systems, ALIFE XV (p. 6).

From a cellular automata perspective:

>Cheney, N., & Lipson, H. (2016). *Topological evolution for embodied cellular automata.* Theoretical Computer Science, 633, 19-27.

A more complex task:

>Cheney, N., Bongard, J., & Lipson, H. (2015). *Evolving Soft Robots in Tight Spaces.* In Proceedings of the 2015 Annual Conference on Genetic and Evolutionary Computation (pp. 935-942). ACM.

A more complex controller:

>Cheney, N., Clune, J., & Lipson, H. (2014). *Evolved electrophysiological soft robots.* In ALIFE (Vol. 14, pp. 222-229).


Voxelyze & VoxCad:

>Hiller, J., & Lipson, H. (2014). *Dynamic simulation of soft multimaterial 3d-printed objects.* Soft Robotics, 1(1), 88-101.


Evolution of soft robots using generative encodings (CPPN):

>Cheney, N., MacCurdy, R., Clune, J., & Lipson, H. (2013). *Unshackling evolution: evolving soft robots with multiple materials and a powerful generative encoding.* In Proceedings of the 15th annual conference on Genetic and evolutionary computation (pp. 167-174). ACM.

Voxelyze & VoxCad plus evolving (and printing) soft robots:

>Hiller, J., & Lipson, H. (2012). *Automatic design and manufacture of soft robots.* IEEE Transactions on Robotics, 28(2), 457-466.


