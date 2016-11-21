evosoro: soft robot simulator
=======================================

<div class="row">
<a href=https://youtu.be/EXuR_soDnFo>
<img src="https://github.com/skriegman/test/blob/master/nick.png" height="135" width="135">
</a>

<a href=https://www.youtube.com/watch?v=HgWQ-gPIvt4>
<img src="https://github.com/skriegman/test/blob/master/electro.png" height="135" width="135">
</a>

<a href=https://youtu.be/4ZqdvYrZ3ro>
<img src="https://github.com/skriegman/test/blob/master/swimming.png" height="135" width="135">
</a>

<a href=https://youtu.be/Cw2SwPNwcfM>
<img src="https://github.com/skriegman/test/blob/master/plant1.png" height="135" width="135">
</a>

<a href=https://www.youtube.com/watch?v=XqIUJcuOgmwl>
<img src="https://github.com/skriegman/test/blob/master/teeth1.png" height="135" width="135">
</a>

<a href=https://www.youtube.com/watch?v=r_SL8VUt-wA>
<img src="https://github.com/skriegman/test/blob/master/cage.png" height="135" width="135">
</a>

</div>

Evosoro is a Python soft robot simulation library based on Voxelyze. It provides a high-level interface for dynamic simulation of soft multimaterial robots.






Citing
------

If using this code for academic purposes please cite the following two papers:

Hiller, J., & Lipson, H. (2014). 
*Dynamic simulation of soft multimaterial 3d-printed objects.*
Soft Robotics, 1(1), 88-101.


Cheney, N., MacCurdy, R., Clune, J., & Lipson, H. (2013). 
*Unshackling evolution: evolving soft robots with multiple materials and a powerful generative encoding.* 
In Proceedings of the 15th annual conference on Genetic and evolutionary computation (pp. 167-174). ACM.


<!--
Dependencies
------------

- Python 2.7

### Mandatory

- [numpy](http://www.numpy.org/)

- [scipy](http://www.scipy.org/)

- [networkx](http://networkx.github.io/)

- [pandas](http://pandas.pydata.org/)

### Recommended

- [matplotlib](http://matplotlib.org/)

- [seaborn](http://seaborn.pydata.org/)
-->

Installation
------------

<!--To install the released version, just do-->
    
<!--    pip install evosoro-->

<!--You may instead want to use the development version from Github, by running-->

<!--    pip install git+git://github.com/skriegman/evosoro.git#egg=evosoro-->

It is recommended that you install [anaconda](https://docs.continuum.io/anaconda/install#) as your Python distribution. Anaconda is a free package manager and Python distribution that includes all of the dependies required for evosoro. However if you instead choose to manually install Python,

    sudo apt-get install python-dev python-pip
    sudo pip install numpy networkx scipy decorator pandas

Install Qt and QMake if you have not already done so, specifically these packages: "libqt4-dev", "qt4-qmake", "libqwt-dev", "freeglut3-dev" and "zlib1g-dev".

    sudo apt-get install libqt4-dev qt4-qmake libqwt-dev freeglut3-dev zlib1g-dev


Install git if you have not already done so.

    sudo apt-get install git

Navigate to your working directory (e.g. your home).

    cd /home

Clone the repo.

    git clone https://github.com/skriegman/evosoro.git

Navigate to the _voxcad directory and compile.

    cd evosoro/evosoro/_voxcad/
    make

Install voxelyze lib:

    cd Voxelyze
    make
    make installusr
    cd ../voxelyzeMain/
    make

Navigate back out to the examples folder and run basic.py
    
    cd ../../examples
    python basic.py
    
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


Examples
--------

After running basic.py for some time, you can start having a look at some of the evolved morphologies and behaviors by opening up some of the generated .vxa files with the graphical interface VoxCAD.

    ./evosoro/evosoro/_voxcad/release/VoxCad

Then selecting the desired .vxa file from 

    "File -> Import -> Simulation"

The .vxa files for the best performing individuals will be saved in evosoro/evosoro/basic_data/bestSoFar/fitOnly.


<!--
---------------------------------------------

The examples:

1. **basic.py** evolve running soft robots in a terrestrial environment. 

2. **swimming_basic.py** a simple mesh-based fluid model, that allows you to observe how soft morphologies evolve in a fluid environment.

3. **swimming_complex.py** evolve soft robots in a simple fluid environment along with additional attribtues that allow for more complexe behaviors.

4. **growth_basic.py** plants that grow towards a light source.

--------------------------------------------
 
-->

Documentation
-------------

Online documentation for Voxelyze is available [here](http://jonhiller.github.io/Voxelyze/annotated.html).


License
-------

Released under a MIT License (MIT)



<!--Other papers-->

<!--Hiller, J., & Lipson, H. (2012). -->
<!--Automatic design and manufacture of soft robots. -->
<!--IEEE Transactions on Robotics, 28(2), 457-466.-->

<!--F. Corucci, N. Cheney, H. Lipson, C. Laschi, and J. Bongard, "Evolving swimming soft-bodied creatures", ALIFE XV, The Fifteenth International Conference on the Synthesis and Simulation of Living Systems, 2016 (late breaking abstract)-->
<!--https://youtu.be/4ZqdvYrZ3ro-->
<!--F. Corucci, N. Cheney, H. Lipson, C. Laschi, and J. Bongard, “Material properties affect evolution’s ability to exploit morphological computation in growing soft-bodied creatures,”  ALIFE XV, The Fifteenth International Conference on the Synthesis and Simulation of Living Systems, 2016-->
<!--https://youtu.be/Cw2SwPNwcfM-->



