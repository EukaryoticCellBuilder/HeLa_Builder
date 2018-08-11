.. HeLa Cell documentation master file, created by
   sphinx-quickstart on Tue Mar 27 21:01:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HeLa Cell
=========

Code to automatically generate HeLa cell input files with a variety of initial conditions for `Lattice Microbes <http://www.scs.illinois.edu/schulten/lm/index.html>`_.

About
-----
The HeLa package is designed to provide the constructed HeLa cell model for spatial simulations for use with Lattice Microbes. It allows customization of the cell geometry and can be adapted to add or remove organelles and compartments. In addition, arbitrary cellular processes in the frame of reaction-diffusion master equations can be studied within the HeLa cell geometry by adding the appropriate reaction and diffusion models. Its chief function is to generate arrays of input files with varying initial conditions in a Pythonic fashion.

The details of the model construction can be found in the publication: Ghaemi et al. **An** *in silico* **Mammalian Wole-Cell Model Reveals the Influence of Spatial Organization on RNA Splicing Efficiency,** *Journal* 2018. Examples include methods for constructing the whole cell HeLa model and the HeLa nucleus model used in the above referenced paper.

NOTE: The code requires Python version 3.6+.

Installation
------------
Install `Lattice Microbes & pyLM <http://www.scs.illinois.edu/schulten/lm/index.html>`_. Then having Lattice Microbes in your path, apply the patch :code:`pylm-builder.patch` in the Main folder by executing:

.. code-block:: bash

    git apply  pylm-builder.patch

No additional dependencies are required by Lattice Microbes and pyLM are needed for HeLa. Once done, the HeLa cell package can be installed simply by executing:

.. code-block:: bash

    python setup.py install

Creating an Input File
----------------------
To create an input file with the geometry from the manuscript, merely execute:

.. code-block:: bash

    python HeLa/hela.py

Documentation
=============

.. toctree::
   :maxdepth: 2

   overview
   examples
   tutorial

