.. HeLa Cell documentation master file, created by
   sphinx-quickstart on Tue Mar 27 21:01:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HeLa Cell Model
===============

The code in this repository can be used to generate input files for simulations of HeLa cell with Lattice Microbes software. The code files in the directory :code:`HeLa` construct the cell geometry, and reaction and diffusion models for simulation. This code can be customized to generate input files with different initial conditions (e.g., varying the nucleus size). This page introduces differet parts of the code. 

The files used to generate the input files include:

  - :code:`hela.py` - Driver code that actually constructs the input file.
  - :code:`hela_geometry.py` - Defines a function that sets up the spatial geometry for the HeLa cell
  - :code:`diffusion.py` - Defines a function that sets up a sample diffusion model 
  - :code:`reaction.py` - Defines a function that sets up a sample reaction model 

The functions are defined as:

.. code-block:: python
    
    # Spatial geometry model
    def geometryModel(sim,
                      cellCenter = lm.point(*micron(9,9,9)),
                      cellRadius = micron(8.9),
                      membraneThickness= micron(0.128),
                      nuclSize= micron(4.15),
                      speckleRadius=micron(0.35),
                      cajalRadius=micron(0.5),
                      poreRadius=micron(0.083),
                      n_cajals=4,
                      n_speckles=20,
                      n_NPCs=1818,
                      n_mito=2000,
                      fb=0.9,
                      fd=0.8,
                      steps=11,
                      lsize=288,
                      limits=[lambda x: x<= 73.0**2,lambda x: x> 139.0**2],
                      mitochondriaDL = lambda s: int(fd*27*(steps-0.25*s)**2/steps**2)
        ):
        ...

    # Diffusion model
    def diffusionModel(sim,
                       d_premRNA = 6.1e-13,
                       d_S = 2.6e-13,
                       d_SpremRNA = 2.6e-13,
                       d_mRNA = 6.1e-13):
        ...


    # Reaction model
    def reactionModel(sim,
                      Ktr=0.0034,
                      Kon=4.66e7,
                      Koff=0.062,
                      Ksplice=0.05
                      ):
        ...

    # Particle model
    def particleModel(sim,
                      n_gene=200,
                      n_premRNA=0,
                      n_S=1000,
                      n_mRNA=0,
                      n_SpremRNA=0):
        ...     


These function "signatures" are importand as they define the parameters that can be varied. Any "named" parameter (e.g., :code:`d_mRNA` or :code:`nuclSize`) can be varied by changig the value.

Driver code (:code:`HeLa/hela.py`):

.. literalinclude:: ../HeLa/hela.py
   :language: python
   :linenos:

Code to define cell geometry (:code:`HeLa/hela_geometry.py`):

.. literalinclude:: ../HeLa/hela_geometry.py
   :language: python
   :linenos:

The geomrtry contains plasma membrane, cytoplasm, mitochondria, Golgi, ER, nucleus, inside which there are nuclear speckles and Cajal bodies, and there are nuclear pore complexes on nuclear envelope. As explained in details in manuscript by Ghaemi et al., 2018, the ER was generated using a cellular automata program as below: 


ER Generation Code
------------------
The ER is generated using a cellular automata program. The code can be found within the :code:`hela.py` file. The function begins with:

.. code-block:: python

    def createERCellularAutomaton(size, 
                                  dimensions=2, 
                                  survivalRate = 0.5, 
                                  deathLimitFxn = 3,
                                  birthNumberFxn = 4,
                                  steps = 10,
                                  limits = []):
        ...

The definition of these parameters can be found within the manuscript. This is called at the end of the geometry generation code. The density and structure of the ER is sensitive to these parameters. These parameters can be accessed via:

    - :code:`fb` - Birth number
    - :code:`fd` - Death number
    - :code:`steps` - Number of iterations of the cellular automata program
    - :code:`lsize` - Lattice size (should equal that of the RDMESimulation)
    - :code:`limits` - Functions that take a distance, measured in lattice sites from the center of the cell, outside of which no automatons can exist
    - :code:`deathLimitFxn` - an iteration dependent modifier for the death number
    - :code:`birthNumberFxn` - an iteration dependent modifier for the birth number

