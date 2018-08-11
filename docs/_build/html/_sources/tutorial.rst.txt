.. HeLa Cell documentation master file, created by
   sphinx-quickstart on Tue Mar 27 21:01:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Advanced Usage
==============

While the HeLa cell example demonstrates the utility of the :py:mod:`lmarray`; the Python class contained within (:py:class:`lmarray.cell.Cell`) is a general class for creating sets of input files. It can be thought of as an "input file factory". The user registers functions specifying various model features like the reaction model, diffusion model, geometry model, initial particle counts, etc. with the :py:class:`lmarray.cell.Cell` class and can then use the cell to generate a variety of input files with varied input parameters. It does this using the *names* of the parameters to these function along user specified ranges for these. This will all become clear shortly. 

It is highly recommended that you familiarize yourself with the use of pyLM prior to reading this tutorial. See the `user guide <http://www.scs.illinois.edu/schulten/lm/download/lm23/InstructionGuide.pdf>`_ and `reference guide <http://www.scs.illinois.edu/schulten/lm/documentation/index.html>`_.

Minimal Example
---------------

First we import the relevant libraries from :py:mod:`pyLM`:

.. code-block:: python

    from pyLM.units import *
    from pyLM.RDME import RDMESimulation

    import lmarray

Next, we create a function to create the simulation volume, and the various parameters that will be common to all the input files such as regions within the volume, types of chemical species, and timestep parameters.

.. code-block:: python

    def createSimulation():
        # This defines a box that is ~1 cubic micron that is discretized
        #  into volumes that are 32nm^3. 
        sim = RDMESimulation(dimension=micron(1.024,1.024,1.024), spacing=nm(32))

        # Next we define a region called "nucleus" which we will use later to 
        #  define spatial geometries within the simulation domain. Note,
        #  there in an implicit region named "default" which is created in the 
        #  code line above.
        sim.addRegion("nucleus")

        # Next we define the chemical species for the simulation
        sim.defineSpecies(["dna", "TF", "dna:TF", "mRNA"])
        
        # Next we will define several simulation parameters
        sim.setTimestep(0.1)             # Discretize time into 0.1 second increments
        sim.setWriteInterval(1.0)        # Save particle counts every second
        sim.setLatticeWriteInterval(1.0) # ... and the lattice too
        sim.setSimulationTime(3600.0)    # Run for one hour

        # Finally, the function returns the skeleton of the simulation,
        #  i.e., a simulation without any interesting parameters specified.
        return sim


So this is not very interesting so far. So far, we have a function that creates a box with discretized grid, that will have a region within it called nucleus". We've also defined several chemical species that will exist, and defined some simulation parameters like the write frequency and timestep size. To spice things up, we will create a simple reaction model.

.. code-block:: python

    def reactionModel(sim,
                      k_bind=10.0,
                      k_unbind=5.0,
                      k_transcribe=2.0):
        # Now we will define a reaction model
        sim.addReaction(("dna","TF"), "dna:TF", k_bind)
        sim.addReaction("dna:TF", ("dna","TF"), k_unbind)
        sim.addReaction("dna:TF", "mRNA", k_transcribe)


This function adds three reactions to the model: 1) A transcription factor binding to a gene, 2) the associated unbinding reaction, and 3) the transcription of the gene when the TF is bound.Notice the named parameters, and default values. Within the :py:class:`lmarray.cell.Cell` class, all parameters must have a default parameter, or be specified when calling :py:meth:`lmarray.cell.Cell.generateLMFiles` (see below).

Next, we create the diffusion model:

.. code-block:: python

    def diffusionModel(sim,
                       d_DNA=0.0,
                       d_TF=1e-12,
                       d_mRNA=1e-13,
                       d_mRNA_nuc_cyt=1e-14):
        # Get handles to the nucleus and the cytoplasm
        cyt = sim.modifyRegion("default")
        nuc = sim.modifyRegion("nucleus")
        
        # Specify the diffusion rates for each type within
        #  each region (the default is 0)
        nuc.setDiffusionRate(species='dna', rate=d_DNA)
        nuc.setDiffusionRate(species='dna:TF', rate=d_DNA)
        nuc.setDiffusionRate(species='TF', rate=d_TF)
        nuc.setDiffusionRate(species='mRNA', rate=d_mRNA)
        cyt.setDiffusionRate(species='mRNA', rate=d_mRNA)

        # Specify diffusion of mRNA between regions
        sim.setTwoWayTransitionRate(species='mRNA', one='nucleus', two='default', rate=d_mRNA_nuc_cyt)
        
This defins the diffusion of species within the various regions and also allows the mRNA to transition between regions.

Next we define a function that creates the nucleus:

.. code-block:: python

    def geometryModel(sim,
                      nucleusRadius):
        # Create the nucleus and add to the simulation
        nucleus = lm.Sphere(micron(0.512,0.512,0.512), nucleusRadius, sim.siteTypes['nucleus'])
        nucleus.thisown = 0
        sim.lm_builder.addRegion(nucleus)

Note that we don't define the nuclear radius here, so it must be specified below in the :py:meth:`lmarray.Cell.generateLMFiles` function. 

Finally, we define the 

.. code-block:: python

    def particleModel(sim,
                      n_DNA=2,
                      n_TF=20,
                      n_mRNA=0):
        # Get handles to the nucleus and the cytoplasm
        cyt = sim.modifyRegion("default")
        nuc = sim.modifyRegion("nucleus")
        
        # Add the actual particles
        nuc.addParticles("dna", n_DNA)
        cyt.addParticles("TF", n_TF)
        cyt.addParticles("mRNA", n_mRNA)

The rest of the script is used to generate a battery of input files. It will be described below.

.. code-block:: python

    # Create a "cells" object with a default parameter for the nucleus radius
    cells = lmarray.Cell(simulationBase=initFunction,
                         defaultParameters={"nucleusRadius":micron(0.5)})

    # Next we add all the models defined above
    cells.setReactionModel(reactionModel)
    cells.setDiffusionModel(diffusionModel)
    cells.setGeometryModel(geometryModel)
    cells.setParticleCounts(particleModel)

    # Create a single file with all default parameters
    filename = cell.generateLMFiles("DefaultParameters")

    # Create an array of simulations (a 3D grid in fact) varying
    #  the transcription factor count, the binding rate, and the 
    #  nucleus radius
    filenames, parameters = cell.generateLMFiles("VariedParameters",
                {"n_TF":[10,20,30,40],
                 "k_bind":[5.0,10.0,15.0,20.0],
                 "nucleusRadius":[micron(0.3), micron(0.4), micron(0.5)]}

So as we can see from the example above, this functionality can easily be used to create a battery of simulations with various different parameters very easily. The last line of the code demonstrates how by passing a dictionary mapping names of parameter to lists of inputs, we can vary the parameters used to create the input files. The :py:class:`lmarray.cell.Cell` class will map these to the names of the parameters within each of the functions defined above and pass the appropriate value. The final example creates a 3D grid of input files with 4 differnt transcription factor counts, 4 different binding rats and 3 different nucleus radii for a total of 48 simulations. When creating arrays of jobs like this, the function will return a list of filenames (which are randomized) and a list of associated parameters (as a dictionary), so the user can reference them later.

Hopefully, you can see the utility of this approach.

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
