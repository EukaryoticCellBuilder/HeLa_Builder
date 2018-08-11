#
# University of Illinois Open Source License
# Copyright 2018 University of Illinois Board of Trustees
# All rights reserved.
# 
# Developed by: Luthey-Schulten Group
# University of Illinois at Urbana-Champaign
# http://www.scs.uiuc.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to
# do so, subject to the following conditions:
# 
#   - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimers.
#   
#   - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimers in the documentation
#   and/or other materials provided with the distribution.
#   
#   - Neither the names of the Luthey-Schulten Group, University of Illinois at
#   Urbana-Champaign, nor the names of its contributors may be used to endorse or
#   promote products derived from this Software without specific prior written
#   permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS WITH THE SOFTWARE.
#

# Standard imports
import copy
import functools
import itertools
import random
import string

# Lattice Microbes imports
import lm
from pyLM.units import *
from pyLM.RDME import RDMESimulation

class Cell(object):
    """A class used as a template for a cell simulation.

    This class simplifies generation of multiple input files. The user
    registers functions that define geometry, reaction, diffusion, and 
    particle models for a Lattice Microbes simulation. It can then be
    used to generate many input files varying one or more of the parameters
    in an automated fashion.

    If specified, the functions will be applied in this order:
        1) :code:`geometryModelFxn` - Specifies regions and creats spatial geometry.
        2) :code:`reactionModelFxn` - Adds reactions among particle types.
        3) :code:`diffusionModelFxn` - Specifies diffusion in and among regions.
        4) :code:`particleModelFxn` - Adds particles of specified counts to different regions.

    Args:
        - simulationBase: A function that creates the initial :py:class:`pyLM.RDME.RDMESimulation` object. This must return the newly created simulation object.
        - defaultParameters (:py:obj:`dict`): A dictionary containing all the parameters necessary to call the reaction model, the diffusion model and the geometry model.

    .. code-block:: python

        def initFxn():
            sim = RDMESimulation(dimensions=micron(1024,1024,1024), spacing=nm(32))

        # Reaction Model
        def rxnModel(sim, ...):
            pass
        # Diffusion Model
        def diffModel(sim, ...):
            pass
        # Geometry Model
        def geomModel(sim, ...):
            pass
        # Particle counts
        def partModel(sim, ...):
            pass

        cell = lmarray.Cell(sim)
        cell.setReactionModel(rxnModel)
        cell.setDiffusionModel(diffModel)
        cell.setGeometryModel(geomModel)
        cell.setParticleCounts(partModel)

        # Create a single simulation file
        cell.generateLMFiles("MySimulation")
        
        # Create an array of simulation files
        filnames, parameters = cell.generateLMFiles("MySimulations", {"k":[1,2,3], "a":[8,9,10], ...})

    Atributes:
        - reactionModelFxn: A function to add reactions among particle types (optional).
        - diffusionModelFxn: A function to specify diffusion in and among regions (optional).
        - geometryModelFxn: A function to specify regions and create the spatial geometry (optional).
        - particleModelFxn: A function to add particles to the various regions (optional).
        - templateRDME (:py:class:`pyLM.RDME.RDMESimbulation`): The template simulation including timingparameters, species defined in the simulation, and optionally any other common simulation specifications.
    """
    def __init__(self, simulationBase, defaultParameters = {}):
        if not callable(simulationBase):
            raise TypeError("'simulationBase' incorrect type.")
        self.templateRDME = simulationBase

        # Set up HeLa.Cell specific features
        self.reactionModelFxn  = None
        self.diffusionModelFxn = None
        self.geometryModelFxn  = None
        self.particleModelFxn  = None

        # Default parameter values
        self.defaultParameters = defaultParameters

    ###################
    # Setup Functions #
    ###################
    def setReactionModel(self, fxn):
        '''Specify the reaction model. This function should generally 
        specify all the reactions associated with the cell. The function
        should take at least one parameter, and the first parameter should
        be the :py:class:`pyLM.RDME.RDMESimulation`. Other parameters
        should be named and may optionally have a default parameter.'''
        if not callable(fxn):
            raise TypeError("Reaction model must be a function or functor.")
        if fxn.__code__.co_argcount < 1:
            raise TypeError("Reaction model must take at least 1 argument.")
        self.reactionModelFxn = fxn

    def setDiffusionModel(self, fxn):
        '''Specify the diffusion model. This function should generally 
        specify all the diffusions associated with the cell. The function
        should take at least one parameter, and the first parameter should
        be the :py:class:`pyLM.RDME.RDMESimulation`. Other parameters
        should be named and may optionally have a default parameter.'''
        if not callable(fxn):
            raise TypeError("Diffusion model must be a function or functor.")
        if fxn.__code__.co_argcount < 1:
            raise TypeError("Diffusion model must take at least 1 argument.")
        self.diffusionModelFxn = fxn

    def setGeometryModel(self, fxn):
        '''Specify the geometry model. This function should generally 
        specify all the shapes associated with the cell. The function
        should take at least one parameter, and the first parameter should
        be the :py:class:`pyLM.RDME.RDMESimulation`. Other parameters
        should be named and may optionally have a default parameter.'''
        if not callable(fxn):
            raise TypeError("Geometry model must be a function or functor.")
        if fxn.__code__.co_argcount < 1:
            raise TypeError("Geometry model must take at least 1 argument.")
        self.geometryModelFxn = fxn

    def setParticleCounts(self, fxn):
        '''Add particles to the cells. This function should generally 
        specify all the parameter counts associated with the cell. The function
        should take at least one parameter, and the first parameter should
        be the :py:class:`pyLM.RDME.RDMESimulation`. Other parameters
        should be named and may optionally have a default parameter.'''
        if not callable(fxn):
            raise TypeError("Particle model must be a function or functor.")
        if fxn.__code__.co_argcount < 1:
            raise TypeError("Particle model must take at least 1 argument.")
        self.particleModelFxn = fxn

    ##############################
    # Model Generation Functions #
    ##############################
    def generateLMFiles(self, filename, parameters={}):
        '''Create one or more LM input files (.lm) to 

        If no parameters are passed to the function, this will create a
        a single file by binding the default parameters associated with the
        :py:class:`Cell`. The filename will be returned.

        If parameters are passed in the dictionary, they should be specified
        with a key and an associated array of values (e.g. "k1" : [1,2,3]).
        This function will create a Cartesian product of all the specified 
        parameters values and produce one LM file for each set of parameters
        by passing them to the associated functions (e.g. reactionModelFxn, 
        diffusionModelFxn, etc.). A randomly generated filename of the form
        will be created (filename_[8 characters].lm). Finally, the function
        will return a list of these filenames and a list of the associated
        parameters.

        Args:
            - filename (:py:obj:`str`): Base for the filename; don't include the ".lm".
            - parameters (:py:obj:`dict`): If specified, will create an array of jobs spanning each parameter value.

        Returns:
            If generating only a single file, will return the filename. Otherwise
            will return a list of file names (randomly generated) and a list of 
            associated parameters (as a dictionary).

        Example:

        .. code-block:: python

            cell = ...
            filenames, parameters = cell.generateLM("inputFileBase_", parameters = {"k":[1,2,3], "d":[4,5,6]})
            with open("InputFiles.txt", "w") as of:
                for fn, params in zip(filenames, parameters):
                    of.write(fn + "\\t" + ";".join(["%s=%f"%(k,v) for k,v in params.items()] + "\\n")
        '''
        def getParameterSet(fxn, paramDict):
            # This function reads the function arguments
            #  and maps from the paramDict to map them.
            # It returns a dictionary that can be passed
            #  to the function using ** to expand all parameters.
            # If not all required parameters are specified,
            #  this will raise an exception.

            # Handles to parameter names
            requiredParams = [x for x in fxn.__code__.co_freevars]
            optionalParams = [x for x in fxn.__code__.co_varnames if x not in fxn.__code__.co_freevars]

            # Determine the mapping
            mapped = {}
            for var, value in paramDict.items():
                # Set the required parameters
                if var in requiredParams:
                    mapped[var] = value
                    requiredParams.remove(var)
                # Set optional parameters
                if var in optionalParams:
                    mapped[var] = value
            
            if len(requiredParams) > 0:
                raise KeyError("Not all required parameters could be mapped for: %s"%(fxn))

            # Return parameter sets
            return mapped
            

        if len(parameters) > 0: # Create an array of files
            keys = sorted(parameters.keys())
            parameterLists = [parameters[x] for x in keys]
            parameterSets = itertools.product(*parameterLists)

            filenames = []
            associatedParameters = []
            for i, params in enumerate(parameterSets):
                print("Generating %d/%d simulation files."%(i+1, functools.reduce(lambda x,y: x*y, map(len, parameterLists))))
                # Create parameter dictionary
                localParameters = {}
                for k, v in zip(keys, params):
                    print("\t%s=%s"%(str(k), str(v)))
                    localParameters[k] = v

                sim = self.templateRDME()
                
                # Apply geometry model
                if self.geometryModelFxn != None:
                    geomParamMap = getParameterSet(self.geometryModelFxn, self.defaultParameters)
                    for k, v in getParameterSet(self.geometryModelFxn, localParameters).items():
                        geomParamMap[k] = v
                    self.geometryModelFxn(sim, **geomParamMap)

                # Apply reaction model
                if self.reactionModelFxn != None:
                    rxnParamMap = getParameterSet(self.reactionModelFxn, self.defaultParameters)
                    for k, v in getParameterSet(self.reactionModelFxn, localParameters).items():
                        rxnParamMap[k] = v
                    self.reactionModelFxn(sim, **rxnParamMap)
                    
                # Apply diffusion model
                if self.diffusionModelFxn != None:
                    diffParamMap = getParameterSet(self.diffusionModelFxn, self.defaultParameters)
                    for k, v in getParameterSet(self.diffusionModelFxn, localParameters).items():
                        diffParamMap[k] = v
                    self.diffusionModelFxn(sim, **diffParamMap)

                # Apply particle model
                if self.particleModelFxn != None:
                    partParamMap = getParameterSet(self.particleModelFxn, self.defaultParameters)
                    for k, v in getParameterSet(self.particleModelFxn, localParameters).items():
                        partParamMap[k] = v
                    self.particleModelFxn(sim, **partParamMap)
           
                # Save the file
                localFilename = filename + "_" + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8)) + ".lm"
                sim.save(localFilename)
                filenames.append(localFilename)
                associatedParameters.append(localParameters)
                print("Done")
                print("")

            return filenames, associatedParameters

        else: # Just create 1 file
            print("Creating simulation file.")
            sim = self.templateRDME()
            
            # Apply geometry model
            if self.geometryModelFxn != None:
                geomParamMap = getParameterSet(self.geometryModelFxn, self.defaultParameters)
                self.geometryModelFxn(sim, **geomParamMap)

            # Apply reaction model
            if self.reactionModelFxn != None:
                rxnParamMap = getParameterSet(self.reactionModelFxn, self.defaultParameters)
                self.reactionModelFxn(sim, **rxnParamMap)
                
            # Apply diffusion model
            if self.diffusionModelFxn != None:
                diffParamMap = getParameterSet(self.diffusionModelFxn, self.defaultParameters)
                self.diffusionModelFxn(sim, **diffParamMap)

            # Apply particle model
            if self.particleModelFxn != None:
                partParamMap = getParameterSet(self.particleModelFxn, self.defaultParameters)
                self.particleModelFxn(sim, **partParamMap)
           
            # Save the file
            sim.save(filename + ".lm")
            print("Done")

            return filename + ".lm"

if __name__ == "__main__":
    import os
    import sys
    import numpy as np

    ######################
    # Test the Cell code #
    ######################a
    print("Testing HeLa.Cell code...")
    def initFxn():
        # Create RDME simulation
        sim = RDMESimulation(dimensions=micron(4.096,4.096,4.096), spacing=nm(32.0), defaultRegion="extracellular")
    
        # Create simple gene expression model
        sim.defineSpecies(["D","pM","M","P","S"])
    
        # Create regions
        sim.addRegion("Nucleus")
        sim.addRegion("Cytoplasm")

        return sim

    #################
    # Define Models #
    #################
    def geomModel(sim,
                  cellCenter=lm.point(*micron(2.0,2.0,2.0)),
                  nucleusRadius=micron(0.4),
                  cellRadius=micron(2.0)):
        # Add locations
        cytoplasm = lm.Sphere(cellCenter, cellRadius,    sim.siteTypes["Cytoplasm"])
        nucleus   = lm.Sphere(cellCenter, nucleusRadius, sim.siteTypes["Nucleus"])
        sim.addShape(cytoplasm)
        sim.addShape(nucleus)

    def rxnModel(sim, k1, k2, k3, k4=np.log(2.0)/6.0, k5=np.log(2.0)/12.0):
        # Get handles for regions
        nucleus   = sim.modifyRegion("Nucleus")
        cytoplasm = sim.modifyRegion("Cytoplasm")

        # Add gene expression model
        nucleus.addReaction("D", ("D","pM"), k1)
        nucleus.addReaction(("pM","S"), ("M","S"), k2)
        cytoplasm.addReaction("M", ("M","P"), k3)
        cytoplasm.addReaction("M","", k4)
        cytoplasm.addReaction("P","", k5)
        
    def diffModel(sim, dM, dS, dP, dD=0.0):
        # Get handles for regions
        nucleus   = sim.modifyRegion("Nucleus")
        cytoplasm = sim.modifyRegion("Cytoplasm")

        # Define within region diffusions
        nucleus.setDiffusionRate("D", rate=dD)
        nucleus.setDiffusionRate("pM", rate=dM)
        nucleus.setDiffusionRate("S", rate=dS)
        nucleus.setDiffusionRate("M", rate=dM)
        cytoplasm.setDiffusionRate("M", rate=dM)
        cytoplasm.setDiffusionRate("P", rate=dM)

        # Define between region diffusions
        sim.setTransitionRate("M", "Nucleus", "Cytoplasm", dM/4.0)

    def partModel(sim):
        sim.addParticles("D", "Nucleus", 1)
        sim.addParticles("S", "Nucleus", 200)
        sim.addParticles("pM", "Nucleus", 5)
        sim.addParticles("M", "Nucleus", 3)
        sim.addParticles("M", "Cytoplasm", 10)
        sim.addParticles("P", "Cytoplasm", 50)

    ####################
    # Create HeLa Cell #
    ####################
    print("Creating HeLa.Cell object...")
    cell = Cell(simulationBase=initFxn,
                defaultParameters={"k1":np.log(2.0)/0.05,"k2":np.log(2.0)/0.001,"k3":np.log(2)/12, "dM":1e-13, "dP":1e-12, "dS":5e-13})

    cell.setGeometryModel(geomModel)
    cell.setDiffusionModel(diffModel)
    cell.setReactionModel(rxnModel)
    cell.setParticleCounts(partModel)

    print("Testing generation of one cell with no parameters...")
    try:
        filename = cell.generateLMFiles(filename="test1")
        os.remove(filename)
    except Exception as e:
        print("Failed creating single cell with no parameters.\nError: %s"%(str(e)))
        raise e
        sys.exit(-1)

    print("Testing generation of multiple cell with some parameters...")
    try:
        filenames, parmaeters = cell.generateLMFiles(filename="test2",
                                        parameters={"k1":[5.0,10.0],
                                                    "k2":[0.05],
                                                    "cellCenter":[lm.point(*micron(1.9,1.9,1.0))],
                                                    "cellRadius":[micron(1.5), micron(1.6), micron(1.7)]})
        for f in filenames:
            os.remove(f)
    except Exception as e:
        print("Failed creating multiple cell with some parameters.\nError: %s"%(str(e)))
        sys.exit(-1)

