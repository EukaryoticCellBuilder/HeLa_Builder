.. HeLa Cell documentation master file, created by
   sphinx-quickstart on Tue Mar 27 21:01:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Customization of The Model
==========================

Different parts of the code can be customized, by changing lines of the driver file (:code:`hela.py`), reaction, diffusion, particle counts or geometry which is our main focus here. Reaction (:code:`reactions.py`) and diffusion (:code:`diffusion.py`) can be customized corresponding to the process of interest. 
As an example, we can change the nucleus size and number of nuclear pore complexes (NPCs) of HeLa cell and generate series of input files.   

The line of the code in (:code:`hela.py`) that should be changed when changing parameters looks like:

.. code-block:: python

    savedFile = cell.generateLMFiles(filename)

Changing the nucleus size
-------------------------

To change the nucleus size of the HeLa cell and count of NPCs we can call:

.. code-block:: python

    files, parameter = cell.generateLMFiles(filename,
                                            {"nuclSize":[micron(3.74), micron(4.67), micron(5.29)],
                                             "n_NPCs":[1230,1515,2461]})

This will create 9 input files (the cross product of the "nuclSize" array and the "n_NPCs" array. One note, when generating multiple files, the :code:`generateLMFiles` function will return two lists:

  - :code:`files` - randomly generated filenames
  - :code:`parameters` - parameters associated with each file stored in a "name":"value" dictionary

Beacuse the filenames are randomly generated, you should save these and the associated parameters to file for further use:

.. code-block:: python

    with open("InputFiles.txt", "w") as f:
        for f, p in zip(files, parameters):
            f.write(f + "\t" + ["%s:%s"%(k, str(v)) for k,v in p.items()] + "\n")


Adapting ER 
------------

Because the ER starts from outside nucleus and expands to the plasma membrane, when the nuclear radius is changed the ER needs to change accordingly. 

The parameters within :code:`createERCellularAutomaton` that should be changed for different nucleus sizes are:

    - 3.74 micron
       - :code:`limits = [lambda x: x<= 58.0**2, lambda x: x> 139.0**2]`
       - :code:`fb = 0.9`
       - :code:`fd = 0.8`
       - :code:`steps = 13`
    - 4.15 micron
       - :code:`limits = [lambda x: x<= 65.0**2, lambda x: x> 139.0**2]`
       - :code:`fb = 0.9`
       - :code:`fd = 0.8`
       - :code:`steps = 11`
    - 4.67 micron
       - :code:`limits = [lambda x: x<= 73.0**2, lambda x: x> 139.0**2]`
       - :code:`fb = 0.9`
       - :code:`fd = 0.8`
       - :code:`steps = 11`
    - 5.29 micron
       - :code:`limits = [lambda x: x<= 82.0**2, lambda x: x> 139.0**2]`
       - :code:`fb = 0.9`
       - :code:`fd = 0.8`
       - :code:`steps = 11`

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
