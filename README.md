# HeLa Cell

Code to automatically generate HeLa cells for [Lattice Microbes](http://www.scs.illinois.edu/schulten/lm/index.html).

## About
The HeLa package is designed to provide the constructed HeLa cell model for spatial simulations for use with Lattice Microbes software. It allows customization of the cell geometry and can be adapted to add or remove organelles and compartments. In addition, arbitrary cellular processes in the frame of reaction-diffusion master equations can be studied within the HeLa cell geometry by adding the appropriate reaction and diffusion models. Its chief function is to generate arrays of input files with varying initial conditions in a Pythonic fashion.

The details of the model construction can be found in the publication: Ghaemi et al. **An** *in silico* **Human Cell Model Reveals the Influence of Spatial Organization on RNA Splicing Efficiency,** *Journal* 2019. Examples include methods for constructing the whole cell HeLa model and the HeLa nucleus model used in the above referenced paper.

NOTE: The code requires Python version 3.6+.

Any questions about using/modifying the code can be addressd to: Zhaleh Ghaemi ghaemi@illinois.edu

## Installation
Install [Lattice Microbes & pyLM](http://www.scs.illinois.edu/schulten/lm/index.html). Then apply the patch code:`pylm-builder.patch` in the Main folder by executing:

``` bash
    git apply pylm-builder.patch
````

No additional dependencies above those required by Lattice Microbes and pyLM are needed for HeLa. Once done, the HeLa cell package can be installed simply by executing:

```bash
  python setup.py install
```

## The HeLa model Documentation

The documentation can then be accessed via:

```bash
  <webbrowser> docs/_build/html/index.html
```


