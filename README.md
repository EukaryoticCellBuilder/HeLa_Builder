# HeLa Cell

Code to automatically generate HeLa cells for [Lattice Microbes](http://faculty.scs.illinois.edu/schulten/Software2.0.html).

## About
The HeLa package is designed to provide the constructed HeLa cell model for spatial simulations for use with Lattice Microbes software. It allows customization of the cell geometry and can be adapted to add or remove organelles and compartments. In addition, arbitrary cellular processes in the frame of reaction-diffusion master equations can be studied within the HeLa cell geometry by adding the appropriate reaction and diffusion models. Its chief function is to generate arrays of input files with varying initial conditions in a Pythonic fashion.

The details of the model construction can be found in the publication: Ghaemi et al. [**An** *in silico* **Human Cell Model Reveals the Influence of Spatial Organization on RNA Splicing Efficiency,**](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007717) PLoS Comput Biol, 2020. Examples include methods for constructing the whole cell HeLa model and the HeLa nucleus model used in the above referenced paper.

NOTE: The code requires Python version 3.6+.

Any questions about using/modifying the code can be addressd to: Zhaleh Ghaemi ghaemi@illinois.edu

## Installation
Install [Lattice Microbes & pyLM](http://faculty.scs.illinois.edu/schulten/Software2.0.html). Then apply the patch code:`pylm-builder.patch` in the Main folder by executing:

``` bash
    git apply pylm-builder.patch
````

No additional dependencies above those required by Lattice Microbes and pyLM are needed for HeLa. Once done, the HeLa cell package can be installed simply by executing:

```bash
  python setup.py install
```
## Quick documentation

To generate the cell model:

```bash
   $python Hela/hela.py
```
To start the reaction-diffusion master equation (RDME) simulations using Lattice Microbes (lm): 

```bash
   $lm -r 1 -sp -sl lm::rdme::MpdRdmeSolver -cr 8 -gr 8 -f cell.lm 
```
This command runs one replica (r) of the RDME simulation on the generated cell (cell.lm) using 8 CPU and 8 GPU cores

## The full HeLa model Documentation

The full documentation can then be accessed via:

```bash
  <webbrowser> docs/_build/html/index.html
```


