Workflow
========

A typical BICePs sampling includes four core steps:`Preparation`,`Restraint`,`Posteriorsampler` and`Analysis`.

Preparation
-----------

The `Preparation` class converts all raw
data to BICePs readable format. It asks for experimental data as well as
correspondingly precomputed experimental observables from simulation. We
recommend users to use [MDTraj](http://mdtraj.org) to compute all
necessary experimental quantities from simulation or use our prepared
functions in the `toolbox`. A tutorial jupyter notebook is availabel [here](https://github.com/vvoelz/biceps/blob/master/BICePs_2.0/tutorials/Preparation/Preparation.ipynb)

Restraint
---------

The `Restraint` class initializes all necessary functions to construct numpy array containing information for BICePs sampling. As a parent class, it also includes child classes for different experimental restraints.

Posteriorsampler
----------------

The `Posteriorsampler` class is closely working with the `Restraint` class. A Markov chain Monte Carlo sampling is performed based on the [Metroplis-Hastings criterion](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm).

Analysis
--------

The `Analysis` is consist of two parts: 
1. using [MBAR](https://pymbar.readthedocs.io/en/master/index.html) algorithm to compute populations and [BICePs scores](https://github.com/vvoelz/biceps/blob/master/BICePs_2.0/markdown/Theory.ipynb) 
2. plot the figures to show posterior distribution of population and nuisance parameters.


## To get started check this [tutorial](https://github.com/vvoelz/biceps/blob/master/BICePs_2.0/tutorials/BICePs_example/BICePs_example.ipynb) first.




