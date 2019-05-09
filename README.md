# Detelectro

This repository contains a software I developed following my MSc in Physics Engineering, where we applied determinant Quantum Monte Carlo (QMC) to a problem of interacting fermions.
If you are interested, my thesis is available at

[Thesis](https://github.com/fmonteir/publications-fmobrito/blob/master/MScThesisFranciscoBrito2018-final.pdf)

The aim of my MSc project was to carry out a theoretical study (with particular emphasis on numerical aspects) of the properties of a TMD (transition metal dichalcogenide) nanoribbon - a graphene-like 2D nanostructure where electron interactions are particularly relevant - using a QMC method.

## Quantum Monte Carlo for Interacting Fermions

Capturing the effects of electron correlations is not an easy task. The difficulty lies in devising a numerical method to solve the many-body Schrödinger equation in a reasonable amount of computer time. The development and application of unbiased methods is a central point in correlated electron systems, particularly in 2D. Naive methods have exponential complexity in the system size, which motivates an approach based on Monte Carlo sampling.

Quantum Monte Carlo (QMC) methods are among the few unbiased methods available to date. In general, they circumvent the exponential complexity hurdle making it algebraic instead. However, for fermionic systems, a sign oscillation deems the algorithm exponential in the system size and inverse temperature, hence not very effective. This is due to the antisymmetric nature of the many-fermion wavefunction. Nowadays, a very relevant question in the field of strongly correlated electrons is whether the sign problem plagues the simulations of a given model, and whether it is possible to circumvent it.
For example, one of the main tasks of my MSc project was to investigate whether this issue would impede the simulation of a minimal model of 2D nanostructures made out of novel graphene-like 2D materials.

Some methods of dealing with the sign problem exist. Of course, they perform differently depending on the problem at hand. However, it might be the case that the sign problem is either absent, as it happens for a class of models (for a given range of parameters), or is not very severe, still allowing an unbiased and accurate simulation of the system at hand. Going beyond the fermion sign problem barrier in quantum simulations is an open topic of research and many possibilities have been put forward recently, namely an approach based on neural networks (see thesis and references therein).

In short, we implement an algorithm to deal with tight-binding problems for 2D interacting electronic models.
An example of an interesting property that we wish to  investigate is magnetism. One might be interested in the different phases that arise within the system and in how do the transitions between them occur.

## Getting started

To run our auxiliary field QMC code simply compile the code using


make


Then, you can run a simulation by running the command


./simulation 1 16 4.65 512 64 2


which corresponds to a edge-antiferromagnetic phase of a TMD nanoribbon.

The arguments of *simulation* are described below.

To maximize the efficiency of the code, some simulation parameters must be known at compile time. To set them, you must provide arguments when running _make_. If you simply run _make_, you will obtain the default parameters that give the aforementioned results.


To change the number of sites, inverse Trotter error, inverse temperature, the frequency of recomputing the Green's functions, the used source file, as well whether informations on the run are printed as the code runs, type:


make clean


make nsites=\<Number of sites\> dt_inv=\<Inverse Trotter Error\>beta=\<Inverse Temperature\> green_afresh_freq=\<Frequency of Recomputing G\> verbose=\<0 - no or 1 - yes> source=<source_file>


To run another simulation, simply type ./simulation followed by its arguments:


./simulation \<tmd\> \<U\> \<mu\> \<Total Number of Sweeps (Space-Time)\> \<Number of Warm-up Sweeps (Space-Time)\>  \< Number of Auto-correlation Sweeps (Space-Time) \>

where the first argument is the TMD (if you choose 1, you get MoS2), the second one is the on-site interaction, followed by the chemical potential.
The last three parameters are related to the Monte Carlo method.

The results of the simulation will be saved in a directory named _tmp_
and will be deleted once you run _make clean_. To save them, instead of using the
makefile explicily, use the _run.sh_ script. You start by setting the parameters.
For example:

TMD=1
NX=10
NY=5
BETA=8
U=16
MU=4.65

Then run!

cd examples

In the _examples/plot-src_ directory there are python scripts that use *matplotlib*
to plot the results of the simulations saved in _examples/data_

## Built with

*C++*

*Python*

## Contributing

*Francisco Monteiro de Oliveira Brito*
*João M. V. P. Lopes*
*Eduardo V. Castro*

## Versioning

v1 - latest update 09.05.2019

## Authors

*Francisco Monteiro de Oliveira Brito*

## Acknowledgements

*João M. V. P. Lopes, Eduardo V. Castro* - supervisors
