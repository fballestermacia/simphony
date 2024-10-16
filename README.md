<a name="readme-top"></a>
![Phonnier_Logo](images/Phonnier_Logo.png)
<div align="center">



<h1 align="center">Welcome to Phonnier!</h3>
 <h3 align="center"><ins>Phon</ins>on Wan<ins>nier</ins></h3>

This code tries to adapt the Wannier Orbital formalism to the calculation of topological propperties of Phonon systems, with emphasis on polar materials, i.e., in the case of LO-TO splitting.

</div>

<!-- Index -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#what-is-phonnier">What is Phonnier?</a>
    </li>
    <li>
      <a href="#how-to-use-phonnier">How to use Phonnier</a>
    </li>
    <li>
	<a href="#capabilities">Capabilities</a>
    </li>
    <li><a href="#examples">Examples</a></li>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
  </ol>
</details>



## What is Phonnier?

Phonnier stants for Phonnon Wannier, a code based on [WannierTools](http://www.wanniertools.com). Phonnier intends to build a Tight-Binding model for phonon systems. 

Since the maximally localized wannier functions for phonons are just delta functions (3 at each atomic position, representing directions x, y and z), the Tight-Binding Hamiltonian is based on the Dynamical Matrix that can be obtained by a previous *ab initio* calculation.

The TB Hamiltonian is read from a phononTB_hr.dat that follows the format of the wannier90_hr.dat file. A simple [python script](utility_scripts/QE2TBDAT) that generates the phononTB_hr.dat file from a QuantumESPRESSO output can be found on the *utility_scripts* directory. The script parses the .dyn files using SSCHA's CellConstructor python library and, if pressent, also removes the long-range part of the dynamical matrix for any dynamical matrix outside of $\Gamma$. This last step is important, since the long range part cannot be interpolated and must be added afterwards. Therefore, if your system has this long-range term, set the 

```
LOTO_correction
```
tag to true.



On the *utility_scripts* folder you can also find some python scripts for plotting and for generating the pn.in input file.

## How to use Phonnier

Once you have installed Phonnier (see <a href="#installation">Installation</a>) and included the [bin](bin) folder in your PATH, you can run Phonnier by executing 
```
> pn.x
```
on the directory where you have your input files.

If you are running on multicore configuration, use
```
> mpirun -np 4 pn.x
```
instead.

The code will read your *pn.in* file and the corresponding *wannier90_hr.dat* file. It will run the calculations designated in your input file and write, in the same folder,a *PN.out* file detailing the execution of the program and the appropiate output files. In case of a failed execution or errors, the *PN.out* file contains all the information about how the code read the input parameters and what calculations executed. An error message will apear in this file in case of failure.

### Input Files

Phonnier only needs two input files: *pn.in* and *wannier90_hr.dat*. 

The latter, as explained above, can be generated automatically using the [QE2TBDAT](utility_scripts/QE2TBDAT) script found on the *utility_scripts* folder. To run it, simply use

```
> QE2TBDAT [name].dyn [# of irr. q points] (ASR) ()=optional
```
on the folder where you hav stored your *.dyn* files. It will substract the long-range part of the dynamical matrix (if pressent) and generate a *[name]_phononTB_hr.dat*. 

By default, Phonnier searches for a file named *wannier90_hr.dat* on the current directory but the filename can be changed by including the following namelist on your *pn.in* file:
```
&TB_FILE
Hrfile = 'NAMEOFYOURhr.datFILE'
/
```

The former, however, tells Phonnier about the calculations you want to make and the parameters for those calculations. It also includes the relevant structure information.

A generic pn.in file can be created using the [createPhonnierInput](utility_scripts/createPhonnierInput) script, that reads the crystal structure from your *.dyn* files and writes a *pn.in* file with the main necessary parameters. It also checks wether your system incudes LO-TO correction or not and writes the dielectric tensor and the born effective charges on the input file. 

For examples on input files, we refer to the *examples* directory.


## Capabilities
Phonnier is intended to be a multipurpose tool for calculating topological quantitites and band structures of phonon systems. What follows is a non-exhaustive list of the main calculations that Phonnier is capable of performing:
<ol>
    <li> Bulk band structure </li>
    <li> Slab band structure </li>
    <li> Wilson loops </li>
    <li> Surface states and Fermi arcs </li>
    <li> Weyl node localization and topological charge </li>

</ol>

For more details on how to perform these calculations, check the <a href="#examples">Examples</a> section.


## Examples

### Surface states on an Obstructed Atomic Band Representation: $\text{AgP}_2$ 

### Weyl nodes, topological charges and Fermi arcs: $\text{Al}_2\text{ZnTe}_4$



## Installation

### Prerequisites

To install Phonnier you need the following
<ol>
    <li> A Fortran compiler: gfortran or ifort </li>
    <li> Lapack and Blas libraries </li>
    <li> OPTIONAL: MPICH v>2.1.5 </li>
    <li> OPTIONAL: ARPACK </li>
</ol>

### Compilation

To compile Phonnier, clone the repository wherever you want your installation to be, using:
```
> git clone https://github.com/fballestermacia/phonnier.git
```
or download the zip file directly from the [github repository](https://github.com/fballestermacia/phonnier).

Then, on the *phonnier/src* directory, copy or rename the Makefile.* that best adjusts to your system and/or needs into a file named simply *Makefile*. Then run:

```
> make
```

Note: the makefiles are prepared for somewhat standard installations and cases, if this didn't work for your system, you might need to modify them slightly. Allthough if your installation wasn't standard, you probably knew what you were doing and don't need help when modifiying makefiles.


Then, inside the *phonnier/bin* directory, an executable file named *pn.x* should appear. Add the *phonnier/bin* folder to your PATH by including the following line

```
export PATH=PATHTOYOURPHONNIERINSTALLATION/phonnier/bin:$PATH
```
in your *.bashrc* file in your home directory.

With that, you should be ready to go! Thank you again for using Phonnier and check the <a href="#how-to-use-phonnier">How to use Phonnier</a> section for mor details.

## Roadmap

[x] Add LO-TO splitting to bulk

[x] Finish Testing Slab Systems

[x] Add the utility scripts to GitHub

[ ] ~~Finish~~ Write documentation

[ ] Add examples

[ ] Continue testing topological quantities with LO-TO


