<p align="center">
   <a href="https://stolllab.github.io/chiLife/">
      <img src="https://github.com/StollLab/chiLife/raw/main/img/chiLife_logo.png" width="300">
   </a>
</p>

<p align="center">
   <a href="https://stolllab.github.io/chiLife/">
      <img src="https://img.shields.io/badge/docs-latest-brightgreen.svg" width="300">
   </a>
</p>


| **UPDATES**                                                                                                                                                                                                                                                                     |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 2025.02.17 : Check out our new examples for trajectory analysis [Example_12](https://stolllab.github.io/chiLife/gallery/12-Analyzing_MD_Simulations.html) and using OpenMM for repacking [Example_13](https://stolllab.github.io/chiLife/gallery/13-OpenMM_Score_Function.html). |
| 2024.12.26 : Changes to the "energy funciton" backend. forcefield, forgive factors, etc should all be parts of the energy function now, and not parts of a rotamer library  Checkout [Example_07](https://stolllab.github.io/chiLife/gallery/13-OpenMM_Score_Function.html)    |                                                                               
| 2024.08.15 : chiLife 1.0 released!                                                                                                                                                                                                                                              |
| 2024.08.15 : chiLife can now be used to create arbitrary peptides with natural and NCAAs. Checkout the `make_peptide` function.                                                                                                                                                 |
| 2024.08.15 : chiLife can make NCAA structures from smiles. Checkout `smiles2residue`. Note: Requires RDKit to be installed.                                                                                                                                                     | 
| chiLife now supports arbitrary backbone attachments including DNA and RNA labels and more!                                                                                                                                                                                      | 



# chiLife
chiLife is a Python package for modeling non-canonical amino acid side chain ensembles and using those ensembles to
predict experimental observables. Currently, it is focused primarily on site-directed spin labels (SDSLs). The goal of
chiLife is to provide a simple, flexible and interoperable Python interface to protein side chain ensemble modeling,
allowing for rapid  development of custom analysis and modeling pipelines. This is facilitated by the use of `RotamerEnsemble`
and  `SpinLabel` objects with standard interfaces for all supported side chain types, side chain modeling methods and 
protein modeling methods. Flexibility is achieved by allowing users to create and use custom `RotamerEnsemble` and 
`SpinLabel` objects as well as custom side chain modeling methods. Interoperability sought by interactions with other 
Python-based molecular modeling packages. This enables the use of experimental data, like double electron-electron 
resonance (DEER), in other standalone protein modeling applications that allow user-defined restraints, such as 
PyRosetta and NIH-Xplor.
 
## Getting Started
Stable distributions of chiLife can be installed using `pip`. 
```bash
pip install chiLife
```

Alternatively, the development version can be installed by downloading and unpacking the GitHub repository, or using 
`git clone` followed by a standard Python setuptools installation.
```bash
git clone https://github.com/StollLab/chiLife.git
cd chiLife
pip install -e .   # Install as editable and update using `git pull origin main`
```  
***
## chiLife Module
The central entity of chiLife is the `SpinLabel` object, which inherits from the more abstract `RotamerEnsemble` 
object. While most people will primarily use `SpinLabel` objects, be aware that most properties and functions 
discussed are also functional on `RotamerLibrary` objects as well. `SpinLabel` objects can be created and "attached" to 
protein models easily and quickly, allowing for scriptable analysis or on-the-fly simulation of distance distributions while modeling.
Attaching a `SpinLabel` to a protein does not alter the protein in any way, allowing the protein model to retain the native amino acid.

### Simple rotamer-library based SpinLabel modeling

```python
import numpy as np
import matplotlib.pyplot as plt
import chilife as xl

# Download protein structure from PDB
MBP = xl.fetch('1omp', save=True)

# Create Spin lables
SL1 = xl.SpinLabel('R1M', site=20, chain='A', protein=MBP)
SL2 = xl.SpinLabel('R1M', site=238, chain='A', protein=MBP)

# Calculate distribution
r = np.linspace(0, 100, 256)
P = xl.distance_distribution(SL1, SL2, r=r)

# Plot distribution
fig, ax = plt.subplots(figsize=(6, 3))
ax.plot(r, P)
ax.set_yticks([])
ax.set_xlabel('Distance ($\AA$)')
for spine in ['left', 'top', 'right']:
    ax.spines[spine].set_visible(False)
plt.show()
```

![MBP L20R1 S238R1](https://github.com/StollLab/chiLife/raw/main/img/L20R1_S238R1_Pr.png)

The side chain ensembles can then be saved using a simple `save` function that accepts an arbitrary number of `
RotamerEnsemble`, `SpinLabel`, `MDAnalyisis.Universe` and `MDAnalyiss.AtomGroup` objects. Because 
`RotamerEnsemble`/`SpinLabel` objects do not mutate the underlying protein, they are saved as separate multi-state 
objects and can be visualized with applications like PyMOL. If you do wish to permanently alter the underlying protein
structure, you can use the [`mutate`](#mutating-protein-structures) function described below.

```
# Save structure
xl.save('MBP_L20R1_S238R1.pdb', SL1, SL2, MBP)
```

![MBP L20R1 S238R1 Structure](https://github.com/StollLab/chiLife/raw/main/img/L20R1_S238R1_Structure.png)

### Off-rotamer sampling and local repacking 
One of the benefits of chiLife is the variety and customizable nature of spin label modeling methods. This includes 
methods to sample side chain conformations that deviate from canonical dihedral angles and fixed rotamer libraries 
(off-rotamer sampling) and methods to repack a `SpinLabel` and its neighboring amino acids.

```python
import chilife as xl

MBP = xl.fetch('1omp')

# Create a SpinLabel object using the MTSSLWizard 'Accessible Volume' Approach
SL1 = xl.SpinLabel.from_wizard('R1M', site=20, chain='A', protein=MBP)

# Create a SpinLabel object by sampling off-rotamer dihedral conformations using the rotamer library as a prior 
SL2 = xl.SpinLabel('R1M', site=238, chain='A', sample=2000, protein=MBP)

# Create a SpinLabel object from a ProEPR.repack trajectory
traj, de = xl.repack(SL1, SL2, protein=MBP)
```

The `repack` function performs a Markov chain Monte Carlo sampling (MCMC) repack of the spin labels, `SL1` and `SL2` and 
neighboring side chains, returning an `MDAnalysis.Universe` object containing all accepted structures of the MCMC 
trajectory, the energy function changes at each acceptance step and new SpinLabel objects attached to the lowest-energy 
structure of the trajectory.

SpinLabel objects and neighboring side chains can be repacked using off-rotamer sampling by using the `off_rotamer=True`
option. In the event off-rotamer sampling is being used for repacking, it is likely that the desired `SpinLabel` object is 
not the default rotamer ensemble attached to the lowest-energy structure, but instead the ensemble of side chains 
created in the MCMC sampling trajectory. This can be done using the `from_trajectory` class method. 

```python
# Create a SpinLabel object from a xl.repack trajectory with off-rotamer sampling
traj, de = xl.repack(SL1, SL2, protein=MBP, off_rotamer=True) 
SL1 = xl.SpinLabel.from_trajectory(traj, site=238)
```

>Note: if you are creating a SpinLabel object from a label that is unknown to chilife you will have to specify which
> atoms the spin density primarily resides on. this is done with the ``spin_atoms`` kwarg, e.g.
> ```python
> SL1 = xl.SpinLabel.from_trajectory(traj, site=238, spin_atoms=['N1', 'O1'])
> ```

When repacking, off-rotamer sampling can be controlled for each dihedral angle separately by passing a list of bools to 
the off_rotamer keyword. For example, passing `off_rotamer = [False, False, False, True, True]` will allow for off-rotamer
sampling of only &chi;<sub>4</sub> and &chi;<sub>5</sub>.


### Mutating protein structures
Sometimes you don't want an entire rotamer ensemble, but just a protein structure mutated at a particular site with 
the most probable spin label conformer. This can be done easily with the `mutate` function.

```python
import chilife as xl

MBP = xl.fetch('1omp')
SL = xl.SpinLabel('R1M', 238, protein=MBP)
MBP_S238R1 = xl.mutate(MBP, SL)
xl.save('MBP_S238R1.pdb', MBP_S238R1)
```

chiLife can mutate several sites at once, and it can mutate canonical amino acids as well.

```python
SL1 = xl.SpinLabel('R1M', 20, protein=MBP)
SL2 = xl.SpinLabel('R1M', 238, protein=MBP)
L284V = xl.RotamerEnsemble('VAL', 284, protein=MBP)
```

Mutating adjacent sites is best done with the `repack` function to avoid clashes between `SpinLabels`/`RotamerEnsembles`. 
This returns a trajectory which can be used to pick the last or lowest-energy frame as the mutated protein.

```python
MBP_L284V_L20R1_S238R1, _, _ = xl.repack(SL1, SL2, L284V, protein=MBP)
```

### Adding user defined spin labels
Site-directed spin labels, and other non-canonical amino acids, are constantly being developed. Additionally, rotamer 
libraries for existing labels continuously undergo incremental improvements or modification to suit particular needs, 
e.g. a rotamer library specifically for transmembrane residues. In fact chiLife itself can be used to develop 
new and improved, or application-specific rotamer libraries. To this end, chiLife makes it easy to create user-defined 
spin labels and custom rotamer libraries. To create a custom rotamer library, all that is needed is (1) a pdb file of 
the spin label, (2) a list of the rotatable dihedral bonds, and (3) a list of the atoms that carry spin density.

```python
xl.create_library(name='TRT_1.0',
                  resname='TRT',
                  pdb='test_data/trt.pdb',
                  dihedral_atoms=[['N', 'CA', 'CB', 'SG'],
                                  ['CA', 'CB', 'SG', 'SD'],
                                  ['CB', 'SG', 'SD', 'CAD'],
                                  ['SG', 'SD', 'CAD', 'CAE'],
                                  ['SD', 'CAD', 'CAE', 'OAC']],
                  spin_atoms='CAQ')
```

This function creates a portable `TRT_1.0_rotlib.npz` file that can be provided when creating a `SpinLabel` object.

```python
xl.SpinLabel('TRT', site=238, protein=MBP, rotlib='TRT_1.0', sample=5000)
```

Thus, the file can be easily shared with coworkers, collaborators or with other chiLife users.

> NOTE: In the above example, the `rotlib` keyword is only used for demonstration purposes. chiLife always searches the 
> current working directory for rotamer library files first. If there is a `XYZ_rotlib.npz` in the working directory 
> and you specify `xl.SpinLabel('XYZ', ...)`, chiLife will assume you want to use the `XYZ_rotlib.npz` rotamer library.

User-defined labels can be constructed from a single-state pdb file or a multi-state PDB file. If constructed from a 
single-state pdb file, a list of dihedral angles and weights can be passed via the `dihedrals` and `weigts` keyword
arguments. For each set of dihedral angles, chiLife create a rotamer and store the whole library using the specified 
name. Alternatively, using a multi-state PDB file can add some additional information, such as isomeric heterogeneity of 
the rotamer library.

For more information on how to use chiLife see [examples](https://stolllab.github.io/chiLife/examples.html) and the 
[workshop](http://github.com/mtessmer/chilife_workshop) repository.

## References

When you are using chiLife in your work, please cite:

>Tessmer, M.H.; Stoll, S. chiLife: An open-source Python package for in silico spin labeling and integrative 
> protein modeling. Plos Comput Biol. 2023, 19:e1010834. https://doi.org/10.1371/journal.pcbi.1010834

If you are using off-rotamer sampling, please cite:

> Tessmer, M.H.; Canarie, E.R.; Stoll, S. Comparative evaluation of spin label modeling methods for protein 
> structural studies. Biophys J. 2022, 121, 3508-3519. https://doi.org/10.1016/j.bpj.2022.08.002

When using bifunctional label modeling, please cite:

> Tessmer, M.H.; Stoll, S. A Rotamer Library Approach to Modeling Side Chain Ensembles of the Bifunctional 
> Spin Label RX. Appl. Magn. Reson. 2023, 55, 127–140. https://doi.org/10.1007/s00723-023-01576-1
>
> Hasanbasri, Z.; Tessmer, M.H.; Stoll, S.; Saxena, S. Modeling of Cu(ii)-based protein spin labels
> using rotamer libraries. Phys. Chem. Chem. Phys. 2024, 26, 6806-6816. https://doi.org/10.1039/D3CP05951K

Note than many rotamer libraries have their own references. Please use the ``chilfe.rotlib_info()`` function 
on the rotamer libraries to check if there are any additional citations that should be referenced.
