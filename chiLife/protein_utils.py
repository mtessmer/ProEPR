import logging, os, urllib, pickle, itertools
from pathlib import Path
from typing import Set, List, Union, Tuple
from numpy.typing import ArrayLike
from dataclasses import dataclass
from collections import Counter, defaultdict
from functools import partial
import MDAnalysis
import numpy as np
from scipy.spatial import cKDTree
from MDAnalysis.core.topologyattrs import Atomindices, Resindices, Segindices, Segids
from MDAnalysis.topology import guessers
import MDAnalysis as mda
import freesasa
from .RotamerLibrary import RotamerLibrary
import chiLife
from .numba_utils import _ic_to_cart

# TODO: Align Atom and ICAtom names with MDA
# TODO: Implement Internal Coord Residue object

@dataclass
class ICAtom:
    """
    Internal coordinate atom class. Used for building the ProteinIC (internal coordinate proteine).

    :param name: str
        Atom name
    :param atype: str
        Atom type
    :param index: int
        Atom number
    :param resn: str
        Name of the residue that the atom belongs to
    :param resi: int
        The residue index/number that the atom belongs to
    :param atom_names: tuple
        The names of the atoms that define the Bond-Angle-Torsion coordinates of the atom
    :param bond_idx: int
        Index of the coordinate atom bonded to this atom
    :param bond: float
        Distance of the bond in angstroms
    :param angle_idx: int
        Index of the atom that creates an angle with this atom and the bonded atom.
    :param angle: float
        The value of the angle between this atom, the bonded atom and the angle atom in radians.
    :param dihedral_idx: int
        The index of the atom that defines the dihedral with this atom, the bonded atom and the angled atom.
    :param dihedral: float
        The value of the dihedral angle defined above in radians.
    :param dihederal_resi: int
        The residues index of the residue that the dihedral angle belongs to. Note that this is not necessarily the same
        residue as the atom, e.g. the location of the nitrogen atom of the i+1 residue often defines the phi dihedral
        angle of the ith residue.
    """
    name: str
    atype: str
    index: int
    resn: str
    resi: int
    atom_names: tuple

    bond_idx: int = np.nan
    bond: float = np.nan
    angle_idx: int = np.nan
    angle: float = np.nan
    dihedral_idx: int = np.nan
    dihedral: float = np.nan
    dihedral_resi: int = np.nan

    def __post_init__(self):
        """If no dihedral_resi is defined default to the atom residue"""
        if np.isnan(self.dihedral_resi):
            self.dihedral_resi: int = self.resi


@dataclass
class Atom:
    """Atom class for atoms in cartesian space."""

    name: str
    atype: str
    index: int
    resn: str
    resi: int
    coords: np.ndarray

    @property
    def x(self):
        return self.coords[0]

    @x.setter
    def x(self, x):
        self.coords[0] = x

    @property
    def y(self):
        return self.coords[1]

    @y.setter
    def y(self, y):
        self.coords[1] = y

    @property
    def z(self):
        return self.coords[1]

    @z.setter
    def z(self, z):
        self.coords[1] = z


class ProteinIC:

    def __init__(self, ICs, **kwargs):
        """
        Object collecting internal coords of atoms making up a protein.
        :param ICs: dict
            dictionary of dicts, of dicts containing ICAtom objects. Top dict specifies chain, middle dict specifies
            residue and bottom dict holds ICAtom objects.
        """

        # Convert to conventional dictionary in case default dict was being used
        self.ICs = {key1: {key2: {key3: ic
                                  for key3, ic in resi.items()}
                           for key2, resi in chain.items()}
                    for key1, chain in ICs.items()}

        # Store important ProteinIC variables
        self.atoms = [ic for chain in ICs.values() for resi in chain.values() for ic in resi.values()]
        self.resis = [key for i in self.ICs for key in self.ICs[i]]
        self.resnames = {j: next(iter(self.ICs[i][j].values())).resn for i in self.ICs for j in self.ICs[i]}
        self.chains = [name for name, chain in ICs.items() for resi in chain.values() for ic in resi.values()]
        self.chain_operators = kwargs.get('chain_operators', None)

        self.bonded_pairs = np.array(kwargs.get('bonded_pairs', None))

        if self.bonded_pairs is None or self.bonded_pairs.any() is None:
            self.nonbonded_pairs = None
        else:
            bonded_pairs = [(a, b) for a, b in self.bonded_pairs]
            possible_bonds = set(itertools.combinations(range(len(self.atoms)), 2))
            self.nonbonded_pairs = np.array(list(possible_bonds - set(bonded_pairs)))

        self.perturbed = False
        self._coords = None
        self.dihedral_defs = self.collect_dih_list()


    def __iter__(self):
        """Iterate over all atoms."""
        return iter(self.atoms)

    def copy(self):
        """Create a deep copy of an ProteinIC instance"""
        return ProteinIC(self.ICs)

    def __len__(self):
        """Return the number of atoms in the structure"""
        return len(self.atoms)

    @property
    def chain_operators(self):
        """Chain_operators are a set of coordinate transformations that can orient multiple chains that are not
        covalently linked. e.g. structures with missing residues or protein complexes"""
        return self._chain_operators

    @chain_operators.setter
    def chain_operators(self, op):
        """
        Assign or calculate operators for all chains.

        :param op: dict
            Dictionary containing an entry for each chain in the ProteinIC molecule. Each entry must contain a
            rotation matrix, 'mx' and translation vector 'ori'.
        """
        if op is None:
            logging.info('No protein chain origins have been provided. All chains will start at [0, 0, 0]')
            op = defaultdict(dict)
            for chain in self.ICs.keys():
                op[chain]['ori'] = np.array([0, 0, 0])
                op[chain]['mx'] = np.identity(3)

        self._chain_operators = op

    @property
    def coords(self):
        if (self._coords is None) or not self.perturbed:
            self._coords = self.to_cartesian()
            self.perturbed = False
        return self._coords 

    def set_dihedral(self, dihedrals, resi, atom_list, chain=None):
        """
        Set one or more dihedral angles of a single residue in internal coordinates for the atoms defined in atom list.

        :param dihedrals: float, ndarray
            Angle or array of angles to set the dihedral(s) to

        :param resi: int
            Residue number of the site being altered

        :param atom_list: ndarray, list, tuple
            Names or array of names of atoms involved in the dihedral(s)

        :return coords: ndarray
            ProteinIC object with new dihedral angle(s)
        """
        self._perturbed = True

        if chain is None and len(self.ICs) == 1:
            chain = list(self.ICs.keys())[0]
        elif chain is None and len(self.ICs) > 1:
            raise ValueError('You must specify the protein chain')

        dihedrals = np.atleast_1d(dihedrals)
        atom_list = np.atleast_2d(atom_list)

        deltas, done = {}, []
        for i, (dihedral, atoms) in enumerate(zip(dihedrals, atom_list)):
            atoms = tuple(atoms)
            ratoms = tuple(reversed(atoms))
            done.append(atoms)

            if ratoms in self.ICs[chain][resi]:
                delta = self.ICs[chain][resi][ratoms].dihedral - dihedral
                self.ICs[chain][resi][ratoms].dihedral = dihedral
                deltas[tuple(atoms[:3])] = delta
            else:
                raise ValueError(
                    f'Dihedral with atoms {atoms}, {ratoms} not found in chain {chain} on resi {resi} internal coordinates:\n' + 
                    '\n'.join([ic for chain in self.ICs.values() for resi in chain.values() for ic in resi]))

        # Check for additional atoms with the same dihedral stem (first three atoms) and also change accordingly
        for atom in self.ICs[chain][resi].values():
            # Skip any atoms that do not have the same stem or have already been set.
            if len(atom.atom_names[:-4:-1]) != 3 or tuple(atom.atom_names) in done or tuple(
                    atom.atom_names[::-1]) in done:
                continue

            # Rotate atoms with the same stem by the same rotation applied to the designated dihedral
            elif tuple(atom.atom_names[:-4:-1]) in deltas:
                atom.dihedral -= deltas[tuple(atom.atom_names[:-4:-1])]

        return self

    def get_dihedral(self, resi, atom_list, chain=None):
        """
        Get the dihedral angle(s) of one or more atom sets at the specified residue. Dihedral angles are returned in
        radians

        :param resi: int
            Residue number of the site being altered

        :param atom_list: ndarray, list, tuple
            Names or array of names of atoms involved in the dihedral(s)

        :return angles: ndarray
            array of dihedral angles corresponding to the atom sets in atom list
        """

        if chain is None and len(self.ICs) == 1:
            chain = list(self.ICs.keys())[0]
        elif chain is None and len(self.ICs) > 1:
            raise ValueError('You must specify the protein chain')

        atom_list = np.atleast_2d(atom_list)
        dihedrals = []
        for atoms in atom_list:
            atoms = tuple(atoms)
            ratoms = tuple(reversed(atoms))

            if ratoms in self.ICs[chain][resi]:
                dihedrals.append(self.ICs[chain][resi][ratoms].dihedral)
            else:
                raise ValueError(
                    f'Dihedral with atoms {atoms}, {ratoms} not found in chain {chain} on resi {resi} internal coordinates:\n' + 
                    '\n'.join([' '.join(list(ic) + [str(resi)]) for chain in self.ICs 
                                                                    for resi in self.ICs[chain] 
                                                                        for ic in self.ICs[chain][resi]]))
        return dihedrals[0] if len(dihedrals) == 1 else np.array(dihedrals)

    def to_cartesian(self):
        """
        Convert internal coordinates into cartesian coordinates.

        :return coord: ndarray
            Array of cartesian coordinates corresponding to ICAtom list atoms
        """
        coord_arrays = []
        for segid in self.chain_operators:
            # Prepare variables for numba compiled function
            ICArray = np.array([[ic.bond_idx, ic.angle_idx, ic.dihedral_idx,
                                 ic.bond, ic.angle, ic.dihedral] for
                                res in self.ICs[segid].values() for ic in res.values()])

            IC_idx_array, ICArray = ICArray[:, :3].astype(int), ICArray[:, 3:]
            cart_coords = _ic_to_cart(IC_idx_array, ICArray)

            # Apply chain operations if any exist
            if ~np.allclose(self.chain_operators[segid]['mx'], np.identity(3)) and \
               ~np.allclose(self.chain_operators[segid]['ori'], np.array([0, 0, 0])):
                cart_coords = cart_coords @ self.chain_operators[segid]['mx'] + \
                            self.chain_operators[segid]['ori']

            coord_arrays.append(cart_coords)

        return np.concatenate(coord_arrays)

    def save_pdb(self, filename: str, mode='w'):
        """
        Save a pdb structure file from a ProteinIC object

        :param filename: str
            Name of file to save

        :param mode: str
            file open mode
        """
        if 'a' in mode:
            with open(str(filename), mode, newline='\n') as f:
                f.write("MODEL\n")
            chiLife.save_pdb(filename, self.atoms, self.to_cartesian(), mode=mode)
            with open(str(filename), mode, newline='\n') as f:
                f.write("ENDMDL\n")
        else:
            chiLife.save_pdb(filename, self.atoms, self.to_cartesian(), mode=mode)

    def has_clashes(self, distance=1.5):
        """
        Checks for an internal clash between nonbonded atoms
        :param distance: float
            Minimum distance allowed between non-bonded atoms

        :return:
        """
        diff = self.coords[self.nonbonded_pairs[:, 0]] - self.coords[self.nonbonded_pairs[:, 1]]
        dist = np.linalg.norm(diff, axis=1)
        has_clashes = np.any(dist < distance)
        return has_clashes

    def get_resi_dihs(self, resi:int):
        """
        Gets the list of heavy atom dihedral definitions for the provided residue.

        :param resi: int
            Index of the residue whose dihedrals will be returned
        """

        if resi==0 or resi==list(self.ICs[1].keys())[-1]: 
            dihedral_defs = []
        elif resi==1:
            dihedral_defs = [["CH3","C", "N", "CA"],
                             ["C", "N", "CA",  "C"],
                             ["N", "CA", "C", "N"]]
        else:
            dihedral_defs = [["C", "N", "CA", "C"],
                             ["N", "CA", "C", "N"]]
        return dihedral_defs

    def collect_dih_list(self) -> list:
        """
        Returns a list of all heavy atom dihedrals

        :return dihs: list
            list of protein heavy atom dihedrals
        """
        dihs = []
    
        # Get backbone and sidechain dihedrals for the provided universe
        for resi in self.ICs[1]:
            resname = self.resnames[resi]
            res_dihs = self.get_resi_dihs(resi)
            res_dihs += chiLife.dihedral_defs.get(resname, [])
            dihs += [(resi, d) for d in res_dihs]

        return dihs


def get_ICAtom(atom: mda.core.groups.Atom, offset: int = 0, preferred_dihedral: List = None) -> ICAtom:
    """
    Construct internal coordinates for an atom given that atom is linked to an MDAnalysis Universe.

    :param atom: MDAnalysis.Atom
        Atom object to obtain internal coordinates for.

    :param offset: int
        Index offset used to construct internal coordinates for separate chains.

    :param preferred_dihedral: list
        Atom names defining the preferred dihedral to be use in the bond-angle-torsion coordinate system.

    :return: ICAtom
        Atom with internal coordinates.
    """
    U = atom.universe
    if (atom.index - offset) == 0:
        return ICAtom(atom.name, atom.type, atom.index - offset, atom.resname, atom.resid, (atom.name,))

    elif (atom.index - offset) == 1:
        atom_names = (atom.name, U.atoms[atom.bonds.indices[0][0]].name)
        return ICAtom(atom.name, atom.type, atom.index - offset, atom.resname, atom.resid, atom_names,
                      atom.bonds.indices[0][0] - offset, atom.bonds.bonds()[0])

    elif (atom.index - offset) == 2:
        atom_names = (atom.name, U.atoms[atom.bonds.indices[0][0]].name, U.atoms[atom.angles.indices[0][0]].name)
        return ICAtom(atom.name, atom.type, atom.index - offset, atom.resname, atom.resid, atom_names,
                      atom.bonds.indices[0][0] - offset, atom.bonds.bonds()[0],
                      atom.angles.indices[0][0] - offset, atom.angles.angles()[0])

    else:

        # Skip dihedrals that are not preferred
        i, j, k = 0, 0, 0
        if preferred_dihedral is not None:
            preferred_stems = [tuple(dh[:3]) for dh in preferred_dihedral]
            for dh_idx, dh in enumerate(atom.dihedrals):
                if any(stem == tuple(dh.atoms.names[:3]) for stem in preferred_stems):
                    if atom.name == dh.atoms.names[-1]:
                        k = dh_idx
                        break

        found = False
        while not found:

            # Check that the bond, angle and dihedral are all defined by the same atoms in the same order
            try:
                condition = atom.index == atom.angles.indices[j][-1]
            except:
                print('hw')
            condition = condition and all(atom.bonds.indices[i] == atom.angles.indices[j][-2:])
            condition = condition and all(atom.angles.indices[j] == atom.dihedrals.indices[k][-3:])
            condition = condition and all(atom.dihedrals.indices[k] >= offset)

            i_idx = atom.bonds.indices[i][0]
            j_idx = atom.angles.indices[j][0]
            k_idx = atom.dihedrals.indices[k][0]

            # If the atom being placed is part of a side chain dihedral that is near the backbone.
            if (atom.name not in ['N', 'CA', 'C', 'O', 'CB', 'CB1', 'CB2', 'CD']) and (atom.type != 'H'):

                # Check if the dihedral definition is coming from the carboxyl end of the residue
                condition = condition and (not any(
                    [x.name in ['C', 'O'] and x.resnum == atom.resnum for x in (atom.angles[j].atoms[0], atom.universe.atoms[k_idx])]))

            # If all checks have been pased use bond i, angle j, and dihedral k
            if condition:
                found = True

            # If all angles and dihedrals containing bond i have been searched increment bond
            elif k == len(atom.dihedrals.indices) - 1 and j == len(atom.angles.indices) - 1:
                k = 0
                j = 0
                i += 1
            # if all dihedrals containing bond i and angle j have been searched increment angle
            elif k == len(atom.dihedrals.indices) - 1:
                k = 0
                j += 1
            # If this dihedral does not contain atoms of the angle then try the next dihedral
            else:
                k += 1

        atom_names = (atom.name, U.atoms[i_idx].name, U.atoms[j_idx].name, U.atoms[k_idx].name)

        if atom_names != ('N', 'C', 'CA', 'N'):
            dihedral_resi = atom.resnum
        else:
            dihedral_resi = atom.resnum - 1

        return ICAtom(atom.name, atom.type, atom.index-offset, atom.resname, atom.resnum, atom_names,
                      i_idx - offset, atom.bonds.bonds()[i],
                      j_idx - offset, atom.angles.angles()[j],
                      k_idx - offset, atom.dihedrals.dihedrals()[k],
                      dihedral_resi=dihedral_resi)


def get_internal_coords(mol: Union[MDAnalysis.Universe, MDAnalysis.AtomGroup],
                        resname: str = None, preferred_dihedrals: List = None) -> ProteinIC:
    """
    Gather a list of Internal

    :param mol: MDAnalysis.Universe, MDAnalysis.AtomGroup
        Molecule to convert into internal coords.

    :param resname: str
        Residue name (3-letter code) of any non-canonical or otherwise unsupported amino acid that should be included in
        the ProteinIC object.

    :param preferred_dihedrals: list
        Atom names of  preffered dihedral definitions to be used in the bond-angle-torsion coordinate system. Often
        used to specify dihedrals of user defined or unsupported amino acids that the user wishes to directly interact
        with.

    :return ICs: ProteinIC
        An ProteinIC object of the supplied molecule
    """
    U = mol.universe
    # Making things up to get the right results but vdwradii should not determine bond length so talk to the MDA devs
    extra_radii = {'S': 2., 'Br': 4.189, 'Cu': 2.75}
    if not hasattr(U, 'bonds'):
        U.add_TopologyAttr('bonds', guessers.guess_bonds(U.atoms, U.atoms.positions, vdwradii=extra_radii))
    elif len(U.bonds) < len(U.atoms):
        U.add_TopologyAttr('bonds', guessers.guess_bonds(U.atoms, U.atoms.positions, vdwradii=extra_radii))
    if not hasattr(U, 'angles'):
        U.add_TopologyAttr('angles', guessers.guess_angles(U.bonds))
    if not hasattr(U, 'dihedrals'):
        U.add_TopologyAttr('dihedrals', guessers.guess_dihedrals(U.angles))
    if not hasattr(U, 'impropers'):
        U.add_TopologyAttr('impropers', guessers.guess_improper_dihedrals(U.angles))

    if resname is not None:
        protein = mol.select_atoms(f'(protein or resname {" ".join(list(chiLife.SUPPORTED_RESIDUES) + [resname])})'
                                 f' and not altloc B')
    else:
        protein = mol.select_atoms(f'(protein or resname {" ".join(list(chiLife.SUPPORTED_RESIDUES))}) and not altloc B')

    all_ICAtoms = defaultdict(partial(defaultdict, dict))
    chain_operators = defaultdict(dict)
    offset = 0
    segid = 0
    for atom in protein.atoms:

        chaind = atom.segid
        if atom.index == protein.atoms[0].index or atom.resid - U.atoms[atom.index - 1].resid not in (0, 1):
            segid += 1
            offset = atom.index
            mx, ori = chiLife.ic_mx(*U.atoms[offset:offset + 3].positions)

            chain_operators[segid]['ori'] = ori
            chain_operators[segid]['mx'] = mx

        ICatom = get_ICAtom(atom, offset=offset, preferred_dihedral=preferred_dihedrals)
        all_ICAtoms[segid][ICatom.dihedral_resi][ICatom.atom_names] = ICatom

    return ProteinIC(all_ICAtoms, chain_operators=chain_operators, bonded_pairs=U.bonds.to_indices())


def save_rotlib(name: str, atoms: ArrayLike, coords: ArrayLike = None) -> None:
    """
    Save a rotamer library as multiple states of the same molecule.

    :param name: str
        file name to save rotamer library to

    :param atoms: list, tuple
        list of Atom objects

    :param coords: np.ndarray
        Array of atom coordinates corresponding to Atom objects
    """

    if not name.endswith('.pdb'):
        name += '.pdb'

    if coords is None and isinstance(atoms[0], list):
        with open(name, 'w', newline='\n') as f:
            for i, model in enumerate(atoms):
                f.write(f'MODEL {i + 1}\n')
                for atom in model:
                    f.write(f"ATOM  {atom.index + 1:5d}  {atom.name:<4s}{atom.resn:3s} {'A':1s}{atom.resi:4d}   "
                            f"{atom.coords[0]:8.3f}{atom.coords[1]:8.3f}{atom.coords[2]:8.3f}{1.0:6.2f}{1.0:6.2f}        "
                            f"  {atom.atype:>2s}\n")
                f.write('ENDMDL\n')

    elif len(coords.shape) > 2:
        with open(name, 'w', newline='\n') as f:
            for i, model in enumerate(coords):
                f.write(f'MODEL {i + 1}\n')
                for atom, coord in zip(atoms, model):
                    f.write(f"ATOM  {atom.index + 1:5d}  {atom.name:<4s}{atom.resn:3s} {'A':1s}{atom.resi:4d}   "
                            f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{1.0:6.2f}{1.0:6.2f}          {atom.atype:>2s}\n")
                f.write('ENDMDL\n')

    else:
        save_pdb(name, atoms, coords)


def save_pdb(name: Union[str, Path], atoms: ArrayLike, coords: ArrayLike, mode: str = 'w') -> None:
    """
    Save a single state pdb structure of the provided atoms and coords

    :param name: str
        Name of file to save

    :param atoms: list, tuple
        List of Atom objects to be saved

    :param coords: np.ndarray
         Array of atom coordinates corresponding to atoms

    :param mode: str
        File open mode. Usually used to specify append ("a") when you want to add structures to a PDB rather than
        overwrite that pdb.
    """
    name = Path(name) if isinstance(name, str) else name
    name = name.with_suffix('.pdb')

    with open(name, mode, newline='\n') as f:
        for atom, coord in zip(atoms, coords):
            f.write(f"ATOM  {atom.index + 1:5d}  {atom.name:<4s}{atom.resn:3s} {'A':1s}{atom.resi:4d}   "
                    f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{1.0:6.2f}{1.0:6.2f}          {atom.atype:>2s}  \n")


def get_missing_residues(protein: Union[MDAnalysis.Universe, MDAnalysis.AtomGroup],
                         ignore: Set[int] = None,
                         use_H: bool = False) -> List:
    """
    Get a list of RotamerLibrary objects corresponding to the residues of the provided protein that are missing heavy
    atoms
    :param protein: MDAnalysis.Universe, MDAnalysis.AtomGroup
        Protein to search for residues with missing atoms.
    :param ignore: ArrayLike
        List of residue numbers to ignore. Usually sites you plan to label or mutate.
    :param use_H: bool
        Whether the new side chain should have hydrogen atoms
    :return: ArrayLike
        List of RotamerLibrary objects corresponding to residues with missing heavy atoms
    """
    ignore = set() if ignore is None else ignore
    missing_residues = []
    cache = {}

    for res in protein.residues:
        # Only consider supported residues because otherwise chiLife wouldn't know what's missing
        if res.resname not in chiLife.SUPPORTED_RESIDUES or res.resnum in ignore or res.resname in ['ALA', 'GLY']:
            continue

        # Check if there are any missing heavy atoms
        heavy_atoms = res.atoms.types[res.atoms.types != 'H']
        if len(heavy_atoms) != cache.get(res.resname, len(RotamerLibrary(res.resname).atom_names)):
            missing_residues.append(RotamerLibrary(res.resname, res.resnum, protein=protein, chain=res.segid, use_H=use_H))

    return missing_residues


def mutate(protein: MDAnalysis.Universe,
           *rotlibs: RotamerLibrary,
           add_missing_atoms: bool = True,
           random_rotamers: bool = False) -> MDAnalysis.Universe:
    """
    Create a new Universe where the native residue is replaced with the highest probability rotamer from a
    RotamerLibrary or SpinLabel object.

    :param protein: MDAnalysis.Universe
        Universe containing protein to be spin labeled

    :param rotlibs: RotamerLibrary, SpinLabel
        Precomputed RotamerLibrary or SpinLabel object to use for selecting and replacing the spin native amino acid

    :param random_rotamers: bool
        Randomize rotamer conformations

    :param add_missing_atoms:
        Remodel side chains missing atoms

    :return U: MDAnalysis.Universe
        New Universe with a copy of the spin labeled protein with the highest probability rotamer
    """
    if add_missing_atoms:
        if len(rotlibs) > 0 and all(not hasattr(lib, 'H_mask') for lib in rotlibs):
            use_H = True
        elif any(not hasattr(lib, 'H_mask') for lib in rotlibs):
            raise AttributeError('User provided some rotlibs with hydrogen atoms and some without. Make sure all '
                                 'rotlibs either do or do not use hydrogen')
        else:
            use_H = False

        missing_residues = get_missing_residues(protein, ignore={res.site for res in rotlibs}, use_H=use_H)
        rotlibs = list(rotlibs) + missing_residues

    label_sites = {int(spin_label.site): spin_label for spin_label in rotlibs}

    protein = protein.select_atoms(f'(protein or resname {" ".join(chiLife.SUPPORTED_RESIDUES)}) and not altloc B')
    label_selstr = " or ".join([f'({label.selstr})' for label in rotlibs])
    other_atoms = protein.select_atoms(f'not ({label_selstr})')

    # Get new universe information
    n_residues = len(other_atoms.residues) + len(rotlibs)
    n_atoms = len(other_atoms) + sum(len(spin_label.atom_names) for spin_label in rotlibs)
    resids = [res.resid for res in protein.residues]

    # Allocate lists for universe information
    atom_info = []
    res_names = []
    segidx = []

    # Loop over residues in old universe
    for i, res in enumerate(protein.residues):

        # If the residue is the spin labeled residue replace it with the highest probability spin label
        if res.resnum in label_sites:
            atom_info += [(i, name, atype) for name, atype in
                          zip(label_sites[res.resnum].atom_names, label_sites[res.resnum].atom_types)]

            # Add missing Oxygen from rotamer libraries
            res_names.append(label_sites[res.resnum].res)
            segidx.append(label_sites[res.resnum].segindex)

        # Else retain the atom information from the parent universe
        else:
            atom_info += [(i, atom.name, atom.type) for atom in res.atoms if atom.altLoc != 'B']
            res_names.append(res.resname)
            segidx.append(res.segindex)

    # Unzip atom information into individual lists
    residx, atom_names, atom_types = zip(*atom_info)
    segids = list(Counter(protein.residues.segids))
    # Allocate a new universe with the appropriate information
    U = mda.Universe.empty(n_atoms, n_residues=n_residues, atom_resindex=residx,
                           residue_segindex=segidx, trajectory=True)

    # Add necessary topology attributes
    U.add_TopologyAttr('name', atom_names)
    U.add_TopologyAttr('type', atom_types)
    U.add_TopologyAttr('resname', res_names)
    U.add_TopologyAttr('resid', resids)
    U.add_TopologyAttr('altLoc', ['A' for atom in range(n_atoms)])
    U.add_TopologyAttr('resnum', resids)
    U.add_TopologyAttr('segid')

    for i, segid in enumerate(segids):
        if i == 0:
            i_segment = U.segments[0]
            i_segment.segid=segid
        else:
            i_segment = U.add_Segment(segid=str(segid))

        mask = np.argwhere(np.asarray(segidx) == i).squeeze()
        U.residues[mask.tolist()].segments = i_segment

    U.add_TopologyAttr(Segids(np.array(segids)))
    U.add_TopologyAttr(Atomindices())
    U.add_TopologyAttr(Resindices())
    U.add_TopologyAttr(Segindices())

    # Apply old coordinates to non-spinlabel atoms
    new_other_atoms = U.select_atoms(f'not ({label_selstr})')
    new_other_atoms.atoms.positions = other_atoms.atoms.positions

    # Apply most probable spin label coordinates to spin label atoms
    for spin_label in label_sites.values():
        sl_atoms = U.select_atoms(spin_label.selstr)
        if random_rotamers:
            sl_atoms.atoms.positions = spin_label.coords[np.random.choice(len(spin_label.coords), p=spin_label.weights)]
        else:
            sl_atoms.atoms.positions = spin_label.coords[np.argmax(spin_label.weights)]

    return U


def randomize_rotamers(protein: Union[mda.Universe, mda.AtomGroup],
                       rotamer_libraries: List[RotamerLibrary],
                       **kwargs) -> None:
    """
    Modify a protein object in place to randomize side chain conformations.

    :param protein: MDAnalysis.Universe, MDAnalysis.AtomGroup
        Protein object to modify.

    :param rotamer_libraries: list
        RotamerLibrary objects attached to the protein corresponding to the residues to be repacked/randomized.
    """
    for rotamer in rotamer_libraries:
        coords, weight = rotamer.sample(off_rotamer=kwargs.get('off_rotamer', False))
        mask = ~np.isin(protein.ix, rotamer.clash_ignore_idx)
        protein.atoms[~mask].positions = coords


def get_sas_res(protein: Union[mda.Universe, mda.AtomGroup], cutoff: float = 30) -> Set[Tuple[int, str]]:
    """
    Run FreeSASA to get solvent accessible surface residues in the provided protein

    :param protein: MDAnalysis.Universe, MDAnalysis.AtomGroup
        Protein object to measure Solvent Accessible Surfaces (SAS) area of and report the SAS residues.

    :param cutoff:
        Exclude residues from list with SASA below cutoff in angstroms squared.

    :return: SAResi
        Set of solvent accessible surface residues
    """
    freesasa.setVerbosity(1)
    # Create FreeSASA Structure from MDA structure
    fSASA_structure = freesasa.Structure()

    [fSASA_structure.addAtom(atom.name, atom.resname, str(atom.resnum), atom.segid,
                             *atom.position) for atom in protein.atoms]

    residues = [(resi.segid, resi.resnum) for resi in protein.residues]

    # Calculate SASA
    SASA = freesasa.calc(fSASA_structure)
    SASA = freesasa.selectArea((f'{resi}_{chain}, resi {resi} and chain {chain}' for chain, resi in residues),
                               fSASA_structure, SASA)

    SAResi = [key.split('_') for key in SASA if SASA[key] >= cutoff]
    SAResi = [(int(R), C) for R, C in SAResi]

    return set(SAResi)


def fetch(pdbid: str, save: bool = False) -> MDAnalysis.Universe:
    """
    Fetch pdb file from the protein data bank and optionally save to disk.

    :param pdbid: str
        4 letter structure ID

    :param save: bool
        If true the fetched PDB will be saved to the disk.

    :return U: MDAnalysis.Universe
        MDAnalysis Universe object of the protein corresponding to the provided PDB ID
    """

    if not pdbid.endswith('.pdb'):
        pdbid += '.pdb'

    urllib.request.urlretrieve(f'http://files.rcsb.org/download/{pdbid}', pdbid)

    U = mda.Universe(pdbid, in_memory=True)

    if not save:
        os.remove(pdbid)

    return U


def presort_key(pdb_line: str) -> Tuple[str, int, int]:
    """
    Assign a base rank to sort atoms of a pdb.

    :param pdb_line: str
        ATOM line from a pdb file as a string.

    :return: chainid, resid, name_order
        ordered ranking of atom for sorting the pdb.
    """
    chainid = pdb_line[21]
    res_name = pdb_line[17:20].strip()
    resid = int(pdb_line[24:27].strip())
    atom_name = pdb_line[12:17].strip()
    atom_type = pdb_line[76:79].strip()
    if res_name == 'ACE':
        name_order = {'CH3':0, 'C': 1, 'O': 2}.get(atom_name, 4) if atom_type != 'H' else 5
    else:
        name_order = atom_order.get(atom_name, 4) if atom_type != 'H' else 5

    return chainid, resid, name_order


def sort_pdb(pdbfile: Union[str, List]) -> List[str]:
    """
    Read ATOM lines of a pdb and sort the atoms according to chain, residue index, backbone atoms and side chain atoms.
    Side chain atoms are sorted by distance to each other/backbone atoms with atoms closest to the backbone coming
    first and atoms furthest from the backbone coming last. This sorting is essential to making internal-coordinates
    with consistent and preferred dihedral definitions.

    :param pdbfile: str, list
        Name of the PDB file or a list of strings containing ATOM lines of a PDB file

    :return lines: list
        Sorted list of strings corresponding to the ATOM entries of a PDB file.
    """

    if isinstance(pdbfile, str):
        with open(pdbfile, 'r') as f:
            lines = f.readlines()

        start_idxs = []
        end_idxs = []

        for i, line in enumerate(lines):
            if 'MODEL' in line:
                start_idxs.append(i)
            elif 'ENDMDL' in line:
                end_idxs.append(i)

        if start_idxs != []:
            return [sort_pdb(lines[s:e]) for s, e in zip(start_idxs, end_idxs)]

    elif isinstance(pdbfile, list):
        lines = pdbfile

    lines = [line for line in lines if line.startswith('ATOM')]

    # Presort
    lines.sort(key=presort_key)

    coords = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in lines])

    # get residue groups
    chain, resi = lines[0][21], int(lines[0][24:27].strip())
    start = 0
    resdict = {}
    for curr, pdb_line in enumerate(lines):

        if chain != pdb_line[21] or resi != int(pdb_line[24:27].strip()):
            resdict[chain, resi] = start, curr
            start = curr
            chain, resi = pdb_line[21], int(pdb_line[24:27].strip())

    resdict[chain, resi] = start, curr + 1
    midsort_key = []
    for key in resdict:
        start, stop = resdict[key]
        kdtree = cKDTree(coords[start:stop])

        # Get all nearest neighbors and sort by distance
        pairs = kdtree.query_pairs(2.2, output_type='ndarray')
        distances = np.linalg.norm(kdtree.data[pairs[:, 0]] - kdtree.data[pairs[:, 1]], axis=1)
        idx_sort = np.argsort(distances)
        pairs = pairs[idx_sort]

        idx = 1
        sorted_args = list(range(np.minimum(4, stop-start)))
        i = 0

        while len(sorted_args) < stop - start:
            if 'search_len' in locals():
                appendid = []
                for idx in sorted_args[-search_len:]:
                    for pair in pairs:
                        if idx in pair:
                            ap = pair[0] if pair[0] != idx else pair[1]
                            if ap not in sorted_args and ap not in appendid:
                                appendid.append(ap)
                if appendid != []:
                    # appendid = list(set(appendid))
                    # tmp = [lines[i + start][15:17].strip() for i in appendid]
                    # tmp = np.argsort(tmp)
                    # appendid = [appendid[i] for i in tmp]
                    sorted_args += appendid
                    search_len = len(appendid)
                else:
                    search_len += 1

            else:
                appendid = []
                for pair in pairs:
                    if idx in pair:
                        ap = pair[0] if pair[0] != idx else pair[1]
                        if ap not in sorted_args:
                            appendid.append(ap)

                if appendid != []:
                    # tmp = [lines[i + start][14:17].strip() for i in appendid]
                    # tmp = np.argsort(tmp)
                    # appendid = [appendid[i] for i in tmp]
                    sorted_args += appendid
                    search_len = len(appendid)
                else:
                    pass
                    # idx += 1
        midsort_key += [x + start for x in sorted_args]

    lines[:] = [lines[i] for i in midsort_key]
    lines.sort(key=presort_key)
    lines = [line[:6] + f"{i + 1:5d}" + line[11:] for i, line in enumerate(lines)]

    return lines

# Define rotamer dihedral angle atoms
with open(os.path.join(os.path.dirname(__file__), 'data/DihedralDefs.pkl'), 'rb') as f:
    dihedral_defs = pickle.load(f)

with open(os.path.join(os.path.dirname(__file__), 'data/rotamer_libraries/RotlibIndexes.pkl'), 'rb') as f:
    rotlib_indexes = pickle.load(f)

atom_order = {'N': 0, 'CA': 1, 'C': 2, 'O': 3}
