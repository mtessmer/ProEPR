import networkx as nx
from joblib import Parallel, delayed
from typing import Union
from itertools import combinations
import logging
import warnings
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
import scipy.optimize as opt
import MDAnalysis as mda

import chilife


class dRotamerEnsemble:
    backbone_atoms = ["H", "N", "CA", "HA", "C", "O"]

    def __init__(self, res, sites, protein=None, chain=None, rotlib=None, **kwargs):
        """ """
        self.res = res
        self.site1, self.site2 = sorted(sites)
        self.site = self.site1
        self.increment = self.site2 - self.site1
        self.kwargs = kwargs

        self.protein = protein
        self.chain = chain if chain is not None else self.guess_chain()
        self.protein_tree = self.kwargs.setdefault("protein_tree", None)

        self.forgive = kwargs.setdefault("forgive", 0.95)
        self.clash_radius = kwargs.setdefault("clash_radius", 14.0)
        self._clash_ori_inp = kwargs.setdefault("clash_ori", "cen")
        self.restraint_weight = kwargs.pop("restraint_weight") if "restraint_weight" in kwargs else 222  # kcal/mol/A^2
        self.alignment_method = kwargs.setdefault("alignment_method", "bisect".lower())
        self.dihedral_sigmas = kwargs.setdefault("dihedral_sigmas", 25)
        self._exclude_nb_interactions = kwargs.setdefault('exclude_nb_interactions', 3)
        self._minimize = kwargs.pop("minimize", True)
        self.min_method = kwargs.pop('min_method', 'L-BFGS-B')
        self.eval_clash = kwargs.pop("eval_clash", True)
        self.energy_func = kwargs.setdefault("energy_func", chilife.get_lj_energy)
        self.temp = kwargs.setdefault("temp", 298)
        self.get_lib(rotlib)
        self.create_ensembles()

        self.RL1.backbone_to_site()
        self.RL2.backbone_to_site()

        self.cst_idx1 = np.where(self.RL1.atom_names[None, :] == self.csts[:, None])[1]
        self.cst_idx2 = np.where(self.RL2.atom_names[None, :] == self.csts[:, None])[1]
        self.rl1mask = np.argwhere(~np.isin(self.RL1.atom_names, self.csts)).flatten()
        self.rl2mask = np.argwhere(~np.isin(self.RL2.atom_names, self.csts)).flatten()

        self.name = self.res
        if self.site1 is not None:
            self.name = f"{self.RL1.nataa}{self.site1}{self.RL1.nataa}{self.site2}{self.res}"
        if self.chain is not None:
            self.name += f"_{self.chain}"

        self.selstr = (
            f"resid {self.site1} {self.site2} and segid {self.chain} and not altloc B"
        )

        self._graph = nx.Graph()
        self._graph.add_edges_from(self.bonds)

        self.protein_setup()
        self.sub_labels = (self.RL1, self.RL2)

    def protein_setup(self):
        if isinstance(self.protein, (mda.AtomGroup, mda.Universe)):
            if not hasattr(self.protein.universe._topology, "altLocs"):
                self.protein.universe.add_TopologyAttr('altLocs', np.full(len(self.protein.universe.atoms), ""))

        self.protein = self.protein.select_atoms("not (byres name OH2 or resname HOH)")
        self.clash_ignore_idx = self.protein.select_atoms(
            f"resid {self.site1} {self.site2} and segid {self.chain}"
        ).ix

        self.resindex = self.protein.select_atoms(self.selstr).residues[0].resindex
        self.segindex = self.protein.select_atoms(self.selstr).residues[0].segindex

        if self.protein_tree is None:
            self.protein_tree = cKDTree(self.protein.atoms.positions)

        protein_clash_idx = self.protein_tree.query_ball_point(
            self.clash_ori, self.clash_radius
        )
        self.protein_clash_idx = [
            idx for idx in protein_clash_idx if idx not in self.clash_ignore_idx
        ]

        _, self.irmin_ij, self.ieps_ij, _ = chilife.prep_internal_clash(self)
        _, self.ermin_ij, self.eeps_ij, _ = chilife.prep_external_clash(self)

        self.aidx, self.bidx = [list(x) for x in zip(*self.non_bonded)]

        if self._minimize:
            self.minimize()

        if self.eval_clash:
            self.evaluate()

    def guess_chain(self):
        if self.protein is None:
            chain = "A"
        elif len(set(self.protein.segments.segids)) == 1:
            chain = self.protein.segments.segids[0]
        elif np.isin(self.protein.residues.resnums, self.site1).sum() == 0:
            raise ValueError(
                f"Residue {self.site1} is not present on the provided protein"
            )
        elif np.isin(self.protein.residues.resnums, self.site1).sum() == 1:
            chain = self.protein.select_atoms(f"resid {self.site1}").segids[0]
        else:
            raise ValueError(
                f"Residue {self.site1} is present on more than one chain. Please specify the desired chain"
            )
        return chain

    def get_lib(self, rotlib):
        if self.protein is not None:
            # get site1 backbone information from protein structure
            sel_txt = f"resnum {self.site1} and segid {self.chain}"


        rotlib = self.res if rotlib is None else rotlib
        if 'ip' not in rotlib:
            rotlib += f'ip{self.increment}'

        # Check if any exist
        rotlib_path = get_possible_rotlibs(rotlib)

        if rotlib_path is None:
            # Check if libraries exist but for different i+n
            rotlib_path = get_possible_rotlibs(rotlib.replace(f'ip{self.increment}', ''), all=True)

            if rotlib_path is None:
                raise NameError(f'There is no rotamer library called {rotlib} in this directory or in chilife')

            else:
                warnings.warn(f'No rotamer library found for the given increment (ip{self.increment}) but rotlibs '
                              f'were found for other increments. chiLife will combine these rotlib to model '
                              f'this site1 pair and but they may not be accurate! Because there is no information '
                              f'about the relative weighting of different rotamer libraries all weights will be '
                              f'set to 1/len(rotlib)')

        if isinstance(rotlib_path, Path):
            libA, libB, csts = chilife.read_library(rotlib_path)

        elif isinstance(rotlib_path, list):
            concatable = ['dihedrals', 'internal_coords']

            cctA, cctB, ccsts = {}, {}, {}
            libA, libB, csts = chilife.read_library(rotlib_path[0])
            unis = []

            for lib in (libA, libB):
                names, types = lib['atom_names'], libA['atom_types']
                residxs =  np.zeros(len(names), dtype=int)
                resnames, resids = np.array([self.res]), np.array([1])
                segidx = np.array([0])

                uni = chilife.make_mda_uni(names, types, resnames, residxs, resids, segidx)
                unis.append(uni)


            for p in rotlib_path:
                tlibA, tlibB, tcsts = chilife.read_library(p)
                for lib, tlib, cct, uni in zip((libA, libB), (tlibA, tlibB), (cctA, cctB), unis):

                    # Libraries must have the same atom order
                    if not np.all(np.isin(tlib['atom_names'], lib['atom_names'])) and \
                           np.all(tlib['dihedral_atoms'] == lib['dihedral_atoms']):

                        raise ValueError(f'Rotlibs {rotlib_path[0].stem} and {p.stem} are not compatable. You may'
                                         f'need to rename one of them.')


                    # Map coordinates
                    ixmap = [np.argwhere(tlib['atom_names'] == aname).flat[0] for aname in lib['atom_names']]
                    cct.setdefault('coords', []).append(tlib['coords'][:, ixmap])

                    # Create new internal coords if they are defined differently
                    lib_ic, tlib_ic =lib['internal_coords'][0], tlib['internal_coords'][0]
                    if np.any(lib_ic.atom_names != tlib_ic.atom_names):
                        uni.load_new(cct['coords'][-1])
                        ics = [chilife.get_internal_coords(uni, lib['dihedral_atoms'], lib_ic.bonded_pairs)
                               for ts in uni.trajectory]
                        tlib['internal_coords'] = ics

                    for field in concatable:
                        cct.setdefault(field, []).append(tlib[field])

                    for field in concatable + ['coords']:
                        lib[field] = np.concatenate(cct[field])

            libA['weights'] = libB['weights'] = np.ones(len(libA['coords'])) / len(libA['coords'])

        self.csts = csts
        self.libA, self.libB = libA, libB
        self.kwargs["eval_clash"] = False


    def create_ensembles(self):

        self.RL1 = chilife.RotamerEnsemble(self.res,
                                           self.site1,
                                           self.protein,
                                           self.chain,
                                           self.libA,
                                           **self.kwargs)

        self.RL2 = chilife.RotamerEnsemble(self.res,
                                           self.site2,
                                           self.protein,
                                           self.chain,
                                           self.libB,
                                           **self.kwargs)


    def save_pdb(self, name=None):
        if name is None:
            name = self.name + ".pdb"
        if not name.endswith(".pdb"):
            name += ".pdb"

        chilife.save(name, self.RL1, self.RL2)

    def minimize(self):
        dummy = self.copy()

        scores =  [self._min_one(i, ic1, ic2, dummy) for i, (ic1, ic2) in
                 enumerate(zip(self.RL1.internal_coords, self.RL2.internal_coords))]

        scores = np.asarray(scores)
        SSEs = np.linalg.norm(self.RL1.coords[:, self.cst_idx1] - self.RL2.coords[:, self.cst_idx2], axis=2).sum(axis=1)
        MSD = SSEs/len(self.csts)
        MSDmin= MSD.min()

        if MSDmin > 0.1:
            warnings.warn(f'The minimum MSD of the cap is {MSD.min()}, this may result in distorted spin label. '
                          f'Check that the structures make sense.')

        if MSDmin > 0.25:

           raise RuntimeError(f'chiLife was unable to connect residues {self.site1} and {self.site2} with {self.res}. '
                              f'Please double check that this is the intended labeling site1. It is likely that these '
                              f'sites are too far apart.')

        self.RL1.backbone_to_site()
        self.RL2.backbone_to_site()

        scores -= scores.min()

        self.weights *= np.exp(-scores / (chilife.GAS_CONST * self.temp) / np.exp(-scores).sum())
        self.weights /= self.weights.sum()

    def _objective(self, dihedrals, ic1, ic2, dummy):

        ic1.set_dihedral(dihedrals[: len(self.RL1.dihedral_atoms)], 1, self.RL1.dihedral_atoms)
        coords1 = ic1.to_cartesian()[self.RL1.ic_mask]

        ic2.set_dihedral(dihedrals[-len(self.RL2.dihedral_atoms):], 1, self.RL2.dihedral_atoms)
        coords2 = ic2.to_cartesian()[self.RL2.ic_mask]

        diff = np.linalg.norm(coords1[self.cst_idx1] - coords2[self.cst_idx2], axis=1)
        ovlp = (coords1[self.cst_idx1] + coords2[self.cst_idx2]) / 2
        coords = np.concatenate([coords1[self.rl1mask], coords2[self.rl2mask], ovlp], axis=0)
        dummy._coords = np.atleast_3d(coords)
        r = np.linalg.norm(coords[self.aidx] - coords[self.bidx], axis=1)

        # Faster to compute lj here
        lj = self.irmin_ij / r
        lj = lj * lj *lj
        lj = lj * lj

        # attractive forces are needed, otherwise this term will perpetually push atoms apart
        internal_energy = self.ieps_ij * (lj * lj - 2 * lj)
        external_energy = self.energy_func(dummy)
        score = (diff @ diff) * self.restraint_weight / len(diff) + internal_energy.sum() + external_energy.sum()

        return score

    def _min_one(self, i, ic1, ic2, dummy):

        d0 = np.concatenate([ic1.get_dihedral(1, self.RL1.dihedral_atoms),
                             ic2.get_dihedral(1, self.RL2.dihedral_atoms)])

        lb = [-np.pi] * len(d0)  # d0 - np.deg2rad(40)
        ub = [np.pi] * len(d0)  # d0 + np.deg2rad(40)  #
        bounds = np.c_[lb, ub]
        xopt = opt.minimize(self._objective, x0=d0, args=(ic1, ic2), bounds=bounds, method=self.min_method)
        self.RL1._coords[i] = ic1.coords[self.RL1.H_mask]
        self.RL2._coords[i] = ic2.coords[self.RL2.H_mask]
        tors = d0 - xopt.x
        tors = np.arctan2(np.sin(tors), np.cos(tors))
        tors = np.sqrt(tors @ tors)

        return xopt.fun + tors

    @property
    def weights(self):
        return self.RL1.weights

    @weights.setter
    def weights(self, value):
        self.RL1.weights = value
        self.RL2.weights = value

    @property
    def coords(self):
        ovlp = (self.RL1.coords[:, self.cst_idx1] + self.RL2.coords[:, self.cst_idx2]) / 2
        return np.concatenate([self.RL1._coords[:, self.rl1mask], self.RL2._coords[:, self.rl2mask], ovlp], axis=1)

    @coords.setter
    def coords(self, value):
        if value.shape[1] != len(self.atom_names):
            raise ValueError(
                f"The provided coordinates do not match the number of atoms of this ensemble ({self.res})"
            )

        self.RL1._coords[:, self.rl1mask] = value[:, :len(self.rl1mask)]
        self.RL2._coords[:, self.rl2mask] = value[:, len(self.rl1mask):len(self.rl1mask) + len(self.rl2mask)]
        self.RL1._coords[:, self.cst_idx1] = value[:, len(self.rl1mask) + len(self.rl2mask):]
        self.RL2._coords[:, self.cst_idx2] = value[:, len(self.rl1mask) + len(self.rl2mask):]


    @property
    def _lib_coords(self):
        ovlp = (self.RL1._lib_coords[:, self.cst_idx1] + self.RL2._lib_coords[:, self.cst_idx2]) / 2
        return np.concatenate([self.RL1._lib_coords[:, self.rl1mask],
                               self.RL2._lib_coords[:, self.rl2mask], ovlp], axis=1)

    @property
    def atom_names(self):
        return np.concatenate((self.RL1.atom_names[self.rl1mask],
                               self.RL2.atom_names[self.rl2mask],
                               self.RL1.atom_names[self.cst_idx1]))

    @property
    def atom_types(self):
        return np.concatenate((self.RL1.atom_types[self.rl1mask],
                               self.RL2.atom_types[self.rl2mask],
                               self.RL1.atom_types[self.cst_idx2]))

    @property
    def centroid(self):
        return self.coords.mean(axis=(0, 1))

    @property
    def clash_ori(self):

        if isinstance(self._clash_ori_inp, (np.ndarray, list)):
            if len(self._clash_ori_inp) == 3:
                return self._clash_ori_inp

        elif isinstance(self._clash_ori_inp, str):
            if self._clash_ori_inp in ["cen", "centroid"]:
                return self.centroid

            elif (ori_name := self._clash_ori_inp.upper()) in self.atom_names:
                return np.squeeze(self.coords[0][ori_name == self.atom_names])

        else:
            raise ValueError(
                f"Unrecognized clash_ori option {self._clash_ori_inp}. Please specify a 3D vector, an "
                f"atom name or `centroid`"
            )

        return self._clash_ori

    @clash_ori.setter
    def clash_ori(self, inp):
        self._clash_ori_inp = inp

    @property
    def side_chain_idx(self):
        side_chain_idx = np.argwhere(
            np.isin(self.atom_names, dRotamerEnsemble.backbone_atoms, invert=True)
        ).flatten()
        return side_chain_idx

    @property
    def rmin2(self):
        return np.concatenate([self.RL1.rmin2, self.RL2.rmin2])

    @property
    def eps(self):
        return np.concatenate([self.RL1.eps, self.RL2.eps])

    @property
    def bonds(self):
        """ """
        if not hasattr(self, "_bonds"):
            bonds = []

            for bond in self.SL1.bonds:
                bndin = np.isin(bond, self.rl1mask)
                if np.all(bndin):
                    bonds.append(bond)
                elif np.any(bndin):
                    bonds.append([bond[0], np.argwhere(self.atom_names == self.RL1.atom_names[bond[1]]).flat[0]])
                else:
                    bonds.append([np.argwhere(self.atom_names == self.RL1.atom_names[bond[0]]).flat[0],
                                 np.argwhere(self.atom_names == self.RL1.atom_names[bond[1]]).flat[0]])

            for bond in self.SL2.bonds:
                bndin = np.isin(bond, self.rl2mask)
                if np.all(bndin):
                    bonds.append([b + len(self.rl1mask) for b in bond])
                elif not bndin[1]:
                    bonds.append([bond[0] + len(self.rl1mask),
                                  np.argwhere(self.atom_names == self.RL2.atom_names[bond[1]]).flat[0]])

            self._bonds = np.array(sorted(set(map(tuple, bonds))), dtype=int)

        return self._bonds

    @bonds.setter
    def bonds(self, inp):
        """

        Parameters
        ----------
        inp :


        Returns
        -------

        """
        self._bonds = set(tuple(i) for i in inp)
        idxs = np.arange(len(self.atom_names))
        all_pairs = set(combinations(idxs, 2))
        self._non_bonded = all_pairs - self._bonds

    @property
    def non_bonded(self):
        """ """
        if not hasattr(self, "_non_bonded"):
            pairs = dict(nx.all_pairs_shortest_path(self._graph, self._exclude_nb_interactions - 1))
            pairs = {(a, b) for a in pairs for b in pairs[a] if a < b}
            all_pairs = set(combinations(range(len(self.atom_names)), 2))
            self._non_bonded = all_pairs - pairs

        return sorted(list(self._non_bonded))

    @non_bonded.setter
    def non_bonded(self, inp):
        """

        Parameters
        ----------
        inp :


        Returns
        -------

        """
        self._non_bonded = set(tuple(i) for i in inp)
        idxs = np.arange(len(self.atom_names))
        all_pairs = set(combinations(idxs, 2))
        self._bonds = all_pairs - self._non_bonded

    def trim_rotamers(self):
        self.RL1.trim_rotamers()
        self.RL2.trim_rotamers()

    def evaluate(self):
        """Place rotamer ensemble on protein site1 and recalculate rotamer weights."""
        # Calculate external energies
        energies = self.energy_func(self)

        # Calculate total weights (combining internal and external)
        self.weights, self.partition = chilife.reweight_rotamers(energies, self.temp, self.weights)
        logging.info(f"Relative partition function: {self.partition:.3}")

        # Remove low-weight rotamers from ensemble
        self.trim_rotamers()

    def __len__(self):
        return len(self.weights)

    # def copy(self):

def get_possible_rotlibs(rotlib: str, all: bool = False) -> Union[Path, None]:
    """

    """
    cwd = Path.cwd()

    # Assemble a list of possible rotlib paths starting in the current directory
    possible_rotlibs = [Path(rotlib),
                        cwd / rotlib,
                        cwd / (rotlib + '.zip'),
                        cwd / (rotlib + '_drotlib.zip')]

    possible_rotlibs += list(cwd.glob(f'{rotlib}*_drotlib.zip'))
    # Then in the user defined rotamer library directory
    for pth in chilife.USER_RL_DIR:
        possible_rotlibs += list(pth.glob(f'{rotlib}*_drotlib.zip'))

    # Then in the chilife directory
    possible_rotlibs += list((chilife.RL_DIR / 'user_rotlibs').glob(f'*{rotlib}*'))

    if all:
        rotlib = []
    for possible_file in possible_rotlibs:
        if possible_file.exists() and all:
                rotlib.append(possible_file)
        elif possible_file.exists():
            rotlib = possible_file
            break
    else:
        if not isinstance(rotlib, list) or rotlib == []:
            rotlib = None

    # rotlib lists need to be sorted to prevent position mismatches for results with tests.
    if isinstance(rotlib, list):
        rotlib = sorted(rotlib)

    return rotlib
