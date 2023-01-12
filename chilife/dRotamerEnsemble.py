from typing import Union
from itertools import combinations
import logging
import warnings
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
import scipy.optimize as opt

import chilife


class dRotamerEnsemble:

    def __init__(self, res, sites, protein=None, chain=None, rotlib=None, **kwargs):
        """ """
        self.res = res
        self.site, self.site2 = sorted(sites)
        self.increment = self.site2 - self.site
        self.kwargs = kwargs

        self.protein = protein
        self.chain = chain if chain is not None else self.guess_chain()
        self.protein_tree = self.kwargs.setdefault("protein_tree", None)

        self.forgive = kwargs.setdefault("forgive", 1.0)
        self.clash_radius = kwargs.setdefault("clash_radius", 14.0)
        self._clash_ori_inp = kwargs.setdefault("clash_ori", "cen")
        self.restraint_weight = kwargs.pop("restraint_weight") if "restraint_weight" in kwargs else 200  # kcal/mol/A^2
        self.alignment_method = kwargs.setdefault("alignment_method", "bisect".lower())
        self.dihedral_sigmas = kwargs.setdefault("dihedral_sigmas", 25)
        self.minimize = kwargs.pop("minimize", True)
        self.eval_clash = kwargs.pop("eval_clash", True)
        self.energy_func = kwargs.setdefault("energy_func", chilife.get_lj_rep)
        self.temp = kwargs.setdefault("temp", 298)
        self.get_lib(rotlib)
        self.create_ensembles()

        self.name = self.res
        if self.site is not None:
            self.name = f"{self.RL1.nataa}{self.site}{self.RL1.nataa}{self.site2}{self.res}"
        if self.chain is not None:
            self.name += f"_{self.chain}"

        self.selstr = (
            f"resid {self.site} {self.site2} and segid {self.chain} and not altloc B"
        )

        self.protein_setup()
        self.sub_labels = (self.RL1, self.RL2)

    def protein_setup(self):
        self.protein = self.protein.select_atoms("not (byres name OH2 or resname HOH)")
        self.clash_ignore_idx = self.protein.select_atoms(
            f"resid {self.site} {self.site2} and segid {self.chain}"
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

        if self.minimize:
            self._minimize()

        if self.eval_clash:
            self.evaluate()

    def guess_chain(self):
        if self.protein is None:
            chain = "A"
        elif len(set(self.protein.segments.segids)) == 1:
            chain = self.protein.segments.segids[0]
        elif np.isin(self.protein.residues.resnums, self.site).sum() == 0:
            raise ValueError(
                f"Residue {self.site} is not present on the provided protein"
            )
        elif np.isin(self.protein.residues.resnums, self.site).sum() == 1:
            chain = self.protein.select_atoms(f"resid {self.site}").segids[0]
        else:
            raise ValueError(
                f"Residue {self.site} is present on more than one chain. Please specify the desired chain"
            )
        return chain

    def get_lib(self, rotlib):
        if self.protein is not None:
            # get site backbone information from protein structure
            sel_txt = f"resnum {self.site} and segid {self.chain}"


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
                              f'this site pair and but they may not be accurate! Because there is no information '
                              f'about the relative weighting of different rotamer libraries all weights will be '
                              f'set to 1/len(rotlib)')

        if isinstance(rotlib_path, Path):
            libA, libB, csts = chilife.read_library(rotlib_path)
        elif isinstance(rotlib_path, list):
            concatable = ('coords', 'dihedrals', 'internal_coords')

            cctA, cctB, ccsts = {}, {}, {}
            libA, libB, csts = chilife.read_library(rotlib_path[0])
            for p in rotlib_path:
                tlibA, tlibB, tcsts = chilife.read_library(p)

                # Libraries must have the same atom order
                if np.any(tlibA['atom_names'] != libA['atom_names']) or \
                    np.any(tlibB['atom_names'] != libB['atom_names']):

                    raise ValueError(f'Rotlibs {rotlib_path[0].stem} and {p.stem} are not compatable. You may'
                                     f'need to rename one of them.')
                for field in concatable:
                    cctA.setdefault(field, []).append(tlibA[field])
                    cctB.setdefault(field, []).append(tlibB[field])

                ccsts.setdefault('cst_pairs', []).append(tcsts['cst_pairs'])
                ccsts.setdefault('cst_distances', []).append(tcsts['cst_distances'])

            for field in concatable:
                libA[field] = np.concatenate(cctA[field])
                libB[field] = np.concatenate(cctB[field])

            libA['weights'] = libB['weights'] = np.ones(len(libA['coords'])) / len(libA['coords'])
            csts = {key: np.concatenate(val) for key, val in ccsts.items()}

        self.cst_idxs, self.csts = tuple(csts.values())
        self.libA, self.libB = libA, libB
        self.kwargs["eval_clash"] = False


    def create_ensembles(self):

        self.RL1 = chilife.RotamerEnsemble(self.res,
                                           self.site,
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

    def _minimize(self):
        def objective(dihedrals, ic1, ic2, opt):
            coords1 = ic1.set_dihedral(
                dihedrals[: len(self.RL1.dihedral_atoms)], 1, self.RL1.dihedral_atoms
            ).to_cartesian()
            coords2 = ic2.set_dihedral(
                dihedrals[-len(self.RL2.dihedral_atoms):], 1, self.RL2.dihedral_atoms
            ).to_cartesian()

            distances = np.linalg.norm(
                coords1[self.cst_idxs[:, 0]] - coords2[self.cst_idxs[:, 1]], axis=1
            )
            diff = distances - opt
            return diff @ diff

        scores = np.empty_like(self.weights)
        for i, (ic1, ic2) in enumerate(
            zip(self.RL1.internal_coords, self.RL2.internal_coords)
        ):
            d0 = np.concatenate(
                [
                    ic1.get_dihedral(1, self.RL1.dihedral_atoms),
                    ic2.get_dihedral(1, self.RL2.dihedral_atoms),
                ]
            )
            lb = [-np.pi] * len(d0)  # d0 - np.deg2rad(40)  #
            ub = [np.pi] * len(d0)  # d0 + np.deg2rad(40) #
            bounds = np.c_[lb, ub]
            xopt = opt.minimize(
                objective, x0=d0, args=(ic1, ic2, self.csts[i]), bounds=bounds
            )
            self.RL1._coords[i] = ic1.coords[self.RL1.H_mask]
            self.RL2._coords[i] = ic2.coords[self.RL2.H_mask]
            scores[i] = xopt.fun

        self.RL1.backbone_to_site()
        self.RL2.backbone_to_site()

        scores /= len(self.cst_idxs)
        scores -= scores.min()

        self.weights *= np.exp(-scores * self.restraint_weight / (chilife.GAS_CONST * self.temp) / np.exp(-scores).sum())
        self.weights /= self.weights.sum()

    @property
    def weights(self):
        return self.RL1.weights

    @weights.setter
    def weights(self, value):
        self.RL1.weights = value
        self.RL2.weights = value

    @property
    def coords(self):
        return np.concatenate([self.RL1._coords, self.RL2._coords], axis=1)

    @coords.setter
    def coords(self, value):
        if value.shape[1] != self.RL1._coords.shape[1] + self.RL2._coords.shape[1]:
            raise ValueError(
                f"The provided coordinates do not match the number of atoms of this ensemble ({self.res})"
            )

        self.RL1._coords = value[:, : self.RL1._coords.shape[1]]
        self.RL2._coords = value[:, -self.RL2._coords.shape[1]:]

    @property
    def atom_names(self):
        return np.concatenate((self.RL1.atom_names, self.RL2.atom_names))

    @property
    def atom_types(self):
        return np.concatenate((self.RL1.atom_types, self.RL2.atom_types))

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
        return np.concatenate(
            [
                self.RL1.side_chain_idx,
                self.RL2.side_chain_idx + len(self.RL1.atom_names),
            ]
        )

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
            self._bonds = chilife.guess_bonds(self.coords[0], self.atom_types)
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
            idxs = np.arange(len(self.atom_names))
            all_pairs = set(combinations(idxs, 2))
            self._non_bonded = all_pairs - set(tuple(bond) for bond in self.bonds)

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
        """Place rotamer ensemble on protein site and recalculate rotamer weights."""

        # Calculate external energies
        energies = self.energy_func(self)

        # Calculate total weights (combining internal and external)
        self.weights, self.partition = chilife.reweight_rotamers(energies, self.temp, self.weights)
        logging.info(f"Relative partition function: {self.partition:.3}")

        # Remove low-weight rotamers from ensemble
        self.trim_rotamers()


def get_possible_rotlibs(rotlib: str, all: bool = False) -> Union[Path, None]:
    """

    """

    cwd = Path.cwd()

    # Assemble a list of possible rotlib paths starting in the current directory
    possible_rotlibs = [Path(rotlib),
                        cwd / rotlib,
                        cwd / (rotlib + '.zip'),
                        cwd / (rotlib + '_drotlib.zip')]

    # Then in the user defined rotamer library directory
    for pth in chilife.USER_RL_DIR:
        possible_rotlibs += list(pth.glob(f'*{rotlib}*_drotlib.zip'))

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

    return rotlib
