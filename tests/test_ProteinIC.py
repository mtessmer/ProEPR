import hashlib, os, pickle
import pytest
import numpy as np
import MDAnalysis as mda
import chilife as xl

pdbids = ["1ubq", "1a2w", '1az5']
ubq = mda.Universe("test_data/1ubq.pdb", in_memory=True)
mbp = mda.Universe("test_data/2klf.pdb", in_memory=True)

ubqIC = xl.newProteinIC.from_protein(ubq)
mbpIC = xl.newProteinIC.from_protein(mbp)

ICs = xl.get_internal_coords(ubq)

gd_kwargs = [
    {"resi": 28, "atom_list": ["C", "N", "CA", "C"]},
    {"resi": 28, "atom_list": [["C", "N", "CA", "C"], ["N", "CA", "C", "N"]]},
]
gd_ans = [-1.1540443794802524, np.array([-1.15404438, -0.66532042])]


def test_from_prot():

    with np.load('test_data/ubq_zmats.npz') as f:
        zmat_ans = f['zmat']
        zmat_idx_ans = f['zmat_idxs']

    np.testing.assert_equal(ubqIC.z_matrix_idxs, zmat_idx_ans)
    np.testing.assert_almost_equal(ubqIC.trajectory.coordinates_array, zmat_ans)

    for i, idxs in enumerate(ubqIC.z_matrix_idxs):
        idxt = [val for val in idxs if val != -2147483648]
        assert i == idxt[-1]


def test_from_prot_traj():

    with np.load('test_data/mbp_zmats.npz') as f:
        zmat_ans = f['zmat']
        zmat_idx_ans = f['zmat_idxs']

    np.testing.assert_equal(mbpIC.z_matrix_idxs, zmat_idx_ans)
    np.testing.assert_almost_equal(mbpIC.trajectory.coordinates_array, zmat_ans)


@pytest.mark.parametrize("pdbid", pdbids)
def test_get_internal_coordinates(pdbid):
    protein = mda.Universe(f"test_data/{pdbid}.pdb", in_memory=True).select_atoms("protein and not altloc B")
    ICs = xl.newProteinIC.from_protein(protein)
    np.testing.assert_almost_equal(ICs.coords, protein.atoms.positions, decimal=4)

    protein = xl.Protein.from_pdb(f"test_data/{pdbid}.pdb").select_atoms("protein and not altloc B")
    ICs = xl.newProteinIC.from_protein(protein)
    np.testing.assert_almost_equal(ICs.coords, protein.atoms.positions, decimal=4)


def test_copy():
    my_copy = ubqIC.copy()

    assert np.all(my_copy.z_matrix_idxs == ubqIC.z_matrix_idxs)
    assert my_copy.z_matrix_idxs is not ubqIC.z_matrix_idxs

    assert np.all(ubqIC.trajectory.coordinates_array == my_copy.trajectory.coordinates_array)
    assert ubqIC.trajectory.coordinates_array is not my_copy.trajectory.coordinates_array


def test_save_pdb():
    protein = mda.Universe("test_data/alphabetical_peptide.pdb").select_atoms("protein")

    uni_ics = xl.newProteinIC.from_protein(protein)
    xl.save("test_data/postwrite_alphabet_peptide.pdb", uni_ics)

    with open("test_data/postwrite_alphabet_peptide.pdb", "r") as f:
        test = hashlib.md5(f.read().encode("utf-8")).hexdigest()

    with open("test_data/alphabetical_peptide.pdb", "r") as f:
        truth = hashlib.md5(f.read().encode("utf-8")).hexdigest()

    os.remove("test_data/postwrite_alphabet_peptide.pdb")
    assert test == truth


def test_has_clashes():
    assert not ICs.has_clashes()
    ubqIC.set_dihedral(np.pi / 2, 35, ["N", "CA", "C", "N"])
    assert ICs.has_clashes()




def test_set_dihedral():
    ICs = xl.get_internal_coords(ubq)
    ICs.set_dihedral(
        [-np.pi / 2, np.pi / 2], 72, [["C", "N", "CA", "C"], ["N", "CA", "C", "N"]]
    )

    ans = mda.Universe("test_data/test_icset.pdb")

    np.testing.assert_almost_equal(ans.atoms.positions, ICs.coords, decimal=3)


def test_set_dihedral2():
    lys = mda.Universe('../chilife/data/rotamer_libraries/residue_pdbs/lys.pdb')

    ICs = xl.get_internal_coords(lys)
    ICs.set_dihedral(np.pi / 2, 1, ["N", "CA", "CB", "CG"])

    print('DSF')


@pytest.mark.parametrize(["inp", "ans"], zip(gd_kwargs, gd_ans))
def test_get_dihedral(inp, ans):
    dihedral = ICs.get_dihedral(**inp)
    np.testing.assert_almost_equal(dihedral, ans)



def test_polypro():
    polypro = mda.Universe("test_data/PPII_Capped.pdb")
    polyproIC = xl.get_internal_coords(polypro)


def test_PRO():
    pro = mda.Universe("../chilife/data/rotamer_libraries/residue_pdbs/pro.pdb")
    pro_ic = xl.get_internal_coords(pro)
    assert ("CD", "CG", "CB", "CA") in pro_ic.ICs[1][1]

def test_PRO2():
    ubq_IC = xl.get_internal_coords(ubq)

    assert ("C", "CA", "N", "C") in ubq_IC.ICs[1][37]
    assert ("C", "CA", "N", "C") in ubq_IC.ICs[1][38]


def test_set_coords():
    R1A = mda.Universe("test_data/R1A.pdb")
    R1A_IC = xl.get_internal_coords(R1A)
    R1A_IC_c = R1A_IC.copy()
    R1A_IC_c.set_dihedral([np.pi/2, -np.pi/2, np.pi/2], 1, [['N', 'CA', 'CB', 'SG' ],
                                                            ['CA', 'CB', 'SG', 'SD'],
                                                            ['CB', 'SG', 'SD', 'CE']])

    coords = R1A_IC_c.coords

    R1A_IC.coords = coords
    np.testing.assert_allclose(R1A_IC.zmats[1], R1A_IC_c.zmats[1], rtol=1e-6)
    np.testing.assert_almost_equal(R1A_IC.coords, R1A_IC_c.coords, decimal=6)


def test_nonbonded():
    with open('test_data/ic_nb.pkl', 'rb') as f:
        nb_ans = pickle.load(f)

    np.testing.assert_equal(ubqIC.nonbonded, nb_ans)


def test_get_zmat_idxs():
    R1A = mda.Universe("test_data/R1A.pdb")
    R1A_IC = xl.get_internal_coords(R1A)
    idxs, stem, idx = R1A_IC.get_zmat_idxs(1, ['CB', 'SG', 'SD', 'CE'], 1)
    assert stem == ('SD', 'SG', 'CB')
    assert idx == -1

    idxs, stem, idx = R1A_IC.get_zmat_idxs(1, ['CE', 'SD', 'SG', 'CB'], 1)
    assert stem == ('SD', 'SG', 'CB')
    assert idx == 0

def test_phi_idxs():

    # Still need to test chain

    idxs = ICs.phi_idxs(range(4, 11,))
    vals = ICs.zmats[1][idxs, -1]
    ans = [ICs.get_dihedral(i, ('C', 'N', 'CA', 'C')) for i in range(4, 11)]

    np.testing.assert_almost_equal(vals, ans)

    idx = ICs.phi_idxs(1)
    assert len(idx) == 0


def test_psi_idxs():

    # Still need to test chain

    idxs = ICs.psi_idxs(range(4, 11,))
    vals = ICs.zmats[1][idxs, -1]
    ans = [ICs.get_dihedral(i, ('N', 'CA', 'C', 'N')) for i in range(4, 11)]

    np.testing.assert_almost_equal(vals, ans)

    idx = ICs.psi_idxs(76)
    assert len(idx) == 0

def test_phi_psi_idxs_multichain():
    prot = xl.fetch('1a2w').select_atoms('protein')
    ICs = xl.get_internal_coords(prot)

    with pytest.raises(ValueError):
        ICs.psi_idxs(10)

    with pytest.raises(ValueError):
        ICs.psi_idxs(10)

    A = ICs.psi_idxs(10, chain=2)
    B = ICs.phi_idxs(10, chain=2)

    assert A[0] == 80
    assert B[0] == 71

def test_chi_idxs():
    idxs = ICs.chi_idxs(range(4, 11, ))
    idxs = np.concatenate(idxs)

    vals = ICs.zmats[1][idxs, -1]
    ans = np.array([-1.0589305, -3.1378433, -3.0641198,  1.338909 , -1.1897472,
                     1.2340594,  1.7094074,  3.0752275, -2.8944212, -3.0027205,
                    -3.040881 ])

    np.testing.assert_almost_equal(vals, ans)

def test_ic_pref_dihe():

    mol = mda.Universe('test_data/test_ic_pref_dihe.pdb', in_memory=True)
    dih = [['N', 'CA', 'CB', 'CB2'],
           ['CA', 'CB', 'CB2', 'NG'],
           ['ND', 'CE3', 'CZ3', 'C31']]

    IC = xl.get_internal_coords(mol, preferred_dihedrals=dih)
    IC.set_dihedral(np.pi / 2, 1, ['ND', 'CE3', 'CZ3', 'C31'])

    sister_dihe_atom_coords = IC.coords[IC.atom_names == 'C36'].flat
    np.testing.assert_almost_equal( sister_dihe_atom_coords, [2.14953486, -3.40647599, -0.93880627])
