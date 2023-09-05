import os
import hashlib
import numpy as np
import pytest
import MDAnalysis as mda
import chilife as xl

protein = mda.Universe("test_data/1ubq.pdb", in_memory=True)
gb1 = mda.Universe("test_data/4wh4.pdb", in_memory=True).select_atoms("protein and segid A")
SL2 = xl.dSpinLabel("DHC", [28, 28+4], gb1)


def test_add_dlabel():
    Energies = np.loadtxt("test_data/DHC.energies")[:, 1]
    P = np.exp(-Energies / (xl.GAS_CONST * 298))
    P /= P.sum()
    xl.create_dlibrary(
        "___",
        "test_data/DHC.pdb",
        sites=(2, 6),
        weights=P,
        dihedral_atoms=[
            [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
            [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
        ],
        spin_atoms=["Cu1"],
    )

    # Test that chain operators were reset
    libA, libB, csts = xl.read_drotlib('___ip4_drotlib.zip')
    for lib in (libA, libB):
        for ic in lib['internal_coords']:
            np.testing.assert_almost_equal(ic.chain_operators[1]['ori'], np.zeros(3))
            np.testing.assert_almost_equal(ic.chain_operators[1]['mx'], np.eye(3))


    os.remove('___ip4_drotlib.zip')

def test_distance_distribution():
    r = np.linspace(15, 50, 256)
    SL1 = xl.SpinLabel("R1M", 6, gb1)
    dd = xl.distance_distribution(SL1, SL2, r)
    d = r[np.argmax(dd)]
    p = np.max(dd)
    assert abs(d - 26.529411764705884) < 1e-7
    assert abs(p - 0.3548311024322747) < 1e-7


def test_centroid():
    np.testing.assert_almost_equal(SL2.centroid, [19.69981126, -13.95217125,  10.8655417], decimal=5)


def test_side_chain_idx():
    SL3 = xl.dSpinLabel("DHC", [28, 32], gb1)
    ans = np.array(['CB', 'CG', 'CD2', 'ND1', 'NE2', 'CE1', 'CB', 'CG', 'CD2', 'ND1',
                    'NE2', 'CE1', 'Cu1', 'O3', 'O1', 'O6', 'N5', 'C11', 'C9', 'C14',
                    'C8', 'C7', 'C10', 'O2', 'O4', 'O5'], dtype='<U3')
    np.testing.assert_equal(SL3.atom_names[SL3.side_chain_idx], ans)



def test_coords_setter():
    SL3 = xl.dSpinLabel("DHC", [28, 32], gb1)
    old = SL3.coords.copy()
    ans = old + 5
    SL3.coords += 5
    np.testing.assert_almost_equal(SL3.coords, ans)


def test_coord_set_error():
    SL3 = xl.dSpinLabel("DHC", [28, 32], gb1)
    ar = np.random.rand(5, 20, 3)
    with pytest.raises(ValueError):
        SL3.coords = ar


def test_mutate():
    SL2 = xl.dSpinLabel('DHC', (28, 32), gb1, min_method='Powell')
    gb1_Cu = xl.mutate(gb1, SL2)
    xl.save("mutate_dSL.pdb", gb1_Cu)

    with open("mutate_dSL.pdb", "r") as f:
        test = hashlib.md5(f.read().encode("utf-8")).hexdigest()

    with open("test_data/mutate_dSL.pdb", "r") as f:
        ans = hashlib.md5(f.read().encode("utf-8")).hexdigest()

    os.remove("mutate_dSL.pdb")

    assert test == ans


def test_single_chain_error():
    with pytest.raises(RuntimeError):
        xl.create_dlibrary(libname='___',
                           pdb='test_data/chain_broken_dlabel.pdb',
                           sites=(15, 17),
                           dihedral_atoms=[[['N', 'CA', 'C13', 'C5'],
                                       ['CA', 'C13', 'C5', 'C6']],
                                      [['N', 'CA', 'C12', 'C2'],
                                       ['CA', 'C12', 'C2', 'C3']]],
                           spin_atoms='Cu1')


def test_restraint_weight():
    SL3 = xl.dSpinLabel("DHC", [28, 32], gb1, restraint_weight=5)
    ans = np.array([0.51175099, 0.48824901])

    np.testing.assert_almost_equal(SL3.weights, ans, decimal=3)
    assert np.any(SL2.weights != SL3.weights)


def test_alternate_increment():
    with pytest.warns():
        with pytest.raises(RuntimeError):
            xl.dSpinLabel("DHC", (15, 44), gb1)

    SL2 = xl.dSpinLabel("DHC", (12, 37), gb1)
    ans = np.array([[25.50056398,  1.27502619,  3.3972164 ],
                    [24.45518418,  2.0136011 ,  3.74825233],
                    [24.57566352,  1.97824049,  3.73004511],
                    [24.63698352,  1.94542836,  3.70611475],
                    [24.52183104,  2.00650002,  3.74023405]])
    np.testing.assert_almost_equal(SL2.spin_centers, ans, decimal=4)


def test_min_method():
    SL2 = xl.dSpinLabel("DHC", (28, 32), gb1, min_method='Powell')
    ans = np.array([[ 18.6062596, -14.705718 ,  12.0624657],
                    [ 18.5973143, -14.7182376,  12.0220758]])

    np.testing.assert_almost_equal(SL2.spin_centers, ans)


def test_no_min():
    SL2 = xl.dSpinLabel("DHC", (28, 32), gb1, minimize=False)

    bb_coords = gb1.select_atoms('resid 28 32 and name N CA C O').positions
    bb_idx = np.argwhere(np.isin(SL2.atom_names, SL2.backbone_atoms)).flatten()

    for conf in SL2.coords:
        # Decimal = 1 because bisect alignment does not place exactly by definition
        np.testing.assert_almost_equal(conf[bb_idx], bb_coords, decimal=1)


def test_trim_false():
    SL1 = xl.dSpinLabel('DHC', (28, 32), gb1, trim=False)
    SL3 = xl.dSpinLabel('DHC', (28, 32), gb1, eval_clash=False)

    assert len(SL1) == len(SL3)
    most_probable = np.sort(SL1.weights)[::-1][:len(SL2)]
    most_probable /= most_probable.sum()
    np.testing.assert_almost_equal(SL2.weights, most_probable)
    assert np.any(np.not_equal(SL1.weights, SL3.weights))

def test_dmin_callback():
    vals = []
    ivals = []
    def my_callback(val, i):
        vals.append(val)
        ivals.append(i)

    SL1 = xl.dRotamerEnsemble('DHC', (28, 32), gb1, minimize=False)
    SL1.minimize(callback=my_callback)

    assert len(vals) > 0
    assert len(ivals) > 0


def test_dihedral_atoms():
    ans = [['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'ND1'],
           ['CD2', 'NE2', 'Cu1', 'O1'],
           ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'ND1'],
           ['CD2', 'NE2', 'Cu1', 'O1']]

    assert np.all(SL2.dihedral_atoms == ans)


def test_dihedrals():
    ans = np.array([[-2.97070555,  2.20229457, -2.5926639 , -1.21838379, -2.07361964, -2.74030998],
                    [-2.94511378,  2.18546192, -2.59115287, -1.2207707 , -2.05741481, -2.73392102]])

    np.testing.assert_almost_equal(SL2.dihedrals, ans)