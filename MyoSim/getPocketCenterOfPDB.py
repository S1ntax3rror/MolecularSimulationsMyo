import os
import numpy as np
import MDAnalysis

pocket_XE5 = ['resname VAL and resid 68 and name CA', 'resname LEU and resid 29 and name CA',
              'resname GLY and resid 65 and name CA', 'resname ILE and resid 107 and name CA',
              'resname HSD and resid 64 and name CA', 'resname LEU and resid 69 and name CA',
              'resname THR and resid 67 and name CA', 'resname ILE and resid 28 and name CA',
              'resname LEU and resid 104 and name CA', 'resname PHE and resid 43 and name CA',
              'resname GLY and resid 25 and name CA', 'resname LEU and resid 72 and name CA',
              'resname LEU and resid 32 and name CA', 'resname THR and resid 39 and name CA',
              'resname ALA and resid 71 and name CA', 'resname ILE and resid 99 and name CA',
              'resname HSD and resid 93 and name CA', 'resname LEU and resid 61 and name CA',
              'resname LYS and resid 42 and name CA', 'resname LEU and resid 89 and name CA',
              'resname TYR and resid 103 and name CA']

path_pdb = "/home/kaeserj/PycharmProjects/CurveFitMorse/MyoSim/frozen_myo_strucutre_data/1K_withCO/dyna1.pdb"
path_dcd = "../Data/dcdForH2BondwiggleV3/dyna1.dcd"
path_psf = "../Data/dcdForH2BondwiggleV3/step3_pbcsetup.psf"


def calc_pocket_pos(pocket_indices, frame):
    positions_pocket = frame._pos[pocket_indices]
    mean_position_pocket = np.mean(positions_pocket, axis=0)
    return mean_position_pocket


def read_pocket_list(res_list, univ):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(univ.select_atoms(res).ix[0])
    return np.array(pocket_indices)


def format_entry(index, group, atom_type, atom_name, x, y, z, mol_type, res_id, charge):
    return f"{index:>10}{group:>10}  {atom_type:<2}        {atom_name:<7}{x:>20.10f} {y:>20.10f} {z:>20.10f}  {mol_type:<2} {res_id:>10}  {charge:>20.10f}"


def add_co_at_center_of_pocket_5(path_pdb, write_to_path, debug=False):
    atom_mda = MDAnalysis.Universe(path_pdb)
    xe5_formatted = read_pocket_list(pocket_XE5, atom_mda)

    pocket_center_pos = calc_pocket_pos(xe5_formatted, atom_mda.trajectory[0])
    if debug:
        print(pocket_center_pos)

    header = f"{2:>10}  EXT"
    entries = [
        format_entry(1, 1, "CO", "O", pocket_center_pos[0], pocket_center_pos[1], pocket_center_pos[2], "CO", 156, 0.0),
        format_entry(2, 1, "CO", "C", pocket_center_pos[0] + 1.12, pocket_center_pos[1], pocket_center_pos[2], "CO",
                     156, 0.0),
    ]

    output = "\n".join([header] + entries)
    if debug:
        print(output)

    with open(write_to_path, "w") as file:
        file.write(output)


def get_coord_of_xe5(path_pdb):
    atom_mda = MDAnalysis.Universe(path_pdb)
    xe5_formatted = read_pocket_list(pocket_XE5, atom_mda)

    pocket_center_pos = calc_pocket_pos(xe5_formatted, atom_mda.trajectory[0])

    return pocket_center_pos


def get_coord_of_CO(path_pdb):
    atom_mda = MDAnalysis.Universe(path_pdb)
    CO_atoms = ["resname CO and type C", "resname CO and type O"]
    pocket_indices = []
    print(set(atom_mda._topology.resnames.values))
    for res in CO_atoms:
        pocket_indices.append(atom_mda.select_atoms(res).ix[0])

    H2_arr = np.array(pocket_indices)

    return atom_mda.trajectory[0]._pos[H2_arr]


def exec_main(write_to_path, debug=False):
    atom_mda = MDAnalysis.Universe(path_pdb)
    xe5_formatted = read_pocket_list(pocket_XE5, atom_mda)

    pocket_center_pos = calc_pocket_pos(xe5_formatted, atom_mda.trajectory[0])
    if debug:
        print(pocket_center_pos)

    header = f"{2:>10}  EXT"
    entries = [
        format_entry(1, 1, "CO", "O", pocket_center_pos[0], pocket_center_pos[1], pocket_center_pos[2], "CO", 156, 0.0),
        format_entry(2, 1, "CO", "C", pocket_center_pos[0] + 1.12, pocket_center_pos[1], pocket_center_pos[2], "CO",156, 0.0),
    ]

    output = "\n".join([header] + entries)
    if debug:
        print(output)

    with open(write_to_path, "w") as file:
        file.write(output)


if __name__ == "__main__":
    # exec_main("../Data/FrozenMyoglobinPDBs/co.crd")
    get_coord_of_CO(path_pdb)