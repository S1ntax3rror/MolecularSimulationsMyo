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

pocket_XE4 = ["resname GLY and resid 25 and name C", "resname ILE and resid 28 and name C",
              "resname LEU and resid 29 and name C", "resname GLY and resid 65 and name C",
              "resname VAL and resid 68 and name C", "resname LEU and resid 72 and name C"]


path_pdb = "../Data/FrozenMyoglobinPDBs/dyna1.pdb"

def calc_pocket_pos(pocket_indices, frame):
    positions_pocket = frame._pos[pocket_indices]
    mean_position_pocket = np.mean(positions_pocket, axis=0)
    return mean_position_pocket

def read_pocket_list(res_list, univ):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(univ.select_atoms(res).ix[0])
    return np.array(pocket_indices)

def get_xe5_xe4_center(path):
    atom_mda = MDAnalysis.Universe(path)
    xe4_formatted = read_pocket_list(pocket_XE4, atom_mda)
    xe5_formatted = read_pocket_list(pocket_XE5, atom_mda)

    pocket_center_pos_XE5 = calc_pocket_pos(xe5_formatted, atom_mda.trajectory[0])
    pocket_center_pos_XE4 = calc_pocket_pos(xe4_formatted, atom_mda.trajectory[0])
    return (pocket_center_pos_XE5+pocket_center_pos_XE4) / 2

def get_xe5_xe4_distance(path):
    atom_mda = MDAnalysis.Universe(path)
    xe4_formatted = read_pocket_list(pocket_XE4, atom_mda)
    xe5_formatted = read_pocket_list(pocket_XE5, atom_mda)

    pocket_center_pos_XE5 = calc_pocket_pos(xe5_formatted, atom_mda.trajectory[0])
    pocket_center_pos_XE4 = calc_pocket_pos(xe4_formatted, atom_mda.trajectory[0])
    return np.linalg.norm(pocket_center_pos_XE4-pocket_center_pos_XE5)



if __name__ == "__main__":
    print(get_xe5_xe4_center(path_pdb))
    print(get_xe5_xe4_distance(path_pdb))

