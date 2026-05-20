from ..utils.pocket_utils import compute_mouth_parameters


def compute_mouth_features_for_all_pockets(ranked_pockets, simplices, protein_coords, ligand_coords=None):
    return compute_mouth_parameters(
        ranked_pockets=ranked_pockets,
        simplices=simplices,
        protein_coords=protein_coords,
        ligand_coords=ligand_coords,
    )