import copy,json,unittest
from pathlib import Path
import numpy as np
from cast_alpha import reconstruct,measures,free_fraction,compare
from scipy.spatial import cKDTree

def one_tetra():
    return dict(ok=True,index_space='input_atoms',alpha_all_dimensions_exported=True,
                cells_idx=[[0,1,2,3]],neighbors=[[-1]*4],exterior=1,
                successors=[[],[]],outside_ancestors=[],tetra_in_alpha=[False],
                face_in_alpha=[[False]*4],alpha_vertex_rows=[[i,False] for i in range(4)],
                alpha_edge_rows=[[i,j,False] for i in range(4) for j in range(i+1,4)],
                alpha_squared=0.,pocket_tetra_components=[[0]])
class Tests(unittest.TestCase):
    def test_single_tetra_algebra(self):
        r=reconstruct(one_tetra());self.assertEqual(len(r['regions']),1)
        self.assertEqual(len(r['regions'][0]['closure_boundary']['triangles']),4)
        self.assertEqual(r['regions'][0]['dual_set']['tetrahedra'],[[0,1,2,3]])
    def test_outside_reachability(self):
        g=one_tetra();g['successors'][0]=[1];g['outside_ancestors']=[0];g['pocket_tetra_components']=[]
        self.assertEqual(reconstruct(g)['regions'],[])
    def test_wrong_reachability_rejected(self):
        g=one_tetra();g['outside_ancestors']=[0]
        with self.assertRaises(ValueError):reconstruct(g)
    def test_alpha_nonclosure_rejected(self):
        g=one_tetra();g['alpha_edge_rows'][0][2]=True
        with self.assertRaises(ValueError):reconstruct(g)
    def test_degenerate_predicate_rejected(self):
        g=one_tetra();g['zero_face_predicates']=[0]
        with self.assertRaises(ValueError):reconstruct(g)
    def test_bad_neighbors_rejected(self):
        g=one_tetra();g['neighbors'][0][0]=0
        with self.assertRaises(ValueError):reconstruct(g)
    def test_support_volume(self):
        g=one_tetra();r=reconstruct(g)
        m=dict(coordinates=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],atomic_radii=[.01]*4,alpha_squared=0.)
        measures(g,m,r,0,1729)
        self.assertAlmostEqual(r['regions'][0]['tetrahedral_support_volume_A3'],1/6)
    def test_sphere_membership(self):
        c=np.array([[0.,0.,0.]]);r=np.array([1.]);p=np.array([[0.,0.,0.],[2.,0.,0.]])
        self.assertEqual(free_fraction(p,c,r,cKDTree(c)).tolist(),[False,True])
    def test_seed_reproducibility(self):
        g=one_tetra();m=dict(coordinates=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],atomic_radii=[.2]*4,alpha_squared=0.)
        a=reconstruct(g);b=reconstruct(g);measures(g,m,a,32,1729);measures(g,m,b,32,1729)
        self.assertEqual(a,b)
if __name__=='__main__':unittest.main()
