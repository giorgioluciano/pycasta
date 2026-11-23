import unittest
from cast_topology import build_topology

FACES=[[1,2,3],[0,2,3],[0,1,3],[0,1,2]]
def fixture(closed):
    from itertools import combinations
    return dict(ok=True,alpha_all_dimensions_exported=True,
        cells_idx=[[0,1,2,3]],tetra_in_alpha=[False],
        face_in_alpha=[[closed]*4],
        alpha_edge_rows=[[a,b,int(closed)] for a,b in combinations(range(4),2)],
        alpha_vertex_rows=[[i,1] for i in range(4)],
        outside_ancestors=[],pocket_tetra_components=[[0]])

class TopologyTests(unittest.TestCase):
    def test_void(self):
        r=build_topology(fixture(True))['regions'][0]
        self.assertEqual(r['n_mouths'],0)
        self.assertEqual(r['classification'],'void')
        self.assertEqual(r['dual_set']['tetrahedra'],[[0,1,2,3]])
        self.assertTrue(r['boundary_edge_manifold_check'])
    def test_connected_mouth(self):
        r=build_topology(fixture(False))['regions'][0]
        self.assertEqual(r['n_mouths'],1)
        self.assertEqual(len(r['mouths'][0]['triangles']),4)
    def test_missing_alpha_data_rejected(self):
        d=fixture(False);d.pop('alpha_all_dimensions_exported')
        with self.assertRaises(ValueError):build_topology(d)
    def test_invalid_alpha_closure_rejected(self):
        d=fixture(False);d['face_in_alpha'][0][0]=True
        with self.assertRaises(ValueError):build_topology(d)

if __name__=='__main__':unittest.main(verbosity=2)
