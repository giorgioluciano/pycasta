import unittest
from cast_validation import surface_checks,closure,component_count
class ValidationTests(unittest.TestCase):
    def test_tetra_boundary(self):
        S=closure([(0,1,2),(0,1,3),(0,2,3),(1,2,3)])
        self.assertTrue(surface_checks(S)['local_surface_checks_passed'])
    def test_pinched_edge_removed(self):
        faces=[(0,1,2),(0,1,3),(0,1,4),(0,1,5)]
        S=closure(faces);raw=surface_checks(S)
        self.assertTrue(any(x['edge']==[0,1] and x['triangle_degree']==4 for x in raw['included_edge_failures']))
        # Remove all lower dimensional simplices: raw failure must not leak into open faces.
        open_faces=set(faces)
        self.assertTrue(surface_checks(open_faces)['local_surface_checks_passed'])
        self.assertEqual(component_count(open_faces),4)
    def test_vertex_two_cycles_detected(self):
        faces=[(0,1,2),(0,2,3),(0,1,3),(0,4,5),(0,5,6),(0,4,6)]
        S=set(faces)|{(0,)}
        bad=surface_checks(S)['included_vertex_link_failures']
        self.assertEqual(bad[0]['link_components'],2)
if __name__=='__main__':unittest.main(verbosity=2)
