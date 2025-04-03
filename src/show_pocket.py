from pymol.cgo import BEGIN, TRIANGLES, VERTEX, END, COLOR
from pymol import cmd
import numpy as np
from plyfile import PlyData

def draw_extruded_pocket(name, vertices, faces, color=(1.0, 0.5, 0.0)):
    obj = [BEGIN, TRIANGLES, COLOR, *color]
    for face in faces:
        for idx in face:
            v = vertices[idx]
            obj.extend([VERTEX, *v])
    obj.append(END)
    cmd.load_cgo(obj, name)

def load_ply_to_pymol(name, ply_path, color=(1.0, 0.5, 0.0)):
    plydata = PlyData.read(ply_path)
    vertices = np.vstack([plydata['vertex'][axis] for axis in ['x', 'y', 'z']]).T
    faces = np.vstack([f[0] for f in plydata['face'].data])
    draw_extruded_pocket(name, vertices, faces, color=color)