import numpy as np
import math
import os
import sys
import gmsh

# gmsh.model.occ.addBox(x, y, z, dx, dy, dz, tag = -1)
# gmsh.model.occ.addBox(x, y, z, dx, dy, dz, tag = -1)
# gmsh.model.occ.addCone(x, y, z, dx, dy, dz, r1, r2, tag = -1, angle = 2*pi)
# gmsh.model.occ.addCylinder(x, y, z, dx, dy, dz, r, tag = -1, angle = 2*pi)
# gmsh.model.occ.fuse(objectDimTags, toolDimTags, tag = -1, removeObject = true, removeTool = true)
# gmsh.model.occ.cut(objectDimTags, toolDimTags, tag = -1, removeObject = true, removeTool = true)
# gmsh.model.occ.fragment(objectDimTags, toolDimTags, tag = -1, removeObject = true, removeTool = true)
# gmsh.model.occ.removeAllDuplicates()
# gmsh.model.occ.synchronize()


def calculer_centre_cercle(point1, point2, point3):
    # Extrait les coordonnées x et y des points
    x1, y1, _ = point1
    x2, y2, _ = point2
    x3, y3, _ = point3

    # Calcul des milieux des segments
    x12 = (x1 + x2) / 2.0
    y12 = (y1 + y2) / 2.0
    x13 = (x1 + x3) / 2.0
    y13 = (y1 + y3) / 2.0

    # Calcul des coefficients des droites perpendiculaires bissectrices
    if y1 != y2:
        a1 = - (x2 - x1) / (y2 - y1)
        b1 = (y2 + y1) / 2.0 - a1 * (x2 + x1) / 2.0
    else:
        a1 = None
        b1 = None

    if y1 != y3:
        a2 = - (x3 - x1) / (y3 - y1)
        b2 = (y3 + y1) / 2.0 - a2 * (x3 + x1) / 2.0
    else:
        a2 = None
        b2 = None

    # Calcul des coordonnées du centre du cercle
    if a1 is not None and a2 is not None:
        x = (b2 - b1) / (a1 - a2)
        y = a1 * x + b1
    elif a1 is None:
        x = x12
        y = a2 * x + b2
    else:
        x = x13
        y = a1 * x + b1

    return x, y


def distance_pts(point1, point2):
    """
    Calcule la distance entre deux points dans l'espace tridimensionnel.
    
    Arguments :
    point1 : tuple, coordonnées du premier point (x, y, z)
    point2 : tuple, coordonnées du deuxième point (x, y, z)
    
    Returns :
    distance : float, distance entre les deux points
    """
    import math
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance

# Finger diameter [m]
Finger_diameter = 2.75e-2
# LDF radius + depth
dy_ldf = 5e-4
r_ldf  = 2e-3
# 
# Extrusion parameters for the cylinder
dx,dy,dz = 2.25e-3,Finger_diameter,2*r_ldf+2*(2*r_ldf)
# 
# Points to compute the position of the center to create the cylinder for the skin
A = [dx,0,0]
B = [1e-2,1.5e-2,0]
C = [dx,Finger_diameter,0]
# Computes the position of the center
xcenter,ycenter = calculer_centre_cercle(A,B,C)
O               = [xcenter,ycenter,0]
radius          = distance_pts(O, A)
# 
# LDF parameters
O_LDF  = [0,Finger_diameter,dz/2]
# 
# Parameters for the "bone" to compute the boolean operation later on
hauteur_centre_bone = 2e-2
Obone  = [0,hauteur_centre_bone,0]
r_bone = 3.5e-3

lc = 5e-4

gmsh.initialize()
# Envisager un rafinement local
# Some kernel parameters of GMSH
gmsh.clear()
gmsh.option.setNumber("General.Terminal",1)
gmsh.option.setNumber("Mesh.Optimize", True)
gmsh.option.setNumber("Mesh.OptimizeNetgen", True)
gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
print("Min mesh size allowed ", gmsh.option.getNumber("Mesh.MeshSizeMin"))
print("Max mesh size allowed ", gmsh.option.getNumber("Mesh.MeshSizeMax"))


# Create the bone
column_v = gmsh.model.occ.addBox(-10*dx, -dy, 0, 10*dx, 3*dy, dz, tag = -1)
bone_v   = gmsh.model.occ.addCylinder(Obone[0], Obone[1], 0, 0, 0, dz, r_bone, tag = -1, angle = 2*np.pi)
gmsh.model.occ.cut([(3, bone_v)],[(3, column_v)])
gmsh.model.occ.synchronize()

# Ass the soft tissue
column_v = gmsh.model.occ.addBox(-10*dx, -dy, 0, 11*dx, 3*dy, dz, tag = -1)
skin_v = gmsh.model.occ.addCylinder(O[0], O[1], O[2], 0, 0, dz, radius, tag = -1, angle = 2*np.pi)
gmsh.model.occ.cut([(3, skin_v)],[(3, column_v)])
gmsh.model.occ.cut([(3, skin_v)],[(3, bone_v)], removeTool=False)
gmsh.model.occ.synchronize()

# Add LDF volume
column_v = gmsh.model.occ.addBox(-10*dx, -dy, 0, 10*dx, 3*dy, dz, tag = -1)
LDF = gmsh.model.occ.addCylinder(O_LDF[0], O_LDF[1], O_LDF[2], 0, -dy_ldf, 0, r_ldf, tag = -1, angle = 2*np.pi)
gmsh.model.occ.cut([(3, LDF)],[(3, column_v)])
gmsh.model.occ.synchronize()

# Add the left box, remove bone and fuse
box_left = gmsh.model.occ.addBox(0, 0, 0, dx, dy, dz, tag = -1)
gmsh.model.occ.fuse([(3, box_left)],[(3, skin_v)])
total = gmsh.model.occ.cut([(3, box_left)],[(3, bone_v)])
gmsh.model.occ.synchronize()
print(total[0])

# Cut at level of the mid bone
cut_box = gmsh.model.occ.addBox(0, 0, 0, 10*dx, hauteur_centre_bone, dz, tag = -1)
gmsh.model.occ.cut(total[0],[(3, cut_box)])
gmsh.model.occ.synchronize()

# Export geometry for check
gmsh.model.occ.synchronize()
gmsh.write('test.geo_unrolled')


gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()
surfaces, volumes = [gmsh.model.getEntities(d) for d in [ 2, 3]]


# print(surfaces, volumes)
# print(volumes)
gmsh.model.occ.synchronize()
# si cylindre
# for volume in volumes[7:]:
#   gmsh.model.occ.remove([volume], recursive = True)
# gmsh.model.occ.remove([(3,6)], recursive = True)
# gmsh.model.occ.remove([(3,12),(3,14),(3,15),(3,17),(3,18),(3,19),(3,20),(3,21),(3,22),(3,23),(3,24),(3,25)], recursive = True)
# si cone
# for volume in volumes[8:]:
#   gmsh.model.occ.remove([volume], recursive = True)
# gmsh.model.occ.remove([(3,5)], recursive = True)
# gmsh.model.occ.remove([(3,7)], recursive = True)

# Ensure not having surprises
for a_volume in volumes:
    center_of_mass = gmsh.model.occ.getCenterOfMass(a_volume[0], a_volume[1])
    if center_of_mass[0]<0:
        gmsh.model.occ.remove([a_volume], recursive = True)
    elif center_of_mass[1]>=C[1]:
        gmsh.model.occ.remove([a_volume], recursive = True)
    elif center_of_mass[1]<0:
        gmsh.model.occ.remove([a_volume], recursive = True)
gmsh.model.occ.synchronize()


surfaces, volumes = [gmsh.model.getEntities(d) for d in [ 2, 3]]

symmetry_bottom, symmetry, load_ldf, top, left, right, up_cyl, bone = [19],[22,16],[17],[23],[25],[18],[24],[20,21]

# 
# tag volumes: si cylindre
tdim =3 
gmsh.model.addPhysicalGroup(tdim, [7], 1)
gmsh.model.setPhysicalName(tdim, 1, 'Tissue')
gmsh.model.addPhysicalGroup(tdim, [6], 2)
gmsh.model.setPhysicalName(tdim, 2, 'LDF')


tdim =2 
gmsh.model.addPhysicalGroup(tdim, symmetry_bottom, 1)
gmsh.model.setPhysicalName(tdim, 1, 'Symmetry Bottom')
gmsh.model.addPhysicalGroup(tdim, symmetry, 2)
gmsh.model.setPhysicalName(tdim, 2, 'Symmetry')
gmsh.model.addPhysicalGroup(tdim, load_ldf, 3)
gmsh.model.setPhysicalName(tdim, 3 , 'Load_LDF')
gmsh.model.addPhysicalGroup(tdim, top, 4)
gmsh.model.setPhysicalName(tdim, 4 , 'top')
gmsh.model.addPhysicalGroup(tdim, left, 5)
gmsh.model.setPhysicalName(tdim, 5 , 'left')
gmsh.model.addPhysicalGroup(tdim, right, 6)
gmsh.model.setPhysicalName(tdim, 6 , 'right')
gmsh.model.addPhysicalGroup(tdim, bone, 7)
gmsh.model.setPhysicalName(tdim, 7 , 'bone')
gmsh.model.addPhysicalGroup(tdim, up_cyl, 8)
gmsh.model.setPhysicalName(tdim, 8, 'up_cyl')

gmsh.model.occ.synchronize()


gmsh.write('Geom_test_2_poro.geo_unrolled')


gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("Mesh_2_poro.msh")

# if 'close' not in sys.argv:
#   gmsh.fltk.run()

gmsh.finalize()