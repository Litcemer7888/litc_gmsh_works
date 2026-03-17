import gmsh
import sys
from math import pi

gmsh.initialize()

gmsh.model.add("circle")
lc = 1

gmsh.model.occ.addPoint(0,0,0, lc, 100)

gmsh.model.occ.addPoint(0,4,0,lc,1)
gmsh.model.occ.addPoint(4,0,0,lc,2)
gmsh.model.occ.addPoint(0,-4,0,lc,3)
gmsh.model.occ.addPoint(-4,0,0,lc,4)

gmsh.model.occ.addPoint(0,6,0,lc,5)
gmsh.model.occ.addPoint(6,0,0,lc,6)
gmsh.model.occ.addPoint(0,-6,0,lc,7)
gmsh.model.occ.addPoint(-6,0,0,lc,8)


gmsh.model.occ.addCircleArc(1, 100, 2, 1)
gmsh.model.occ.addCircleArc(2, 100, 3, 2)
gmsh.model.occ.addCircleArc(3, 100, 4, 3)
gmsh.model.occ.addCircleArc(4, 100, 1, 4)

gmsh.model.occ.addCircleArc(5, 100, 6, 5)
gmsh.model.occ.addCircleArc(6, 100, 7, 6)
gmsh.model.occ.addCircleArc(7, 100, 8, 7)
gmsh.model.occ.addCircleArc(8, 100, 5, 8)

gmsh.model.occ.addCurveLoop([1,2,3,4], 1)
gmsh.model.occ.addCurveLoop([5,6,7,8], 2)

gmsh.model.occ.addPlaneSurface([1, -2], 1)

angle = 2 * pi #/ 20 # to see crossection
revolved_volumes = gmsh.model.occ.revolve([(2, 1)], 10, 10, 0, 0, 1, 0, angle)
gmsh.model.occ.synchronize()


gmsh.model.mesh.generate(3)

gmsh.write("torus.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
