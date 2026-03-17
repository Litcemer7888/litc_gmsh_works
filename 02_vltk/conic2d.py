import numpy as np
import vtk


class ScewEllipse:

    def __init__(self, num_points, p, e, a0, freq, r_scale=2.0):
        self.size = num_points
        self.p = p
        self.e = e
        self.a0 = min(a0, self.p / (1 + 1.1 * e))
        self.freq = freq
        self.r_scale = r_scale
        self.phase = 0.0
        
        self.phi = np.linspace(0, 2*np.pi, num_points)
        
        det = 1 + self.e * np.cos(self.phi)
        det = np.where(np.abs(det) < 0.01, 0.01, det)
        self.f_phi = self.p / det

        self.nodes = np.zeros(shape=(3, num_points, 1), dtype=np.double)

        self.update_positions_from_phase()

        self.smth = self.r_phi.reshape(num_points, 1)
    
    def update_positions_from_phase(self):
        self.r_phi = self.f_phi * ( 1 + self.a0 * np.cos(self.freq * self.phi + self.phase) )
        self.r_phi *= self.r_scale
        
        self.nodes[0, :, 0] = self.r_phi * np.cos(self.phi)
        self.nodes[1, :, 0] = self.r_phi * np.sin(self.phi)
        self.nodes[2, :, 0] = self.a0 * np.sin(self.freq * self.phi) * self.r_phi
        self.smth = self.r_phi.reshape(self.size, 1)

    def move(self, tau):
        self.phase += tau * 2.0
        self.update_positions_from_phase()

    def snapshot(self, snap_number):
        structuredGrid = vtk.vtkStructuredGrid()
        points = vtk.vtkPoints()

        smth = vtk.vtkDoubleArray()
        smth.SetName("radius")

        number = self.size
        for i in range(0, number):
            for j in range(0, 1):
                points.InsertNextPoint(self.nodes[0][i,j], self.nodes[1][i,j], self.nodes[2][i,j])
                smth.InsertNextValue(self.smth[i,j])

        structuredGrid.SetDimensions(number, 1, 1)
        structuredGrid.SetPoints(points)
        structuredGrid.GetPointData().AddArray(smth)

        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetInputDataObject(structuredGrid)
        writer.SetFileName("ellipse-step-" + str(snap_number).zfill(3) + ".vts")
        writer.SetDataModeToAscii()
        writer.Write()
        


num_points = 50
p = 1.0
e = 0.5
a0 = 0.2
freq = 5
r_scale = 2.0
T = 120

tau = 0.05


m = ScewEllipse(num_points, p, e, a0, freq, r_scale)
m.snapshot(0)

for i in range(1, T):
    m.move(tau)
    m.snapshot(i)

