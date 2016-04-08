from __future__ import print_function
from gpaw.transport.analysor import Transport_Plotter
import numpy as np
import sys
from gpaw import GPAW
from pylab import *
plotter=Transport_Plotter(0)
plotter.read_overhead()

import vtk
from ase.visualize.vtk.atoms import vtkAtoms
usewx = False
try:
    import wx
    usewx = True
except ImportError:
    pass
if usewx:
    from vtk.wx.wxVTKRenderWindow import wxVTKRenderWindow
    app = wx.PySimpleApp()
    frame = wx.Frame(None, -1, 'wxVTKRenderWindow', size=(800,600))
    widget = wxVTKRenderWindow(frame, -1)
    win = widget.GetRenderWindow()
    ren = vtk.vtkRenderer()
    win.AddRenderer(ren)
else:
    ren = vtk.vtkRenderer()
    win = vtk.vtkRenderWindow()
    win.AddRenderer(ren)
    win.SetSize(800,600)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(win)
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
atoms = plotter.atoms.copy()
calc = GPAW()
atoms.set_calculator(calc)
calc.initialize(atoms)
calc.scf.converged = True
if len(sys.argv) <= 2:
    if len(sys.argv[1]) <= 2:
        calc.forces.F_av = plotter.get_info('force', sys.argv[1])
    else:
        tmp = sys.argv[1].split('-')
        ref = int(tmp[1])
        sample = int(tmp[0])
        calc.forces.F_av = plotter.get_info('force', sample) - plotter.get_info('force', ref)
else:
    calc.forces.F_av = plotter.get_info('force', sys.argv[1], sys.argv[2])

print('maximum', np.max(abs(calc.forces.F_av)))
va = vtkAtoms(atoms)
va.add_cell()
va.add_axes()
va.add_forces()
va.add_actors_to_renderer(ren)
if usewx:
    frame.Show()
    app.MainLoop()
else:
    iren.Initialize()
    win.OffScreenRenderingOff()
    win.Render()
    iren.Start()
w2i=vtk.vtkWindowToImageFilter()
w2i.SetInput(win)

