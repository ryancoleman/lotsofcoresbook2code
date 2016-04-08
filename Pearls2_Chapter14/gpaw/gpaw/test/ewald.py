from __future__ import print_function
import numpy as np
from gpaw.utilities.ewald import Ewald

verbose = True

# References for Madelung constants with homogenous background charges:
# M.P. Marder, "Condensed Matter Physics", J. Wiley, (2000),
# scaled to the Wigner-Seitz radius.
# References for Madelung constants of ionic crystals:
# Y. Sakamoto, The J. of Chem. Phys., vol. 28, (1958), p. 164,
# scaled to the nearest neighbor distance.

if 1: #fcc
    cell = np.array([[0., .5, .5],
                      [.5, .0, .5],
                      [.5, .5, .0]])
    basis = np.array([[0.,0.,0.]])
    charges =  np.array([1.])
    r = np.array([0.0,0.0,0.0])
    charges = np.array([1])
    ewald = Ewald(cell)
    e = ewald.get_electrostatic_potential(r, basis, charges, excludefroml0=0)
    a_this = -e * ewald.get_wigner_seitz_radius(sum(charges))
    a_ref = 1.79175
    if verbose:
        print('Madelung energy, fcc:', a_this, a_this-a_ref)
    assert abs(a_this-a_ref) < 1e-5

if 1: #hcp
    a = 1.
    c = np.sqrt(8./3.)*a
    cell =  np.array([[.5*a, -0.5*np.sqrt(3)*a,0],
                      [.5*a, 0.5*np.sqrt(3)*a,0],
                      [0., 0., c]])
    basis = np.array([[.5*a, .5/np.sqrt(3)*a, 0.],
                      [.5*a, -.5/np.sqrt(3)*a, 0.5*c]])
    r = np.array([.5*a, .5/np.sqrt(3)*a, 0.])
    charges =  np.array([1.,1.])
    ewald = Ewald(cell)
    e = ewald.get_electrostatic_potential(r, basis, charges, excludefroml0=0)
    a_this = -e * ewald.get_wigner_seitz_radius(sum(charges))
    a_ref = 1.79168
    if verbose:
        print('Madelung energy, hcp:', a_this, a_this-a_ref)
    assert abs(a_this-a_ref) < 1e-5
    
if 1: #Simple Cubic
    cell = np.diag([1.,1.,1.])
    basis = np.array([[0.,0.,0.]])
    charges =  np.array([1.])
    r = np.array([0.0,0.0,0.0])
    charges = np.array([1])
    ewald = Ewald(cell)
    e = ewald.get_electrostatic_potential(r, basis, charges, excludefroml0=0)
    a_this = -e*ewald.get_wigner_seitz_radius(1)
    a_ref = 1.76012
    if verbose:
        print('Madelung energy, sc:',  a_this, a_this-a_ref)
    assert abs(a_this-a_ref) < 1e-5

if 1: # NaCl
    cell = np.array([[0., .5, .5],
                      [.5, .0, .5],
                      [.5, .5, .0]])
    basis = np.array([[0.,0.,0.],
                      [.5,.5,.5]])
    r = np.array([0.0,0.0,0.0])
    charges =  np.array([1,-1])
    ewald = Ewald(cell)
    e = ewald.get_electrostatic_potential(r, basis, charges, excludefroml0=0)
    a_this = - 0.5 * e
    a_ref = 1.7475645946331822
    if verbose:
        print('Madelung constant, NaCl (B1):', a_this, a_this-a_ref)
    assert abs(a_this-a_ref) < 1e-13


if 0: # NaCl
    cell =np.diag((1.,1.,1.))
    basis = np.array([[0.,0.,0.],
                      [.5,.5,.0],
                      [.5,.0,.5],
                      [.0,.5,.5],
                      [.5,0.,0.],
                      [.0,.5,.0],
                      [.0,.0,.5],
                      [.5,.5,.5]]
                     )
    r = np.array([0.0,0.0,0.0])
    charges = np.array([1,1,1,1,-1,-1,-1,-1])
    ewald = Ewald(cell)
    e = ewald.get_electrostatic_potential(r, basis, charges, excludefroml0=0)
    a_this = - 0.5 * e
    a_ref = 1.7475645946331822
    if verbose:
        print('Madelung constant, NaCl (B1):', a_this, a_this-a_ref)
    assert abs(a_this-a_ref) < 1e-13

if 1: # CsCl
    cell = np.diag((1.,1.,1.))
    basis = np.array([[0.,0.,0.],
                      [.5, .5, .5]])
    r = np.array([0.0,0.0,0.0])
    charges =  np.array([1,-1])
    ewald = Ewald(cell)
    e = ewald.get_electrostatic_potential(r, basis, charges, excludefroml0=0)
    a_this = - .5 * np.sqrt(3.) * e
    a_ref = 1.76267477307099
    if verbose:
        print('Madelung constant, CsCl (B2):', a_this, a_this-a_ref)
    assert abs(a_this-a_ref) < 1e-13
    
if 1: # Zincblende, cubic ZnS
    cell =  np.array([[.0,.5,.5],
                      [.5,.0,.5],
                      [.5,.5,.0]])
    basis = np.array([[0.,0.,0.],
                      [.25,.25,.25]])
    r = np.array([0.0,0.0,0.0])
    charges =  np.array([1,-1])
    ewald = Ewald(cell)
    e = ewald.get_electrostatic_potential(r, basis, charges, excludefroml0=0)
    a_this = - 0.25 * np.sqrt(3) * e
    a_ref = 1.63805505338879
    if verbose:
        print('Madelung constant, Zincblende (B3):', a_this, a_this-a_ref)
    assert abs(a_this-a_ref) < 1e-13

if 1: # Ideal Wurtzite, ZnS, (B4)
    a = 1.
    c = np.sqrt(8./3.)*a
    u = 3./8.
    cell =  np.array([[.5*a, -0.5*np.sqrt(3)*a,0],
                      [.5*a, 0.5*np.sqrt(3)*a,0],
                      [0., 0., c]])
    basis = np.array([[.5*a, .5/np.sqrt(3)*a, 0.],
                       [.5*a, -.5/np.sqrt(3)*a, 0.5*c],
                       [.5*a, .5/np.sqrt(3)*a, u*c],
                       [.5*a, -.5/np.sqrt(3)*a, (.5+u)*c]])
    r = np.array([.5*a, .5/np.sqrt(3)*a, 0.])
    charges =  np.array([1.,1.,-1,-1])
    ewald = Ewald(cell)
    e = ewald.get_electrostatic_potential(r, basis, charges, excludefroml0=0)
    a_this = - 1. * u * c * e
    a_ref = 1.64132162737
    if verbose:
        print('Madelung constant, Ideal Wurtzite (B4):', a_this, a_this-a_ref)
    assert abs(a_this-a_ref) < 1e-11
