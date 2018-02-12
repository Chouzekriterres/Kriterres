import numpy as np
import os
#import sys, h5py, binascii

import pygimli as pg
from pygimli.meshtools import appendTriangleBoundary, merge2Meshes, mesh
from pygimli.mplviewer import drawMesh
from pygimli.viewer import showMesh
# from pygimli.meshtools.mesh import exportHDF5Mesh
# from pygimli.meshtools import convertMesh


# Def des paramètres de géométrie du modéle
xmin, xmax = 0, 40  # en cm
emin, emax = 0, 80  # elevation en cm
dtrou = 30  # prof du trou en cm
etrou = emax - dtrou  # elevation du fond du trou
r = 2  # rayon du trou en cm
dx = 0.1  # On définit le pas de la maille
zaff = 20  # profondeur en cm jusqu'où on souhaite un maillage affiné.
eaff = etrou-zaff


assert dtrou + zaff < emax

xtrou_reg = np.arange(xmin, r + dx, dx, 'float')
etrou_reg = np.arange(etrou, emax + dx, dx, 'float')


efin_reg = np.arange(eaff, etrou + dx, dx, 'float')
h_0 = -95  # charge initiale en cm
h_1 = 10  # hauteur d'eau au fond du trou en cm

etrou_reg[len(etrou_reg)-2:0:-1]

# A présent on crée une zone grâce à un polygone

poly = pg.Mesh(2)  # empty 2d mesh
nStart = poly.createNode(r, emax, 0.0) #On crée un noeud de départ, on travaille en 2D donc le dernier terme vaut 0.0

nA = nStart #On définit le noeud de départ
for e in etrou_reg[len(etrou_reg)-2:0:-1]: #On démarre du haut et on se balade sur l'axe des z en créant un noeud à chaque fois
    nB = poly.createNode(r, e, 0.0)
    poly.createEdge(nA, nB) #On définit un côté entre le noeud précédemment crée et le nouveau
    nA = nB #On remplace le noeud de départ par le noeud nouvellement crée

nB = poly.createNode(r, etrou, 0.0) #On crée un noeud    
poly.createEdge(nA, nB)
nA = nB
    
for x in xtrou_reg[len(xtrou_reg)-2:0:-1]:
    nB = poly.createNode(x, etrou, 0.0)
    poly.createEdge(nA, nB) #On définit un côté entre le noeud précédemment crée et le nouveau
    nA = nB #On remplace le noeud de départ par le noeud nouvellement crée
    
nB = poly.createNode(0, etrou, 0.0) #On crée un noeud    
poly.createEdge(nA, nB)
nA = nB

for e in efin_reg[len(efin_reg)-2:0:-1]: #On démarre de 1 et on se balade sur l'axe des x en créant un noeud à chaque fois
    nB = poly.createNode(0, e, 0.0)
    poly.createEdge(nA, nB) #On définit un côté entre le noeud précédemment crée et le nouveau
    nA = nB #On remplace le noeud de départ par le noeud nouvellement crée
    
nC = poly.createNode(0.0, 0.0, 0.0)
poly.createEdge(nB, nC)
nD = poly.createNode(xmax, 0.0, 0.0)
poly.createEdge(nC, nD)
nE = poly.createNode(xmax, emax, 0.0)
poly.createEdge(nD, nE)
poly.createEdge(nE, nStart) #On ferme le polygone!

tri = pg.TriangleWrapper(poly) #On appelle la fonction triangle
tri.setSwitches('-pzeAfaq31')

mesh = pg.Mesh(2) 
tri.generate(mesh) 

for cell in mesh.cells():
    cell.setMarker(2)

showMesh(mesh)
