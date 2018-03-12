

import numpy as np
import pygimli as pg

def definition_maillage(xmin, xmax, emin, emax, dtrou, etrou, r, dx, zaff, eaff):
    """Définition du maillage de base à base triangulaire pour l'entrée dans SWMS_2D"""
    
    assert dtrou + zaff < emax

    xtrou_reg = np.arange(xmin, r + dx, dx, 'float')
    etrou_reg = np.arange(etrou, emax + dx, dx, 'float')


    efin_reg = np.arange(eaff, etrou + dx, dx, 'float')
    
    #A présent on crée une zone grâce à un polygone

    poly = pg.Mesh(2)  # empty 2d mesh
    nStart = poly.createNode(0.0, 0.0, 0.0) #On crée un noeud de départ, on travaille en 2D donc le dernier terme vaut 0.0

    nA = nStart #On définit le noeud de départ

    #nB = poly.createNode(xmin, eaff, 0.0)
    #poly.createEdge(nA, nB)
    #nB = nA

    for e in efin_reg[0:len(efin_reg)]: #On démarre de 1 et on se balade sur l'axe des x en créant un noeud à chaque fois
        nB = poly.createNode(0, e, 0.0)
        poly.createEdge(nA, nB) #On définit un côté entre le noeud précédemment crée et le nouveau
        nA = nB #On remplace le noeud de départ par le noeud nouvellement crée
    
    #nB = poly.createNode(0, etrou, 0.0) #On crée un noeud    
    #poly.createEdge(nA, nB)
    #nA = nB

    #Définir les noeuds au fond du trou    
    for x in xtrou_reg[1:len(xtrou_reg)]:
        nB = poly.createNode(x, etrou, 0.0)
        poly.createEdge(nA, nB) #On définit un côté entre le noeud précédemment crée et le nouveau
        nA = nB #On remplace le noeud de départ par le noeud nouvellement crée

    #nB = poly.createNode(r, etrou, 0.0) #On crée un noeud    
    #poly.createEdge(nA, nB)
    #nA = nB    
    
    #Définir les noeuds le long du trou
    #nA = nStart #On définit le noeud de départ
    for f in etrou_reg[1:len(etrou_reg)]: #On démarre du haut et on se balade sur l'axe des z en créant un noeud à chaque fois
        nB = poly.createNode(r, f, 0.0)
        poly.createEdge(nA, nB) #On définit un côté entre le noeud précédemment crée et le nouveau
        nA = nB #On remplace le noeud de départ par le noeud nouvellement crée

    nC=nB
    #nC = poly.createNode(r, emax, 0.0)
    #poly.createEdge(nB, nC)
    nD = poly.createNode(xmax, emax, 0.0)
    poly.createEdge(nC, nD)
    nE = poly.createNode(xmax, emin, 0.0)
    poly.createEdge(nD, nE)
    poly.createEdge(nE, nStart) #On ferme le polygone!
    
    mesh=pg.meshtools.createMesh(poly, quality=33, area=5, smooth=[1,10])

    #tri = pg.TriangleWrapper(poly) #On appelle la fonction triangle
    #tri.setSwitches('-pzeAfaq31')
    # Ici on a :
    # p : planar straight line graph ==> fichier poly
    # z : on démarre le comptage à 0
    # A : assigne un attribut à chaque triangle qui indique à quel segment il appartient et est lié
    # f : algorithme de triangulation (?)
    # a : impose une surface contrainte pour chaque triangle on peut ajouter un nombre si on veut préciser
    # q31 : impose que les triangles générés aient au minimun des angles de 20° O
    
    
    #A présent on génère le maillage hétérogène

    #mesh = pg.Mesh(2) #On appelle le maillage
    #tri.generate(mesh) #On génère les triangles au sein du polygone précédemment crée
    
    #from pygimli.meshtools import polytools as plc
    #import math

    #c1 = plc.createCircle(pos=[r, 50], segments=12, radius=20.0, start=-math.pi/2, end=(math.pi)/2, isClosed=True)
    #c2=pg.meshtools.createMesh(c1)
    #nc1=c1.createNode(0.0, 50, 0.0)
    #nc2=c1.createNode(2.0, 50, 0.0)
    #c1.createEdge(nc1, nc2)
    #showMesh(c1)
    #cercle=pg.meshtools.createMesh([mesh,c1])
    #cercle=merge2Meshes(mesh,c1)
    #showMesh(cercle)
    
    pg_pos = mesh.positions()
    mesh_pos = np.array((np.array(pg.x(pg_pos)), np.array(pg.y(pg_pos)), np.array(pg.z(pg_pos)))).T #On crée une matrice contenant la position des noeuds
    mesh_cells = np.zeros((mesh.cellCount(), 3)) #Matrice vide de la taille du nombre de cellules
    for i, cell in enumerate(mesh.cells()): #On rentre les cellules das une matrice
        mesh_cells[i] = cell.ids()
    
    return mesh, pg_pos, mesh_pos, mesh_cells

