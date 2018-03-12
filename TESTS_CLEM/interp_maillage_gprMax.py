
def interp_maillage_gprMax(xmin, xmax, zmin, zmax, dx, param, mesh, mesh_pos, fich_thetas, temps):    
    
    import numpy as np
    import pygimli as pg
    from pygimli.meshtools import appendTriangleBoundary, merge2Meshes
    from pygimli.mplviewer import drawMesh
    from pygimli.viewer import showMesh
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    from pygimli.mplviewer import drawMesh, drawModel
    from pygimli.meshtools import interpolate
    from pygimli.meshtools import nodeDataToCellData

    
    xreg = np.arange(xmin, xmax + dx, dx, 'float')
    zreg = np.arange(zmin, zmax + dx, dx, 'float')

    mesh2 = pg.Mesh(3)
    mesh2.createGrid(xreg, zreg)
    for c in mesh2.cells():
        c.setMarker(3)
        
    pg_pos2 = mesh2.positions()
    mesh2_pos2 = np.array((np.array(pg.x(pg_pos2)), np.array(pg.y(pg_pos2)), np.array(pg.z(pg_pos2)))).T #On crée une matrice contenant la position des noeuds
    mesh2_cells2 = np.zeros((mesh2.cellCount(), 4)) #Matrice vide de la taille du nombre de cellules
    for i, cell in enumerate(mesh2.cells()): #On rentre les cellules dans une matrice
        mesh2_cells2[i] = cell.ids()
        
    mx = pg.x(mesh2.cellCenter())
    my = pg.y(mesh2.cellCenter())
    mesh_pos2=mesh2_pos2[:,0:2]
    
    xv, yv = np.meshgrid(xreg, zreg, sparse=False, indexing='ij')
    
    T=temps
    
    maillage = mesh_pos #maillage triangulaire que l'on a défini pour SWMS_2D
    #maillage = np.loadtxt("TEST_AFF/mesh-40-80-2-30.dat")
    (x,z)=np.shape(maillage)
    
    eps_w=80.1
    eps_s=2.5
    p=param[1]
    
    def CRIM(x,eps_w,eps_s,p):
        import math
        y=round(( math.sqrt(eps_w)*x+(1-p)*math.sqrt(eps_s)+(p-x) )**2,3)
        return(y)
    
    theta=np.loadtxt(fich_thetas) #Ouvrir le fichier
    
    eps=np.zeros(len(theta))

    for i in range(0,len(theta)):
        eps[i]=CRIM(theta[i], eps_w, eps_s, p)
    
    eps_mat=np.zeros([x,int((len(eps)/x))])
    for i in range(0,len(T)+1) : #
        xi=i*x
        eps_mat[:,i]=eps[xi:(xi+x)]
    
    fig, ax = plt.subplots(len(T)+1, figsize=(80, 500))
    grid_lin=np.zeros([len(mx),len(T)+1])
    grid_mat={}
    for i in range(0, len(T)+1) : #
        grid=np.zeros([len(xv[:,0]), len(xv[0,:])])
        outdata=interpolate(mesh2,mesh,eps_mat[:,i], fill_value=min(eps))
        outdata2=nodeDataToCellData(mesh2,outdata)
        for j in range(0,len(xv[0,:])):
            k=j*len(xv[:,0])
            kk=len(xv[:,0])
            grid[:,j]=np.around(outdata[k:(k+kk)], decimals=1)
        grid_mat[i]=grid.T
        a=np.where(grid_mat[i]==0.0)
        grid_mat[i][a]=min(eps)
        #drawModel(ax[i], mesh2 , outdata2)
    
    return xv, yv, T, mx, my, mesh2, grid, grid_mat, eps_mat

