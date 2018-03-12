import numpy as np

def interp_charge_trou(mesh_pos, d, radius, h_0, h_1):
    """Interpolation de la valeur de la charge sur les parois du trou"""
    
    p = np.zeros((len(mesh_pos), 4), dtype = 'float') #matrice vide
    p[:,:2] = np.array(mesh_pos[:,:2]) #coordonnées xy des noeuds du maillage
    
    def charge_trou(p,d,radius,h_0,h_1):
    
        for i in range(0,len(p)) :
            if p[i,1]==d and p[i,0]<=radius :
                p[i,2]=1
                p[i,3]=h_1
            elif p[i,1]==0 :
                p[i,2]=-3
                p[i,3]=h_0
            else :
                p[i,3]=h_0
           
    charge_trou(p,d,radius,h_0,h_1)        
        
    #Interpolation pour la valeur de charge au fond du trou
    r = np.array(np.where((p[:,1]<=d+h_1) & (p[:,0]==radius))).T
    for i in range(1,len(r)):
        p[r[i,0],2] = 1
        p[r[i,0],3] = h_1-(p[r[i,0],1]-d)
     
    n_charge_imp=len(np.array(np.where(p[:,2]==1)).T)
    n_free_drainage=len(np.array(np.where(p[:,2]==-3)).T)
    n_total=n_charge_imp+n_free_drainage
        
    return p, n_charge_imp, n_free_drainage, n_total

