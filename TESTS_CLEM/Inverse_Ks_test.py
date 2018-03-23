
# coding: utf-8

# # Setup file python 

# In[ ]:





# In[17]:


from definition_maillage import *
from ecriture_fich_gprMax import *
from ecriture_in import *
from interp_charge_trou import *
from interp_maillage_gprMax import *
    
#get_ipython().run_line_magic('matplotlib', 'inline')

import numpy as np
import os
import matplotlib.pyplot as plt
import sys, h5py, binascii

import pygimli as pg
from pygimli.meshtools import appendTriangleBoundary, merge2Meshes, mesh
from pygimli.mplviewer import drawMesh
from pygimli.viewer import showMesh

#Def des paramètres de géométrie du modéle
xmin,xmax = 0, 40 # en cm
emin,emax = 0, 80 #  elevation en cm
dtrou = 30 # prof du trou en cm
etrou = emax - dtrou # elevation du fond du trou
r=2 # rayon du trou en cm
dx = .1 #On définit le pas de la maille
zaff= 20 #profondeur en cm jusqu'où on souhaite un maillage affiné. 
eaff=etrou-zaff

h_0=-95 #charge initiale en cm, soit l'état initial du sol (teneur en eau exprimée en charge)
h_1=10 #hauteur d'eau au fond du trou en cm

#Définition du maillage triangulaire
[mesh, pg_pos, mesh_pos, mesh_cells]=definition_maillage(xmin, xmax, emin, emax, dtrou, etrou, r, dx, zaff, eaff)

def run(Ks):
#Début de la boucle de calcul
#for i in range(0,1):
    tr = 0.06
    ts = 0.3
    alpha = 0.016
    n = 8.52
    Ks = Ks[0]
    param = [tr, ts, tr, ts, alpha, n, Ks, Ks, ts] #Paramètres d'entrée tr, ts, alpha, n, Ks
    temps=[0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 1.17, 1.33, 1.50, 1.67, 1.83, 2.00, 2.17, 2.33, 2.50, 2.67, 2.83, 3.00, 3.17, 3.33, 3.50, 3.67, 3.83, 4.00, 4.17, 4.33, 4.50, 4.67, 4.83, 5.00, 5.17, 5.33, 5.50, 5.67, 5.83, 6.00]#Temps de calcul en minutes
    #temps=[1.00, 2.00, 3.00, 4.00, 5.00, 6.00]#Temps de calcul en minutes
    print(param)

    #Interpolation de la charge sur les bords du trou
    [p, charge_imp, n_free_drainage, n_total]=interp_charge_trou(mesh_pos, etrou, r, h_0, h_1)

    #Création des fichiers .in nécessaires pour SMWS_2D
    ecriture_in(mesh, param, temps, p, n_total)

    #Lancement SWMS_2D
    get_ipython().system('cp Grid.in SWMS_2D.IN/')
    get_ipython().system('cp Selector.in SWMS_2D.IN/')
    get_ipython().system('~/Codes/HD2/H2D')

    #Interpolation du maillage sur une grille rectangulaire 2D pour gprMax
    #On crée un maillage rectangulaire avec les dimensions du modèle
    xmin, xmax = 0.0, 40.
    zmin, zmax = 0.0, 80.

    dx = 1.0 #en cm

    #fich_thetas = "TEST_AFF/th.out" #Fichier contenant les thetas
    fich_thetas = "SWMS_2D.OUT/th.out" #Fichier contenant les thetas

    [xv, yv, T, mx, my, mesh2, grid, grid_mat, eps_mat] = interp_maillage_gprMax(xmin, xmax, zmin, zmax, dx, param, mesh, mesh_pos[:,:2], fich_thetas, temps)

    import math
    #Lancement gprMax
    nom='Forward'

    sigma=0.0000
    #eps_w=80.1 #Epsilon eau
    eps_w=80.1
    eps_pvc=3 #Epsilon pvc
    eps_s=2.5 #Epsilon sable
    p=param[1] #Porosité = theta s
    #Maxieps=( math.sqrt(eps_w)*0.9*p + (1-p)*math.sqrt(eps_s) +0.1*p )**2
    Maxieps=( math.sqrt(eps_w)*p + (1-p)*math.sqrt(eps_s) )**2
    Minv=0.3/math.sqrt(Maxieps)
    #Minv=0.3/eps_w
    Minlambda=Minv*(10**9)/(2800*(10**6))
    dl=Minlambda/10

    time=0.000000030
    dx=abs((xv.T*0.01)[0,1]-(xv.T*0.01)[0,0])
    dy=abs((yv.T*0.01)[1,0]-(yv.T*0.01)[0,0])
    #dl=dx/2
    #dl=dx
    fac_dt = 0.2

    #Définition de ce qu'on a besoin pour gprMax
    import os
    import math
    from scipy.interpolate import griddata

    #for i in range(0,len(T)) :
    #    grid_z0 = griddata((mx, my), grid[:,i], (xv, yv), method='linear', fill_value=0.0).T

    materiaux = {}

    def materiau(x):
        if x in materiaux:
            return materiaux[x]
        valeur = "sand{}".format(len(materiaux))
        materiaux[x] = valeur
        return valeur
    for i in range(0,len(T)+1):
        for j in range(0,81):
            for k in range(0,41):
                #grid_mat[i][j,k]=round(grid_mat[i][j,k],3)
                #grid_mat[i][j,k]=grid_mat[i][j,k]*2
                #f[i][j,k]=CRIM(grid_mat[i][j,k], eps_w,eps_s,p)
                materiau(grid_mat[i][j,k])
        #nmedia=len((yv.T*0.01)[:,0])*len(grid_mat[i][0,:])
        nmedia=len(materiaux)+2
        ecriture_fich_gprMax(param, xv.T*0.01, yv.T*0.01, T, grid_mat[i], i, param[6], nom, h_1*0.01, etrou*0.01, r*0.01, sigma, eps_w, eps_pvc, eps_s, p, nmedia, time, dx, dy, dl, fac_dt, materiaux)
        fichier=nom+str(int(param[6]*1000))+'_'+str(i+1)+'.in'
        command="../gprMax "+fichier
        #print(command)
        os.popen(command).readlines()

    merge=nom+str(int(param[6]*1000))+'_'
    copy=nom+str(int(param[6]*1000))+'__merged.out'
    command2="../gprMaxMerge "+merge
    command3="cp "+copy+" RESULTS/."
    os.popen(command2).readlines()
    os.popen(command3).readlines()

    import h5py
    import math

    def picking(filename) :
        f = h5py.File(filename, 'r')
        path = '/rxs/rx1/'
        #modelruns = f.attrs['Modelruns']
        modelruns=len(temps)+1
        samples = f.attrs['Iterations']
        dt = f.attrs['dt']*1e9
        #positions = f.attrs['Positions'][:,0,0]
        dx = 1
        #dx = np.diff(positions)[0]
        data = np.ones((samples, modelruns))
        t_max = np.zeros(modelruns)
        tt=np.zeros(modelruns)
        tps = np.zeros(modelruns)

        #Calcul du temps d'arrivé de la première réflexion:
        h = math.sqrt(0.3**2 + 0.2**2)
        v_init=0.3/(math.sqrt(3.854))
        t_init = (2*h)/v_init
        itmin=(t_init/dt,0)
        #print(itmin)
        for model in range(0,modelruns):
            data[:,model] = f['%s%s' % (path, 'Ez')][:,model]
            if model==0 :
                itmin_i=itmin
            else :
                itmin_i=(tt[model-1],0)
            #print(itmin_i)
            t_max[model] = np.max(data[int(itmin_i[0]):,model])
            tt[model] = np.where(data[:,model]==t_max[model])[0]
            tps[model] = tt[model]*dt   
        #print(t_max)
        #print(tt)
        tps=tps-tps[0]
        dx_dt=(dx,dt)
        #print(tps)
        return tps

    os.popen("rm -rf *.in")
    os.popen("rm -rf *.out")

    tps=picking('RESULTS/'+copy)
    print(tps)
    return(tps)
    
import spotpy
import numpy as np

class spotpy_setup(object):
    
    def __init__(self):
        self.params=[spotpy.parameter.Uniform('Ks', 0.05 , 0.1 , 0.01 , 0.2 , 0.05 , 0.1)]
            
    def parameters(self):
        return spotpy.parameter.generate(self.params)
    
    def simulation(self, vector):
        x=np.array(vector)
        #print(x)
        simulations=run(x)
        #print(simulations)
        return simulations
    
    def evaluation(self):
        observations=np.loadtxt('Tps_arriv_TEST.dat') #Ouvrir le fichier
        #print(observations)
        return observations
    
    def objectivefunction(self, simulation, evaluation):
        objectivefunction=-spotpy.objectivefunctions.rmse(evaluation, simulation)
        return objectivefunction
    
spotpy_setup=spotpy_setup()

sampler = spotpy.algorithms.sceua(spotpy_setup,dbname='SCEUA_CMF',dbformat='csv', parallel='mpi')
results=[]
algorithms=['SCE-UA']

for algorithm in algorithms:
    sampler.sample(1,ngs=2,kstop=100,pcento=0.001,peps=0.01)
    results=sampler.getdata()


# In[18]:


spotpy_setup=spotpy_setup()

sampler = spotpy.algorithms.sceua(spotpy_setup,dbname='SCEUA_CMF',dbformat='csv')
results=[]
algorithms=['SCE-UA']

for algorithm in algorithms:
    sampler.sample(1,ngs=2,kstop=100,pcento=0.001,peps=0.01)
    results=sampler.getdata()


# In[22]:


get_ipython().system('mpirun -c 4 function_parallel.py')


# In[ ]:


def picking(filename) :
    f = h5py.File(filename, 'r')
    path = '/rxs/rx1/'
    #modelruns = f.attrs['Modelruns']
    modelruns=len(temps)+1
    samples = f.attrs['Iterations']
    dt = f.attrs['dt']*1e9
    #positions = f.attrs['Positions'][:,0,0]
    dx = 1
    #dx = np.diff(positions)[0]
    data = np.ones((samples, modelruns))
    t_max = np.zeros(modelruns)
    tt=np.zeros(modelruns)
    tps = np.zeros(modelruns)

    #Calcul du temps d'arrivé de la première réflexion:
    h = math.sqrt(0.3**2 + 0.2**2)
    v_init=0.3/(math.sqrt(3.854))
    t_init = (2*h)/v_init
    itmin=(t_init/dt,0)
    #print(itmin)
    for model in range(0,modelruns):
        data[:,model] = f['%s%s' % (path, 'Ez')][:,model]
        if model==0 :
            itmin_i=itmin
        else :
            itmin_i=(tt[model-1],0)
        #print(itmin_i)
        t_max[model] = np.max(data[int(itmin_i[0]):,model])
        tt[model] = np.where(data[:,model]==t_max[model])[0]
        tps[model] = tt[model]*dt   
    #print(t_max)
    #print(tt)
    tps=tps-tps[0]
    dx_dt=(dx,dt)
    #print(tps)
    return tps


# In[ ]:


copy='Forward87__merged.out'
temps=[0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 1.17, 1.33, 1.50, 1.67, 1.83, 2.00, 2.17, 2.33, 2.50, 2.67, 2.83, 3.00, 3.17, 3.33, 3.50, 3.67, 3.83, 4.00, 4.17, 4.33, 4.50, 4.67, 4.83, 5.00, 5.17, 5.33, 5.50, 5.67, 5.83, 6.00]#Temps de calcul en minutes
tps=picking('RESULTS/'+copy)
observations=np.loadtxt('Tps_arriv_TEST.dat') #Ouvrir le fichier


# In[11]:


observations


# In[12]:


simulations=tps


# In[ ]:


def objectivefunction(self, simulation, evaluation):
    objectivefunction=-spotpy.objectivefunctions.rmse(evaluation, simulation)
    return objectivefunction


# In[20]:


spotpy.analyser.plot_parametertrace(results)
#spotpy.analyser.plot_parametertrace_algorithms(results,algorithmnames=algorithms,parameternames=['Ks'])

