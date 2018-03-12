import numpy as np
import math

def ecriture_fich_gprMax(param, X, Y, h, grid_z0, name, ks, nom, h_eau, etrou, radius, sigma, eps_w, eps_pvc, eps_s, p, nmedia, time, dx, dy, dl, fac_dt, materiaux) :
    
    fgrid=open(nom+str(int(ks*1000))+'_'+str(name+1)+'.in',"w")
    #fd=fopen([name '_' nom num2str(int16(Ks*1000)) '.in'],'w'); ce qu'on veut en matlab
    
    fgrid.write("""------------------------------------------------\n""")
    #fgrid.write("""#number_of_media: {}\n""".format(nmedia+4))
    #fgrid.write("""#domain: {} {} {}\n""".format(round(2*max(X[0,:])+dx,2), round(max(Y[:,0])*2,2), dl))
    fgrid.write("""#domain: {} {} {}\n""".format(round(2*max(X[0,:])+dx,2), round(max(Y[:,0])*2,2), dl))
    fgrid.write("""#dx_dy_dz: {} {} {}\n""".format(dl, dl, dl))
    fgrid.write("""#time_window: {}\n""".format(time))
    #fgrid.write("""#time_step_stability_factor: {}\n""".format(fac_dt))
    #fgrid.write("""#abc_type: pml\n""")
    #fgrid.write("""#pml_cells: {}\n""".format(10))
    fgrid.write("""#waveform: ricker 1.0 1000e6 Mydipole\n""")
    w=max(X[0,:])+dx/2 #largeur
    fgrid.write("""#hertzian_dipole: z {} {} {} Mydipole\n""".format(round(w+0.18,2), round(max(Y[:,0])+3/2*dy,2), 0.0))#TX
      
    #fgrid.write("""#analysis: 1 {}_{}{}.out b\n""".format(name, nom, int(ks*1000)))
      
    #fgrid.write("""#tx: {} {} Mydipole 0.0 {}\n""".format(w+0.3, max(Y[:,0])+3/2*dy, time))
    fgrid.write("""#rx: {} {} {} \n""".format(round(w+0.22,2), round(max(Y[:,0]+3/2*dy),2), 0.0))#RX
    
    fgrid.write("""#num_threads: 4 \n""")
    #fgrid.write("""#end_analysis:\n""")
    fgrid.write("""------------------------------------------------\n""")
    
    #k=1
    #VWC=np.zeros([len(Y[:,0])*len(grid_z0[0,:])+1, 1])
    #for i in range(0, len(Y[:,0])) :
        #for j in range(0, len(grid_z0[0,:])) :
            #grid_z0[i,j]=( math.sqrt(eps_w)*grid_z0[i,j]+(1-p)*math.sqrt(eps_s)+(p-grid_z0[i,j]) )**2 #Formule de CRIM
            #fgrid.write("""#material: {} {} 1.0 0.0 sand{}\n""".format(grid_z0[i,j], sigma, k))
            #VWC[k-1]=grid_z0[i,j]
            #k = k+1
            
    for i in materiaux:
        fgrid.write("""#material: {} {} 1.0 0.0 {}\n""".format(i, sigma, materiaux[i]))
            
    #fgrid.write("""#material: {} {} 1.0 0.0 eau\n""".format(eps_w,sigma)) #Donne le epsilon de l'eau
    fgrid.write("""#material: {} {} 1.0 0.0 pvc\n""".format(eps_pvc,sigma)) #Donne le epsilon du pvc
    fgrid.write("""------------------------------------------------\n""")
    #Creation du modele en entier par symetrie des resultats SWMS
    #Entier
    w=max(X[0,:])+dx/2 #largeur

    #Demi modèle
    #w=dx/2

    d=dy/2  #hauteur
    k=1 #compteur
    count=1      
    #A=np.zeros([len(mesh_cells[:,0])*10, 3]) #Taille de la matrice = nombre de cellules du nouveau maillage*4
      
    for ii in range(0,len(Y[:,0])) :
          for jj in range(0, len(grid_z0[0,:])) :
              fgrid.write("""#box: {} {} {} {} {} {} {}\n""".format(round(w+X[ii,jj]-dx/2,2), round(d+Y[ii,jj]-dy/2,2), 0.0, round(w+X[ii,jj]+dx/2,2),round(d+Y[ii,jj]+dy/2,2), dl, materiaux[grid_z0[ii,jj]]))
              #A[count,0]=w+X[ii,jj]-dx/2
              #A[count,1]=d+Y[ii,jj]-dy/2
              #A[count,2]=VWC[k-1]
      
              count=count+1
              #A[count,0]=w+X[ii,jj]+dx/2
              #A[count,1]=d+Y[ii,jj]+dy/2
              #A[count,2]=VWC[k-1]
      
              count=count+1
              fgrid.write("""#box: {} {} {} {} {} {} {}\n""".format(round(w-X[ii,jj]-dx/2,2), round(d+Y[ii,jj]-dy/2,2), 0.0, round(w-X[ii,jj]+dx/2,2),round(d+Y[ii,jj]+dy/2,2), dl, materiaux[grid_z0[ii,jj]]))
              #A[count,0]=w-X[ii,jj]-dx/2
              #A[count,1]=d+Y[ii,jj]-dy/2
              #A[count,2]=VWC[k-1]

              k=k+1
              count=count+1
              #A[count,0]=w-X[ii,jj]+dx/2
              #A[count,1]=d+Y[ii,jj]+dy/2
              #A[count,2]=VWC[k-1]
              count=count+1
                
    #Ajout du tuyau pvc dans le trou      
    #fgrid.write("""#box: {} {} {} {} {} {} pvc\n""".format(round(w-dx/2-radius,2), round(etrou+d,2), 0.0, round(w+dx/2+radius,2), round(max(Y[:,0])*1.1,2), dl)) #ajout du tuyau tout le long
    #fgrid.write("""#box: {} {} {} {} {} {} pvc\n""".format(0.38, 0.505, 0.0, 0.43, 0.88, 0.6*dl)) #ajout du tuyau tout le long
    #o=np.where((A[:,0]<=w+dx/2+radius) & (A[:,0]>=w-dx/2-radius) & (A[:,1]<=max(Y[:,0])*1.1) & (A[:,2]>=etrou+d+h_eau))
    #A[o,2]=3

    #Ajout d'eau au fond du trou
    #fgrid.write("""#box: {} {} {} {} {} {} eau\n""".format(w-dx/2-radius+0.007, etrou+d, 0.0, w+dx/2+radius-0.007, etrou+d+h_eau, dl))
    #q=np.where((A[:,0]<=w+dx/2+radius-0.007) & (A[:,0]>=w-dx/2-radius+0.007) & (A[:,1]<=etrou+d+h_eau) & (A[:,1]>=etrou+d))
    #A[q,2]=81

    #Ajout de l'air dans le trou
    fgrid.write("""#box: {} {} {} {} {} {} free_space\n""".format(w-dx/2-radius+0.007, etrou, 0.0, w+dx/2+radius-0.007, max(Y[:,0])*2, dl))
    #s=np.where((A[:,0]<=w+dx/2+radius-0.007) & (A[:,0]>=w-dx/2-radius+0.007) & (A[:,1]<=max(Y[:,0])*1.1) & (A[:,1]>=etrou+d+h_eau))
    #A[s,2]=1
    
    fgrid.write("""------------------------------------------------\n""")
    fgrid.write("""#messages: n\n""")
    fgrid.write("""#geometry_view: 0 0 0 {} {} {} {} {} {} modele_vue n""".format(2*max(X[0,:])+dx, max(Y[:,0])*2, dl, dl, dl, dl))

    fgrid.close()
    
    return
                     

