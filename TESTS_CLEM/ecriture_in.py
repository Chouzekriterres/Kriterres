import numpy as np

def ecriture_in(mesh, param, temps, p, n_total):
    """Fonction qui permet d'écrire tous les fichiers nécessaires pour lancer SWMS_2D"""
    
    
    paramlist=" ".join([str(i) for i in param])
    tempslist=" ".join([str(i) for i in temps])

    s = """*** BLOCK A: BASIC INFORMATION *****************************************
    'Heading'
    Example 1 - Column Test
    LUnit  TUnit  MUnit  BUnit     (units are obligatory for all input data)
     'cm'   'min'  '-'    '-'
    Kat (0:horizontal plane, 1:axisymmetric vertical flow, 2:vertical plane
    1
    MaxIt   TolTh   TolH       (maximum number of iterations and tolerances)
      21    .001   0.5
    lWat    lChem   ChecF   ShortF  FluxF   AtmInF  SeepF  FreeD  DrainF
     t  f  f      t       t       f       f      t      f
    *** BLOCK B: MATERIAL INFORMATION **************************************
    NMat    NLay    hTab1   hTabN   NPar
      1      1      .001    200.     9
    thr     ths     tha     thm     Alfa    n       Ks      Kk      thk
     {}
    *** BLOCK C: TIME INFORMATION ******************************************
    dt      dtMin   dtMax   DMul    DMul2   MPL
     .1    .01     10.     1.1     .7     {}
    TPrint(1),TPrint(2),...,TPrint(MPL)                   (print-time array)
       {}
    *** END OF INPUT FILE SELECTOR.IN************************************
    """.format(paramlist,len(temps),tempslist)
    
    fselector=open("Selector.in","w")

    fselector.write(s)

    fselector.close()
    
    Ncells = mesh.cellCount() #nombre de cellules du maillage
    Nnodes = mesh.nodeCount() #nombre de noeuds du maillage
    Nbc = n_total #nombre de noeuds avec Boundary conditions

    t = np.zeros((Ncells, 3), dtype = 'float')
    for i in range(0,Ncells) :
        c=mesh.cell(i)
        for j in range(0,3) :
            l=c.node(j)
            f=l.id()
            t[i,j]=f
    
    a = (np.array(np.where(p[:,2]==1))).T #Noeuds ayant une charge imposée constante
    b = (np.array(np.where(p[:,2]==-3))).T #Noeuds ayant un drainage
    
    def width(p,a,b) :
        import math
        e = np.zeros((len(a), 3), dtype = 'float')
        f = np.zeros((len(b), 3), dtype = 'float')
        e[:,0] = a[:,0]
        f[:,0] = b[:,0]
        e[:,1] = p[a[:,0],0]
        f[:,1] = p[b[:,0],0]
    
        for i in range(0,len(a[:,0])-1) :
            if i==0 : #Premier noeud
                e[i,2] = ((math.pi)/3) * ( (p[a[i+1,0],0] + 2* p[a[i,0],0]) * (p[a[i+1,0],0]-p[a[i,0],0]))
            elif i==(len(a[:,0])-1) : #Dernier noeud
                e[i,2] = ((math.pi)/3) * ( p[a[i-1,0],0] + 2* p[a[i,0],0] * (p[a[i,0],0]-p[a[i-1,0],0]))
            else : #Autres noeuds
                e[i,2] = ((math.pi)/3) * ( (p[a[i-1,0],0] + 2* p[a[i,0],0]) * (p[a[i,0],0]-p[a[i-1,0],0]) + (p[a[i+1,0],0] + 2* p[a[i,0],0]) * (p[a[i+1,0],0]-p[a[i,0],0]))
            
        for i in range (0,len(b[:,0])-1) :
            if i==0 :
                f[i,2] = ((math.pi)/3) * ((p[b[i+1,0],0] + 2* p[b[i,0],0]) * (p[b[i+1,0],0]-p[b[i,0],0]))
            elif i==(len(b[:,0])-1) :
                f[i,2] = ((math.pi)/3) * ((p[b[i-1,0],0] + 2* p[b[i,0],0]) * (p[b[i,0],0]-p[b[i-1,0],0]))
            else :
                f[i,2] = ((math.pi)/3) * ((p[b[i-1,0],0] + 2* p[b[i,0],0]) * (p[b[i,0],0]-p[b[i-1,0],0]) + (p[b[i+1,0],0] + 2* p[b[i,0],0]) * (p[b[i+1,0],0]-p[b[i,0],0]))
    
        return e , f #Permet de retourner e et f pour s'en servir après
    
    [e, f] = width(p,a,b)
    
    dim = [Nnodes, Ncells, 2, Nbc, 0]

    dim_list=" ".join([str(i) for i in dim])

    s1 = """*** BLOCK H: NODAL INFORMATION **************************************************
          NumNP     NumEl       IJ      NumBP     NObs
       {}
       n  Code    x      z          h       Conc      Q     M   B    Axz   Bxz   Dxz
    """.format(dim_list)

    s2 = """*** BLOCK I: ELEMENT INFORMATION ************************************************
       e   i   j   k   l   Angle  Aniz1  Aniz2  LayNum
    """

    s3="""*** BLOCK J: BOUNDARY GEOMETRY INFORMATION *************************************
        Node number array:
    """
    
    fgrid=open("Grid.in","w")

    fgrid.write(s1)

    #code = p[:,2] #A définir par clemence pour chaque node
    #h = p[:,3]  #A definir par clemence pour chaque node

    # for i in np.arange(Nnodes):
    #   fgrid.write("""{} {} {} {} {} .00E+00 .00E+00 1 0 1 1 1 \n""".format(i + 1, code, node.x(i), node.y(i), h)) 

    for i in range(0,len(p)) :
        fgrid.write("""{} {} {} {} {} .00E+00 .00E+00 1 0 1 1 1 \n""".format(i + 1, int(p[i,2]), round(p[i,0],2), round(p[i,1],2), round(p[i,3],2)))
    
    #Attention, il faut enlever une ligne mais je ne vois pas comment faire...
    fgrid.write(s2)

    for i in range(0,len(t)) :
        fgrid.write("""{} {} {} {} {} 0 1 1 1 \n""".format(i + 1, int(t[i,0]+1), int(t[i,1]+1), int(t[i,2]+1), int(t[i,2]+1)))

    fgrid.write(s3)

    k = 1

    for i in range(0,len(e[:,0])) :
        fgrid.write("""{} \t""".format(int(e[i,0]+1)))
        k = k + 1
        if k==7 :
            fgrid.write("""\n""")
            k = 1
        
    k = 1

    fgrid.write("""\n""")
    for i in range(0,len(f[:,0])) :
        fgrid.write("""{} \t""".format(int(f[i,0])))
        k = k + 1
        if k==7 :
            fgrid.write("""\n""")
            k = 1

    k = 1 

    fgrid.write("""\n""")
    fgrid.write("""Width array:\n""")
    for i in range(0,len(e[:,0])) :
        fgrid.write("""{} \t""".format(round(e[i,2],3)))
        k = k + 1
        if k==7 :
            fgrid.write("""\n""")
            k = 1
    
    fgrid.write("""\n""")

    k = 1 

    for i in range(0,len(f[:,0])) :
        fgrid.write("""{} \t""".format(round(f[i,2],3)))
        k = k + 1
    if k==7 :
        fgrid.write("""\n""")
        k=1

    fgrid.write("""\n""")
    fgrid.write("""Length:\n""")
    fgrid.write("""0\n""")

    fgrid.write("""*** END OF INPUT FILE GRID.IN *************************************************\n""")

    fgrid.close()
    
    return 
    





