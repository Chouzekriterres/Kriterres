# GPR-processing-tools

Ce dossier contient des fonctions de traitement de données radar de sol. Les fonctions sont commentées.

# Pour le faire marcher démarrer

Install miniconda (comment?)

Installer les dépendences:

    > conda install h5py matplotlib scipy

Activer l'environnement conda

    > conda env list  
    # conda environments:
    #
    gprMax                   /home/sainteno/miniconda3/envs/gprMax
    root                  *  /home/sainteno/miniconda3

    > source activate root

Lancer le traitement:
      
    > python3 run.py Talik_Bscan_2D_merged.out




Dans Jupyter:

    conda install jupyter

    jupyter notebook demo_user_gain.ipynb
