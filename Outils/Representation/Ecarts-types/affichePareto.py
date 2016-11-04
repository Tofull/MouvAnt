#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Import des modules
import os #pour gerer les fichiers
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import gridspec






# Préparation de la fenetre d'affichage
gs = gridspec.GridSpec(3,1)
fig = plt.figure()
ax = fig.add_subplot(gs[0])
ay = fig.add_subplot(gs[1])
az = fig.add_subplot(gs[2])
ax.set_xlabel('sigmaX')
ax.set_ylabel('sigmaY')
ay.set_xlabel('sigmaX')
ay.set_ylabel('sigmaZ')
az.set_xlabel('sigmaY')
az.set_ylabel('sigmaZ')





# Chargement du fichier contenant les datums et leurs ecarts-types
fichier = open('paretoBisJoin.txt', "r")
condition=True #Pour afficher la legende qu'une seule fois
#Parcours du fichier
for ligne in fichier:
    calculSigma=ligne.split(' ')
    sigmaX=float(calculSigma[1])
    sigmaY=float(calculSigma[2])
    sigmaZ=float(calculSigma[3])
    if (condition):
        ax.scatter(sigmaX,sigmaY, c="blue", marker="o",label='Non Dominees')
        ay.scatter(sigmaX,sigmaZ, c="blue", marker="o",label='Non Dominees')
        az.scatter(sigmaY,sigmaZ, c="blue", marker="o",label='Non Dominees')
    else:
        ax.scatter(sigmaX,sigmaY, c="blue", marker="o")
        ay.scatter(sigmaX,sigmaZ, c="blue", marker="o")
        az.scatter(sigmaY,sigmaZ, c="blue", marker="o")
    condition=False
fichier.close()
            
            
            
# Chargement du fichier contenant les datums et leurs ecarts-types
fichier = open('fichierFauxPositifsMOGA.txt', "r")
condition=True #Pour afficher la legende qu'une seule fois
#Parcours du fichier
for ligne in fichier:
    calculSigma=ligne.split(' ')
    sigmaX=float(calculSigma[1])
    sigmaY=float(calculSigma[2])
    sigmaZ=float(calculSigma[3])
    if (condition):
        ax.scatter(sigmaX,sigmaY, c="red", marker="o",label='Faux Positif')
        ay.scatter(sigmaX,sigmaZ, c="red", marker="o",label='Faux Positif')
        az.scatter(sigmaY,sigmaZ, c="red", marker="o",label='Faux Positif')
    else:
        ax.scatter(sigmaX,sigmaY, c="red", marker="o")
        ay.scatter(sigmaX,sigmaZ, c="red", marker="o")
        az.scatter(sigmaY,sigmaZ, c="red", marker="o")
    condition=False
fichier.close()



# Chargement du fichier contenant les datums et leurs ecarts-types
fichier = open('fichierVraiPositifsMOGA.txt', "r")
condition=True #Pour afficher la legende qu'une seule fois
#Parcours du fichier
for ligne in fichier:
    calculSigma=ligne.split(' ')
    sigmaX=float(calculSigma[1])
    sigmaY=float(calculSigma[2])
    sigmaZ=float(calculSigma[3])
    if (condition):
        ax.scatter(sigmaX,sigmaY, c="green", marker="o",label='Vrai Positif')
        ay.scatter(sigmaX,sigmaZ, c="green", marker="o",label='Vrai Positif')
        az.scatter(sigmaY,sigmaZ, c="green", marker="o",label='Vrai Positif')
    else:
        ax.scatter(sigmaX,sigmaY, c="green", marker="o")
        ay.scatter(sigmaX,sigmaZ, c="green", marker="o")
        az.scatter(sigmaY,sigmaZ, c="green", marker="o")
    condition=False
fichier.close()





# Affichage de la fenetre
plt.legend()
plt.grid(True)
plt.show()
