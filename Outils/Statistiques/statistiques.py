#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Import des modules
from math import *
import os





def statistiques(nom_fichier):
    """
    Fonction permettant d'effectuer des statistiques sur les solutions obtenues avec l'ACO dans le but de qualifier la convergence de l'algorithme
    """
    # Les statistiques sont comparées aux solutions du front de pareto de toutes les solutions
    try:
        ensembleCompare=set()
        fichierCompare= open('paretobis.txt','r')
        for ligne in fichierCompare:
            tableau=ligne.split()
            ensembleCompare.add(tableau[0])
    except:
        raise
    finally:
        fichierCompare.close()

    # On va écrire tous les resultats dans le fichier ResultatStats.txt
    try:
        fichierResultat=open("ResultatStats.txt","a")
        fichierResultat.write(nom_fichier+'\n\n')

        # Calculs statistiques
        try:
            print(nom_fichier)
            fichier = open(nom_fichier, "r")
            liste_sigmax=[]
            liste_sigmay=[]
            liste_sigmaz=[]
            liste_rac_sigmas=[]
            ensembleACO=set()

            for ligne in fichier:
                ligne_texte=ligne.split()
                ensembleACO.add(ligne_texte[0])
                sigma_x=float(ligne_texte[1])*1000
                sigma_y=float(ligne_texte[2])*1000
                sigma_z=float(ligne_texte[3])*1000
                rac_sigmas=sqrt(sigma_x**2+sigma_y**2+sigma_z**2)
                liste_sigmax.append(sigma_x)
                liste_sigmay.append(sigma_y)
                liste_sigmaz.append(sigma_z)
                liste_rac_sigmas.append(rac_sigmas)

            max_sigmax = max(liste_sigmax)
            fichierResultat.write("max_sigmax : %f" % max_sigmax)
            fichierResultat.write("\n")

            min_sigmax = min(liste_sigmax)
            fichierResultat.write("min_sigmax : %f" % min_sigmax)
            fichierResultat.write("\n")

            max_sigmay = max(liste_sigmay)
            fichierResultat.write("max_sigmay : %f" % max_sigmay)
            fichierResultat.write("\n")

            min_sigmay = min(liste_sigmay)
            fichierResultat.write("min_sigmay : %f" % min_sigmay)
            fichierResultat.write("\n")

            max_sigmaz = max(liste_sigmaz)
            fichierResultat.write("max_sigmaz : %f" % max_sigmaz)
            fichierResultat.write("\n")

            min_sigmaz = min(liste_sigmaz)
            fichierResultat.write("min_sigmaz : %f" % min_sigmaz)
            fichierResultat.write("\n")

            max_rac_sigmas = max(liste_rac_sigmas)
            fichierResultat.write("max_rac_sigmas : %f" % max_rac_sigmas)
            fichierResultat.write("\n")

            min_rac_sigmas = min(liste_rac_sigmas)
            fichierResultat.write("min_rac_sigmas : %f" % min_rac_sigmas)
            fichierResultat.write("\n")

            moyenne = sum(liste_rac_sigmas)/len(liste_rac_sigmas)
            fichierResultat.write("moyenne : %f" % moyenne)
            fichierResultat.write("\n")

            if len(liste_rac_sigmas)%2 == 0:
                mediane = (liste_rac_sigmas[int((len(liste_rac_sigmas)/2))]+liste_rac_sigmas[int(len(liste_rac_sigmas)/2+1)])/2
            else:
                mediane = liste_rac_sigmas[floor(len(liste_rac_sigmas)/2)]
            fichierResultat.write("mediane : %f" % mediane)
            fichierResultat.write("\n")

            barycentre = []
            XG = sum(liste_sigmax)/len(liste_sigmax)
            YG = sum(liste_sigmay)/len(liste_sigmay)
            ZG = sum(liste_sigmaz)/len(liste_sigmaz)
            barycentre.append(XG)
            barycentre.append(YG)
            barycentre.append(ZG)
            fichierResultat.write("barycentre : %s" % barycentre)
            fichierResultat.write("\n")

            dist_G_origine = sqrt(XG**2+YG**2+ZG**2)
            fichierResultat.write("dist_G_origine : %f" % dist_G_origine)
            fichierResultat.write("\n")

            L=[]
            for i in range(len(liste_rac_sigmas)):
                L.append((liste_rac_sigmas[i]-moyenne)**2)
            somme=sum(L)

            ecart_type = sqrt(1/len(liste_rac_sigmas)*somme)
            fichierResultat.write("ecart_type : %f" % ecart_type)
            fichierResultat.write("\n")

            coeff_var = ecart_type/moyenne
            fichierResultat.write("coeff_var : %f" % coeff_var)
            fichierResultat.write("\n")

            ensembleVraiPositif = ensembleACO & ensembleCompare
            ensembleFauxPositif = ensembleACO - ensembleVraiPositif

            taux_faux_positif=((len(ensembleACO)-len(ensembleVraiPositif))/(len(ensembleACO)))*100
            fichierResultat.write("Taux de Faux positifs : %f %%" % taux_faux_positif)
            fichierResultat.write("\n")

            taux_vrai_positif = (len(ensembleVraiPositif)/len(ensembleCompare))*100
            fichierResultat.write("Taux de Vrais positifs : %f %%" % taux_vrai_positif)
            fichierResultat.write("\n")


        except:
            raise
        finally:
            fichier.close()
    except:
        raise
    finally:
        fichierResultat.write("\n\n\n")
        fichierResultat.close()





if __name__=="__main__":
    # Pour chaque fichier du dossier Result, on calcule les statistiques
    for filename in os.listdir('Result'):
        statistiques('Result/'+filename)

