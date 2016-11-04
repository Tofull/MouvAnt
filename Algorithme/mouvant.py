#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Import de modules
from numpy.random import randint #pour gerer le tirage aléatoire des criteres et du nombre de stations
from numpy.random import choice #pour gerer le tirage aléatoire pondere
from ctypes import * #pour gerer l'appel C
import os #pour gerer les fichiers
from numpy import Inf
from math import *

import copy
import sys

#Appel à la bibliothèque C créée avec
#gcc -shared -Wl,-soname,mouvant -o mouvant.so -fPIC mouvant.c
zelib = CDLL(os.getcwd()+"/mouvant.so")


def lst_objet(nom_fichier):
    """
    Lecture du fichier pour le multiobjectif
    """
    try:  
        fichier = open(nom_fichier, "r")
        liste=[]
        lst_critere=[]
        lst_contrainte=[]
        lst_champs=[]
        for ligne in fichier:
            # Les lignes commencant par un % sont considérées comme des commentaires
            if ligne[0]!='%':
                # Les lignes commencant par un # sont les en tete des colonnes
                if ligne[0]=='#':
                    ligne=ligne.split()
                    for k in range(len(ligne)-1):
                        #les entetes commencant par un - sont des contraintes
                        if ligne[k+1][0]=='-':
                            lst_contrainte.append(ligne[k+1][1:len(ligne[k+1])])
                            lst_champs.append(ligne[k+1][1:len(ligne[k+1])])
                        #les entetes commencant par un + sont des criteres
                        elif ligne[k+1][0]=='+':
                            lst_critere.append(ligne[k+1][1:len(ligne[k+1])])
                            lst_champs.append(ligne[k+1][1:len(ligne[k+1])])
                        else:
                            lst_champs.append(ligne[k+1])
                else:
                    item={}
                    ligne=ligne.split()
                    for k in range(len(ligne)):
                        item[lst_champs[k]]=ligne[k]
                    # Ajout des champs phéromones et probabilités                    
                    item['pheromone']={}
                    item['probabilite']={}
                    for k in range(len(ligne)-1):
                        item['pheromone'][lst_champs[k+1]]=0
                        item['probabilite'][lst_champs[k+1]]=0
                    liste.append(item)
            
        return liste,lst_critere,lst_contrainte
    
    except:
        raise
    finally:
        fichier.close()


def check_dominance(a, b,liste_objectif,sigmaxMin,sigmayMin,sigmazMin): 
    """ 
    RELATION DE PARETO
    check_dominance compare deux individus a et b par rapport à une liste d'objectifs fixés
    Si a domine b -> renvoie 1
    Si a est dominé par b -> renvoie -1
    renvoie 0 si a et b sont non dominés (on ne peut pas les comparer entre eux)

    """
    flag1=0
    flag2=0
    for objectif_i in liste_objectif: #Pour tous les objectifs
        if a.fonctionPareto(objectif_i) < b.fonctionPareto(objectif_i): #Si f_i(a)<f_i(b) : a domine au moins une fois b
            flag1=1
        else:
            if a.fonctionPareto(objectif_i) > b.fonctionPareto(objectif_i):#Si f_i(a)>f_i(b) : a est dominé au moins une fois par b
                flag2=1
    if flag1==1 and flag2==0: #a domine au moins une fois b et n'est jamais dominé par b -> renvoie 1
        return 1
    else:
        if flag1==0 and flag2==1: #a est dominé au moins une fois par b et ne domine jamais b -> renvoie -1
            return -1
        else: #a domine au moins une fois b et b domine au moins une fois b -> non dominé -> renvoie 0
            return 0
    
    
def print_matrice(matrice):
    """
    Fonction qui affiche une matrice dans la console
    """
    print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrice]))
    
    

class Ant(object):
    """
    Creation d une classe de fourmis qui possedent certaines caracteristiques
    """
    def __init__(self,**kwargs):
        '''
        Constructeur de la classe.
        '''
        self.Objet=[]
        self.Critere=""
        self.NbStation=0
        self.sigmax=0
        self.sigmay=0
        self.sigmaz=0
        self.datum=""
        self.cout=0

    def fonctionCout(self,critere,sigmaxMin,sigmayMin,sigmazMin):
        """
        Fonction cout définie pour l'ACO var mACO-6 
        La fonction cout est associée à un objectif précis : critere.
        Le resultat est compris entre 0 et 1.
        Plus self."critere" se rapproche du minimum trouvé, plus le resultat de fonctionCout est proche de 1.
        """
        if critere=="sigmax":
            return 1/(1+self.sigmax-sigmaxMin)
        elif critere=="sigmay":
            return 1/(1+self.sigmay-sigmayMin)
        elif critere=="sigmaz":
            return 1/(1+self.sigmaz-sigmazMin)
            
    def fonctionPareto(self,critere):
        if critere=="sigmax":
            return self.sigmax
        elif critere=="sigmay":
            return self.sigmay
        elif critere=="sigmaz":
            return self.sigmaz
    
    
def calculNombreStation(nomFichier):
    """
    Fonction pour déterminer le nombre de station à partir du fichier d'entree au format .dat du lareg
    """
    fichier=open(nomFichier, "r")
    
    count=0
    for ligne in fichier:
        if count==1:
           nbr=int(ligne)
           break
        count+=1    
    
    fichier.close()
    
    return nbr

    
    
def creerFichierStation(nomSemaine,nombreStation):
    """
    Fonction qui permet de generer le fichier .stations necessaire à Mouvant
    """
    fichier=open(nomSemaine+".stations","w")
    fichier.write("% les entetes avec des - sont les contraintes, les + les critères à optimiser\n")
    fichier.write("# station    -poids    +sigmax    +sigmay	+sigmaz\n")
    for k in range(1,nombreStation+1):
        fichier.write(str(k)+"\t1\t0\t0\t0\n")
    fichier.close()
    return nomSemaine+".stations"
    
    
    
def statistiques(nomSemaine):
	"""
	Fonction permettant d'effectuer des statistiques sur les solutions obtenues avec l'ACO dans le but de qualifier la convergence de l'algorithme,
    et d'en extraire le meilleur datum
	"""
	try:  
		fichier = open(nomSemaine+".solutionNonDominees", "r")
		liste_datums=[]
		liste_rac_sigmas=[]

		for ligne in fichier:
			ligne_texte=ligne.split()
			datum=str(ligne_texte[0])
			sigma_x=float(ligne_texte[1])*1000
			sigma_y=float(ligne_texte[2])*1000
			sigma_z=float(ligne_texte[3])*1000
			rac_sigmas=sqrt(sigma_x**2+sigma_y**2+sigma_z**2)
			liste_datums.append(datum)
			liste_rac_sigmas.append(rac_sigmas)

		
		min_rac_sigmas = min(liste_rac_sigmas)
		meilleurDatum=liste_datums[liste_rac_sigmas.index(min_rac_sigmas)]
		

	except:
		raise
	finally:
		fichier.close()    
	try:
	    fichierMeilleurDatum = open(nomSemaine+".meilleurDatum", "w")
	    fichierMeilleurDatum.write(meilleurDatum)
	except:
		raise
	finally:
		fichierMeilleurDatum.close() 
		
    
def main(lst,lst_critere,nomSemaine,nombreStation):
    #===========================================================================
    # #PRINCIPE DE L'ALGORITHME
    #===========================================================================
    
              
    #parametrage du probleme
    #initialiser les traces de phéromones 
    #repeter jusquà nombre maximal d'iteration
        #pour chaque fourmi k dans nombre de fourmi:
            #tirage au sort du nombre de station et de l'objectif considéré 
            #tant que nombreStationsChoisies != nombreStationsTirees 
                #choisir un objet avec une probabilité elevee
                #Solution <- Solution U {nouvel objet}
            #stockage de ces solutions dans un fichier
        #Calcul des ecarts types
        #Depot de pheromone sur les meilleures fourmis de l'iteration pour chaque critere
        #Pareto - Depot general de pheromone pour les fourmis qui concilie au mieux la minimisation des criteres
        #Evaporation
        #Recalculer les probabilites en fonction des phéromones
        #plafonner les phéromones
    #generation du fichier de resultat   



    #===========================================================================
    # # Parametres:
    #===========================================================================
    #Nombre d'iterations
    if nombreStation<=12:
        max_it=14
    elif nombreStation==13:
        max_it=13
    elif nombreStation==14:
        max_it=10
    elif nombreStation==15:
        max_it=16
    elif nombreStation==16:
        max_it=32
    elif nombreStation==17:
        max_it=65
    elif nombreStation==18:
        max_it=130
    elif nombreStation==19:
        max_it=250
    elif nombreStation==20:
        max_it=510
    else: 
        max_it = 1000
        
    #Nombre de fourmis
    if nombreStation<=12:
        nb_fourmis=28
    elif nombreStation==13:
        nb_fourmis=60
    elif nombreStation==14:
        nb_fourmis=160
    else:
        nb_fourmis=200
    
    
    #------------------
    # Nombre de pheromones deposées pour les solutions non dominees du front de pareto
    k = 1 
    # Evaporation des phéromones à chaque tour
    p = k*0.1 
    #------------------
    # Contraintes
    NbStation_min=3
    #------------------
    # Bornes des probas
    proba_max=0.99
    proba_min=0.01
    pheromone_initiale=k
    
    #Initialisation des minima des criteres (pour la fonction cout)
    sigmaxMin=Inf
    sigmayMin=Inf
    sigmazMin=Inf
    #------------------
    # Initialisation des phéromones
    #
    # station en question
    #    pheromone
    #        sigmax
    #        sigmay
    #        sigmaz
    #    probabilite
    #        sigmax
    #        sigmay
    #        sigmaz
    sum_proba=0
    for objet in lst:
        for _ in lst_critere:
            objet['pheromone'][_]=pheromone_initiale 
            objet['probabilite'][_]=round(1/(len(lst)),5) # Equiprobabilite
            if _ == lst_critere[0]:
                sum_proba+=objet['probabilite'][_] # Equiprobabilite
           
    # Ajustement des probabilités pour que la somme soit égale à 1    
    if sum_proba != 1:
        diff=1-sum_proba
        for _ in lst_critere:
            lst[0]['probabilite'][_]=lst[0]['probabilite'][_]+diff
    #------------------
    # Creation des fourmis
    lst_fourmis=[]
    anciennes_non_dominees=[]
    for _ in range(nb_fourmis):
        lst_fourmis.append(Ant())
    #------------------  
    #Allocation des matrices (APPEL C)
    pgmn=zelib.lectureFichier(byref(c_char_p(nomSemaine.encode('utf8'))))
    

    #===========================================================================
    # # Debut ACO
    #===========================================================================
    
    
    for _ in range(max_it):
        #on ouvre le fichier datum
        fichier=open(nomSemaine+".datum","w")
        
        for antk in lst_fourmis:
            # On restart la fourmi
            antk.Objet.clear()
            # Tirage au sort de son nombre de station
            antk.NbStation=randint(NbStation_min,nombreStation+1)
            # Tirage au sort de son critere d'etude
            antk.Critere=lst_critere[randint(0,len(lst_critere))]
            #------------------  
            # Copie des stations à traiter
            lst_station_restante=lst.copy()
             
            # Condition pour les boucles. La condition devient fausse quand la fourmis
            # a traité ses antk.NbStation stations
            condition=True
            compte=0
            # Recuperation de la liste des probabilités pour utiliser la fonction choice
            lst_indice=[]
            lst_probabilite=[]
            indice=0
            for objet in lst:
                lst_indice.append(indice)
                lst_probabilite.append(objet['probabilite'][antk.Critere])
                indice+=1
            del(indice)
                
            while condition:
                # On sort de cette boucle dès que l'objet à ajouter n'a pas déjà été traité
                while True:
                    # Choisir une station avec la probabilité la plus forte
                    indice_objet_choisi=choice(lst_indice, p=lst_probabilite)
                    objet_suivant=lst[indice_objet_choisi]
                    if lst_station_restante[indice_objet_choisi]!=0:
                        break
                    
                # On ajoute la station à la fourmi
                antk.Objet.append(objet_suivant)
                del(objet_suivant)
                # On modifie la liste des stations restantes à traiter
                lst_station_restante[indice_objet_choisi]=0    
                
                compte+=1
                condition = compte!=antk.NbStation #Sera forcement faux a un moment donné car 3<=antk.NbStation<=nombreStation
            
            del(condition)
            del(compte)
            del(lst_indice)
            del(lst_probabilite)   
            del(lst_station_restante)
            
             
            #print("\n ___________ \n Fourmis : \n")
            #print("Critere choisi : "+antk.Critere)
            #print("Stations trouvées : "+str(len(antk.Objet))+"/ "+str(antk.NbStation) +" à calculer")
            #print("Liste des stations selectionnees :")
            
            
            #Creation du datum pour la fourmi
            datum=[0 for variable in range(nombreStation)]
            for station in antk.Objet:
                datum[int(station["station"])-1]=1
            #print(datum)
            
            #Conversion du datum en chaine de caracteres
            chDatum=""
            for selection in datum:
                chDatum+=str(selection)
            del(datum)
            antk.datum=chDatum
            chDatum+="\n"
            #print(chDatum)
                
            #On ecrit le datum dans le fichier datum
            fichier.write(chDatum)
            del(chDatum)
            
        
        #Fermeture du fichier de datum
        fichier.close()
        del(fichier)
        
        #Calcul des ecarts types à partir du datum
        zelib.main(pgmn,byref(c_char_p(nomSemaine.encode('utf8'))))
        
        #Recuperation du fichier de resultat généré par le C
        try:
            fichier=open(nomSemaine+".ecartType","r")
            numeroLigne=0
            for ligne in fichier:
                calculSigma=ligne.split(" ")
                sigmax=calculSigma[0]
                sigmay=calculSigma[1]
                sigmaz=calculSigma[2]
                lst_fourmis[numeroLigne].sigmax=float(sigmax)
                lst_fourmis[numeroLigne].sigmay=float(sigmay)
                lst_fourmis[numeroLigne].sigmaz=float(sigmaz)
                numeroLigne+=1
            del(calculSigma)
            del(numeroLigne)
            del(sigmax)
            del(sigmay)
            del(sigmaz)
        except:
            raise
        finally:
            fichier.close()
            del(fichier)
            
        #print("\nResultat pour chaque fourmis\n")
        for antk in lst_fourmis:
            #print(antk.Critere,antk.datum,antk.sigmax,antk.sigmay,antk.sigmaz)
            if antk.Critere==lst_critere[0]:
                sigmaxMin=min(sigmaxMin,antk.sigmax)
            elif antk.Critere==lst_critere[1]:
                sigmayMin=min(sigmayMin,antk.sigmay)
            elif antk.Critere==lst_critere[2]:
                sigmazMin=min(sigmazMin,antk.sigmaz)
        #print("\nSigma minima")
        #print(sigmaxMin,sigmayMin,sigmazMin)
        
        #Calcul du cout pour chaque fourmi en fonction de son critere et determination du meilleur cout
        meilleursCouts=[0 for variable in lst_critere]
        for antk in lst_fourmis:
            if antk.Critere==lst_critere[0]:
                antk.cout=1/(1+antk.sigmax-sigmaxMin)
                if antk.cout>meilleursCouts[0]:
                    meilleursCouts[0]=antk.cout
            elif antk.Critere==lst_critere[1]:
                antk.cout=1/(1+antk.sigmay-sigmayMin)
                if antk.cout>meilleursCouts[1]:
                    meilleursCouts[1]=antk.cout
            elif antk.Critere==lst_critere[2]:
                antk.cout=1/(1+antk.sigmaz-sigmazMin)
                if antk.cout>meilleursCouts[2]:
                    meilleursCouts[2]=antk.cout

        
        
        #Ajout de pheromone pour les meilleurs de l'iteration pour chaque critere
        for antk in lst_fourmis:
            if antk.Critere==lst_critere[0]:
                if antk.cout==meilleursCouts[0]:
                    for station in antk.Objet:
                        obj=lst[int(station["station"])-1]#numero de la station dans lst
                        obj['pheromone'][antk.Critere]+=antk.cout
            elif antk.Critere==lst_critere[1]:
                if antk.cout==meilleursCouts[1]:
                    for station in antk.Objet:
                        obj=lst[int(station["station"])-1]#numero de la station dans lst
                        obj['pheromone'][antk.Critere]+=antk.cout
            elif antk.Critere==lst_critere[2]:
                if antk.cout==meilleursCouts[2]:
                    for station in antk.Objet:
                        obj=lst[int(station["station"])-1]#numero de la station dans lst
                        obj['pheromone'][antk.Critere]+=antk.cout
        del(meilleursCouts)  
        del(obj)
    
        #=======================================================================
        # #Front de Pareto
        #=======================================================================
        # Pour determiner les tirages qui concilie au mieux la minimisation des trois criteres
        indice_fourmisA=0
        tableau_dominance=[[[] for ___ in lst_fourmis] for __ in lst_fourmis]
        for fourmisA in lst_fourmis:
            indice_fourmisB=0
            for fourmisB in lst_fourmis:
                flag= check_dominance(fourmisA, fourmisB, lst_critere,sigmaxMin,sigmayMin,sigmazMin)
                tableau_dominance[indice_fourmisA][indice_fourmisB]=flag
                indice_fourmisB+=1
            indice_fourmisA+=1
        del(indice_fourmisA)
        del(indice_fourmisB)
        del(flag)
        #
        #             fourmisA
        #fourmis B    | 0|-1|
        #             | 1| 0|  -> A est dominée par B
        #
        
        non_dominees=[]
        compteur=0
        for i in tableau_dominance:
            elimine=False
            for j in i:
                if j==-1:
                    elimine=True
            if elimine==False:
                non_dominees.append(compteur)
            compteur+=1
        #print("\nFourmis non dominees")
        #print(non_dominees)
        #print("\nTableau de dominance")
        #print_matrice(tableau_dominance)
        del(tableau_dominance)
        
        # Ajout des pheromones pour les fourmis dans la liste des non dominees
        for compteur in non_dominees:
            anciennes_non_dominees.append(copy.copy(lst_fourmis[compteur]))
            for station in lst_fourmis[compteur].Objet:
                obj=lst[int(station["station"])-1]#numero de la station dans lst
                for critere in lst_critere:                
                    obj['pheromone'][critere]+=1
                    
        del(obj)
        del(non_dominees)
                 
                 
        # On met à jour la liste des solutions non dominées en incorporant les nouvelles à toutes celles trouvées jusqu'à present
        indice_fourmisA=0
        tableau_dominance=[[[] for ___ in anciennes_non_dominees] for __ in anciennes_non_dominees]
        for fourmisA in anciennes_non_dominees:
            indice_fourmisB=0
            for fourmisB in anciennes_non_dominees:
                flag= check_dominance(fourmisA, fourmisB, lst_critere,sigmaxMin,sigmayMin,sigmazMin)
                tableau_dominance[indice_fourmisA][indice_fourmisB]=flag
                indice_fourmisB+=1
            indice_fourmisA+=1
        
        del(indice_fourmisA)
        del(indice_fourmisB)
        del(flag)
        
        
        non_dominees_archive=[]
        compteur=0
        for i in tableau_dominance:
            elimine=False
            for j in i:
                if j==-1:
                    elimine=True
            if elimine==False:
                non_dominees_archive.append(compteur)
            compteur+=1
        
        del(elimine)
        del(tableau_dominance)
        nouvelles_non_dominee=[]
        for compteur in non_dominees_archive:
            nouvelles_non_dominee.append(copy.copy(anciennes_non_dominees[compteur]))
            
        del(non_dominees_archive)
        del(compteur)
        
        anciennes_non_dominees.clear()
        anciennes_non_dominees+=nouvelles_non_dominee
        del(nouvelles_non_dominee)
        
                    
        #=======================================================================
        # # Evaporation des pheromones       
        #=======================================================================
        #print("\nPheromones de la premiere station avant evaporation")
        #print(lst[0]["pheromone"])
        for objet in lst:
            # Evaporation pour chaque critere
            for critere in lst_critere:
                objet['pheromone'][critere]=(1-p)*objet['pheromone'][critere]
        #print("Pheromones de la premiere station apres evaporation")
        #print(lst[0]["pheromone"])
        
        
        
        
        #=======================================================================
        # # Calcul des probabilites
        #=======================================================================
        # Initialisation
        somme_probabilite={}
        somme_probabilite_apres_bornage={}
        somme_probabilite_apres_normalisation={}
        proba_maximale={}
            
            
        for critere in lst_critere:
            somme_probabilite[critere]=0
            somme_probabilite_apres_bornage[critere]=0
            proba_maximale[critere]=0
            somme_probabilite_apres_normalisation[critere]=0
                
        #Calcul des probabilites avant bornage
        for objet in lst:
            for critere in lst_critere:
                objet['probabilite'][critere]=float(objet['pheromone'][critere])
                somme_probabilite[critere]=somme_probabilite[critere]+objet['probabilite'][critere]
             
        # Bornage des probabilites et calcul de la nouvelle somme
        for objet in lst:
            for critere in lst_critere:
                nouvelleProbabilite=round(objet['probabilite'][critere]/somme_probabilite[critere],5)
                if nouvelleProbabilite > proba_max:
                    objet['probabilite'][critere]=proba_max
                    
                elif nouvelleProbabilite < proba_min:
                    objet['probabilite'][critere]=proba_min
                else:
                    objet['probabilite'][critere]=nouvelleProbabilite
                somme_probabilite_apres_bornage[critere]=somme_probabilite_apres_bornage[critere]+objet['probabilite'][critere]
        del(nouvelleProbabilite)
        del(somme_probabilite)      
        
        # Normalisation des probabilites
        for objet in lst:
            for critere in lst_critere:
                objet['probabilite'][critere]=round(objet['probabilite'][critere]/somme_probabilite_apres_bornage[critere],5)
        del(somme_probabilite_apres_bornage)         
         
                 
        # Pour que la somme des probabilites fasse 1:
        for objet in lst:
            for critere in lst_critere:
                proba_maximale[critere]=max(objet['probabilite'][critere],proba_maximale[critere])
                somme_probabilite_apres_normalisation[critere]+=objet['probabilite'][critere]        
            
        # On corrige la meilleure probabilité     
        for critere in lst_critere:
            if somme_probabilite_apres_normalisation[critere]!=1.:
                diff=1-somme_probabilite_apres_normalisation[critere]
                for objet in lst:
                    if objet['probabilite'][critere]==proba_maximale[critere]:
                        objet['probabilite'][critere]=objet['probabilite'][critere]+diff
                        break
                del(diff)
        
        del(somme_probabilite_apres_normalisation)
        del(proba_maximale)
        #print("{} iterations calculées sur {} ({}%)".format(_+1,max_it,round((_+1)/max_it*100,2)))
    
    zelib.liberergmn(pgmn)
    """ 
    for compteur in non_dominees:
        print(str(lst_fourmis[compteur].datum)+" "+str(lst_fourmis[compteur].sigmax)+" "+str(lst_fourmis[compteur].sigmay)+" "+str(lst_fourmis[compteur].sigmaz))
    """      
    # Creation du fichier .solutionNonDominees
    fichierNonDominees=open(nomSemaine+".solutionNonDominees","w")
    for compteur in anciennes_non_dominees:
        fichierNonDominees.write(str(compteur.datum)+" "+str(compteur.sigmax)+" "+str(compteur.sigmay)+" "+str(compteur.sigmaz)+"\n")
    fichierNonDominees.close()
    return 0
    
            

    
    
    


 
if __name__=="__main__":
    if (sys.argv[0]=="mouvant.py" and len(sys.argv)==1):
        print("Please enter : python3 aco_lib.py file_para.dat file_tran.dat")
        
    elif (sys.argv[0]=="mouvant.py" and len(sys.argv)==3):
        fichierPara=sys.argv[1]
        fichierTran=sys.argv[2]
        nomSemaine=fichierPara.split("_para")[0]
        nombreStation=calculNombreStation(fichierTran)
        nomFichierStation=creerFichierStation(nomSemaine,nombreStation)
        
        listes=lst_objet(nomFichierStation)
        
        liste_objets=listes[0]
        liste_criteres=listes[1]
        
        k=main(liste_objets,liste_criteres,nomSemaine,nombreStation)
        statistiques(nomSemaine)
    else:
        print("Please enter : python3 aco_lib.py file_para.dat file_tran.dat")
        
      
