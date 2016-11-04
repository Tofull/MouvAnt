#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Import des modules
import numpy as np
from numpy.random import choice
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pylab import *
from time import sleep





def cree_matrice(taille_matrice,defaut=0):
    """
    Fonction qui retourne une matrice de taille taille_matrice[0]×taille_matrice[1] avec une valeur par defaut
    """
    matrice=[]
    for k in range(taille_matrice[0]):
        matrice.append([])
        for _ in range(taille_matrice[1]):
            matrice[k].append(defaut)
    return matrice





def print_matrice(matrice):
    """
    Fonction qui affiche une matrice dans la console
    """
    print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in matrice]))
    
    
    
    
    
def modifie_matrice(matrice,ligne,colonne,valeur):
    """
    Fonction qui modifie la valeur de la cellule de la matrice indicée à la ligne et colonne entrées par la valeur valeur 
    """
    matrice[ligne][colonne]=valeur
    
    
    
    
    
def affiche_matrice(matrice,**kwargs):
    """
    Fonction qui affiche une matrice dans une fenetre graphique (**kwargs permet de personnaliser l'affichage)
    """
    ax.matshow(matrice,**kwargs)





def affiche_trajet(lst,matrice,**kwargs):
    """
    Fonction qui permet de creer une matrice représentant un trajet à partir d'une matrice support et de la liste des coordonnees
    """
    ligne=len(matrice)
    colonne=len(matrice[0])
    new_matrice=cree_matrice([ligne,colonne])
    for ville in lst:
        modifie_matrice(new_matrice, ville[1], ville[0], 1)
    ax.matshow(new_matrice,**kwargs)
            
            
            
            
            
def recherche_valeur_matrice(matrice,valeur):
    """
    Permet de trouver l'indice et la colonne de la premiere occurence de valeur trouvée dans la matrice. Utile pour trouver des valeurs uniques, comme le point de départ et le point d'arrivée
    """
    indice_ligne=0
    for ligne in matrice:
        indice_colonne=0
        for _ in ligne:
            if matrice[indice_ligne][indice_colonne]==valeur:
                return [indice_ligne,indice_colonne]
            indice_colonne+=1
        indice_ligne+=1
        
        
        
        
        
def enleve_boucle(lst,valeur):
    """
    Fonction qui enleve une boucle dans une liste
    >>> enleve_boucle([1,2,3,4,5,3,7,8],3)
    >>> [1,2,3,7,8]
    
    """
    trouve1=-1
    trouve2=-1
    indice=0
    for elem in lst:
        if elem==valeur:
            if trouve1!=-1 and trouve2==-1: # Trouve l'indice de la deuxième occurence de valeur
                trouve2=indice
                break
            elif trouve1==-1: # Trouve l'indice de la premiere occurence de valeur
                trouve1=indice
        indice=indice+1
    lst = lst[:trouve1] + lst[trouve2 :] # Enleve la boucle
    return lst

    
    
    

class Ant(object):
    def __init__(self,**kwargs):
        '''
        Constructeur de la classe.
        '''
        self.Objet=[]
        self.Position=[0,0]
        self.Longueur=0
        self.x_min=0
        self.x_max=0
        self.y_min=0
        self.y_max=0
        
    def villesVoisines(self):
        """
        Fonction qui genere la liste des villes voisines à la fourmi
        """
        lst=[]
        if self.Position[1]-1>=self.y_min:
            villeGauche=[self.Position[0],self.Position[1]-1]
            lst.append(villeGauche)
        if self.Position[1]+1<=self.y_max:
            villeDroite=[self.Position[0],self.Position[1]+1]
            lst.append(villeDroite)
        if self.Position[0]-1>=self.x_min:
            villeHaut=[self.Position[0]-1,self.Position[1]]
            lst.append(villeHaut)
        if self.Position[0]+1<=self.x_max:
            villeBas=[self.Position[0]+1,self.Position[1]]
            lst.append(villeBas)
        return lst
    
    def enleve_boucle(self,valeur):
        """
        Fonction qui enleve une boucle dans la liste d'objet de la fourmi
        >>> print(self.Objet)
        >>> [1,2,3,4,5,3,7,8]
        >>> enleve_boucle(3)
        >>> [1,2,3,7,8]
        
        """
        trouve1=-1
        trouve2=-1
        indice=0
        for elem in self.Objet:
            if elem==valeur:
                if trouve1!=-1 and trouve2==-1:
                    trouve2=indice
                    break
                elif trouve1==-1:
                    trouve1=indice
            indice=indice+1 
        self.Objet=self.Objet[:trouve1] + self.Objet[trouve2:]
        return 0
    
    def fonctionCout(self):
        """
        Fonction qui calcule la distance parcourue par une fourmi
        """
        x0=self.Objet[0][0]
        y0=self.Objet[0][1]
        
        distance=0
        for k in range(len(self.Objet)-1):
            x=self.Objet[k+1][0]
            y=self.Objet[k+1][1]
            distance=distance+sqrt((x-x0)**2+(y-y0)**2)
            x0=self.Objet[k][0]
            y0=self.Objet[k][1]
        return distance
            
            
            
            
            
def main(matrice,obstacle=[]):
    #===========================================================================
    # PARAMETRES D'AFFICHAGE
    #===========================================================================
    plt.ion() # Active le mode interactif 
    CM = mpl.colors.ListedColormap(['red','gray', 'green','white','black']) # Map de couleur
        
    #===========================================================================
    # Paramètrage du problème
    #===========================================================================
    villeInitiale=recherche_valeur_matrice(matrice, -1)
    villeFinale=recherche_valeur_matrice(matrice, -2)
    # Bornes du problème
    x_min=0
    x_max=len(matrice[0])-1
    y_min=0
    y_max=len(matrice)-1
    # Paramètres de l'algorithme
    pheromone_initiale=5
    max_it = int(1.5*len(matrice))
    nb_fourmis = len(matrice)
    k = 1 # Nombre de pheromones deposées
    p = 0.1 # Evaporation des phéromones à chaque tour
    beta=1 # Importance du dépot de phéromones dans le calcul des probabilités
    # Bornes des probas
    proba_max=0.99
    proba_min=0.01
    
    # Création des matrices de phéromones et de probabilité
    pheromone=cree_matrice([len(matrice),len(matrice[0])],pheromone_initiale)
    probabilite=cree_matrice([len(matrice),len(matrice[0])],1/(len(matrice)*len(matrice[0])))
    
    # Creation des fourmis
    lst_fourmis=[]
    for _ in range(nb_fourmis):
        lst_fourmis.append(Ant())
    
    # Paramétrage des fourmis
    for fourmis in lst_fourmis:
        fourmis.Position=villeInitiale
        fourmis.x_min=x_min
        fourmis.x_max=x_max
        fourmis.y_min=y_min
        fourmis.y_max=y_max
        
    # Résultats
    resultat=[[[] for __ in range(nb_fourmis)] for _ in range(max_it)]
    best_L=-1
    best_Set=set()
    
    # Pour afficher les graphes dans une meme figure
    fig, ax = plt.subplots()
        
    # Début ACO
    for _ in range(max_it):
        index_fourmi=0
        for antk in lst_fourmis:
            # On restart la fourmi
            antk.Objet.clear()
            # Elles démarrent toutes à la ville initiale
            antk.Objet.append(villeInitiale)
            antk.Position=villeInitiale
             
            # Condition pour les boucles. La condition devient fausse quand la fourmi atteint la ville finale
            condition=True
            while condition:
                # Recuperation de la liste des probabilités des villes voisines (mises à jour à chaque boucle)
                lst_indice=[]
                lst_probabilite=[]
                lst=antk.villesVoisines()
                indice=0
                max_proba=0
                for objet in lst:
                    lst_indice.append(indice)
                    max_proba+=probabilite[objet[0]][objet[1]]
                    indice=indice+1
                    
                for objet in lst:
                    lst_probabilite.append(probabilite[objet[0]][objet[1]]/max_proba)
                    
                # On sort de cette boucle dès que l'objet à ajouter respecte les contraintes (ne sort pas de la grille et n'est pas un obstacle)
                while True:
                    indice_objet_choisi=choice(lst_indice, p=lst_probabilite)
                    x_objet_suivant=lst[indice_objet_choisi][0]
                    y_objet_suivant=lst[indice_objet_choisi][1]
                    objet_suivant=[x_objet_suivant,y_objet_suivant]
                    if x_objet_suivant>=0 and x_objet_suivant<taille_matrice[0] and y_objet_suivant>=0 and y_objet_suivant<taille_matrice[1] and (objet_suivant not in obstacle):
                        break
                    
                # On l'ajoute dans le tour de la fourmi 
                antk.Objet.append(objet_suivant)
                # On met à jour la position de la fourmi
                antk.Position=objet_suivant
                
                # On sort dès qu'on a atteint la ville finale
                condition=objet_suivant!=villeFinale
                
            # On elimine les boucles dans la liste des villes parcourues
            villes_parcourues=[]
            while len(villes_parcourues)!=len(antk.Objet):
                villes_parcourues=[]
                for ville in antk.Objet:
                    if tuple(ville) in villes_parcourues:
                        antk.enleve_boucle(ville)
                        villes_parcourues=[]
                    else:
                        villes_parcourues.append(tuple(ville))
            
            # Depot des phéromones dans les villes parcourues
            for objet in antk.Objet:
                pheromone[objet[0]][objet[1]]+=k
                
                
            # Calcul du cout de la fourmi et ajout au résultat
            resultat[_][index_fourmi].append(antk.fonctionCout()) 
            
            index_fourmi=index_fourmi+1
            
       
        #=======================================================================
        # Affichage des fourmis
        #
        #http://matplotlib.org/examples/color/colormaps_reference.html
        #cmap pour la couleur, alpha pour la transparence
        max_nb_ville=0
        for fourmis in lst_fourmis:
            max_nb_ville=max(max_nb_ville,len(fourmis.Objet))

        
        FOURMIS_X=np.array([[-1 for __ in range(len(lst_fourmis))] for _ in range(max_nb_ville+1)])
        FOURMIS_Y=np.array([[-1 for __ in range(len(lst_fourmis))] for _ in range(max_nb_ville+1)])
        indice=0
        for fourmis in lst_fourmis:
            num_objet=0
            for objet in fourmis.Objet:
                FOURMIS_X[num_objet,indice]=objet[0]
                FOURMIS_Y[num_objet,indice]=objet[1]
                num_objet=num_objet+1
            indice=indice+1
        
        tableau=np.array([[[0 for __ in range(len(matrice[0]))] for _ in range(len(matrice))] for prout in range(max_nb_ville+1)] )
        tableau[:,villeInitiale[1],villeInitiale[0]]=-2
        tableau[:,villeFinale[1],villeFinale[0]]=-5
        if obstacle!=[]:
            for obst in obstacle:
                tableau[:,obst[1],obst[0]]=-3
        
        
        tableauP=np.array([[[0 if _==0 else tableauP[max_ville_prec,pm,__]*0.9 for __ in range(len(matrice[0]))] for pm in range(len(matrice))] for prout in range(max_nb_ville+1)] )
        max_ville_prec=max_nb_ville
        for k in range(max_nb_ville):
            for d in range(len(FOURMIS_X[k])):
                if FOURMIS_Y[k,d]==-1 or FOURMIS_X[k,d]==-1:
                    pass
                else:
                    tableau[k,FOURMIS_Y[k,d],FOURMIS_X[k,d]]=2
                    tableauP[k+1:,FOURMIS_Y[k,d],FOURMIS_X[k,d]]=tableauP[k+1,FOURMIS_Y[k,d],FOURMIS_X[k,d]]+25
       
        
        image2 = plt.imshow(tableauP[1, :, :],cmap='PuBu',alpha=1)
        image = plt.imshow(tableau[0, :, :],cmap=CM,alpha=0.4, interpolation='none')
        
        for k in np.arange(max_nb_ville+1):
            image.set_data(tableau[k, :, :])
            image2.set_data(tableauP[k, :, :])
            plt.draw()
            sleep(0.005)
        
        #
        # 
        #=======================================================================
                
        
        
        # Evaporation des pheromones        
        for ligne in range(len(pheromone)):
            for colonne in range(len(pheromone[0])):
                modifie_matrice(pheromone, ligne, colonne, (1-p)*pheromone[ligne][colonne])

        
        # Calcul des probabilites
        somme_probabilite=0    
        somme_probabilite_apres_bornage=0

        for ligne in range(len(probabilite)):
            for colonne in range(len(probabilite[0])):
                probabilite[ligne][colonne]=pheromone[ligne][colonne]**beta
                somme_probabilite=somme_probabilite+probabilite[ligne][colonne]
                
        # Bornage des probabilites
        for ligne in range(len(probabilite)):
            for colonne in range(len(probabilite[0])):
                nouvelleProbabilite=probabilite[ligne][colonne]/somme_probabilite
                if nouvelleProbabilite > proba_max:
                    probabilite[ligne][colonne]=proba_max
                    
                elif nouvelleProbabilite < proba_min:
                    probabilite[ligne][colonne]=proba_min
                else:
                    probabilite[ligne][colonne]=nouvelleProbabilite
                somme_probabilite_apres_bornage=somme_probabilite_apres_bornage+probabilite[ligne][colonne]
    
        # Normalisation des probabilites après bornage
        for ligne in range(len(probabilite)):
            for colonne in range(len(probabilite[0])):
                probabilite[ligne][colonne]=probabilite[ligne][colonne]/somme_probabilite_apres_bornage
        
        # Recuperation des meilleurs résultats
        if best_L==-1:
            best_L=min(resultat[_])
            for fourmis in lst_fourmis:
                # On recherche les fourmis les plus performantes de la premiere iteration
                if fourmis.fonctionCout()==best_L:
                    trajet=[]
                    for obj in fourmis.Objet:
                        trajet.append(tuple(obj))
                    best_Set.add(tuple(trajet))
        else:
            best_L_precedent=best_L
            best_L=min(min(resultat[_]),best_L_precedent)
            # Si best_L a changé, alors une fourmis a trouvé un meilleure chemin, donc on efface le meilleur chemin enregistré
            if best_L_precedent!=best_L: 
                best_Set=set()
                
            for fourmis in lst_fourmis:
                # On recherche les fourmis les plus performantes de la kieme iteration
                if fourmis.fonctionCout()==best_L:
                    trajet=[]
                    for obj in fourmis.Objet:
                        trajet.append(tuple(obj))
                    best_Set.add(tuple(trajet))
                    
        print("{} iterations calculées sur {} ({}%)".format(_+1,max_it,round((_+1)/max_it*100,2)))   
    
   
    return best_Set



    
    
if __name__=="__main__":
    # Configuration du probleme
    taille_matrice=(10,10)
    matrice=cree_matrice(taille_matrice)
    # Placement de la ville initiale
    modifie_matrice(matrice, 1, 1, -1)
    # Placement de la ville finale
    modifie_matrice(matrice, taille_matrice[0]-2, taille_matrice[0]-2, -2)

    # Obstacle
    obstacle=[[5,1],[5,2],[5,3],[5,4],[5,5],[5,6],[5,8],[4,8],[4,6],[7,8],[7,7],[7,9],[6,1],[7,1],[7,3],[9,3],[8,3],[6,5],[7,5]]
    #obstacle=[]
    
    # Resultat
    resultat=[]
    meilleur_trajet=main(matrice,obstacle)
    for subset in meilleur_trajet:
        lst_ville=[]
        for ville in subset:
            lst_ville.append(list(ville))
        resultat.append(lst_ville)
        print(lst_ville)
        
    # Pour afficher les graphes dans une meme figure
    fig, ax = plt.subplots()    
    # Affichage d'un des meilleurs trajets trouvés
    affiche_trajet(obstacle,matrice,alpha=0.5,cmap='hot')
    affiche_trajet(resultat[0],matrice,alpha=0.5,cmap='summer')
    plt.ioff()
    
    # Affichage global
    plt.show()
    
    
    
    
    
    
    
