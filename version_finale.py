import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation, PillowWriter
import random as rng

class Simulation :
    
    def __init__(self, liste_de_poissons, liste_de_predateurs, liste_obstacles, N, dt, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, rayon_cohesion, rayon_separation, rayon_alignement, rayon_predation, rayon_proies):
        self.liste_de_poissons = liste_de_poissons
        self.liste_de_predateurs = liste_de_predateurs
        self.liste_obstacles = liste_obstacles
        self.N = N
        self.dt = dt
        
        self.alpha_cohesion = alpha_cohesion
        self.alpha_separation = alpha_separation
        self.alpha_alignement = alpha_alignement
        self.a_rng = a_rng
        self.rayon_cohesion = rayon_cohesion
        self.rayon_separation = rayon_separation
        self.rayon_alignement = rayon_alignement
        self.rayon_predation = rayon_predation
        self.rayon_proies = rayon_proies
        
        self.initialiser_matrices_boids()
        
    def initialiser_matrices_boids(self):
        """
        Cree la matrice des poissons
        """
        for poisson in self.liste_de_poissons:
            poisson.ajouter_simulation(self)
        for predateur in self.liste_de_predateurs:
            predateur.ajouter_simulation(self)
        
    def voisin_le_plus_proche(self, p, i): 			
        """
        Renvoie le boid le plus proche
        Args:
            i: indice temporel concerné de la matrice positions
            p: indice du poisson concerné dans liste_de_poissons

        Returns:
            poisson_voisin : le poisson le plus proche 
        """
        poisson_1 = self.liste_de_poissons[p]		
        poisson_voisin = None
        distance_min = np.inf
        for n in range(len(self.liste_de_poissons)):
            if n != p :
                poisson_2 = self.liste_de_poissons[n]
                if poisson_2.vivant[i] == True:
                    distance = poisson_1.distance(poisson_2, i)
                    if distance < distance_min :
                        distance_min = distance
                        poisson_voisin = poisson_2
        return poisson_voisin
                    
    def centre_masse_voisins(self, poisson, i):
        """
        Renvoie le centre de masse des poissons
        Args:
            poisson : objet poisson duquel on calcule le centre de masse des voisins
            i: indice temporel concerné de la matrice positions

        Returns:
            point_centre : les coordonnees du centre de masse des poissons
        """
        liste_voisins = poisson.poissons_dans_rayon_cohesion
        point_centre = np.zeros((2))
        nb_voisins = len(liste_voisins)
        if nb_voisins != 0:
            for voisin in liste_voisins :
                point_centre += voisin.positions[i]
            point_centre = point_centre / nb_voisins
            return 1, point_centre
        else:
            return None, point_centre
        
    def vitesse_voisins(self, poisson, i): # i est l'indice temporel concerné de la matrice positions
        """
        Renvoie la vitesse du banc de poissons
        Args:
            i: indice temporel concerné de la matrice positions

        Returns:
            vitesse_moy : les coordonnées du vecteur vitesse du banc de poissons
        """
        vitesse_moy = np.zeros((2))
        liste_voisins = poisson.poissons_dans_rayon_alignement
        nb_voisins = len(liste_voisins)

        if nb_voisins != 0:
            for voisin in liste_voisins :
                vitesse_moy += voisin.vitesses[i]
            vitesse_moy = vitesse_moy / nb_voisins
            return 1, vitesse_moy
        else:
            
            return None, vitesse_moy

    def voisin_visible(self,boid,voisin,i):
        visible = True
        j = 0
        while visible and j < len(self.liste_obstacles): #while car non visible si visible selon un obstacle et non visible selon un autre
            obstacle = self.liste_obstacles[j]
            x1, y1, x2, y2 = obstacle.liste_limites #coordonées du point 1 et 2 de l'obstacle
                
            #on calcule les équations de droite (si elles existent) de: l'obstacle, la droite coupant le poisson et le point 1, la droite coupant le poisson et le point 2
            #y = ax + b
            
            if x1!=boid.positions[i,0]:
                a1 = (y1-boid.positions[i,1])/(x1-boid.positions[i,0]); 
                b1 = y1 - a1*x1
                cond1 = (voisin.positions[i,1] > a1*voisin.positions[i,0] + b1) == (y2>a1*x2 + b1)
            else:
                cond1 = (voisin.positions[i,0] > x1) == (x2>x1)
                
            if x2!=boid.positions[i,0]:
                a2 = (y2-boid.positions[i,1])/(x2-boid.positions[i,0])
                b2 = y2 - a2*x2
                cond2 = (voisin.positions[i,1] > a2*voisin.positions[i,0] + b2) == (y1>a2*x1 + b2)
            else:     
                cond2 = (voisin.positions[i,0] > x2) == (x1>x2)     
                   
            if x2!=x1:
                a_ob = (y2-y1)/(x2-x1)
                b_ob = y1 - a_ob*x1
                cond3 = (boid.positions[i,1] > a_ob*boid.positions[i,0] + b_ob) != (voisin.positions[i,1] > a_ob*voisin.positions[i,0] + b_ob) 
            else:
                cond3 = (boid.positions[i,0] > x1) != (voisin.positions[i,0] > x1)       
            
            #Si un poisson n'est pas visible, c'est qu'il est dans l'intersection de l'espace de l'autre côté de l'obstacle et de celui "entre" les droites 1 et 2
            if cond1 and cond2 and cond3:
                visible = False
            j += 1
        return visible
        
    def mise_a_jour_liste_poissons(self, poisson, i, p, angle_de_vision = 150):

        poisson.poissons_dans_rayon_cohesion   = []
        poisson.poissons_dans_rayon_alignement = []
        poisson.poissons_dans_rayon_separation = []
        poisson.predateur_dans_vision = []

        dico_voisin = {'cohesion' : [],
                       'separation' : [],
                       'alignement' : []}

        dico_rayons = {'cohesion' : self.rayon_cohesion,
                       'separation' : self.rayon_separation,
                       'alignement' : self.rayon_alignement}
        
        liste_rayon_ordonnes = sorted(dico_rayons.items(), key=lambda x: x[1])
        
        cos_angle_vision = np.cos(np.deg2rad(angle_de_vision))

        for i_poisson in range (0, len(self.liste_de_poissons)):

            if i_poisson != p :
                voisin = self.liste_de_poissons[i_poisson]
                distance_poisson_voisin_i = poisson.distance(voisin, i)

            
                if distance_poisson_voisin_i < max(self.rayon_alignement, self.rayon_cohesion, self.rayon_separation):
                    vecteur_vitesse = poisson.vitesses[i, :]
                    vecteur_poisson_voisin = voisin.positions[i, :] - poisson.positions[i, :]
                    norme_vitesse = np.linalg.norm(vecteur_vitesse)
            
                    cos_angle = np.dot(vecteur_poisson_voisin, vecteur_vitesse) / (norme_vitesse *distance_poisson_voisin_i)
                    
                    if cos_angle > cos_angle_vision :

                        if distance_poisson_voisin_i < self.rayon_separation:

                            poisson.poissons_dans_rayon_separation.append(voisin)
                        if distance_poisson_voisin_i > self.rayon_separation and distance_poisson_voisin_i < self.rayon_alignement:

                            poisson.poissons_dans_rayon_alignement.append(voisin)
                        if distance_poisson_voisin_i > self.rayon_alignement and distance_poisson_voisin_i < self.rayon_cohesion:

                            poisson.poissons_dans_rayon_cohesion.append(voisin)
                            
        poisson.poissons_dans_rayon_cohesion   = dico_voisin['cohesion']
        poisson.poissons_dans_rayon_alignement = dico_voisin['alignement']
        poisson.poissons_dans_rayon_separation = dico_voisin['separation']
            
        for predateur in self.liste_de_predateurs:
            distance_predateur = poisson.distance(predateur, i)
            
            if distance_predateur < self.rayon_predation : 
                vecteur_vitesse = poisson.vitesses[i, :]
                vecteur_predateur = predateur.positions[i, :] - poisson.positions[i, :]
                norme_vitesse = np.linalg.norm(vecteur_vitesse)
                cos_angle = np.dot(vecteur_predateur, vecteur_vitesse) / (norme_vitesse *distance_predateur)
                                            
                if cos_angle > cos_angle_vision and self.voisin_visible(poisson,predateur,i):
                    poisson.predateur_dans_vision.append(predateur)
      
    def mise_a_jour_liste_predateurs(self, predateur, i, p, angle_de_vision = 160):
        predateur.poisson_le_plus_proche = None
        predateur.predateur_dans_rayon_separation = []
                
        cos_angle_vision = np.cos(np.deg2rad(angle_de_vision))
        
        chasse = False
        liste_proies_en_visu = []
        liste_distance_proies = []
        
        for i_predateur in range(0, len(self.liste_de_predateurs)):
            if i_predateur != p :
                voisin = self.liste_de_predateurs[i_predateur]
                distance_voisin = predateur.distance(voisin, i)
                
                if distance_voisin < (1.5*self.rayon_separation) : 
                    vecteur_vitesse = predateur.vitesses[i, :]
                    vecteur_voisin = voisin.positions[i, :] - predateur.positions[i, :]
                    norme_vitesse = np.linalg.norm(vecteur_vitesse)
                    cos_angle = np.dot(vecteur_voisin, vecteur_vitesse) / (norme_vitesse *distance_voisin)
                            
                    if cos_angle > cos_angle_vision and self.voisin_visible(predateur,voisin,i): 
                        predateur.predateur_dans_rayon_separation.append(voisin)
        
        for proie in self.liste_de_poissons : 
            distance_proie = predateur.distance(proie, i)
            if proie.vivant[i] == True :
                if distance_proie < self.rayon_proies :
                    vecteur_vitesse = predateur.vitesses[i, :]
                    vecteur_proie = proie.positions[i, :] - predateur.positions[i, :]
                    norme_vitesse = np.linalg.norm(vecteur_vitesse)
                    cos_angle = np.dot(vecteur_proie, vecteur_vitesse) / (norme_vitesse *distance_proie)
                                         
                    if cos_angle > cos_angle_vision and self.voisin_visible(predateur,proie,i):
                        if chasse != True:
                            if distance_proie < 25 : 
                                chasse == True
                            liste_proies_en_visu.append(proie)
                            liste_distance_proies.append(distance_proie)
    
        if np.size(liste_distance_proies) > 0:
            minimum = liste_distance_proies[0]
            indice_min = 0 
            for i in range(1,len(liste_distance_proies)):
                if liste_distance_proies[i]<minimum:
                    minimum = liste_distance_proies[i]
                    indice_min = i
            predateur.poisson_le_plus_proche = liste_proies_en_visu[indice_min]                

    def composante_acceleration_separation(self, poisson, i):
        res = np.zeros((2))
        liste_voisins = poisson.poissons_dans_rayon_separation
        compteur = len(liste_voisins)

        if compteur == 0:
            return res
        else:
            for voisin in liste_voisins:
                position_voisin = voisin.positions[i, :] 
                position_poisson = poisson.positions[i, :]
                distance_poisson_voisin = poisson.distance(voisin, i)
                res = res + self.alpha_separation * (position_poisson - position_voisin) / distance_poisson_voisin**2
            return res/compteur
        
    def composante_acceleration_cohesion(self, poisson, i):
        res = np.zeros((2))
        valeur, centre_masse_voisins = self.centre_masse_voisins(poisson, i)   ## centre de masse des voisins dans le rayon de vision de cohesion du poisson

        if valeur == None:
            return res
        else:
            res = self.alpha_cohesion * (centre_masse_voisins - poisson.positions[i, :])
            return res
        
    def composante_acceleration_alignement(self, poisson, i):
        res = np.zeros((2))
        valeur, vitesse_voisins = self.vitesse_voisins(poisson, i)   ## vecteur vitesse des voisins

        if valeur == None:
            return res
        else:
            res = self.alpha_alignement * ( vitesse_voisins - poisson.vitesses[i, :])
            return res
        
    def composante_acceleration_predation(self, poisson,i):
        res = np.zeros(2)
        liste_predateurs_vus = poisson.predateur_dans_vision
        compteur = len(liste_predateurs_vus)
        if compteur == 0:
            return res
        else:
            for predateur in liste_predateurs_vus:
                position_predateur = predateur.positions[i, :] 
                position_poisson = poisson.positions[i, :]
                distance_poisson_predateur = poisson.distance(predateur, i)
                res = res + 40*self.alpha_separation *(position_poisson - position_predateur)/distance_poisson_predateur**2
        return res/compteur
            
    def composante_acceleration_proie(self, predateur,i):
        res = np.zeros(2)
        proie =predateur.poisson_le_plus_proche
        if proie != None:
            position_predateur = predateur.positions[i, :]
            position_proie = proie.positions[i, :]
            vecteur_direction = position_proie - position_predateur
            vec_dir =vecteur_direction/np.linalg.norm(vecteur_direction)
            res = 150*self.alpha_alignement * vec_dir
        return res
    
    def composante_acceleration_separation_predateurs(self, predateur, i):
        res = np.zeros((2))
        liste_voisins = predateur.predateur_dans_rayon_separation
        compteur = len(liste_voisins)
        if compteur == 0:
            return res
        else:
            for voisin in liste_voisins:
                position_voisin = voisin.positions[i, :] 
                position_predateur = predateur.positions[i, :]
                distance_predateur_voisin = predateur.distance(voisin, i)
                res = res + 15*self.alpha_separation * (position_predateur - position_voisin) / distance_predateur_voisin**2
            return res/compteur
    
    def composante_acceleration_obstacles(self, poisson, i):
        res = np.zeros((2))
        for obstacle in self.liste_obstacles:
            vecteur_normal, projete_tan = poisson.infos_distance(obstacle, i)
            projete_normal = np.linalg.norm(vecteur_normal)
            if projete_normal <= obstacle.distance_repulsion and abs(projete_tan)<obstacle.longueur/2 and np.linalg.norm(vecteur_normal) != 0:
                res += obstacle.coefficient_repulsion*vecteur_normal/(np.linalg.norm(vecteur_normal))**2
        return res 
        
    def test_manger(self, i):
        """
        Regarde si un poisson et un prédateur sont superposés et si c'est le cas, supprime le poisson
        Entrée :
            i : l'indice temporel de la simulation
        """
        for predateur in self.liste_de_predateurs:
            if len(self.liste_de_poissons) != 0:
                for poisson in self.liste_de_poissons:
                    if poisson.vivant[i] == True :
                        if abs(predateur.distance(poisson,i)) < 25:
                            #print("je mange")
                            poisson.vivant[i] = False
                         
    
    
    def calcul_parametre_ordre(self, indice):
        """
        Calcul la valeur du paramètre d'ordre
        Args:
            indice : l'indice temporel de la simulation
        Returns:
            parametre_ordre_simulation : le parametre d'ordre de la simulation à un temps donné

        """
        
        vitesse_totale = np.array([0.,0.])
        vitesse_moyenne = 0
        for element in self.liste_de_poissons: 
            if element.vivant[indice] == True:
                vitesse_totale += np.array(element.vitesses[indice])
                vitesse_moyenne += np.linalg.norm(element.vitesses[indice])
                parametre_ordre_simulation = np.linalg.norm(vitesse_totale)/vitesse_moyenne
        return parametre_ordre_simulation
    
    def moyennage_parametre_ordre(self, indice, niter):
        sum = 0
        for i in range(indice, indice+niter):
            sum += self.calcul_parametre_ordre(i)
        return sum/niter

    def calcul_tableaux(self):
        for i in range(0, self.N):
 
            for p in range(0, len(self.liste_de_poissons)):
                poisson = self.liste_de_poissons[p]
                
                if poisson.vivant[i] == True:
                    self.mise_a_jour_liste_poissons(poisson, i, p)
                    
    
                    a_cohesion = self.composante_acceleration_cohesion(poisson, i)
    
                    a_separation = self.composante_acceleration_separation(poisson, i)
    
                    a_alignement = self.composante_acceleration_alignement(poisson, i)
                    
                    a_predation = self.composante_acceleration_predation(poisson, i)
                    
                    a_obstacle = self.composante_acceleration_obstacles(poisson, i)
                    
                    poisson.angle_accel_alea = poisson.angle_accel_alea + rng.uniform(-0.5, 0.5)   ## le vecteur vitesse aleatoire comprend des valeurs allant de -0,4 x a_rng à 0,6 x a_rng
                
                    a_aleatoire = np.array([self.a_rng * np.cos(poisson.angle_accel_alea), self.a_rng * np.sin(poisson.angle_accel_alea)])
    
                    poisson.accelerations[i+1, :] = a_cohesion + a_separation + a_alignement + a_predation + a_obstacle + a_aleatoire
    
    
                    
                    ## Calcul du vecteur vitesse au rang n+1
                    vecteur_vitesse = poisson.vitesses[i, :] + self.dt * poisson.accelerations[i+1, :]
    
                    ## Ajustement du vecteur vitesse pour avoir une vitesse inférieur à v_max
                    norme_vitesse = np.linalg.norm(vecteur_vitesse)
                    if norme_vitesse > poisson.v_max :
                        poisson.vitesses[i+1, :] = (vecteur_vitesse * poisson.v_max) / norme_vitesse
                    else:
                        poisson.vitesses[i+1, :] = vecteur_vitesse
    
                    ## Calcul du vecteur position au rang n+1
                    poisson.positions[i+1, :] = poisson.positions[i, :] + self.dt * poisson.vitesses[i+1, :]
                    self.test_manger(i)
                    poisson.vivant[i+1] = poisson.vivant[i]
                else :
                    poisson.vitesses[i+1, :] = poisson.vitesses[i, :]
                    poisson.positions[i+1, :] = poisson.positions[i, :]
                
            
       
            for p in range(0, len(self.liste_de_predateurs)):
                predateur = self.liste_de_predateurs[p]
            
                self.mise_a_jour_liste_predateurs(predateur, i, p)
            
                a_proies = self.composante_acceleration_proie(predateur, i)
                a_separation_predateurs = self.composante_acceleration_separation_predateurs(predateur, i)
                a_obstacle = self.composante_acceleration_obstacles(predateur, i)
                
                predateur.angle_accel_alea = predateur.angle_accel_alea + rng.uniform(-0.5, 0.5)   ## le vecteur vitesse aleatoire comprend des valeurs allant de -0,4 x a_rng à 0,6 x a_rng
               
                a_aleatoire = np.array([self.a_rng * np.cos(predateur.angle_accel_alea), self.a_rng * np.sin(predateur.angle_accel_alea)])

                predateur.accelerations[i+1, :] = a_proies + a_separation_predateurs + a_obstacle + a_aleatoire


            
                ## Calcul du vecteur vitesse au rang n+1
                vecteur_vitesse = predateur.vitesses[i, :] + self.dt * predateur.accelerations[i+1, :]
                    
                ## Ajustement du vecteur vitesse pour avoir une vitesse inférieur à v_max
                norme_vitesse = np.linalg.norm(vecteur_vitesse)
                if norme_vitesse > predateur.v_max :
                    predateur.vitesses[i+1, :] = (vecteur_vitesse * predateur.v_max) / norme_vitesse
                else:
                    predateur.vitesses[i+1, :] = vecteur_vitesse

                ## Calcul du vecteur position au rang n+1
                predateur.positions[i+1, :] = predateur.positions[i, :] + self.dt * predateur.vitesses[i+1, :]
                
   

class GUI:
    def __init__(self, simulation, vitesse_lecture = 1.0, coord_lim = 300, suivi = True, gif = False):
        self.simulation = simulation
            
        self.coord_lim = coord_lim #utile pour update_suivi(frame)
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(-coord_lim, coord_lim)
        self.ax.set_ylim(-coord_lim, coord_lim)
        self.ax.set_aspect('equal')
    
        self.obstacles = []
        self.init_obstacles()
        self.triangles = []
        self.init_triangles()

        if suivi == True:
            self.ax.grid()
            self.ani = FuncAnimation(self.fig, self.update_suivi, frames=np.arange(0,self.simulation.N-1,1), interval=self.simulation.dt*1000/vitesse_lecture, blit=False)
        else:
            self.ani = FuncAnimation(self.fig, self.update_non_suivi, frames=np.arange(0,self.simulation.N-1,1), interval=self.simulation.dt*1000/vitesse_lecture, blit=True) 
        if gif:
            writer = PillowWriter(fps=15,
                                metadata=dict(artist='Maxence'),
                                bitrate=1800)
            self.ani.save('test.gif', writer=writer)
        else:
            plt.show()
    
    def init_triangles(self):
        """
        Cree les triangles des boids
        
    """
        for poisson in self.simulation.liste_de_poissons:
            pos, vit = poisson.position_initiale, poisson.vitesse_initiale
            coords = self.coords_triangle(pos, vit)
            triangle = Polygon(coords, closed=True, color='skyblue')
            self.ax.add_patch(triangle)
            self.triangles.append(triangle)
        for predateur in self.simulation.liste_de_predateurs:
            pos, vit = predateur.position_initiale, predateur.vitesse_initiale
            coords = self.coords_triangle(pos, vit)
            triangle = Polygon(coords, closed=True, color='red')
            self.ax.add_patch(triangle)
            self.triangles.append(triangle)
    
    def init_obstacles(self):
        """
        crée les lignes représentant les obstacles
        """
        for obstacle in self.simulation.liste_obstacles:
            x1, y1, x2, y2 = obstacle.liste_limites
            ligne = Line2D([x1,x2],[y1,y2], color='crimson')
            self.ax.add_line(ligne)
            self.obstacles.append(ligne)
    
    def update_suivi(self,frame):
        """
        Met  à jour les triangles des boids
    Args:
        frame : frame actuelle de l'animation'

    Returns:
        self.triangles : les coordonnées des triangles représentant les boids
    """
        for i in range(len(self.simulation.liste_de_poissons)):
            pos,vit = self.simulation.liste_de_poissons[i].positions[frame], self.simulation.liste_de_poissons[i].vitesses[frame]
            new_coords = self.coords_triangle(pos, vit)
            self.triangles[i].set_xy(new_coords)
            self.triangles[i].set_fill(self.simulation.liste_de_poissons[i].vivant[frame])
        for j in range(len(self.simulation.liste_de_predateurs)):
            pos,vit = self.simulation.liste_de_predateurs[j].positions[frame], self.simulation.liste_de_predateurs[j].vitesses[frame]
            new_coords = self.coords_triangle(pos, vit)
            self.triangles[j+len(self.simulation.liste_de_poissons)].set_xy(new_coords)
        if len(self.simulation.liste_de_predateurs) != 0:
            self.ax.set_xlim(-self.coord_lim+self.simulation.liste_de_predateurs[0].positions[frame][0], self.coord_lim+self.simulation.liste_de_predateurs[0].positions[frame][0])
            self.ax.set_ylim(-self.coord_lim+self.simulation.liste_de_predateurs[0].positions[frame][1], self.coord_lim+self.simulation.liste_de_predateurs[0].positions[frame][1])
        else : 
            self.ax.set_xlim(-self.coord_lim+self.simulation.liste_de_poissons[0].positions[frame][0], self.coord_lim+self.simulation.liste_de_poissons[0].positions[frame][0])
            self.ax.set_ylim(-self.coord_lim+self.simulation.liste_de_poissons[0].positions[frame][1], self.coord_lim+self.simulation.liste_de_poissons[0].positions[frame][1])
        return self.triangles
    
    def update_non_suivi(self,frame):
        """
        Met  à jour les triangles des boids
    Args:
        frame : frame actuelle de l'animation'

    Returns:
        self.triangles : les coordonnées des triangles représentant les boids
    """
        for i in range(len(self.simulation.liste_de_poissons)):
            pos,vit = self.simulation.liste_de_poissons[i].positions[frame], self.simulation.liste_de_poissons[i].vitesses[frame]
            new_coords = self.coords_triangle(pos, vit)
            self.triangles[i].set_xy(new_coords)
            self.triangles[i].set_fill(self.simulation.liste_de_poissons[i].vivant[frame])
        for j in range(len(self.simulation.liste_de_predateurs)):
            pos,vit = self.simulation.liste_de_predateurs[j].positions[frame], self.simulation.liste_de_predateurs[j].vitesses[frame]
            new_coords = self.coords_triangle(pos, vit)
            self.triangles[j+len(self.simulation.liste_de_poissons)].set_xy(new_coords)
        return self.triangles
    
    
    def coords_triangle(self,pos,vit,size=0.5):
        """
        Calcule les coordonnées du triangle représentant le boid et l'aligne avec sa vitesse'
    Args:
        pos (tuple): tuple contenant x,y les coordonnées du boid
        vit (tuple): tuple contenant vx,vy correspondant au vecteur vitesse du boid
        size (float, optional): facteur de taille du boid

    Returns:
        triangle : coordonnées des faces d'un triangle représentant le boid
    """
        x, y = pos
        xvit, yvit = vit
        hauteur_tetraedre = 50*size             #hauteur du tetraedre
        Largeur_base = 1/2*hauteur_tetraedre      #largeur d'un coté de la base

        #si vitesse nulle attention à ne pas diviser par une norme nulle -> direction par défaut = (1,0,0)
        if xvit==0 and yvit==0 : 
            xvit = 1
            yvit = 0
        norme_vit = np.sqrt(xvit**2+yvit**2)

        #on cherche un vecteur orthonormale avec la direction du boid:
        xa, ya = xvit/norme_vit, yvit/norme_vit                  #selon la direction/vitesse
        xb, yb = -ya, xa                                         #vecteur orthogonal choisi arbitrairement

        
        #on calcule les points limites de la base et celui de l'apex du tetraedre: 
        pointe_haut = [x + 2/3*hauteur_tetraedre*xa,y+2/3*hauteur_tetraedre*ya ]
        pointe_gauche = [x - 1/3*hauteur_tetraedre*xa -Largeur_base/2*xb,y- 1/3*hauteur_tetraedre*ya - Largeur_base/2*yb]
        pointe_droite = [x - 1/3*hauteur_tetraedre*xa +Largeur_base/2*xb,y- 1/3*hauteur_tetraedre*ya + Largeur_base/2*yb]
        

        #on renseigne les faces du tetraedre
        triangle =[pointe_haut, pointe_gauche, pointe_droite]
        
        return triangle

class Boid :
    
    def __init__(self, position_initiale, vitesse_initiale,v_max):
        self.position_initiale = position_initiale
        self.vitesse_initiale = vitesse_initiale
        self.angle_accel_alea = rng.uniform(0, 2 * np.pi)
        self.v_max = v_max
        
    def ajouter_simulation(self, simulation):
        self.simulation = simulation
        self.initialisation_matrices()
        
    def initialisation_matrices(self):
        self.positions = np.zeros((self.simulation.N+1, 2))
        self.vitesses = np.zeros((self.simulation.N+1, 2))
        self.accelerations = np.zeros((self.simulation.N+1, 2))
        self.vivant = np.zeros(self.simulation.N+1)
        
        self.vivant[0] = True
        self.positions[0] = self.position_initiale
        self.vitesses[0] = self.vitesse_initiale

    def __str__(self):
        text = "Positions :\n" + str(self.positions) + "\nVitesses :\n" + str(self.vitesses) + "\nAccelerations :\n" + str(self.accelerations)
        return text

    def distance(self, boid, i): 	#distance avec un boid à l'indice i
        diff_position = self.positions[i] - boid.positions[i]
        norme = np.linalg.norm(diff_position)
        return norme
    
    def infos_distance(self, obstacle, i): 	#distance avec un boid à l'indice i
        vecteur_arrete = np.array([obstacle.liste_limites[2]-obstacle.liste_limites[0],obstacle.liste_limites[3]-obstacle.liste_limites[1]]/np.linalg.norm([obstacle.liste_limites[2]-obstacle.liste_limites[0],obstacle.liste_limites[3]-obstacle.liste_limites[1]]))
        vecteur_normal = np.array([-vecteur_arrete[1],vecteur_arrete[0]])
        difference_position = self.positions[i]-[(obstacle.liste_limites[2]+obstacle.liste_limites[0])/2,(obstacle.liste_limites[3]+obstacle.liste_limites[1])/2]
        projete_normal = np.dot(difference_position, vecteur_normal)
        projete_tang = np.dot(difference_position, vecteur_arrete)
        return projete_normal*vecteur_normal, projete_tang


class Poisson(Boid):
    
    def __init__(self, position_initiale, vitesse_initiale, v_max):
        super().__init__(position_initiale, vitesse_initiale,v_max)
        
    
    def initialisation_matrices(self):
        super().initialisation_matrices()
        self.poissons_dans_rayon_cohesion   = []
        self.poissons_dans_rayon_separation = []
        self.poissons_dans_rayon_alignement = []
        self.predateur_dans_vision = []
        
        
class Predateur(Boid):
    
    def __init__(self, position_initiale, vitesse_initiale, v_max):
        super().__init__(position_initiale, vitesse_initiale,v_max)
    
    
    def initialisation_matrices(self):
        super().initialisation_matrices()
        self.poisson_le_plus_proche = []
        self.predateur_dans_rayon_separation = []
        
class Obstacle():
    def __init__(self, liste_limites):        #####liste_limites de la forme [x1,y1,x2,y2]
        self.liste_limites=np.array(liste_limites)
        self.distance_repulsion = 200
        self.coefficient_repulsion = 100000
        self.longueur = np.linalg.norm([liste_limites[2]-liste_limites[0],liste_limites[3]-liste_limites[1]])    





def generate_poissons():
    return [Poisson([rng.random()*500, rng.random()*500-250],
                    [rng.random()*100-250, rng.random()*100-50],
                    500)
            for _ in range(25)]
    
def copier_poissons(liste_poissons):
    return [Poisson(poisson.position_initiale,poisson.vitesse_initiale,poisson.v_max) for poisson in liste_poissons]

def meshgrid():
    distance_seuil = 100; alpha_separation = 10000; a_rng = 60
    r_cohesion = 400; r_separation = 60; r_alignement = 5
    N = 500
    poissons = generate_poissons()
    
    alpha_cohesion = np.arange(0,1001,200)
    alpha_alignement = np.arange(0,51,10)
    
    X, Y = np.meshgrid(alpha_cohesion, alpha_alignement)
    
    Z = np.zeros((6,6),dtype=float)
    for i in range(6):
        for j in range(6):
            liste = copier_poissons(poissons)
            sim = Simulation(liste, [], N, 0.01, distance_seuil, alpha_cohesion[i], alpha_separation, alpha_alignement[j], a_rng, r_cohesion, r_separation, r_alignement)
            sim.calcul_tableaux()
            print('/',end='')
            Z[i,j] = sim.moyennage_parametre_ordre(400,100)
        print('   ',end='')
    # Affichage avec imshow
    plt.imshow(Z, extent=(alpha_cohesion.min(), alpha_cohesion.max(), alpha_alignement.min(), alpha_alignement.max()), origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar(label="paramètre d'ordre")
    plt.xlabel('alpha_cohesion')
    plt.ylabel('alpha_alignement')
    plt.title("Carte de paramètre d'ordre")
    plt.show()
    
def graphe_pos(simulation,lim=[0,0,0,0]): #xmin,xmax,ymin,ymax
    if lim == [0,0,0,0]:
        xmax = np.max(simulation.liste_de_poissons[0].positions[:,0])
        ymax = np.max(simulation.liste_de_poissons[0].positions[:,1])
        xmin = np.min(simulation.liste_de_poissons[0].positions[:,0])
        ymin = np.min(simulation.liste_de_poissons[0].positions[:,1])
        for poisson in simulation.liste_de_poissons:
            xma = np.max(poisson.positions[:,0])
            yma = np.max(poisson.positions[:,1])
            xmi = np.min(poisson.positions[:,0])
            ymi = np.min(poisson.positions[:,1])
            if xma>xmax:
                xmax = xma
            if yma>ymax:
                ymax = yma
            if xmi<xmin:
                xmin = xmi
            if ymi<ymin:
                ymin=ymi
        lim = [xmin,xmax,ymin,ymax]
    fig, ax = plt.subplots(1,1,figsize=(15,5))
    couleurs = ['b', 'g', 'c', 'm', 'y']
    for i in range(len(simulation.liste_de_poissons)):
        ax.plot(simulation.liste_de_poissons[i].positions[:,0],simulation.liste_de_poissons[i].positions[:,1],couleurs[i%len(couleurs)]+':')
    for predateur in simulation.liste_de_predateurs:
        ax.plot(predateur.positions[:,0],predateur.positions[:,1],'r--')
    for obstacle in simulation.liste_obstacles:
        xo1,yo1,xo2,yo2 = obstacle.liste_limites
        ax.plot([xo1,xo2],[yo1,yo2],'k')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_xlim([lim[0],lim[1]])
    ax.set_ylim([lim[2],lim[3]])
    ax.grid()
    plt.title("positions des boïds (Y en fonction de X)")    
    plt.axis("equal")
    plt.show()

def recherche_alignement():
    distance_seuil = 100; alpha_cohesion = 3; alpha_separation = 10000; a_rng = 60
    r_cohesion = 500; r_separation = 60; r_alignement = 5
    N = 3000
    poissons = []
    for i in range(20):
        a = rng.random()*500
        b = rng.random()*500-250
        c = rng.random()*100-250
        d = rng.random()*100-50
        poissons.append(Poisson([a,b],[c,d],500))
    for alpha_alignement in range(0, 10, 1):
        nouvelle_simu = Simulation(poissons, [], N, 0.01, distance_seuil, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, temps_calcul_ordre=20.0)
        nouvelle_simu.calcul_tableaux()

        print(f"le parametre d'ordre est de {nouvelle_simu.parametre_ordre}")


