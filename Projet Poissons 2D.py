import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.animation import FuncAnimation
N = 500
poisson = Poisson([0, 0], [0, 1], 10, 3, 20)
simu = Simulation([poisson],[],N, 0.1)
simu.initialiser_matrices_poissons()
gui = GUI(simu)

class Simulation :
    
    def __init__(self, liste_de_poissons, liste_de_predateurs, N, dt):
        self.liste_de_poissons = liste_de_poissons
        self.liste_de_predateurs = liste_de_predateurs
        self.N = N
        self.dt = dt
        self.initialiser_matrices_poissons()
        
    def initialiser_matrices_poissons(self):
        """
        Cree la matrice des poissons
        """
        for poisson in self.liste_de_poissons:
            poisson.ajouter_simulation(self)
        
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
                distance = poisson_1.distance(poisson_2, i)
                if distance < distance_min :
                    distance_min = distance
                    poisson_voisin = poisson_2
        return poisson_voisin
                    
    def centre_masse_poissons(self, i):
        """
        Renvoie le centre de masse des poissons
        Args:
            i: indice temporel concerné de la matrice positions

        Returns:
            point_centre : les coordonnees du centre de masse des poissons
        """
        point_centre = np.zeros((2))
        for poisson in self.liste_de_poissons :
            point_centre += poisson.positions[i]
        point_centre = point_centre / len(self.liste_de_poissons)
        return point_centre
        
    def vitesse_banc_poissons(self, i): # i est l'indice temporel concerné de la matrice positions
        """
        Renvoie la vitesse du banc de poissons
        Args:
            i: indice temporel concerné de la matrice positions

        Returns:
            vitesse_moy : les coordonnées du vecteur vitesse du banc de poissons
        """
        vitesse_moy = np.zeros((2))
        for poisson in self.liste_de_poissons :
            vitesse_moy += poisson.vitesses[i]
        vitesse_moy = vitesse_moy / len(self.liste_de_poissons)
        return vitesse_moy

        def calcul_tableaux(self):       

        for i in range(0, self.N):

            for p in range(0, len(self.liste_de_poissons)):
                
                poisson = self.liste_de_poissons[p]
                
                ## Calcul accélérations au rang n+1
                voisin = self.voisin_le_plus_proche(p, i)  ## objet du poisson le plus proche
                dist_voisin = poisson.distance(voisin, i)     ## distance avec le poisson le plus proche 

                if dist_voisin <= self.distance_seuil:         ## distance seuil à définir
        
                    a_cohesion = np.zeros((2))
                    a_separation = np.zeros((2))
                    a_alignement = np.zeros((2))
                    
                    ## récupération centre masse banc et vitesse banc
                    centre_masse_banc = self.centre_masse_poissons(i)   ## centre de masse du banc de poisson 
                    vitesse_banc = self.vitesse_banc_poissons(i)             ## vecteur vitesse du banc de poisson (je pense moyenne de la vitesse sur x puis sur y puis sur z)

                    a_cohesion = self.alpha_cohesion * ( centre_masse_banc - poisson.positions[i, :])

                    a_separation = self.alpha_separation * (poisson.positions[i, :] - voisin.positions[i, :]) / dist_voisin

                    a_alignement = self.alpha_alignement * ( vitesse_banc - poisson.vitesses[i, :])

                    poisson.accelerations[i+1, :] = a_cohesion + a_separation + a_alignement

#                 else:   
#                     n_rand = randint()
#                     poisson.accelerations[i, :] = n_rand * self.a_max
                
                ## Calcul du vecteur vitesse au rang n+1
                poisson.vitesses[i+1, :] = poisson.vitesses[i, :] + self.dt * poisson.accelerations[i+1, :]

                ## Ajustement du vecteur vitesse pour avoir une vitesse inférieur à v_max
#                 vecteur_vitesse = poisson.vitesses[i+1, :]
#                 poisson.vitesses[i+1, :] = (poisson.vitesses[i+1, :] * poisson.v_max) / np.linalg.norm(vecteur_vitesse)

                ## Calcul du vecteur position au rang n+1
                poisson.positions[i+1, :] = poisson.positions[i, :] + self.dt * poisson.vitesses[i+1, :]




class GUI:
    def __init__(self, simulation):
        self.simulation = simulation

        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(-200, 200)
        self.ax.set_ylim(-200, 200)
        self.ax.set_aspect('equal')

        self.triangles = []
        self.init_triangles()

        self.ani = FuncAnimation(self.fig, self.update, frames=np.arange(0,self.simulation.N-1,1), interval=10, blit=True)
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
            triangle = Polygon(coords, closed=True, color='skyblue')
            self.ax.add_patch(triangle)
            self.triangles.append(triangle)
    
    def update(self,frame):
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
        for j in range(len(self.simulation.liste_de_predateurs)):
            pos,vit = self.simulation.liste_de_predateurs[j].positions[frame], self.simulation.liste_de_predateurs[j].vitesses[frame]
            new_coords = self.coords_triangle(pos, vit)
            self.triangles[j+len(self.simulation.liste_de_poissons)].set_xy(new_coords)
        return self.triangles
    
    
    def coords_triangle(self,pos,vit,size=1.0):
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
    
    def __init__(self, position_initiale, vitesse_initiale):
        self.position_initiale = position_initiale
        self.vitesse_initiale = vitesse_initiale
        
    def ajouter_simulation(self, simulation):
        self.simulation = simulation
        self.initialisation_matrices()
        
    def initialisation_matrices(self):
        self.positions = np.zeros((self.simulation.N+1, 2))
        self.vitesses = np.zeros((self.simulation.N+1, 2))
        self.accelerations = np.zeros((self.simulation.N+1, 2))
        self.positions[0] = self.position_initiale
        self.vitesses[0] = self.vitesse_initiale

    def __str__(self):
        text = "Positions :\n" + str(self.positions) + "\nVitesses :\n" + str(self.vitesses) + "\nAccelerations :\n" + str(self.accelerations)
        return text

    def distance(self, boid, i): 	#distance avec un boid à l'indice i
        diff_position = self.positions[i] - boid.positions[i]
        norme = np.linalg.norm(diff_position)
        return norme


class Poisson(Boid):
    
    def __init__(self, position_initiale, vitesse_initiale, R_attraction_poisson, R_repulsion_poisson, v_max):
        super().__init__(position_initiale, vitesse_initiale)
        self.R_attraction_poisson = R_attraction_poisson
        self.R_repulsion_poisson = R_repulsion_poisson
        self.v_max = v_max



def jeu_de_test():
    poisson1 = Poisson([-1, -2], [0, 0], 10, 3, 20)
    poisson2 = Poisson([2, 2], [5, 0], 10, 3, 20)
    poisson3 = Poisson([0, 0], [1, -3], 10, 3, 20)
    
    la_simu = Simulation([poisson1, poisson2, poisson3], [], 5, 0.1)
    
    print(poisson1)
    print()
    print(poisson2)
    print()
    
    distance = poisson1.distance(poisson2, 0)
    print(distance)
    print()
    
    poisson_voisin = la_simu.voisin_le_plus_proche(0, 0)
    print(poisson_voisin)
    print()
    
    centre = la_simu.centre_masse_poissons(0)
    print(centre)
    v_moy = la_simu.vitesse_banc_poissons(0)
    print(v_moy)
