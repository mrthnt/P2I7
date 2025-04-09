import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Simulation :
    
    def __init__(self, liste_de_poissons, liste_de_predateurs, N, dt):
        self.liste_de_poissons = liste_de_poissons
        self.liste_de_predateurs = liste_de_predateurs
        self.N = N
        self.dt = dt
        self.initialiser_matrices_poissons()
        
    def initialiser_matrices_poissons(self):
        for poisson in self.liste_de_poissons:
            poisson.ajouter_simulation(self)
        
    def voisin_le_plus_proche(self, p, i): 			# i est l'indice temporel concerné de la matrice positions
        poisson_1 = self.liste_de_poissons[p]		# p est l'indice du poisson concerné dans liste_de_poissons
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
                    
    def centre_masse_poissons(self, i): # i est l'indice temporel concerné de la matrice positions
        point_centre = np.zeros((2))
        for poisson in self.liste_de_poissons :
            point_centre += poisson.positions[i]
        point_centre = point_centre / len(self.liste_de_poissons)
        return point_centre
        
    def vitesse_banc_poissons(self, i): # i est l'indice temporel concerné de la matrice positions
        vitesse_moy = np.zeros((2))
        for poisson in self.liste_de_poissons :
            vitesse_moy += poisson.vitesses[i]
        vitesse_moy = vitesse_moy / len(self.liste_de_poissons)
        return vitesse_moy

class GUI:
    def __init__(self, simulation):
        self.simulation = simulation

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(11, projection='2d')
        self.ax.set_xlim([-2, 2])
        self.ax.set_ylim([-2, 2])
        self.ax.set_box_aspect([1, 1])

        self.poly = []
        self.init_poly()

        self.ani = FuncAnimation(self.fig, self.update, frames=self.simulation.N, interval=10, blit=False)
        plt.show()
    #####
    MODIFIER LA SUITE POUR AFFICHER EN 2D
    ####
    
    def init_tetras(self):
        for poisson in self.simulation.liste_de_poissons:
            pos, vit = poisson.position_initiale, poisson.vitesse_initiale
            faces = self.faces_tetra(pos, vit)
            tetra = PolyCollection(faces, facecolors='orange', edgecolors = 'k', alpha = 0.8)
            self.ax.add_collection3d(tetra)
            self.tetras.append(tetra)
        for predateur in self.simulation.liste_de_predateurs:
            pos, vit = predateur.position_initiale, predateur.vitesse_initiale
            faces = self.faces_tetra(pos, vit)
            tetra = PolyCollection(faces, facecolors='orange', edgecolors = 'k', alpha = 0.8)
            self.ax.add_collection3d(tetra)
            self.tetras.append(tetra)
    
    def update(self,frame):
        for i in range(len(self.simulation.liste_de_poissons)):
            pos,vit = self.simulation.liste_de_poissons[i].positions[frame], self.simulation.liste_de_poissons[i].vitesses[frame]
            nouveau_triangle = self.triangle(pos, vit)
            self.tetras[i].set_verts(nouveau_triangle)
        for j in range(len(self.simulation.liste_de_predateurs)):
            pos,vit = self.simulation.liste_de_predateurs[j].positions[frame], self.simulation.liste_de_predateurs[j].vitesses[frame]
            nouvelles_faces = self.faces_tetra(pos, vit)
            self.tetras[j+len(self.simulation.liste_de_poissons)].set_verts(nouvelles_faces)
        return self.tetras #rajouter le bail de size
        
    def triangle(pos,vit,size=1.0):
        """

    Args:
        pos (tuple): tuple contenant x,y,z les coordonées du boid
        vit (tuple): tuple contenant vx,vy,vz correspondant au vecteur vitesse du boid
        size (float, optional): facteur de taille du boid

    Returns:
        coordonées des faces d'un tetraedre représentant le boid
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
