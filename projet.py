import numpy as np
import random as rng
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.animation import FuncAnimation

class Simulation :
    
    def __init__(self, liste_de_poissons, liste_de_predateurs, N, dt, distance_seuil, alpha_cohesion, alpha_separation, alpha_alignement, a_rng):
        self.liste_de_poissons = liste_de_poissons
        self.liste_de_predateurs = liste_de_predateurs
        self.N = N
        self.dt = dt
        self.initialiser_matrices_poissons()
        self.distance_seuil = distance_seuil
        self.alpha_cohesion = alpha_cohesion
        self.alpha_separation = alpha_separation
        self.alpha_alignement = alpha_alignement
        self.a_rng = a_rng
        
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
        point_centre = np.zeros((3))
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
        vitesse_moy = np.zeros((3))
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

        
                a_cohesion = np.zeros((3))
                a_separation = np.zeros((3))
                a_alignement = np.zeros((3))
                    
                ## récupération centre masse banc et vitesse banc
                centre_masse_banc = self.centre_masse_poissons(i)   ## centre de masse du banc de poisson 
                vitesse_banc = self.vitesse_banc_poissons(i)             ## vecteur vitesse du banc de poisson (je pense moyenne de la vitesse sur x puis sur y puis sur z)

                a_cohesion = self.alpha_cohesion * ( centre_masse_banc - poisson.positions[i, :])

                a_separation = self.alpha_separation * (poisson.positions[i, :] - voisin.positions[i, :]) / dist_voisin

                a_alignement = self.alpha_alignement * ( vitesse_banc - poisson.vitesses[i, :])
                
                n1_rand = rng.random()*2 - 1 # sur l'intervalle [-1 ; 1]
                n2_rand = rng.random()*2 - 1
                n3_rand = rng.random()*2 - 1
                a_aleatoire = [n1_rand * self.a_rng, n2_rand * self.a_rng, n3_rand * self.a_rng]

                poisson.accelerations[i+1, :] = a_cohesion + a_separation + a_alignement + a_aleatoire


                
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

class GUI:
    def __init__(self, simulation, vitesse_lecture = 1.0, coord_lim = 200):
        self.simulation = simulation #assigne la simulation
                    #initialisation figure
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlim([-coord_lim, coord_lim])
        self.ax.set_ylim([-coord_lim, coord_lim])
        self.ax.set_zlim([-coord_lim, coord_lim])
        self.ax.set_box_aspect([1, 1, 1])
                    #réglages d'affichage
        self.ax.grid(False)
        self.ax.xaxis.pane.fill = False
        self.ax.yaxis.pane.fill = False
        self.ax.zaxis.pane.fill = False
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_zticks([])
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
                    #initialisation des tetraèdres représentant les boids
        self.tetras = []
        self.init_tetras()
                    #animation des tetraèdres représentant les boids
        self.ani = FuncAnimation(self.fig, self.update, frames=self.simulation.N, interval=self.simulation.dt*1000/vitesse_lecture, blit=False)
        plt.show()
    
    
    def init_tetras(self):
        """
        Intialise les tetraedres
        """
        for poisson in self.simulation.liste_de_poissons:
            pos, vit = poisson.position_initiale, poisson.vitesse_initiale
            faces = self.faces_tetra(pos, vit)
            tetra = Poly3DCollection(faces, facecolors='orange', edgecolors = 'k', alpha = 0.8)
            self.ax.add_collection3d(tetra)
            self.tetras.append(tetra)
        for predateur in self.simulation.liste_de_predateurs:
            pos, vit = predateur.position_initiale, predateur.vitesse_initiale
            faces = self.faces_tetra(pos, vit)
            tetra = Poly3DCollection(faces, facecolors='orange', edgecolors = 'k', alpha = 0.8)
            self.ax.add_collection3d(tetra)
            self.tetras.append(tetra)
    
    def update(self,frame):
        """
        Met à jour les coordonnées des tetraedres
        Args:
            frame: frame renseignée de l'animation au temps i
        Returns:
            self.tetras : liste des coordonnées des tetraedres représentant les boids
        """
        for i in range(len(self.simulation.liste_de_poissons)):
            pos,vit = self.simulation.liste_de_poissons[i].positions[frame], self.simulation.liste_de_poissons[i].vitesses[frame]
            nouvelles_faces = self.faces_tetra(pos, vit)
            self.tetras[i].set_verts(nouvelles_faces)
        for j in range(len(self.simulation.liste_de_predateurs)):
            pos,vit = self.simulation.liste_de_predateurs[j].positions[frame], self.simulation.liste_de_predateurs[j].vitesses[frame]
            nouvelles_faces = self.faces_tetra(pos, vit)
            self.tetras[j+len(self.simulation.liste_de_poissons)].set_verts(nouvelles_faces)
        return self.tetras #rajouter le bail de size
        
    def faces_tetra(self,pos,vit,size=1.0):
        """
        Args:
            pos (tuple): tuple contenant x,y,z les coordonnées du boid
            vit (tuple): tuple contenant vx,vy,vz correspondant au vecteur vitesse du boid
            size (float, optional): facteur de taille du boid

        Returns:
            coordonnées des faces d'un tetraedre représentant le boid
        """
        x, y, z = pos
        xvit, yvit, zvit = vit
        h = 50*size             #hauteur du tetraedre
        H = 30*size             #hauteur de la base
        L = 2/np.sqrt(3)*H      #largeur d'un coté de la base
    
        #si vitesse nulle attention à ne pas diviser par une norme nulle -> direction par défaut = (1,0,0)
        if xvit==0 and yvit==0 and zvit==0: 
            xvit = 1
            yvit = 0
            zvit = 0
        norme_vit = np.sqrt(xvit**2+yvit**2+zvit**2)
    
        #on cherche trois vecteurs orthonormaux avec 1 selon la direction du boid:
        xa, ya, za = xvit/norme_vit, yvit/norme_vit, zvit/norme_vit     #selon la direction/vitesse
        xb, yb, zb = -ya, xa, 0                                         #vecteur orthogonal choisi arbitrairement
        xc, yc, zc = -za*xa, -ya*za, xa**2 + ya**2                      #produit vectoriel simplifié des deux précédents

        xbase, ybase, zbase = x - xa*h/3, y - ya*h/3, z - za*h/3        #coordonées du centre de la base du tetraedre
    
        #on calcule les points limites de la base et celui de l'apex du tetraedre:
        base = [[xbase + xb*2/3*H,ybase + yb*2/3*H,zbase + zb*2/3*H],[xbase - xb*H/3 + L/2*xc,ybase -yb*H/3 + L/2*yc,zbase -zb*H/3 + L/2*zc],[xbase - xb*H/3 - L/2*xc,ybase -yb*H/3 - L/2*yc,zbase -zb*H/3 -L/2*zc]]
        apex = [x + 2/3*h*xa, y + 2/3*h*ya, z + 2/3*h*za]
    
        #on renseigne les faces du tetraedre
        faces = [
            [base[0], base[1], base[2]],  # face de la base
            [base[0], base[1], apex],     # faces qui coincident avec l'apex
            [base[1], base[2], apex],
            [base[2], base[0], apex]
        ]
        return faces

class Boid :
    
    def __init__(self, position_initiale, vitesse_initiale):
        self.position_initiale = position_initiale
        self.vitesse_initiale = vitesse_initiale
        
    def ajouter_simulation(self, simulation):
        self.simulation = simulation
        self.initialisation_matrices()
        
    def initialisation_matrices(self):
        self.positions = np.zeros((self.simulation.N+1, 3))
        self.vitesses = np.zeros((self.simulation.N+1, 3))
        self.accelerations = np.zeros((self.simulation.N+1, 3))
        self.positions[0] = self.position_initiale
        self.vitesses[0] = self.vitesse_initiale
        
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
    
       def __init__(self, position_initiale, vitesse_initiale, v_max):
        super().__init__(position_initiale, vitesse_initiale)
        self.v_max = v_max



  # test de classes
def test_1():
    distance_seuil = 100; alpha_cohesion = 1; alpha_separation = 1; alpha_alignement = 1; a_rng = 0.2
    
    poisson1 = Poisson([-1, -2, -3], [0, 0, 0], 20)
    poisson2 = Poisson([2, 2, 2], [5, 0, -5], 20)
    poisson3 = Poisson([0, 0, 0], [1, -3, -7], 20)
    
    la_simu = Simulation([poisson1, poisson2, poisson3], [], 5, 0.1, distance_seuil, alpha_cohesion, alpha_separation, alpha_alignement, a_rng)
    
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
    print("\n")
    
def test_2():
    distance_seuil = 100; alpha_cohesion = 5; alpha_separation = 10000; alpha_alignement = 50; a_rng = 0 
    N = 1000
    poisson4 = Poisson([-100, -100, 0], [100, 100, 50], 500)
    poisson5 = Poisson([100, 100, 100], [-50, -100, 0], 500)
    poisson6 = Poisson([50, 0, -100], [-20, 10, 0], 500)
    

    
    nouvelle_simu = Simulation([poisson4, poisson5, poisson6], [], N, 0.004, distance_seuil, alpha_cohesion, alpha_separation, alpha_alignement, a_rng)
    nouvelle_simu.calcul_tableaux()
    
    fenetre = GUI(nouvelle_simu, 1)
