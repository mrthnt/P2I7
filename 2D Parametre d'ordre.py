import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.animation import FuncAnimation
import random as rng
from scipy.spatial import cKDTree

class Simulation :
    
    def __init__(self, liste_de_poissons, liste_de_predateurs, N, dt, distance_seuil, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, temps_calcul_ordre=3.0):
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
        self.voisins_plus_proches = np.zeros(len(self.liste_de_poissons))
        self.temps_calcul_ordre = temps_calcul_ordre
        
    def initialiser_matrices_poissons(self):
        """
        Cree la matrice des poissons
        """
        for poisson in self.liste_de_poissons:
            poisson.ajouter_simulation(self)
        
    def mettre_a_jour_voisins(self, i):
        positions = np.array([poisson.positions[i] for poisson in self.liste_de_poissons])
        tree = cKDTree(positions)
        # Renvoie l'indice du plus proche voisin pour chaque poisson
        _, indices_voisins = tree.query(positions, k=2)  # k=2 pour exclure soi-même
        return indices_voisins[:, 1]
                    
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
        taille_banc = len(self.liste_de_poissons)
        for i in range(0, self.N):

            ## récupération centre masse banc et vitesse banc
            centre_masse_banc = self.centre_masse_poissons(i)   ## centre de masse du banc de poisson 
            vitesse_banc = self.vitesse_banc_poissons(i)             ## vecteur vitesse du banc de poisson (je pense moyenne de la vitesse sur x puis sur y puis sur z)
            self.voisins_plus_proches = self.mettre_a_jour_voisins(i)
            for p in range(0, len(self.liste_de_poissons)):
                
                poisson = self.liste_de_poissons[p]
                
                ## Calcul accélérations au rang n+1
                voisin = self.liste_de_poissons[self.voisins_plus_proches[p]]  ## objet du poisson le plus proche
                dist_voisin = poisson.distance(voisin, i)     ## distance avec le poisson le plus proche 

        
                a_cohesion = np.zeros((2))
                a_separation = np.zeros((2))
                a_alignement = np.zeros((2))

                a_cohesion = self.alpha_cohesion * ( centre_masse_banc - poisson.positions[i, :])

                a_separation = self.alpha_separation * (poisson.positions[i, :] - voisin.positions[i, :]) / dist_voisin**2

                a_alignement = self.alpha_alignement * ( vitesse_banc - poisson.vitesses[i, :])
                
                n1_rand = rng.random()*2 - 1 # sur l'intervalle [-1 ; 1]
                n2_rand = rng.random()*2 - 1
                a_aleatoire = [n1_rand * self.a_rng, n2_rand * self.a_rng]

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
            
            temps_simu = i*self.dt
            liste_parametre_ordre = []
            if temps_simu >= self.temps_calcul_ordre and temps_simu <= (self.temps_calcul_ordre + 50*self.dt) :
                vitesse_totale = [0.,0.]
                for element in self.liste_de_poissons : 
                    vitesse_totale += np.array(element.vitesses[i])
                liste_parametre_ordre.append(np.linalg.norm(vitesse_totale)/(taille_banc * np.linalg.norm(vitesse_banc)))
        liste_parametre_ordre = np.array(liste_parametre_ordre)
        self.parametre_ordre=np.mean(liste_parametre_ordre)
    
    




class GUI:
    def __init__(self, simulation, vitesse_lecture = 1.0, coord_lim = 200):
        self.simulation = simulation
        
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(-coord_lim, coord_lim)
        self.ax.set_ylim(-coord_lim, coord_lim)
        self.ax.set_aspect('equal')

        self.triangles = []
        self.init_triangles()

        self.ani = FuncAnimation(self.fig, self.update, frames=np.arange(0,self.simulation.N-1,1), interval=self.simulation.dt*1000/vitesse_lecture, blit=True)
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
    
    def __init__(self, position_initiale, vitesse_initiale, v_max):
        super().__init__(position_initiale, vitesse_initiale)
        self.v_max = v_max



def test_1():
    distance_seuil = 100; alpha_cohesion = 1; alpha_separation = 1; alpha_alignement = 1; a_rng = 0.2
    
    poisson1 = Poisson([-1, -2], [0, 0], 20)
    poisson2 = Poisson([2, 2], [5, 0], 20)
    poisson3 = Poisson([0, 0], [1, -3], 20)
    
    la_simu = Simulation([poisson1, poisson2, poisson3], [], 5, 0.1, distance_seuil, alpha_cohesion, alpha_separation, alpha_alignement, a_rng)
    
    print(poisson1)
    print()
    print(poisson2)
    print()
    
    distance = poisson1.distance(poisson2, 0)
    print(distance)
    print()
    
    centre = la_simu.centre_masse_poissons(0)
    print(centre)
    v_moy = la_simu.vitesse_banc_poissons(0)
    print(v_moy)
    print("\n")

def test_2():
    distance_seuil = 100; alpha_cohesion = 5; alpha_separation = 100000; alpha_alignement = 50; a_rng = 1000
    N = 5000
    poisson4 = Poisson([-100, -100], [10, 10], 500)
    poisson5 = Poisson([100, 100], [-10, -10], 500)
    poisson6 = Poisson([50, 0], [0, 0], 500)
    

    
    nouvelle_simu = Simulation([poisson4, poisson5, poisson6], [], N, 0.01, distance_seuil, alpha_cohesion, alpha_separation, alpha_alignement, a_rng)
    nouvelle_simu.calcul_tableaux()
    
    fenetre = GUI(nouvelle_simu, vitesse_lecture = 1.0, coord_lim=500)
    
def test_3():
    distance_seuil = 100; alpha_cohesion = 5; alpha_separation = 10000; alpha_alignement =5; a_rng = 1000
    N = 5000
    poissons = []
    for i in range(20):
        a = rng.random()*500
        b = rng.random()*500-250
        c = rng.random()*100-250
        d = rng.random()*100-50
        poissons.append(Poisson([a,b],[c,d],500))
    nouvelle_simu = Simulation(poissons, [], N, 0.01, distance_seuil, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, temps_calcul_ordre=2.0)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu,1,500)
    print(f"le parametre d'ordre est de {nouvelle_simu.parametre_ordre}")
    
test_3()
