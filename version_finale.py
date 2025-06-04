import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.animation import FuncAnimation
import random as rng

class Simulation :
    
    def __init__(self, liste_de_poissons, liste_de_predateurs, N, dt, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, rayon_cohesion, rayon_separation, rayon_alignement):
        self.liste_de_poissons = liste_de_poissons
        self.liste_de_predateurs = liste_de_predateurs
        self.N = N
        self.dt = dt
        self.initialiser_matrices_poissons()
        
        self.alpha_cohesion = alpha_cohesion
        self.alpha_separation = alpha_separation
        self.alpha_alignement = alpha_alignement
        self.a_rng = a_rng
        self.rayon_cohesion = rayon_cohesion
        self.rayon_separation = rayon_separation
        self.rayon_alignement = rayon_alignement
        
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

    def mise_a_jour_liste_voisin(self, poisson, i, p, angle_de_vision = 140):

        poisson.poissons_dans_rayon_cohesion   = []
        poisson.poissons_dans_rayon_alignement = []
        poisson.poissons_dans_rayon_separation = []
        
        cos_angle_vision = np.cos(np.deg2rad(angle_de_vision))

        for i_poisson in range (0, len(self.liste_de_poissons)):

            if i_poisson != p :
                voisin = self.liste_de_poissons[i_poisson]
                distance_poisson_voisin_i = poisson.distance(voisin, i)

            
                if distance_poisson_voisin_i < max(self.rayon_alignement, self.rayon_cohesion, self.rayon_separation):
                    vecteur_vitesse = poisson.vitesses[i, :]
                    vecteur_poisson_voisin = voisin.positions[i, :] - poisson.positions[i, :]
                    norme_vitesse = np.linalg.norm(vecteur_vitesse)
            
                    cos_angle = np.dot(vecteur_poisson_voisin, vecteur_vitesse) / (norme_vitesse * distance_poisson_voisin_i)
                    
                    if cos_angle > cos_angle_vision :

                        if distance_poisson_voisin_i < self.rayon_separation:

                            poisson.poissons_dans_rayon_separation.append(voisin)
                            
                        if distance_poisson_voisin_i > self.rayon_separation and distance_poisson_voisin_i < self.rayon_alignement:

                            poisson.poissons_dans_rayon_alignement.append(voisin)
                        if distance_poisson_voisin_i > self.rayon_alignement and distance_poisson_voisin_i < self.rayon_cohesion:

                            poisson.poissons_dans_rayon_cohesion.append(voisin)

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
        liste_voisins = poisson.poissons_dans_rayon_alignement
        valeur, vitesse_voisins = self.vitesse_voisins(poisson, i)   ## vecteur vitesse des voisins

        if valeur == None:
            return res
        else:
            res = self.alpha_alignement * ( vitesse_voisins - poisson.vitesses[i, :])
            return res
    
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
        for element in self.liste_de_poissons : 
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
        somme_parametre_ordre = 0.
        for i in range(0, self.N):
 
            for p in range(0, len(self.liste_de_poissons)):
                poisson = self.liste_de_poissons[p]
                
                self.mise_a_jour_liste_voisin(poisson, i, p)

                a_cohesion = self.composante_acceleration_cohesion(poisson, i)

                a_separation = self.composante_acceleration_separation(poisson, i)

                a_alignement = self.composante_acceleration_alignement(poisson, i)
                
                poisson.angle_accel_alea = poisson.angle_accel_alea + rng.uniform(-0.5, 0.5)   ## le vecteur vitesse aleatoire comprend des valeurs allant de -0,4 x a_rng à 0,6 x a_rng
                
                a_aleatoire = np.array([self.a_rng * np.cos(poisson.angle_accel_alea), self.a_rng * np.sin(poisson.angle_accel_alea)])

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
    def __init__(self, simulation, vitesse_lecture = 1.0, coord_lim = 300, suivi = False):
        self.simulation = simulation
        
        self.coord_lim = coord_lim #utile pour update_suivi(frame)
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(-coord_lim, coord_lim)
        self.ax.set_ylim(-coord_lim, coord_lim)
        self.ax.set_aspect('equal')

        self.triangles = []
        self.init_triangles()

        if suivi == True:
            self.ax.grid()
            self.ani = FuncAnimation(self.fig, self.update_suivi, frames=np.arange(0,self.simulation.N-1,1), interval=self.simulation.dt*1000/vitesse_lecture, blit=False)
        else:
            self.ani = FuncAnimation(self.fig, self.update_non_suivi, frames=np.arange(0,self.simulation.N-1,1), interval=self.simulation.dt*1000/vitesse_lecture, blit=True) 
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
        for j in range(len(self.simulation.liste_de_predateurs)):
            pos,vit = self.simulation.liste_de_predateurs[j].positions[frame], self.simulation.liste_de_predateurs[j].vitesses[frame]
            new_coords = self.coords_triangle(pos, vit)
            self.triangles[j+len(self.simulation.liste_de_poissons)].set_xy(new_coords)
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
        for j in range(len(self.simulation.liste_de_predateurs)):
            pos,vit = self.simulation.liste_de_predateurs[j].positions[frame], self.simulation.liste_de_predateurs[j].vitesses[frame]
            new_coords = self.coords_triangle(pos, vit)
            self.triangles[j+len(self.simulation.liste_de_poissons)].set_xy(new_coords)
        return self.triangles
    
    
    def coords_triangle(self,pos,vit,size=0.1):
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
        self.angle_accel_alea = rng.uniform(0, 2 * np.pi)
        
    def ajouter_simulation(self, simulation):
        self.simulation = simulation
        self.initialisation_matrices()
        
    def initialisation_matrices(self):
        self.positions = np.zeros((self.simulation.N+1, 2))
        self.vitesses = np.zeros((self.simulation.N+1, 2))
        self.accelerations = np.zeros((self.simulation.N+1, 2))

        self.positions[0] = self.position_initiale
        self.vitesses[0] = self.vitesse_initiale

        self.poissons_dans_rayon_cohesion   = []
        self.poissons_dans_rayon_separation = []
        self.poissons_dans_rayon_alignement = []

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



def test_3():
    alpha_cohesion = 0; alpha_separation = 2500; alpha_alignement = 0; a_rng = 0
    r_cohesion = 400; r_separation = 100; r_alignement = 200
    N = 2500
    
    #poissons = generate_poissons()
    poissons = [Poisson([100, -100],
                    [50, 50],
                    100),
                Poisson([100, 100],
                    [50, -50],
                    100)]
    nouvelle_simu = Simulation(poissons, [], N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu,1,500)
    print(f"le parametre d'ordre à 2 secondes est de {nouvelle_simu.moyennage_parametre_ordre(1800,700)}")

def generate_poissons():
    return [Poisson([rng.random()*200-100, rng.random()*200-100],
                    [rng.random()*200-100, rng.random()*200-100],
                    100)
            for _ in range(25)]
    
def copier_poissons(liste_poissons):
    return [Poisson(poisson.position_initiale,poisson.vitesse_initiale,poisson.v_max) for poisson in liste_poissons]

def meshgrid():
    alpha_cohesion = 2; alpha_separation = 10000; alpha_alignement = 3; a_rng = 100
    r_separation = 60
    N = 1200
    poissons = generate_poissons()
    
    r_cohesion = np.linspace(0, 200, 15)
    r_alignement = np.linspace(0, 200, 15)
    
    X, Y = np.meshgrid(alpha_cohesion, alpha_alignement)
    
    Z = np.zeros((15,15),dtype=float)
    for i in range(15):
        for j in range(15):
            liste = copier_poissons(poissons)
            sim = Simulation(liste, [], N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion[i], r_separation, r_alignement[j])
            sim.calcul_tableaux()
            print('/',end='')
            Z[i,j] = sim.moyennage_parametre_ordre(900,300)
        print(' ',end='')
    # Affichage avec imshow
    print(Z)
    plt.imshow(Z, extent=(0, 200, 0, 200), origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar(label="paramètre d'ordre")
    plt.xlabel('r_alignement')
    plt.ylabel('r_cohesion')
    plt.title("Carte de paramètre d'ordre")
    plt.show()

test_3()
