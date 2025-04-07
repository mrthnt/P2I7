import numpy as np



class Simulation :
    
    def __init__(self, liste_de_poissons, liste_de_predateurs, N):
        self.liste_de_poissons = liste_de_poissons
        self.liste_de_predateurs = liste_de_predateurs
        self.N = N

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

class GUI:
    def __init__(self, simulation):
        self.simulation = simulation

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlim([-2, 2])
        self.ax.set_ylim([-2, 2])
        self.ax.set_zlim([0, 2])
        self.ax.set_box_aspect([1, 1, 1])

        self.poly = []
        self.init_poly()

        self.ani = FuncAnimation(self.fig, self.update, frames=self.simulation.N, interval=10, blit=False)
        plt.show()
    
    
    def init_tetra(self):
        pass
    
    def update(self,frame):
        for i in range(len(simulation.liste_de_poissons)):
            pos,vit = self.simulation.liste_de_poissons[i].positions[frame], self.simulation.liste_de_poissons[i].vitesses[frame]
            nouvelles_faces = faces_tetra(pos, vit)
            self.poly[i].set_verts(nouvelles_faces)
        return self.poly
        
    def faces_tetra(pos,vit,size=1.0):
    """

    Args:
        pos (tuple): tuple contenant x,y,z les coordonées du boid
        vit (tuple): tuple contenant vx,vy,vz correspondant au vecteur vitesse du boid
        size (float, optional): facteur de taille du boid

    Returns:
        coordonées des faces d'un tetraedre représentant le boid
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
    
    def __init__(self, position_initiale, vitesse_initiale, N, dt):
        self.position_initiale = position_initiale
        self.vitesse_initiale = vitesse_initiale
        self.N = N
        self.dt = dt
        self.initialisation_matrices()
        
    def initialisation_matrices(self):
        self.positions = np.zeros((self.N+1, 3))
        self.vitesses = np.zeros((self.N+1, 3))
        self.accelerations = np.zeros((self.N+1, 3))
        
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
    
    def __init__(self, position_initiale, vitesse_initiale, R_attraction_poisson, R_repulsion_poisson, v_max, N, dt):
        super().__init__(position_initiale, vitesse_initiale, N, dt)
        self.R_attraction_poisson = R_attraction_poisson
        self.R_repulsion_poisson = R_repulsion_poisson
        self.v_max = v_max



  # test de classes
def test_de_classes():
    poisson1 = Poisson([0, 0, 1], [0, 0, 0], 10, 3, 20, 5, 0.1)
    poisson2 = Poisson([0, 2, 0], [0, 0, 0], 10, 3, 20, 5, 0.1)
    print(poisson1)
    print()
    print(poisson2)
    print()

    distance = poisson1.distance(poisson2, 0)
    print(distance) 
    print()

    poisson_voisin = la_simu.voisin_le_plus_proche(0, 0)
    print(poisson_voisin)




