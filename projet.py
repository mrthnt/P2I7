import numpy as np



class Simulation :
    
    def __init__(self, liste_de_poissons, liste_de_predateurs, N):
        self.liste_de_poissons = liste_de_poissons
        self.liste_de_predateurs = liste_de_predateurs

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




