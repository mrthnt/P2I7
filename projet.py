import numpy as np



class Simulation :
    
    def __init__(self, liste_de_poissons, liste_de_predateurs, N):
        self.liste_de_poissons = liste_de_poissons
        self.liste_de_predateurs = liste_de_predateurs



class Boid :
    
    def __init__(self, position_initiale, vitesse_initiale, N, dt):
        self.position_initiale = position_initiale
        self.vitesse_initiale = vitesse_initiale
        self.N = N
        self.initialisation_matrices()
        self.dt = dt
        
    def initialisation_matrices(self):
        self.positions = np.zeros((self.N+1, 3))
        self.vitesses = np.zeros((self.N+1, 3))
        self.accelerations = np.zeros((self.N+1, 3))
        
        self.positions[0] = self.position_initiale
        self.vitesses[0] = self.vitesse_initiale

    def __str__(self):
        text = "Positions :\n" + str(self.positions) + "\nVitesses :\n" + str(self.vitesses) + "\nAccelerations :\n" + str(self.accelerations)
        return text



class Poisson(Boid):
    
    def __init__(self, position_initiale, vitesse_initiale, R_attraction_poisson, R_repulsion_poisson, v_max, N, dt):
        super().__init__(position_initiale, vitesse_initiale, N, dt)
        self.R_attraction_poisson = R_attraction_poisson
        self.R_repulsion_poisson = R_repulsion_poisson
        self.v_max = v_max

  
