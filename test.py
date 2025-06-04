###ca fait une boite ou les poissons et les prédateurs sont enfermés
def boite():
    alpha_cohesion = 20; alpha_separation = 10000; alpha_alignement = 10; a_rng = 60
    r_cohesion = 400; r_separation = 60; r_alignement = 200; r_predation = 600; r_proies = 700; 
    N = 500
    predateur = Predateur([0,0], [0,15], 600)
    predateur2 = Predateur([6,0], [15,0], 600)
    obstacle1 = Obstacle([-201,400,-201,-400])
    obstacle2 = Obstacle([-201,400,601,400])
    obstacle3 = Obstacle([-200,-401,600,-401])
    obstacle4 = Obstacle([603,401,600,-401])
    liste_obstacle = [obstacle1,obstacle2,obstacle3,obstacle4]
    predateurs = [predateur, predateur2]
    poissons = generate_poissons()
    nouvelle_simu = Simulation(poissons,predateurs,liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu,1,500)

###ca fait une batiment avec 2 assaillants et 17 civils
def batiment1():
    alpha_cohesion = 20; alpha_separation = 10000; alpha_alignement = 10; a_rng = 60
    r_cohesion = 400; r_separation = 60; r_alignement = 200; r_predation = 600; r_proies = 700; 
    N = 500
    predateur = Predateur([0,0], [0,15], 600)
    predateur2 = Predateur([6,0], [15,0], 600)
    bord_gauche = Obstacle([-200,400,-200,-400])
    bord_haut = Obstacle([-200,400,600,400])
    bord_bas_gauche = Obstacle([-200,-400,100,-400])
    bord_bas_droit = Obstacle([600,-400,300,-400])
    bord_droit = Obstacle([600,400,600,-400])
    centre_gauche_haut = Obstacle([0,50,100,150])
    centre_gauche_bas = Obstacle([0,50,100,-50])
    centre_droit_haut = Obstacle([200,50,100,150])
    centre_droit_bas = Obstacle([200,50,100,-50])
    ligne_bas = Obstacle([200,-200,600,-200])
    liste_obstacle = [bord_bas_droit,bord_droit,bord_haut,bord_gauche,bord_bas_gauche,
                      centre_droit_bas,centre_droit_haut,centre_gauche_haut,centre_gauche_bas,
                      ligne_bas
                      ]
    poisson1 = Poisson([200,300], [5,0], 500)
    poisson2 = Poisson([250,300], [5,0], 500)
    poisson3 = Poisson([300,300], [5,0], 500)
    poisson4 = Poisson([400,300], [5,0], 500)
    poisson5 = Poisson([450,300], [5,0], 500)
    poisson6 = Poisson([200,10], [5,0], 500)
    poisson7 = Poisson([-200,0], [5,0], 500)
    poisson12 = Poisson([-200,0], [5,0], 500)
    poisson8 = Poisson([250,-300], [5,0], 500)
    poisson9 = Poisson([400,-325], [5,0], 500)
    poisson10 = Poisson([200,400], [5,0], 500)
    poisson11 = Poisson([200,400], [5,0], 500)
    poisson13 = Poisson([-150,-200],[3,1], 500)
    poisson14 = Poisson([-100,-200],[-2,0], 500)
    poisson15 = Poisson([-50,250],[3,0],500)
    poisson16 = Poisson([-75,150], [0,-100], 500)
    poisson17 = Poisson([0,5],[9,0],500)
    predateurs = [predateur,predateur2]
    poissons = [poisson1,poisson2,poisson3,poisson4,poisson5,poisson6,poisson7,poisson8,poisson9,poisson10,poisson11,poisson12,poisson13,poisson14,poisson15,poisson16]
    nouvelle_simu = Simulation(poissons, predateurs,liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu,1,500,False)


### verifiaction que le mur coupe la vision
def test_avec_mur():
    alpha_cohesion = 20; alpha_separation = 10000; alpha_alignement = 10; a_rng = 60
    r_cohesion = 400; r_separation = 60; r_alignement = 5; r_predation = 600; r_proies = 700; 
    N = 500
    ligne = Obstacle([0,300,0,-300])
    liste_obstacle = [ligne]
    poisson1 = Poisson([100,300], [5,0], 500)
    poisson2 = Poisson([150,300], [5,0], 500)
    poisson3 = Poisson([100,300], [5,0], 500)
    poisson4 = Poisson([100,300], [5,0], 500)
    poisson5 = Poisson([150,300], [5,0], 500)
    poisson6 = Poisson([100,10], [5,0], 500)
    poisson7 = Poisson([-100,0], [5,0], 500)
    poisson12 = Poisson([-100,0], [5,0], 500)
    poisson8 = Poisson([-150,-300], [5,0], 500)
    poisson9 = Poisson([100,-325], [5,0], 500)
    poisson10 = Poisson([100,400], [5,0], 500)
    poisson11 = Poisson([100,400], [5,0], 500)
    poisson13 = Poisson([-150,-200],[3,1], 500)
    poisson14 = Poisson([-100,-200],[-2,0], 500)
    poisson15 = Poisson([-50,250],[3,0],500)
    poisson16 = Poisson([-75,150], [0,-100], 500)
    predateurs = []
    poissons = [poisson1,poisson2,poisson3,poisson4,poisson5,poisson6,poisson7,poisson8,poisson9,poisson10,poisson11,poisson12,poisson13,poisson14,poisson15,poisson16]
    nouvelle_simu = Simulation(poissons, predateurs,liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu)

def test_sans_mur():
    distance_seuil = 100; alpha_cohesion = 20; alpha_separation = 10000; alpha_alignement = 10; a_rng = 60
    r_cohesion = 400; r_separation = 60; r_alignement = 5; r_predation = 600; r_proies = 700; 
    N = 500
    liste_obstacle = []
    poisson1 = Poisson([100,300], [5,0], 500)
    poisson2 = Poisson([150,300], [5,0], 500)
    poisson3 = Poisson([100,300], [5,0], 500)
    poisson4 = Poisson([100,300], [5,0], 500)
    poisson5 = Poisson([150,300], [5,0], 500)
    poisson6 = Poisson([100,10], [5,0], 500)
    poisson7 = Poisson([-100,0], [5,0], 500)
    poisson12 = Poisson([-100,0], [5,0], 500)
    poisson8 = Poisson([-150,-300], [5,0], 500)
    poisson9 = Poisson([100,-325], [5,0], 500)
    poisson10 = Poisson([100,400], [5,0], 500)
    poisson11 = Poisson([100,400], [5,0], 500)
    poisson13 = Poisson([-150,-200],[3,1], 500)
    poisson14 = Poisson([-100,-200],[-2,0], 500)
    poisson15 = Poisson([-50,250],[3,0],500)
    poisson16 = Poisson([-75,150], [0,-100], 500)
    predateurs = []
    poissons = [poisson1,poisson2,poisson3,poisson4,poisson5,poisson6,poisson7,poisson8,poisson9,poisson10,poisson11,poisson12,poisson13,poisson14,poisson15,poisson16]
    nouvelle_simu = Simulation(poissons, predateurs,liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu)

def cas_test1():
    alpha_cohesion = 0; alpha_separation = 0; alpha_alignement = 0; a_rng = 0
    r_cohesion = 400; r_separation = 60; r_alignement = 200; r_predation = 600; r_proies = 700; 
    N = 250
    
    poisson1 = Poisson([-100,100], [100,-100], 500)
    poisson2 = Poisson([-100,-100], [100,100], 500)
    predateurs = []
    liste_obstacle = []
    poissons = [poisson1,poisson2]
    nouvelle_simu = Simulation(poissons, predateurs,liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu)

def cas_test2():
    alpha_cohesion = 0; alpha_separation = 10000; alpha_alignement = 0; a_rng = 0
    r_cohesion = 400; r_separation = 60; r_alignement = 200; r_predation = 600; r_proies = 700; 
    N = 250
    
    poisson1 = Poisson([-100,100], [100,-100], 500)
    poisson2 = Poisson([-100,-100], [100,100], 500)
    predateurs = []
    liste_obstacle = []
    poissons = [poisson1,poisson2]
    nouvelle_simu = Simulation(poissons, predateurs,liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu)

def cas_test3():
    alpha_cohesion = 20; alpha_separation = 0; alpha_alignement = 0; a_rng = 0
    r_cohesion = 400; r_separation = 60; r_alignement = 200; r_predation = 600; r_proies = 700; 
    N = 500
    
    poisson1 = Poisson([-100,100], [100,-100], 500)
    poisson2 = Poisson([-100,-100], [100,100], 500)
    predateurs = []
    liste_obstacle = []
    poissons = [poisson1,poisson2]
    nouvelle_simu = Simulation(poissons, predateurs,liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu,1,500)
    
def cas_test4():
    alpha_cohesion = 0; alpha_separation = 0; alpha_alignement = 10; a_rng = 0
    r_cohesion = 400; r_separation = 60; r_alignement = 200; r_predation = 600; r_proies = 700; 
    N = 250
    
    poisson1 = Poisson([-100,100], [100,-100], 500)
    poisson2 = Poisson([-100,-100], [100,100], 500)
    predateurs = []
    liste_obstacle = []
    poissons = [poisson1,poisson2]
    nouvelle_simu = Simulation(poissons, predateurs,liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu,1,500)

#### PENSER A MODIF LES VALEURS D'ANGLE VISION
def cas_test5():
    alpha_cohesion = 20; alpha_separation = 10000; alpha_alignement = 10; a_rng = 60
    r_cohesion = 400; r_separation = 60; r_alignement = 200; r_predation = 600; r_proies = 700; 
    N = 250
    
    poisson1 = Poisson([-50,50], [0,100], 500)
    poisson2 = Poisson([0,0], [100,0], 500)
    predateurs = []
    liste_obstacle = []
    poissons = [poisson1,poisson2]
    nouvelle_simu = Simulation(poissons, predateurs,liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu,1,500)
