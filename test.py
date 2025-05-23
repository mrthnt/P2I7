###ca fait une boite ou les poissons et les prédateurs sont enfermés
def test_3():
    distance_seuil = 100; alpha_cohesion = 20; alpha_separation = 10000; alpha_alignement = 10; a_rng = 60
    r_cohesion = 400; r_separation = 60; r_alignement = 5; r_predation = 600; r_proies = 700; 
    N = 500
    predateur = Predateur([0,0], [0,15], 600)
    predateur2 = Predateur([6,0], [15,0], 600)
    obstacle1 = Obstacle([-201,400,-201,-400])
    obstacle2 = Obstacle([-201,400,601,400])
    obstacle3 = Obstacle([-200,-401,600,-401])
    obstacle4 = Obstacle([603,401,600,-401])
    liste_obstacle = [obstacle1,obstacle2,obstacle3,obstacle4]
    poissons = generate_poissons()
    nouvelle_simu = Simulation(poissons, [predateur,predateur2],liste_obstacle, N, 0.01, alpha_cohesion, alpha_separation, alpha_alignement, a_rng, r_cohesion, r_separation, r_alignement, r_predation, r_proies)
    nouvelle_simu.calcul_tableaux()
    fenetre = GUI(nouvelle_simu,1,500)
    print(f"le parametre d'ordre à 2 secondes est de {nouvelle_simu.moyennage_parametre_ordre(int(2.0/nouvelle_simu.dt),200)}")

















































