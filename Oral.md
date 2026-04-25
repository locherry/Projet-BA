# Oral BA
## Dimensionnement à l'effort tranchant


Pourquoi l'effort tranchant est-il dangereux ?

Une poutre fléchie génère des contraintes tangentes τ dans la section. Pour une section rectangulaire homogène, la distribution est parabolique — nulle aux fibres extrêmes, maximale à l'axe neutre.

Dans une poutre en béton armé, cette contrainte de cisaillement crée des fissures inclinées à ~45°. Sans armature transversale, la rupture est brutale. Avec des cadres, les fissures sont « cousues » et la ruine devient ductile.

### L'analogie du treillis de Ritter-Mörsch

Pour calculer les armatures, on remplace la poutre réelle par un treillis fictif dont les barreaux représentent les flux de forces internes :
Membrure supérieure
Béton comprimé — reprend la résultante F_c de la zone comprimée.
Membrure inférieure
Acier tendu — reprend la résultante F_s des armatures longitudinales.
Bielles de béton
Diagonales comprimées d'angle θ. La limite est V_Rd,max.
Montants (cadres)
Acier transversal en traction — cousent les fissures à l'angle θ.


### Les deux calculs fondamentaux

1. Vérif. 1Résistance des bielles — V_Rd,max
Les diagonales de béton sont des bielles comprimées. Si V_Ed dépasse leur résistance, le béton s'écrase quelle que soit la quantité d'acier.
V_Rd,max = α_cw · b_w · z · ν₁ · f_cd / (cot θ + tan θ)

2. Vérif. 2Dimensionnement des cadres — A_sw/s_t
Chaque cadre « coud » plusieurs fissures inclinées. On calcule la densité d'armatures nécessaire :
densité d'armatures transversales
A_sw / s_t = V_Ed / (z · f_ywd · cot θ)


### Choisir l'angle de bielle θ

L'EC2 autorise un angle θ tel que 1 ≤ cot θ ≤ 2.5, soit θ entre 21.8° et 45°.
Stratégie du code : on part de cot θ = 2.5 (cot_theta_optimal()). Si V_Ed,max dépasse V_Rd,max, on cherche le cot θ minimal. Si même cot θ = 1 ne suffit pas, la section doit être augmentée.

## Épure d'arrêt des barres
L'idée de base : on ne met que ce qu'il faut, là où il le faut
Une poutre fléchie ne sollicite pas ses armatures de la même façon sur toute sa longueur. Le moment est maximal en travée et nul aux appuis — inutile de faire courir toutes les barres d'un bout à l'autre. L'épure d'arrêt répond à la question : à partir d'où peut-on se passer d'une couche de barres ?

1. Le décalage de l'épure — pourquoi on ne coupe pas là où M_Ed = M_Rd
En présence d'armatures transversales, la traction dans les barres longitudinales est augmentée par l'effet de l'analogie du treillis. EC2 §9.2.1.3 impose de décaler l'épure des moments d'une longueur :
a_l = z · cot θ / 2

2. Les couches et leur moment résistant cumulé
Le ferraillage est organisé en couches empilées (bas → haut). Chaque couche apporte une contribution supplémentaire au moment résistant. On calcule M_Rd pour les sous-ensembles cumulatifs :

3. Les points d'arrêt théoriques — où couper ?
Pour chaque couche k, le point d'arrêt théorique est l'abscisse où le moment décalé M_Ed,shifted(x) croise le moment résistant M_Rd de la couche k−1 (c'est-à-dire sans la couche k). En dehors de cet intervalle, la couche k n'apporte rien.

4. L'ancrage — on ne coupe pas net
Couper une barre au point théorique serait dangereux : il faut que la contrainte dans la barre puisse être transmise au béton sur une certaine longueur. EC2 §8.4 impose une longueur d'ancrage l_bd calculée à partir de la résistance d'adhérence f_bd :
f_bd = 2,25 · η₁ · η₂ · f_ctd
puis :
l_b,rqd = (φ/4) · (σ_sd / f_bd)

5. La vérification d'ancrage — est-ce qu'il y a assez de place ?
La dernière vérification consiste à s'assurer que l_bd "rentre" entre l'appui et le point d'arrêt.
Si cette condition n'était pas vérifiée, il faudrait soit réduire le diamètre des barres (l_bd ∝ φ), soit prévoir un anchorage spécial (crochet, platine d'about).


s_t(x) = n · A_t,leg · z · f_ywd · cot θ / V_Ed(x)
s_t est inversement proportionnel à V_Ed. Là où V_Ed est grand → s_t est petit. 
Là où V_Ed → 0 → s_t → ∞, borné par s_t,max.