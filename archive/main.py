import numpy as np
from plot_diagramms import plot_diagramms
from longitudinal_reinforcement import plot_neutral_axis

# Paramètres 
L = 15.25  # m

# Charge à l'ELU
q = 91e3  # N/m

Ra = q * L / 2
Rb = q * L / 2

x = np.linspace(0, L, 1000)
T = Ra - q * x
M = Ra * x - q * x**2 / 2

T_max = Ra
T_min = -Rb
M_max = Ra**2 / (2 * q)
x_Mmax = Ra / q

# Section (T)
h_tot = 1.25  # m
b_eff = 2  # m total width
h_f = 0.15  # m height of the flanges
b_w = 0.5  # m width of the web
h_w = h_tot - h_f

d = 0.9 * h_tot  # Distance to the longitudinal reinforcement

# Areas
A_f = b_eff * h_f
A_w = b_w * h_w
A_tot = A_f + A_w

# Centroid y (measured from bottom)
y_f = h_w + h_f / 2  # centre of flanges from bottom
y_w = h_w / 2  # centre of web from bottom

y_G = (A_f * y_f + A_w * y_w) / A_tot

# Distances from centroid
d_f_y = y_f - y_G
d_web_y = y_w - y_G
d_f_z = 0
d_web_z = 0

# Inertia
I_f_z = (b_eff * h_f**3) / 12 + A_f * d_f_y**2
I_web_z = (b_w * h_w**3) / 12 + A_w * d_web_y**2
I_tot_z = I_f_z + I_web_z

I_f_y = (b_eff**3 * h_f) / 12 + A_f * d_f_z**2
I_web_y = (b_w**3 * h_w) / 12 + A_w * d_web_z**2
I_tot_y = I_f_y + I_web_y


# Matériaux
## Béton C30/37
# slide 20 https://moodle.insa-lyon.fr/pluginfile.php/438491/mod_resource/content/1/3%20-%20Caract%C3%A9ristiques%20du%20b%C3%A9ton%20et%20de%20lacier%20darmature.pdf
f_ck = 30.0e6  # Pa
f_ctm = 2.9e6  # Pa
gamma_c = 1.5
f_cd = f_ck / gamma_c  # Pa

## Acier FeE500B
gamma_s = 1.15
f_yk = 500e6  # Limite d'élasticité caractéristique fyk (Pa)
f_yd = f_yk / gamma_s

M_Ed = M_max  # N·m — moment de calcul à l'ELU (mi-travée)


plot_diagramms(L, q, Ra, Rb, x, T, M, T_max, T_min, M_max, x_Mmax)
plot_neutral_axis(M_Ed, f_cd, f_ctm, f_yk, f_yd, h_tot, b_eff, h_f, b_w, d)
