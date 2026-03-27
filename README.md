# Projet BA — TNP Villeurbanne

Dimensionnement d'une poutre en béton armé en section en T selon l'Eurocode 2.

## Structure

```
├── main.py
├── classes/
│   ├── TSection.py               # Géométrie de la section en T
│   ├── Materials.py              # Béton (Concrete) et acier (Steel)
│   ├── Beam.py                   # Poutre isostatique + charge répartie
│   ├── LoadCalculator.py         # Combinaisons ELU/ELS
│   ├── ReinforcementLayout.py    # Positionnement des armatures
│   └── FlexuralDesign.py         # Calcul AN, A_s
└── plot/
    ├── BeamPlotter.py            # Diagrammes T(x) et M(x)
    └── SectionPlotter.py         # Section + profil de déformation
```

## Données
*Toutes les données doivent être définies dans main.py*
| Paramètre | Valeur |
|---|---|
| Portée | L = 15.25 m |
| Section | T : b_eff = 2.0 m, h = 1.25 m, h_f = 0.15 m, b_w = 0.5 m |
| Béton | C30/37 — f_ck = 30e6 Pa |
| Acier | HA B500B — f_yk = 500e6 Pa |

## Méthode

Flexion simple à l'ELU sans charge axiale (EC2, chapitre 5).
L'axe neutre est localisé dans la table de compression (M_Ed < M_ct),
la section est donc traitée comme une section rectangulaire b_eff × h.

## Lancer le projet

```bash
python main.py
```
