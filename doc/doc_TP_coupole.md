---
title: Codes pour l'ASTronomie à ORsay
subtitle: Installation et utilisation
thanks: <https://pypi.org/project/castor-orsay/>
date: Septembre 2021
shortauthor: A. Stcherbinine & G. Pelouze
author:
    - Aurélien Stcherbinine (<aurelien.stcherbinine@ias.u-psud.fr>)
    - Gabriel Pelouze (<gabriel.pelouze@ias.u-psud.fr>)

lang: fr
toc: False
numbersections: true
toc-depth: 2
papersize: a4
geometry:
    - top=30mm
    - bottom=32mm
    - left=20mm
    - right=30mm
linkcolor: blue
urlcolor: blue
toccolor: blue
header-includes: |
    \usepackage[mono=false]{libertine}

    \let\oldtableofcontents\tableofcontents
    \renewcommand{\tableofcontents}{\oldtableofcontents\clearpage}

    \setlength{\parindent}{0pt}
    \setlength{\parskip}{1ex}
    \renewcommand{\labelitemi}{–}

    \usepackage{fancyhdr}
    \usepackage{titling}
    \fancyhead{}
    \fancyhead[c]{\textit{A. Stcherbinine \& G. Pelouze :} \thetitle}
    \renewcommand{\headrulewidth}{0.4pt}
    \pagestyle{fancy}

    \usepackage[linecolor=lightgray,linewidth=1pt,skipabove=12pt,skipbelow=15pt]{mdframed}
    \surroundwithmdframed{Shaded}
...

# Installation de CASTOR sur les sessions UPSaclay
En raison des restrictions utilisateurs sur les machines du département, il nous
est impossible d'installer le module via `pip` (même en local).
Toutefois, il s'avère que cette commande est autorisée dans les sessions émulées
lors d'une connexion ssh. 
L'installation devra donc se faire depuis une session ssh avant de pouvoir utiliser
le module en local.

## Identifiants
Ce manuel vise à détailler l'installation et l'utilisation du module `castor-orsay` sur les
machines des salles informatiques du bâtiment 625 ($\hbar$) de l'Université Paris-Saclay.
Vos identifiants de connexion aux postes sont rappelés ci-dessous (à condition de posséder un 
compte Adonis de l'Université Paris-Saclay).

* **login** : 
  - entier : prenom.nom *(ex: gerard.dupont)*
  - ou court : 1ère lettre prénom + 6 nom *(ex: gdupont)*
* **password** : Votre mot de passe Adonis.

## Étape 1 : en ssh
Ouvrir un terminal et se connecter en ssh à sa session en tapant (remplacer *login* par votre identifiant
personnel) :

~~~bash
$ ssh login@ssh1.pgip.universite-paris-saclay.fr
~~~

ou

~~~bash
$ ssh login@ssh2.pgip.universite-paris-saclay.fr
~~~

Puis installer le module python `castor-orsay` dans votre répertoire :

~~~bash
$ pip3 install castor-orsay --user
~~~

Se déconnecter ensuite de la session ssh en tapant :

~~~bash
$ exit
~~~


## Étape 2 : en local
Le module `castor-orsay` est désormais installé dans votre répertoire `/home`
(option `--user` de `pip3`), mais il faut maintenant indiquer à python où aller le chercher.

Pour cela, ouvrir les fichiers `.bashrc` et `.bash_profile` (si le fichier n'existe pas, le créer) via :

~~~bash
$ gedit .bashrc .bash_profile &
~~~

Puis, ajouter la ligne suivante à la fin de chacun d'entre eux.

~~~bash
export PATH=$PATH:$HOME/.local/bin
~~~

Recharger enfin le fichier `.bashrc` (via la commande `source ~/.bashrc`) ou relancer un terminal.

**Bravo vous avez installé CASTOR !**

**Note :** Les fichiers `.bashrc` et `.bash_profile` sont des fichiers cachés
situés dans votre `/home`. 
Vous pouvez les afficher via la commande :

~~~bash
$ ls -a
~~~


# Utilisation

Pour voir les différentes option disponibles, taper dans un terminal :

~~~bash
$ castor_ <Tab><Tab>
~~~

Il vous sera alors affiché le résultat suivant, à savoir la liste des complétions
possibles. Ici les différentes fonctions de CASTOR.

~~~bash
$ castor_
castor_align                    castor_align_spectra    castor_exoplanet_analysis       
castor_pointing_analysis        castor_prepare          castor_rotate_spectra
castor_wavelength_calibration
$ castor_
~~~

Puis vous pouvez afficher l'aide de chaque fonction avec `-h`, tapez par exemple :

~~~bash
$ castor_prepare -h

usage: castor_prepare [-h] [--sci-path SCI_PATH]
                      [--sci-dark-path SCI_DARK_PATH] [--flat-path FLAT_PATH]
                      [--flat-dark-path FLAT_DARK_PATH] [-o OUTPUT] [-O]
                      target_name

Prepare a series of images.

positional arguments:
  target_name           Name of the target

optional arguments:
  -h, --help            show this help message and exit
  --sci-path SCI_PATH   Directory containing the science FITS. Default:
                        {target_name}/sci
  --sci-dark-path SCI_DARK_PATH
                        Directory containing the science dark FITS. Default:
                        {target_name}/sci_dark
  --flat-path FLAT_PATH
                        Directory containing the flat FITS. Default:
                        {target_name}/flat
  --flat-dark-path FLAT_DARK_PATH
                        Directory containing the flat dark FITS. Default:
                        {target_name}/flat_dark
  -o OUTPUT, --output OUTPUT
                        Output cube FITS. Default:
                        {target_name}/cube_prepared.fits
  -O, --overwrite       Overwrite output if it already exists.
~~~

## Liste des fonctions CASTOR
* `castor_prepare` : Réduction des données et stockage dans un unique fichier FITS.

* `castor_align` : Alignement d'un cube d'images.

* `castor_rotate_spectra` : Rotation automatique d'un cube de spectres réduits.
L'angle de rotation peut être déterminé manuellement, ou automatiquement à l'aide d'une transformée de Radon sur les
spectres d'étalonnage (lampe ArNe).

* `castor_align_spectra` : Alignement vertical d'un cube de spectres retournés.

* `castor_wavelength_calibration` : Étalonnage en longueur d'onde des pixels à partir d'au moins deux points de référence.
Génère un fichier texte à deux colonnes contenant la correspondance de chaque pixel de l'axe horizontal en longueur d'onde.

* `castor_exoplanet_analysis` : Outil pour générer une courbe de transit d'exoplanète à partir de données brutes.

* `castor_pointing_analysis` : *Outils pour la mise en station de la monture, non utilisé pour les TP.*


**Pour plus d'information sur chaque fonction, consulter leur documentation comme indiqué ci-dessus.**

## Ordre d'utilisation
L'ordre des étapes du traitement étant important, la plupart des fonctions CASTOR prennent en entrée
un cube préalablement généré par une autre.
L'ordre d'appel des différentes fonctions doit donc être respecté.

**Note :** Remplacer `target_name` par l'adresse de votre répertoire contenant les données 
(`.` pour le répertoire courant).

### Imagerie (caméra CCD seule)
1. `castor_prepare [options] target_name`
2. `castor_align [options] target_name`

### Spectroscopie
1. `castor_prepare [options] target_name`
2. `castor_rotate_spectra [options] target_name`
3. `castor_align_spectra [options] target_name` **(Pas pour la Lune !)**
4. `castor_wavelength_calibration [options] target_name`

### Transit d'exoplanète
1. `castor_exoplanet_analysis [options] target_name`


## Lecture des données en Python
### Ouverture des cubes finaux
~~~python
from astropy.io import fits
f = fits.open('cube_prepared.fits')
f.info() # affiche une synthèse du contenu du fichier
hdu = f[0] # sélectionne le 1er HDU du FITS
data = hdu.data # tableau contenant les données
header = hdu.header # dictionnaire contenant l’en-tête
~~~

### Lecture du fichier texte d'étalonnage
~~~python
import numpy as np
px, wvl = np.loadtxt('wavelength_array.txt', unpack=True)
# Importation des deux colonnes du fichier texte en deux tableaux distincts
~~~

 * `px` : les indices de chaque pixel de l'image le long de l'axe spectral
 * `wvl` : les longueurs d'ondes associées

**Note :** Pour plus d'informations sur l'utilisation de Python, voir [*tuto_python_astro.pdf*](https://github.com/gpelouze/tuto_python_astro/releases/latest/download/tuto_python_astro.pdf) par G. Pelouze.


<!-- Documentation de cette documentation

Ce document a pour vocation à être converti en PDF. Le document est converti
par pandoc en passant par xelatex. À chaque fois qu’un nouveau tag est pushé
sur le repo, une action Github est déclanchée. Cette action (configurée dans
`.github./workflows/main.yml`) compile le PDF, puis l’enregistre avec la
release associée au tag. Le README (première page du repo) contient un lien
vers la dernière version du PDF.

La compilation du PDF est effectuée par `.github/workflows/entrypoint.sh`,
appelé par l’action Github. L’action échoue s’il y a un problème de
compilation, et il faut résoudre le problème (sinon le PDF accessible depuis le
README n’est pas à jour!). Il y a deux options:

1. Bidouiller `doc_TP_coupole.md` ou `entrypoint.sh`, commiter, faire un tag
   temporaire, puis pusher. Ça déclenche les actions Github, donc la
   compilation du PDF. Puis recommencer jusqu’à ce que le problème soit résolu…
   Je recommande très peu cette méthode, car elle génère plein de tags (et donc
   de releases) inutiles, sur Github ET sur PyPI. En plus, vous perdez pas mal
   de temps à commiter, tagger, pusher, et attendre le déclenchement de
   l’action.

2. Exécuter Docker en local, de sorte à faire *exactement* ce que fait Github,
   mais sans les étapes commit/tag/push/attente. Il suffit d’exécuter en local
   depuis la racine du repo:

   sudo docker run --volume "$(pwd):/data" --user $(id -u):$(id -g) --entrypoint "/data/.github/workflows/entrypoint.sh" pandoc/latex:2.9

   Les erreurs de compilation s’affichent dans le terminal. Le PDF généré (s’il
   y en a un) est dans `doc/doc_TP_coupole.pdf`. (Et le dossier `texmfhome`
   peut être supprimé après la compilation.)

   Une fois que la compilation se passe bien, commiter, tagger, et pusher pour
   générer et enregistrer le PDF sur Github.

-->
