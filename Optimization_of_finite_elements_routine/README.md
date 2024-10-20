# ENSTA Paris - SIM203 - Projet _"Éléments finis discontinus"_

Consignes : [https://sim203.pages.math.cnrs.fr/](https://sim203.pages.math.cnrs.fr/)

## Objectifs pour le premier TD

* Utiliser les commandes expliquées dans ce fichier README
* Comprendre les grandes étapes du fichier `main.cpp`
* Faire le lien entre les formules du cours théorique et l’étape numéro 2 du fichier `main.cpp`
* Comprendre la signification des paramètres du fichier `setup`

## Commandes de compilation

1. Compiler le programme avec la commande `make`. La compilation crée l’exécutable _'dg'_ dans le répertoire principal. Une re-compilation avec `make` ne re-compilera que les fichiers modifiés.
2. Pour effacer les fichiers de compilation et l’exécutable : `make clean`

## Exécution d’un benchmark

1. Une fois le programme compilé, allez dans un dossier de benchmark : `cd benchmark/cube`. Ce dossier contient une géométrie au format _'gmsh'_ _('mesh.geo')_, un fichier de configuration de benchmark pour le programme _'dg'_ _('setup')_ et un script qui peut être utiliser pour visualiser la solution avec gmsh _('visu.script')_.
2. Avant l’exécution, générez le maillage : ``gmsh mesh.geo -3``
3. Exécutez le programme à partir du dossier du benchmark : `../../dg`
4. Pendant l'exécution, lisez ce qui s’affiche.
Attention, l’exécution doit se faire dans le dossier du benchmark !
Le programme _'dg'_ lit automatiquement le maillage _'mesh.msh'_ et le fichier _'setup'_.
5. Regardez le contenu du dossier _'output'_ à la fin de l’exécution.

## Visualisation des résultats

1. A la fin de l’exécution, téléchargez le maillage _'mesh.msh'_ et les fichiers de solution qui se trouvent dans le dossier _'benchmark/cube/output'_
2. Sur votre ordinateur, vous pouvez visualiser la solution avec le logiciel _gmsh_ : `gmsh mesh.msh output/*`
2. Pour une visualisation plus jolie, vous pouvez utiliser le script qui se trouve dans le dossier _'benchmark/cube'_ : `gmsh vizu.script`
