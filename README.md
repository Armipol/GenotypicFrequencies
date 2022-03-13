# GenotypicFrequencies

Les différents codes ont été développés de sorte à fonctionner ainsi :

- Pour générer les données réelles dans un fichier :
  - Modifier les 4 filepaths du fichier ExcelGeneratorReal.py afin de mettre vos chemins d'accès aux données
  - Créer un document Excel dans le dossier qui comprend les codes avec le nom souhaité
  - Modifier le paramètre "excelName" de ExcelGeneratorReal.py en bas de page pour donner le nom du fichier excel créé
  - Indiquer dans le paramètre "numeroMelanges" les indices des mélanges à traiter (ex : [12] pour traiter uniquement le mélange Tm1012)
  - Lancer le programme

  
- Pour traiter ces données :
  - Dans TraitementExcelReal.py, modifier les 4 filepaths
  - Modifier le paramètre "excelName" vers le document souhaité
  - Indiquer dans "numeroMelanges" les indices des mélanges à traiter
  - Lancer le programme
  - Ouvrir ensuite le fichier excel pour visualiser les résultats


