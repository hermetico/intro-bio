### Setup
f_score.py and slilhouett_value.py use both python 2.7 and numpy. Intro-bio-kmeans.jar uses java and maven to manage the dependencies, they can be found in the pom.xml file.

### Running intro-bio-kmeans.jar
```
java -jar intro-bio-kmeans.jar [options...] arguments...
 -f (--file) VAL     : The input File
 -i (--iterations) N : The numberof iterations (default: 20)
 -k N                : The number of k centroids (default: 3)
 -o (--output) VAL   : The output file name
```
### Running f_score.py
```
python f_score.py [-h] -f FILE -g GOLD

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  File containing results
  -g GOLD, --gold GOLD  File containing gold standards
```
### Running silhouette_value.py
```
python silhouette_value.py [-h] -f FILE -g GENES

optional arguments:
  -h, --help                show this help message and exit
  -f FILE, --file FILE      File containing results
  -g GENES, --genes GENES   File containing gene expressions

```

Made with love by

* Katharina Krabel
* Malte Steinke
* Juan Marin