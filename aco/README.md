### Setup
aco.py uses python 2.7 and two packages that are not included in the standard library. These are numpy and matplotlib
and can be installed using PIP just running:

pip install matplotlib numpy

Depending on your python distribution you might need to install python TKinter. If you need so just run

sudo apt install python-tk


### Running aco
usage: aco.py [-h] -f FILE [-a ALPHA] [-b BETA] [-e EVAPORATION] [-n ANTS]
              [-v]

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  TSP file containing the data
  -a ALPHA, --alpha ALPHA
                        Alpha value should be >= 0
  -b BETA, --beta BETA  Beta value should be >= 1
  -e EVAPORATION, --evaporation EVAPORATION
                        Evaporation value 0-1
  -n ANTS, --ants ANTS  Num of ants
  -v, --visualization   Visualize routes on real time


Made with love by

* Katharina Krabel
* Malte Steinke
* Juan Marin