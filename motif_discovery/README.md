### Setup
seqmotif.py uses python 2.7 and two packages that are not included in the standard library. These are numpy and matplotlib
and can be installed using PIP just running:

pip install matplotlib numpy

Depending on your python distribution you might need to install python TKinter. If you need so just run

sudo apt install python-tk


### Running seqmotif
In order to run the code just run

python seqmotif.py -f input_filename [ -l motif_length  -th threshold -o seqlogo.png]

If no -l or -o are added, seqmotif will assume length 19 and output a file called 'seqlogo.png'.
-th is used to accept an starting position in Z, and also to check if the results have converged.
If no threshold is added, seqmotif will use 0.01.

Made with love by

* Katharina Krabel
* Malte Steinke
* Juan Marin