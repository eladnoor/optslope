optslope
========

A toolbox for running OptSlope, an algorithm for metabolic engineering

Dependencies:
* cobrapy (latest from https://github.com/opencobra/cobrapy)
* libSBML 5
* numpy 1.7.1
* matplotlib 1.2.1
* pysvg 0.2.2
* PuLP 1.5.4
* IBM ILOG CPLEX 12.6

Specific instructions for Ubuntu 13.10:
* Run the following commands in the terminal:
  - sudo apt-get install libsbml5-python
  - sudo pip install pulp
  - sudo pip install pysvg
  - git clone https://github.com/opencobra/cobrapy ~/git/cobrapy
  - git clone https://github.com/eladnoor/optslope ~/git/optslope

* Make sure that the current directory is in PYTHONPATH.
  This can be achieved by adding this line to ~/.bashrc:
  export PYTHONPATH=.:$PYTHONPATH

* Install CPLEX (we recommend applying to IBM Academic Initiative)
