
aNKflag v1.0
=============

Flagger is written by Apurba Bera, NCRA-TIFR, Pune, India

python3 wrapper is developed by Devojyoti Kansabanik, NCRA-TIFR, Pune, India

Requirements
============

aNKflag will only work in python3.6

aNKflag requires GNU gsl library, which is provided locally

aNKflag uses casa6, casa6 works only in python3.6 and need gfortran3.

Install aNKflag
================

python3 setup.py install --user

Usuage
=======

from aNKflag import runank   (Import runank from aNKflag module)

ankobj=runank.ANKFLAG() (Create aNKflag object)

ankobg.runank(params)  (Function to execute aNKflag; type "ankobj.runank?" in IPython shell to get details of the parameters)
 
