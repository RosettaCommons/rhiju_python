#!/usr/bin/python

from os import system
from sys import argv

for pdb in argv:
    system( 'java -cp ~rhiju/src/dangle/chiropraxis.jar chiropraxis.dangle.Dangle "alpha, beta, gamma, delta, epsilon, zeta" %s | suitename -chart ' % pdb )
