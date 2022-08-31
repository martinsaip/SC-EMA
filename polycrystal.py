__author__ = ['titrian', 'aydin', 'Martin Saip']
# Conversion to Python >= 3.10, corrections, improvements & refactoring done by Martin Saip in 2020 & 2021.

import copy, sys
from sympy import * #Symbol, nsolve, var, subs

from mpmath import *
from sympy.matrices import *
from numpy import linalg

#import crystalStack as st

EPS = 1.e-8

import numpy as np
import tools.utilities as u

#global crystallist, varlist

#crystallist = []
#varlist = []



def __call__(self, *args, **kwargs):
    return self.subs(kwargs)

class polycrystal:
    def __init__(self):

        self.f1 = None
        self.f2 = None
        self.conc = 1.

        self.K0 = None
        self.mue0 = None

#    def addCrystal(self, name, params):
#        crystallist.append(str(name) + str(params))

#    def removeCrystal(self, name, params):
#        for k in crystallist:
#            if k == str(name) + str(params):
#                crystallist.remove(k)

    def setYoungsMod(self):
#        """ This module calculates the youngs modulus E.
#            K0 and mue0 (bulk modulus and shear modulus) is
#            assumed to be known
#        """
        if self.mue0 is None or self.K0 is None:
            print("bulk modulus or shear modulus unknown")
            exit()
        self.E = 2 * self.mue0 * (1 + (3 * self.K0 - 2 * self.mue0)/(6 * self.K0 + 2 * self.mue0))

    def setpoissonratio(self):
    #        """ This module calculates the poisson's ratio v.
    #            K0 and mue0 (bulk modulus and shear modulus) is
    #            assumed to be known
    #        """
        if self.mue0 is None or self.K0 is None:
            print("poission's ratio unknown")
            exit()
        self.v = (3 * self.K0 - 2 * self.mue0) / (2 * (3 * self.K0 + self.mue0 ))

    def areamoduli(self):
        if self.EE1 is None:
            print("area moduli unknown")
            exit()
        if self.EE2 is None:
            print("area moduli unknown")
            exit()
        if self.EE3 is None:
            print("area moduli unknown")
            exit()
        if self.EE1 == 0:
            self.A1 = "not calculated"
        if self.EE2 == 0:
            self.A2 = "not calculated"
        if self.EE3 == 0:
            self.A3 = "not calculated"
        else:
            self.A1 = 1 / (self.EE1 + (3*self.K0)**-1)
            self.A2 = 1 / (self.EE2 + (3*self.K0)**-1)
            self.A3 = 1 / (self.EE3 + (3*self.K0)**-1)

    def __add__(self, other):
        newFunction = copy.copy(self)
        newFunction.f1 = self.f1 + other.f1
        newFunction.f2 = self.f2 + other.f2
        return newFunction

    def setK0andMue0(self):
        f = [lambda a, b: self.f1.subs(dict(K0=a, mue0=b)),
             lambda a, b: self.f2.subs(dict(K0=a, mue0=b))]

        k = 3; a = -1.; b = -1.;
        while a <= EPS or b <= EPS:
            a = 0.5; b = 0.5
            try:
                a, b = findroot(f, (k *0.1 , k *0.1 ))# TODO: we need an accurate initial guess
                #print (a, b)
            except:
                a, b = -1.,-1.
            k += 1

            self.K0 = a#a*100.
            self.mue0 = b#b*100.

    def setConc(self, conc):
        if conc > 1 or conc < 0:
            print("Concentration not in the range <0, 1> => exiting!")
            exit()
        self.conc = conc
        self.f1 = self.conc * self.local_1
        self.f2 = self.conc * self.local_2


class cubic(polycrystal):
    def __init__(self, C11 = None, C12 = None, C44 = None):
        """
        This is the cubic polycrystal class. Inherits from polycrystal class... Containing the mathematical functions for
        """
        if (C11 is None or C12 is None or C44 is None):
            #u.inputError("C parameters not set")
            exit()

        self.C11 = C11/100    # Thanks to Python 3, we do not need to explicitly state "float" in these 'C's
        self.C12 = C12/100
        self.C44 = C44/100

        correctInput, msg = self.checkCond()

        if correctInput == True:

            polycrystal.__init__(self)
            self.crystalname = "cubic"
            self.Cparamlist = (C11, C12, C44)
#            self.addCrystal(self.crystalname, self.Cparamlist)


            K0   = Symbol('K0')
            mue0 = Symbol('mue0')

            m_cubic = Matrix (([self.C11,self.C12,self.C12,      0,        0,       0],
                               [self.C12,self.C11,self.C12,      0,        0,       0],
                               [self.C12,self.C12,self.C11,      0,        0,       0],
                               [       0,       0,       0,self.C44,       0,       0],
                               [       0,       0,       0,       0,self.C44,       0],
                               [       0,       0,       0,       0,       0,self.C44]))

            S_cubic = linalg.inv(m_cubic)

            L11 = 0
            L12 = 0
            L13 = 1.

            L21 = 1 / sqrt(2)    # Thanks to Python 3, we do not need to explicitly state "float"
            L22 = L21
            L23 = 0

            L31 = 1 / sqrt(3)    # Thanks to Python 3, we do not need to explicitly state "float"
            L32 = L31
            L33 = L31

            self.EE1 = ( L11**4 + 2 * (L11**2) * (L12**2) * S_cubic [0,1] +  2 * (L11**2) * (L13**2) * S_cubic [0,2] + L12**4 * S_cubic [1,1]\
                         + 2 * (L12**2) * (L13**2) * S_cubic [1,2] +  L13**4 * S_cubic [2,2] + (L12**2) * (L13**2) * S_cubic [3,3]\
                         + (L11**2) * (L13**2) * S_cubic [4,4] + (L11**2) * (L12**2) * S_cubic [5,5] )

            self.EE2 = ( L21**4 + 2 * (L21**2) * (L22**2) * S_cubic [0,1] +  2 * (L21**2) * (L23**2) * S_cubic [0,2] + L22**4 * S_cubic [1,1]\
                         + 2 * (L22**2) * (L23**2) * S_cubic [1,2] +  L23**4 * S_cubic [2,2] + (L22**2) * (L23**2) * S_cubic [3,3]\
                         + (L21**2) * (L23**2) * S_cubic [4,4] + (L21**2) * (L22**2) * S_cubic [5,5] )

            self.EE3 = ( L31**4 + 2 * (L31**2) * (L32**2) * S_cubic [0,1] +  2 * (L31**2) * (L33**2) * S_cubic [0,2] + L32**4 * S_cubic [1,1]\
                         + 2 * (L32**2) * (L33**2) * S_cubic [1,2] +  L33**4 * S_cubic [2,2] + (L32**2) * (L33**2) * S_cubic [3,3]\
                         + (L31**2) * (L33**2) * S_cubic [4,4] + (L31**2) * (L32**2) * S_cubic [5,5] )


            mue = self.C44
            nue = 0.5 * (self.C11 - self.C12)
    #        K   = K0
            beta  = -3 * (K0 + 2 * mue0)/(5 * mue0 * (3 * K0 + 4 * mue0))
            coeff = 1/5    # Thanks to Python 3, we do not need to explicitly state "float"
            denom = (3 - ( self.C11 + 2 * self.C12 - 3 * K0))


            # these are the main functions
            self.local_1 = coeff*((2*nue - 2*mue0)/(1 - 2*beta*(nue - mue0)) + (3*mue - 3*mue0)/(1 - 2*beta*(mue - mue0)))
            #coeff * (1/(self.C11 - self.C12 - 2* mue0) - beta)**(-1)\
             #       + 3 * (1/(self.C44 - mue0) - 2*beta)**(-1)

            self.local_2 = (3 * (self.C11 + 2 * self.C12) - 9. * K0)/denom

            self.f1 = self.local_1
            self.f2 = self.local_2

            # -------------------------------  Upper & Lower Bound  ------------------------------- #
            #[NOT IN GPa yet]
            # BULK MODULUS
            # 1/ Voigt
            Kvoigt = (1/3. * (self.C11 + 2 * self.C12))   #Bulk Modulus

            # 2/ Reuss = Voigt
            Kreuss = Kvoigt

            # 3/ HS-
            KHSn = Kvoigt

            # 4/ HS+
            KHSp = Kvoigt

            # SHEAR MODULUS
    #        S11 = (self.C11 + self.C12) / ((self.C11 - self.C12) * (self.C11 + 2 * self.C12))
    #        S12 = -(self.C12)/ ((self.C11 - self.C12) * (self.C11 + 2 * self.C12))
    #        S44 = 1/self.C44
            G1 = 0.5*(self.C11 - self.C12)
            G2 = self.C44
            K = 1/3 *(self.C11 + 2*self.C12)
            # 1/ Voigt
    #        Gvoigt = (1/5.*(self.C11 - self.C12 + 3 * self.C44))*100.
            Gvoigt = 1/5.*( 2 * G1 + 3 * G2)

            # 2/ Reuss
    #        Greuss = 5 / (4 * (S11 - S12) + 3 * S44)
    #        Greuss = (5 * (self.C11 - self.C12) * self.C44 / 4 * self.C44 + 3 * (self.C11 - self.C12))*100    # Thanks to Python 3, we do not need to explicitly state "float"
            Greuss = (5 * G1 * G2) / (2 * G2 + 3 * G1)
            # 3/ HS-
            beta1 = - (3 * (K + 2 * G1)) / (5 * G1 * (3 * K + 4 * G1))
            GHSn = G1 + 3*((5/(G2-G1))- 4*beta1)**(-1)
            # 4/ HS+
            beta2 = - (3 * (K + 2 * G2)) / (5 * G2 * (3 * K + 4 * G2))
            GHSp = G2 + 2*((5/(G1-G2))- 6*beta2)**(-1)

            self.voigt_bulk   = Kvoigt
            self.reuss_bulk   = Kreuss
            self.HS_min_bulk  = KHSn
            self.HS_post_bulk = KHSp

            self.voigt_shear   = Gvoigt
            self.reuss_shear   = Greuss
            self.HS_min_shear  = GHSn
            self.HS_post_shear = GHSp


        else:
            print('warning' + '<br>')
            print(msg)
            exit()

    def getElasticDict(self):
        mydict = {'C11' : self.C11*100,    # Thanks to Python 3, we do not need to explicitly state "float" in this dict
                  'C12' : self.C12*100,
                  'C44' : self.C44*100}
        return mydict

    def checkCond(self):
        # Check condition
            m_cubic = Matrix (([self.C11,self.C12,self.C12,      0,        0,       0],
                               [self.C12,self.C11,self.C12,      0,        0,       0],
                               [self.C12,self.C12,self.C11,      0,        0,       0],
                               [       0,       0,       0,self.C44,       0,       0],
                               [       0,       0,       0,       0,self.C44,       0],
                               [       0,       0,       0,       0,       0,self.C44]))

            cubic_det         = det(m_cubic)
            cubic_condition_1 = self.C11 + 2 * self.C12
            cubic_condition_2 = self.C44
            cubic_condition_3 = self.C11 - self.C12

            if cubic_condition_1 > 0:
                pass
            else:
                #u.inputError('Please enter right value: C11 + 2 C12 > 0')
                return False, 'Please enter right value: C11 + 2 C12 > 0 (cubic phase).'
            if cubic_condition_2 > 0:
                pass
            else:
                #u.inputError('Please enter right value: C44 > 0')
                return False, 'Please enter right value: C44 > 0 (cubic phase).'

            if cubic_condition_3 > 0:
                pass
            else:
                #u.inputError('Please enter right value: C11 - C12 > 0')
                return False, 'Please enter right value: C11 - C12 > 0 (cubic phase).'

            return True, ''


class tetragonal(polycrystal):
    def __init__(self, C11 = None, C12 = None, C13 = None, C33 = None, C44 = None, C66 = None):

        if (C11 is None or C12 is None or C13 is None or C33 is None or C44 is None or C66 is None):
            u.inputError("C parameters not set!")
            exit()

        self.C11 = C11/100    # Thanks to Python 3, we do not need to explicitly state "float" in these 'C's
        self.C12 = C12/100
        self.C13 = C13/100
        self.C33 = C33/100
        self.C44 = C44/100
        self.C66 = C66/100

        correctInput, msg = self.checkCond()

        if correctInput == True:

            polycrystal.__init__(self)
    #       self.conc = 1
            self.crystalname = "tetragonal"
            self.Cparamlist = (C11, C12, C13, C33, C44, C66)
#            self.addCrystal(self.crystalname, self.Cparamlist)


            K0   = Symbol('K0')
            mue0 = Symbol('mue0')

            m_tetragonal = Matrix (([self.C11,self.C12,self.C13,       0,       0,       0],
                                    [self.C12,self.C11,self.C13,       0,       0,       0],
                                    [self.C13,self.C13,self.C33,       0,       0,       0],
                                    [       0,       0,       0,self.C44,       0,       0],
                                    [       0,       0,       0,       0,self.C44,       0],
                                    [       0,       0,       0,       0,       0,self.C66]))

            S_tetragonal = linalg.inv(m_tetragonal)
            #        print m_tetragonal
            #        print S_tetragonal
            #        print S_tetragonal [0,5]

            L11 = 0
            L12 = 0
            L13 = 1.

            L21 = 1 / sqrt(2)    # Thanks to Python 3, we do not need to explicitly state "float"
            L22 = L21
            L23 = 0

            L31 = 1 / sqrt(3)    # Thanks to Python 3, we do not need to explicitly state "float"
            L32 = L31
            L33 = L31

            self.EE1 = ( L11**4 + 2 * (L11**2) * (L12**2) * S_tetragonal [0,1] +  2 * (L11**2) * (L13**2) * S_tetragonal [0,2] + L12**4 * S_tetragonal [1,1]\
                         + 2 * (L12**2) * (L13**2) * S_tetragonal [1,2] +  L13**4 * S_tetragonal [2,2] + (L12**2) * (L13**2) * S_tetragonal [3,3]\
                         + (L11**2) * (L13**2) * S_tetragonal [4,4] + (L11**2) * (L12**2) * S_tetragonal [5,5] )

            self.EE2 = ( L21**4 + 2 * (L21**2) * (L22**2) * S_tetragonal [0,1] +  2 * (L21**2) * (L23**2) * S_tetragonal [0,2] + L22**4 * S_tetragonal [1,1]\
                         + 2 * (L22**2) * (L23**2) * S_tetragonal [1,2] +  L23**4 * S_tetragonal [2,2] + (L22**2) * (L23**2) * S_tetragonal [3,3]\
                         + (L21**2) * (L23**2) * S_tetragonal [4,4] + (L21**2) * (L12**2) * S_tetragonal [5,5] )

            self.EE3 = ( L31**4 + 2 * (L31**2) * (L32**2) * S_tetragonal [0,1] +  2 * (L31**2) * (L33**2) * S_tetragonal [0,2] + L32**4 * S_tetragonal [1,1]\
                         + 2 * (L32**2) * (L33**2) * S_tetragonal [1,2] +  L33**4 * S_tetragonal [2,2] + (L32**2) * (L33**2) * S_tetragonal [3,3]\
                         + (L31**2) * (L33**2) * S_tetragonal [4,4] + (L31**2) * (L32**2) * S_tetragonal [5,5] )

            Knue     = 1./9 * ( self.C33 + 2 * (self.C11 + self.C12) + 4 * self.C13)
            M        = self.C11 + self.C12 + 2 * self.C33 - 4 * self.C13
            C2       = self.C33 * (self.C11 + self.C12) - 2 * self.C13**2
            psi      = self.C11 + self.C12 + self.C33 - 3 * K0 - 2 * mue0
            deltakk  = C2 - K0 * ( M - 6 * mue0 ) - 6 * mue0 * Knue
            sigma    = 9 * Knue- 9./2 * K0 - 6 * mue0 + 3./2 * self.C33 - 6 * self.C13
            beta     = -3 * (K0 + 2 * mue0)/(5 * mue0 * (3 * K0 + 4 * mue0))
            eta      = -3./(3 * K0 + 4 * mue0)
            gamma    = 1./9 * (eta - 3 * beta)

            denom1 = (1 - beta * psi - 9 * gamma *(Knue - K0) + 1./3 * eta * beta * deltakk)
            denom2 = 2 *(1 - beta *(self.C11 - self.C12 -2 * mue0))
            coeff = 1 / 15    # Thanks to Python 3, we do not need to explicitly state "float"

            fone = M - 6 * mue0
            ftwo = (self.C11 - self.C12 - 2 * mue0)*(3 - 2 * beta * sigma - 27 * gamma * (Knue - K0) - 2 * eta * beta * deltakk) - eta * deltakk

            fa = (fone + ftwo) / (denom2 * denom1)
            fb = 6 * (self.C44 - mue0) / (1 - 2 * beta * (self.C44 - mue0))
            fc = 3 * (self.C66 - mue0) / (1 - 2 * beta * (self.C66 - mue0))

            self.local_1 = coeff *( fa + fb + fc )
            self.local_2 = (3 * (Knue - K0) - beta * deltakk)/denom1 #/(1 - beta * psi - 9 * gamma *(Knue -K0)+ 1./3 * eta * beta * deltakk)

            self.f1 = self.local_1
            self.f2 = self.local_2

            # -------------------------------  Upper & Lower Bound  ------------------------------- #

            # Bulk Modulus
            G3 = 0.5 * (self.C11 - self.C12)

            G_low  = [self.C44, self.C66, G3, C2 / (6* Knue)]
            G1 = max(G_low)

            G_up = [self.C44, self.C66, G3, (1/6. * M)]
            G2 = min(G_up)
            if abs(M - 6 * G2) < EPS: # Avoid division by zero. Take the next larger number
                G_up.sort()
                G2 = G_up[1]

            K1 = (C2 - 6 * G1 * Knue) / (M * 6 * G1)
            K2 = (C2 - 6 * G2 * Knue) / (M * 6 * G2)

            beta1 = -3 * (K1 + 2 * G1) / (5 * G1 * (3 * K1 + 4 * G1))
            beta2 = -3 * (K2 + 2 * G2) / (5 * G2 * (3 * K2 + 4 * G2))

            delta1 = C2 - K1*(M - 6*G1) - 6*G1*Knue
            delta2 = C2 - K2*(M - 6*G2) - 6*G2*Knue


            # 1/ Voigt
            Kvoigt = 1/9. * (2 * (self.C11 + self.C12) + self.C33 + 4 * self.C13)

            # 2/ Reuss
            Kreuss = C2 / M

            # 3/ HS-
            KHSn = K1 + ((Knue - K1) -1/3. * beta1 * delta1 )  / ( 1 - 1/3.*beta1 * (M - 6*G1))

            # 4/ HS+
            KHSp = K2 + ((Knue - K2) -1/3. * beta2 * delta2 )  / ( 1 - 1/3.*beta2 * (M - 6*G2))


            # Shear Modulus

            eta1      = -3/(3 * K1 + 4 * G1)    # Thanks to Python 3, we do not need to explicitly state "float" in these
            eta2      = -3/(3 * K2 + 4 * G2)

            gamma1    = 1/9 * (eta1 - 3 * beta1)
            gamma2    = 1/9 * (eta2 - 3 * beta2)


            sigma1    = 9 * Knue- 9./2 * K1 - 6 * G1 + 3./2 * self.C33 - 6 * self.C13
            sigma2    = 9 * Knue- 9./2 * K2 - 6 * G2 + 3./2 * self.C33 - 6 * self.C13

            psi1      = self.C11 + self.C12 + self.C33 - 3 * K1 - 2 * G1
            psi2      = self.C11 + self.C12 + self.C33 - 3 * K2 - 2 * G2

            a1 = (M - 6 * G1 + (self.C11 - self.C12 - 2*G1) * (3 - 2 * beta1 * sigma1 - 27 * gamma1 * (Knue - K1) - 2 * eta1 * beta1 * delta1) - eta1 * delta1) /\
                 (2 * (1 - beta1 * (self.C11 - self.C12 - 2 * G1)) * (1 - beta1 * psi1 - 9 * gamma1 * (Knue - K1) + 1/3. * eta1 * beta1 * delta1))

            a2 = (M - 6 * G2 + (self.C11 - self.C12 - 2*G2) * (3 - 2 * beta2 * sigma2 - 27 * gamma2 * (Knue - K2) - 2 * eta2 * beta2 * delta2) - eta2 * delta2) /\
                 (2 * (1 - beta2 * (self.C11 - self.C12 - 2 * G2)) * (1 - beta2 * psi2 - 9 * gamma2 * (Knue - K2) + 1/3. * eta2 * beta2 * delta2))

            b1 = 6 * (C44 - G1) / (1 - 2 * beta1 * (C44 - G1))
            b2 = 6 * (C44 - G2) / (1 - 2 * beta2 * (C44 - G2))

            c1 = 3 * (C66 - G1) / (1 - 2 * beta1 * (C66 - G1))
            c2 = 3 * (C66 - G2) / (1 - 2 * beta2 * (C66 - G2))

            B21 = 1/15 * (a1 + b1 + c1)    # Thanks to Python 3, we do not need to explicitly state "float" in these 'B's
            B22 = 1/15 * (a2 + b2 + c2)

            # 1/ Voigt
            Gvoigt = 1/30 * (M + 3 * self.C11 - 3 * self.C12 + 12 * self.C44 + 6 * self.C66)    # Thanks to Python 3, we do not need to explicitly state "float"

            # 2/ Reuss
            Greuss = 15 * (18 * Knue /C2 + 6 / (self.C11 - self.C12) + 6 / self.C44 + 3 / self.C66)**-1

            # 3/ HS-
            GHSp = G1 + B21 / (1 + 2 * beta1 * B21)

            # 4/ HS+
            GHSn = G2 + B22 / (1 + 2 * beta2 * B22)

            self.voigt_bulk   = Kvoigt
            self.reuss_bulk   = Kvoigt
            self.HS_min_bulk  = KHSn
            self.HS_post_bulk = KHSp

            self.voigt_shear   = Gvoigt
            self.reuss_shear   = Greuss
            self.HS_min_shear  = KHSn
            self.HS_post_shear = KHSp

        else:
            print('warning' + '<br>')
            print(msg)
            exit()

    def getElasticDict(self):
        mydict = {'C11' : self.C11*100,    # Thanks to Python 3, we do not need to explicitly state "float" in these 'C's
                  'C12' : self.C12*100,
                  'C13' : self.C13*100,
                  'C33' : self.C33*100,
                  'C44' : self.C44*100,
                  'C66' : self.C66*100}
        return mydict

    def checkCond(self):
        # Check condition
        m_tetragonal = Matrix (([self.C11,self.C12,self.C13,       0,       0,       0],
                                [self.C12,self.C11,self.C13,       0,       0,       0],
                                [self.C13,self.C13,self.C33,       0,       0,       0],
                                [       0,       0,       0,self.C44,       0,       0],
                                [       0,       0,       0,       0,self.C44,       0],
                                [       0,       0,       0,       0,       0,self.C66]))

        tetragonal_det         = det(m_tetragonal)

        tetragonal_condition_1 = self.C11

        tetragonal_condition_2 = det(Matrix(([self.C11,self.C12],
                                             [self.C12,self.C11])))

        tetragonal_condition_3 = det(Matrix(([self.C11,self.C12,self.C13],
                                             [self.C12,self.C11,self.C13],
                                             [self.C13,self.C13,self.C33])))

        tetragonal_condition_4 = self.C44
        tetragonal_condition_5 = self.C66

        if tetragonal_condition_1 > 0:
            pass
        else:
    #            u.inputError('Please enter right value: C11 > 0')
            return False, 'Please enter such a value that C11 > 0 (tetragonal phase).'
        if tetragonal_condition_2 > 0:
            pass
        else:
    #            u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12), (C12, C11) )')
            return False, 'Please enter such a value that det (matrix A) > 0!A =  [(C11, C12),(C12, C11)] (tetragonal phase).'
        if tetragonal_condition_3 > 0:
            pass
        else:
    #            u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13), (C12, C11, C13), (C13, C13, C33) )')
            return False, 'Please enter such a value that det (matrix A) > 0! A = [(C11, C12, C13), (C12, C11, C13), (C13, C13, C33)] (tetragonal phase).'
        if tetragonal_condition_4 > 0:
            pass
        else:
    #            u.inputError('Please enter right value: C44 > 0')
            return False, 'Please enter such a value that C44 > 0 (tetragonal phase).'
        if tetragonal_condition_5 > 0:
            pass
        else:
    #            u.inputError('Please enter right value: C66 > 0')
            return False, 'Please enter such a value that C66 > 0 (tetragonal phase).'

        return True, ''

class trigonal(polycrystal):
    def __init__(self, C11 = None, C12 = None, C13 = None, C14 = None, C33 = None, C44 = None):

        if (C11 is None or C12 is None or C13 is None or C14 is None or C33 is None or C44 is None):
#            u.inputError("C parameters not set")
            exit()

        self.C11 = C11/100    # Thanks to Python 3, we do not need to explicitly state "float" in these 'C's
        self.C12 = C12/100
        self.C13 = C13/100
        self.C14 = C14/100
        self.C33 = C33/100
        self.C44 = C44/100

        correctInput, msg = self.checkCond()

        if correctInput == True:

            polycrystal.__init__(self)
            self.crystalname = "trigonal"
            self.Cparamlist = (C11, C12, C13, C14, C33, C44)
#            self.addCrystal(self.crystalname, self.Cparamlist)

            self.EE1 = 0
            self.EE2 = 0
            self.EE3 = 0

            K0   = Symbol('K0')
            mue0 = Symbol('mue0')

            Knue     = 1./9 * ( self.C33 + 2 * (self.C11 + self.C12) + 4 * self.C13)
            M        = self.C11 + self.C12 + 2 * self.C33 - 4 * self.C13
            C2       = self.C33 * (self.C11 + self.C12) - 2 * self.C13**2
            psi      = self.C11 + self.C12 + self.C33 - 3 * K0 - 2 * mue0
            deltakk  = C2 - K0 * ( M - 6 * mue0 ) - 6 * mue0 * Knue
            C66      = 0.5 * (self.C11 - self.C12)
            chi      = 4 * (self.C14**2 - (self.C44 - mue0) - (C66 - mue0))
            beta     = -3 * (K0 + 2 * mue0) / (5 * mue0 * (3 * K0 + 4 * mue0))
            eta      = -3./(3 * K0 + 4 * mue0)
            gamma    = 1./9 * (eta - 3 * beta)

            denom1 = (1 - beta * psi - 9 * gamma *(Knue -K0) + 1./3 * eta * beta * deltakk)
            denom2 = (1 - 2 * beta * (self.C44 + C66 - 2 * mue0) - beta**2 * chi)
            coeff = 1 / 30.

            fa = (M - 6 * mue0 - eta * deltakk)/denom1
            fb = 12 * (self.C44 + C66 - 2 * mue0 + beta * chi)/denom2

            self.local_1 = coeff *( fa + fb )
            self.local_2 = (3 * (Knue - K0) - beta * deltakk)/denom1

            self.f1 = self.local_1
            self.f2 = self.local_2

            # --------------------  Upper and Lower Bound  ----------------------------- #

            # Bulk Modulus
            i_G1 = 0.5 * (self.C44 + C66) - (0.25 * (self.C44 - C66)**2 + self.C14**2)**0.5
            i_G2 = 0.5 * (self.C44 + C66) + (0.25 * (self.C44 - C66)**2 + self.C14**2)**0.5

            G_low  = [i_G1, C2 / (6* Knue)]
            G1 = max(G_low)

            G_up = [i_G2,(1/6. * M)]
            G2 = min(G_up)
            if abs(M - 6 * G2) < EPS: # Avoid division by zero. Take the next larger number
                G_up.sort()
                G2 = G_up[1]

            K1 = (C2 - 6 * G1 * Knue) / (M * 6 * G1)
            K2 = (C2 - 6 * G2 * Knue) / (M * 6 * G2)

            beta1 = -3 * (K1 + 2 * G1) / (5 * G1 * (3 * K1 + 4 * G1))
            beta2 = -3 * (K2 + 2 * G2) / (5 * G2 * (3 * K2 + 4 * G2))

            delta1 = C2 - K1 * (M - 6 * G1) - 6 * G1 * Knue
            delta2 = C2 - K2 * (M - 6 * G2) - 6 * G2 * Knue

            # 1/ Voigt
            Kvoigt = 1/9. * (2 * (self.C11 + self.C12) + self.C33 + 4 * self.C13)

            # 2/ Reuss
            Kreuss = C2 / M

            # 3/ HS-
            KHSn = K1 + ((Knue - K1) - 1/3. * beta1 * delta1 )  / ( 1 - 1/3. * beta1 * (M - 6 * G1))

            # 4/ HS+
            KHSp = K2 + ((Knue - K2) - 1/3. * beta2 * delta2 )  / ( 1 - 1/3. * beta2 * (M - 6 * G2))

            # Shear Modulus

            eta1      = -3/(3 * K1 + 4 * G1)    # Thanks to Python 3, we do not need to explicitly state "float" in these
            eta2      = -3/(3 * K2 + 4 * G2)

            gamma1    = 1/9 * (eta1 - 3 * beta1)
            gamma2    = 1/9 * (eta2 - 3 * beta2)

            omega1    = 4 * (self.C14**2 - (self.C44 - G1) * (C66 - G1))
            omega2    = 4 * (self.C14**2 - (self.C44 - G2) * (C66 - G2))

            a1 = (M - 6 * G1 - eta1 * delta1) / (1 - beta1 * (self.C11 + self.C12 + self.C33 - 3*K1 - 2*G1) - 9 * gamma1 * (Knue - K1) + 1/3. * eta1 * beta1 * delta1)

            a2 = (M - 6 * G2 - eta2 * delta2) / (1 - beta2 * (self.C11 + self.C12 + self.C33 - 3*K2 - 2*G2) - 9 * gamma2 * (Knue - K2) + 1/3. * eta2 * beta2 * delta2)
            b1  = 12 * (self.C44 + C66 - 2 * G1 + beta1 * omega1) / (1- 2 * beta1 * (self.C44 + C66 - 2 * G1) - (beta1**2 * omega1))
            b2  = 12 * (self.C44 + C66 - 2 * G2 + beta2 * omega2) / (1- 2 * beta2 * (self.C44 + C66 - 2 * G2) - (beta2**2 * omega2))

            B21 = 1/30 * (a1 + b1)    # Thanks to Python 3, we do not need to explicitly state "float" in these
            B22 = 1/30 * (a2 + b2)

            # 1/ Voigt
            Gvoigt = 1/30* (M + 12 * self.C44 + 12 * C66)

            # 2/ Reuss
            Greuss = 5/2 *  (C2 * (self.C44 * C66 - self.C14**2) / (3* Knue * (self.C44 * C66 - self.C14**2) + C2 * (self.C44 + C66)))

            # 3/ HS-
            GHSp = G1 + B21 / (1 + 2 * beta1 * B21)

            # 4/ HS+
            GHSn = G2 + B22 / (1 + 2 * beta2 * B22)

            D = [self.C44 , C66]
            D1 = min (D)
            D2 = max (D)

            if G1 < D2:
#                print 'D1 accepted'
                pass
            else:
                print('D1 is False --> Please enter a value corresponding to a mechanically stable system!')
            if G2 > D1:
#                print 'D2 accepted'
                pass
            else:
                print('D2 is False --> Please enter a value corresponding to a mechanically stable system!')

            self.voigt_bulk   = Kvoigt
            self.reuss_bulk   = Kvoigt
            self.HS_min_bulk  = KHSn
            self.HS_pos_bulk  = KHSp

            self.voigt_shear   = Gvoigt
            self.reuss_shear   = Greuss
            self.HS_min_shear  = KHSn
            self.HS_pos_shear  = KHSp
        else:
            print('warning' + '<br>')
            print(msg)
            exit()

    def getElasticDict(self):
        mydict = {'C11' : self.C11*100,    # Thanks to Python 3, we do not need to explicitly state "float" in these 'C's
                  'C12' : self.C12*100,
                  'C13' : self.C13*100,
                  'C14' : self.C14*100,
                  'C33' : self.C33*100,
                  'C44' : self.C44*100}
        return mydict

    def checkCond(self):
            # Check condition
            C66      = 0.5 * (self.C11 - self.C12)

            m_trigonal = Matrix (([self.C11, self.C12,self.C13, self.C14,       0,  0],
                                  [self.C12, self.C11,self.C13,-self.C14,       0,  0],
                                  [self.C13, self.C13,self.C33,        0,       0,  0],
                                  [self.C14,-self.C14,       0, self.C44,       0,  0],
                                  [       0,       0,       0,         0,self.C44,  self.C14/2.],
                                  [       0,       0,       0,         self.C14/2.,       0,C66]))

            trigonal_det         = det(m_trigonal)

            trigonal_condition_1 = self.C11

            trigonal_condition_2 = det(Matrix(([self.C11, self.C12],
                                               [self.C12, self.C11])))

            trigonal_condition_3 = det(Matrix(([self.C11, self.C12,self.C13],
                                               [self.C12, self.C11,self.C13],
                                               [self.C13, self.C13,self.C33])))

            trigonal_condition_4 = det(Matrix(([self.C11, self.C12,self.C13, self.C14],
                                               [self.C12, self.C11,self.C13,-self.C14],
                                               [self.C13, self.C13,self.C33,        0],
                                               [self.C14,-self.C14,       0, self.C44])))

            trigonal_condition_5 = det(Matrix(([self.C11, self.C12,self.C13, self.C14,       0],
                                               [self.C12, self.C11,self.C13,-self.C14,       0],
                                               [self.C13, self.C13,self.C33,        0,       0],
                                               [self.C14,-self.C14,       0, self.C44,       0],
                                               [       0,        0,       0,        0,self.C44])))

            trigonal_condition_6 = det(m_trigonal)

            if trigonal_condition_1 > 0:
                pass
            else:
#                u.inputError('Please enter right value: C11 > 0')
                return False, 'Please enter such a value that C11 > 0 (trigonal phase).'
            if trigonal_condition_2 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12), (C12, C11) )')
                return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12), (C12, C11)] (trigonal phase).'
            if trigonal_condition_3 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13), (C12, C11, C13), (C13, C13, C33) )')
                return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13), (C12, C11, C13), (C13, C13, C33)] (trigonal phase).'
            if trigonal_condition_4 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13, C14), (C12, C11, C13, -C14), (C13, C13, C33, 0), (C14, -C14, 0, C44) )')
                return False, 'Please enter such a value that det (matrix A) > 0! A = [(C11, C12, C13, C14), (C12, C11, C13, -C14), (C13, C13, C33, 0), (C14, -C14, 0, C44)] (trigonal phase).'
            if trigonal_condition_5 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13, C14, 0), (C12, C11, C13, -C14, 0), (C13, C13, C33, 0, 0), (C14, -C14, 0, C44, 0), (0, 0, 0, 0, C44) )')
                return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13, C14, 0), (C12, C11, C13, -C14, 0), (C13, C13, C33, 0, 0), (C14, -C14, 0, C44, 0), (0, 0, 0, 0, C44)] (trigonal phase).'
            if trigonal_condition_6 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13, C14, 0, 0), (C12, C11, C13, -C14, 0, 0), (C13, C13, C33, 0, 0, 0), (C14, -C14, 0, C44, 0, 0), (0, 0, 0, 0, C44, 0), (0, 0, 0, 0, 0, C66) )')
                return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13, C14, 0, 0), (C12, C11, C13, -C14, 0, 0), (C13, C13, C33, 0, 0, 0), (C14, -C14, 0, C44, 0, 0), (0, 0, 0, 0, C44, 0), (0, 0, 0, 0, 0, C66)] (trigonal phase).'

            return True, ''


class hexagonal(polycrystal):
    def __init__(self, C11 = None, C12 = None, C13 = None, C33 = None, C55 = None):

        if (C11 is None or C12 is None or C13 is None or C33 is None or C55 is None):
#            u.inputError("C parameters not set")
            exit()

        self.C11 = C11/100    # Thanks to Python 3, we do not need to explicitly state "float" in these 'C's
        self.C12 = C12/100
        self.C13 = C13/100
        self.C33 = C33/100
        self.C55 = C55/100

        correctInput, msg = self.checkCond()

        if correctInput == True:

            polycrystal.__init__(self)
            self.crystalname = "hexagonal"
            self.Cparamlist = (C11, C12, C13, C33, C55)
#            self.addCrystal(self.crystalname, self.Cparamlist)

            K0   = Symbol('K0')
            mue0 = Symbol('mue0')

            self.EE1 = 0
            self.EE2 = 0
            self.EE3 = 0

            Knue     = 1./9 * ( self.C33 + 2 * (self.C11 + self.C12) + 4 * self.C13)
            M        = self.C11 + self.C12 + 2 * self.C33 - 4 * self.C13
            C2       = self.C33 * (self.C11 + self.C12) - 2 * self.C13**2
            psi      = self.C11 + self.C12 + self.C33 - 3 * K0 - 2 * mue0
            deltakk  = C2 - K0 * ( M - 6 * mue0 ) - 6 * mue0 * Knue
            C66      = 0.5 * (self.C11 - self.C12)
            C44      = self.C55
            beta     = -3 * (K0 + 2 * mue0) / (5 * mue0 * (3 * K0 + 4 * mue0))
            eta      = -3./(3 * K0 + 4 * mue0)
            gamma    = 1./9 * (eta - 3 * beta)

            denom = (1 - beta * psi - 9 * gamma *(Knue -K0) + 1./3 * eta * beta * deltakk)
            coeff = 1 / 30    # Thanks to Python 3, we do not need to explicitly state "float"

            fa = (M - 6 * mue0 - eta * deltakk)/(1 - beta * psi - 9 * gamma * (Knue - K0) + 1./3 * eta * beta * deltakk)
            fb = 12 * ((C44 - mue0) / (1 - 2 * beta * (C44 - mue0))
                     + (C66 - mue0) / (1 - 2 * beta * (C66 - mue0)))

            self.local_1 = self.conc * coeff *( fa + fb )
            self.local_2 = self.conc * (3 * (Knue - K0) - beta * deltakk)/denom

            self.f1 = self.local_1
            self.f2 = self.local_2


            # -------------------------------  Upper & Lower Bound  ------------------------------- #
            # Bulk Modulus
            G_low  = [C44, C66, C2 / (6* Knue)]
            G1 = max(G_low)

            G_up = [C44, C66, (1/6 * M)]    # Thanks to Python 3, we do not need to explicitly state "float"
            G2 = min(G_up)
            if abs(M - 6 * G2) < EPS: # Avoid division by zero. Take the next larger number
                G_up.sort()
                G2 = G_up[1]

            K1 = (C2 - 6 * G1 * Knue) / (M - 6 * G1)
            K2 = (C2 - 6 * G2 * Knue) / (M - 6 * G2)

            beta1 = -3 * (K1 + 2 * G1) / (5 * G1 * (3 * K1 + 4 * G1))
            beta2 = -3 * (K2 + 2 * G2) / (5 * G2 * (3 * K2 + 4 * G2))

            delta1 = C2 - K1 * (M - 6 * G1) - 6 * G1 * Knue
            delta2 = C2 - K2 * (M - 6 * G2) - 6 * G2 * Knue

            # 1/ Voigt
            Kvoigt = (1/9. * (2 * (self.C11 + self.C12) + self.C33 + 4 * self.C13))

            # 2/ Reuss
            Kreuss = (C2 / M)

            # 3/ HS-
            KHSn = (K1 + ((Knue - K1) - 1/3. * beta1 * delta1 )  / ( 1 - 1/3.*beta1 * (M - 6*G1)))

            # 4/ HS+
            KHSp = (K2 + ((Knue - K2) - 1/3. * beta2 * delta2 )  / ( 1 - 1/3.*beta2 * (M - 6*G2)))

           # Shear Modulus

            eta1      = -3/(3 * K1 + 4 * G1)    # Thanks to Python 3, we do not need to explicitly state "float" in these
            eta2      = -3/(3 * K2 + 4 * G2)

            gamma1    = 1/9 * (eta1 - 3 * beta1)
            gamma2    = 1/9 * (eta2 - 3 * beta2)

            a1 = M - 6 * G1 - eta1 * delta1 / (1 - (beta1 * (self.C11 + self.C12 + self.C33 - 3*K1 - 2*G1)) - (9*gamma1*(Knue - K1)) + (1/3. * eta1 * beta1 * delta1))

            a2 = M - 6 * G2 - eta2 * delta2 / (1 - (beta2 * (self.C11 + self.C12 + self.C33 - 3*K2 - 2*G2)) - (9*gamma2*(Knue - K2)) + (1/3. * eta2 * beta2 * delta2))

            b1 = 12 * (C44 - G1) / (1 - 2 * beta1 * (C44 - G1))
            b2 = 12 * (C44 - G2) / (1 - 2 * beta2 * (C44 - G2))

            c1 = 12 * (C66 - G1) / (1 - 2 * beta1 * (C66 - G1))
            c2 = 12 * (C66 - G2) / (1 - 2 * beta2 * (C66 - G2))

            B21 = 1/30 * (a1 + b1 + c1)    # Thanks to Python 3, we do not need to explicitly state "float" in these 'B's
            B22 = 1/30 * (a2 + b2 + c2)

            # 1/ Voigt
            Gvoigt = (1/30.* (M + 12 * C44 + 12 * C66))

            # 2/ Reuss
            Greuss = (5/2. *  (C2 * C44 * C66) / (3* Knue * C44 * C66 + C2 * (C44 + C66)))

            # 3/ HS-
            GHSn = (G1 + B21 / (1 + 2 * beta1 * B21))

            # 4/ HS+
            GHSp = (G2 + B22 / (1 + 2 * beta2 * B22))

            self.voigt_bulk   = Kvoigt
            self.reuss_bulk   = Kvoigt
            self.HS_min_bulk  = KHSn
            self.HS_pos_bulk  = KHSp

            self.voigt_shear   = Gvoigt
            self.reuss_shear   = Greuss
            self.HS_min_shear  = KHSn
            self.HS_pos_shear  = KHSp

        else:
            print('warning' + '<br>')
            print(msg)
            exit()

    def getElasticDict(self):
        mydict = {'C11' : self.C11*100,    # Thanks to Python 3, we do not need to explicitly state "float" in these 'C's
                  'C12' : self.C12*100,
                  'C13' : self.C13*100,
                  'C33' : self.C33*100,
                  'C55' : self.C55*100}
        return mydict

    def checkCond(self):
            # Check condition
            C66      = 0.5 * (self.C11 - self.C12)

            m_hexagonal = Matrix (([self.C11,self.C12,self.C13,       0,       0,  0],
                                   [self.C12,self.C11,self.C13,       0,       0,  0],
                                   [self.C13,self.C13,self.C33,       0,       0,  0],
                                   [       0,       0,       0,self.C55,       0,  0],
                                   [       0,       0,       0,       0,self.C55,  0],
                                   [       0,       0,       0,       0,       0,C66]))

            hexagonal_det         = det(m_hexagonal)

            hexagonal_condition_1 = self.C11

            hexagonal_condition_2 = det(Matrix(([self.C11,self.C12],
                                                [self.C12,self.C11])))

            hexagonal_condition_3 = det(Matrix(([self.C11,self.C12,self.C13],
                                                [self.C12,self.C11,self.C13],
                                                [self.C13,self.C13,self.C33])))

            hexagonal_condition_4 = det(Matrix(([self.C11,self.C12,self.C13,       0],
                                                [self.C12,self.C11,self.C13,       0],
                                                [self.C13,self.C13,self.C33,       0],
                                                [       0,       0,       0,self.C55])))

            hexagonal_condition_5= det(Matrix(([self.C11,self.C12,self.C13,       0,       0],
                                               [self.C12,self.C11,self.C13,       0,       0],
                                               [self.C13,self.C13,self.C33,       0,       0],
                                               [       0,       0,       0,self.C55,       0],
                                               [       0,       0,       0,       0,self.C55])))

            hexagonal_condition_6 = det(m_hexagonal)

            if hexagonal_condition_1 > 0:
                pass
            else:
#                u.inputError('Please enter right value: C11 > 0')
                return False, 'Please enter such a value that C11 > 0 (hexagonal phase).'
            if hexagonal_condition_2 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12), (C12, C11) )')
                return False, 'Please enter such a value that det (matrix A) > 0! A = [(C11, C12), (C12, C11)] (hexagonal phase).'
            if hexagonal_condition_3 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13), (C12, C11, C13), (C13, C13, C33) )')
                return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13), (C12, C11, C13), (C13, C13, C33)] (hexagonal phase).'
            if hexagonal_condition_4 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13, 0), (C12, C11, C13, 0), (C13, C13, C33, 0), (0, 0, 0, C55) )')
                return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13, 0), (C12, C11, C13, 0), (C13, C13, C33, 0), (0, 0, 0, C55)] (hexagonal phase).'
            if hexagonal_condition_5 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13, 0, 0), (C12, C11, C13, 0, 0), (C13, C13, C33, 0, 0), (0, 0, 0, C55, 0), (0, 0, 0, 0, C55) )')
                return False, 'Please enter such a value that det (matrix A) > 0! A = [(C11, C12, C13, 0, 0), (C12, C11, C13, 0, 0), (C13, C13, C33, 0, 0), (0, 0, 0, C55, 0), (0, 0, 0, 0, C55)] (hexagonal phase).'
            if hexagonal_condition_6 > 0:
                pass
            else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13, 0, 0, 0), (C12, C11, C13, 0, 0, 0), (C13, C13, C33, 0, 0, 0), (0, 0, 0, C55, 0, 0), (0, 0, 0, 0, C55, 0), (0, 0, 0, 0, 0, C66) )')
                return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13, 0, 0, 0), (C12, C11, C13, 0, 0, 0), (C13, C13, C33, 0, 0, 0), (0, 0, 0, C55, 0, 0), (0, 0, 0, 0, C55, 0), (0, 0, 0, 0, 0, C66)] (hexagonal phase).'

            return True, ''


class orthorombic(polycrystal):
    def __init__(self, C11 = None, C12 = None, C13 = None, C22 = None,
                 C23 = None, C33 = None, C44 = None, C55 = None, C66 = None):

        if (C11 is None or C12 is None or C13 is None or C22 is None or C23 is None or C33 is None or C44 is None or C55 is None or C66 is None):
            u.inputError("C parameters not set!")
            exit()

        self.C11 = C11/100    # Thanks to Python 3, we do not need to explicitly state "float" in these 'C's
        self.C12 = C12/100
        self.C13 = C13/100
        self.C22 = C22/100
        self.C23 = C23/100
        self.C33 = C33/100
        self.C44 = C44/100
        self.C55 = C55/100
        self.C66 = C66/100

        correctInput, msg = self.checkCond()

        if correctInput == True:

            polycrystal.__init__(self)
            self.crystalname = "orthorombic"
            self.Cparamlist = (C11, C12, C13, C22, C23, C33, C44, C55, C66)
            #            self.addCrystal(self.crystalname, self.Cparamlist)

            K0   = Symbol('K0')
            mue0 = Symbol('mue0')

            m_orthorombic = Matrix (([self.C11,self.C12,self.C13,       0,       0,       0],
                                     [self.C12,self.C22,self.C23,       0,       0,       0],
                                     [self.C13,self.C23,self.C33,       0,       0,       0],
                                     [       0,       0,       0,self.C44,       0,       0],
                                     [       0,       0,       0,       0,self.C55,       0],
                                     [       0,       0,       0,       0,       0,self.C66]))

            S_orthorombic = linalg.inv(m_orthorombic)
            #        print m_orthorombic
            #        print S_orthorombic

            L11 = 0
            L12 = 0
            L13 = 1    # Thanks to Python 3, we do not need to explicitly state "float"

            L21 = 1 / sqrt(2)    # Thanks to Python 3, we do not need to explicitly state "float"
            L22 = L21
            L23 = 0

            L31 = 1 / sqrt(3)    # Thanks to Python 3, we do not need to explicitly state "float"
            L32 = L31
            L33 = L31

            self.EE1 = ( L11**4 + 2 * (L11**2) * (L12**2) * S_orthorombic [0,1] +  2 * (L11**2) * (L13**2) * S_orthorombic [0,2] + L12**4 * S_orthorombic [1,1]\
                         + 2 * (L12**2) * (L13**2) * S_orthorombic [1,2] +  L13**4 * S_orthorombic [2,2] + (L12**2) * (L13**2) * S_orthorombic [3,3]\
                         + (L11**2) * (L13**2) * S_orthorombic [4,4] + (L11**2) * (L12**2) * S_orthorombic [5,5] )

            self.EE2 = ( L21**4 + 2 * (L21**2) * (L22**2) * S_orthorombic [0,1] +  2 * (L21**2) * (L23**2) * S_orthorombic [0,2] + L22**4 * S_orthorombic [1,1]\
                         + 2 * (L22**2) * (L23**2) * S_orthorombic [1,2] +  L23**4 * S_orthorombic [2,2] + (L22**2) * (L23**2) * S_orthorombic [3,3]\
                         + (L21**2) * (L23**2) * S_orthorombic [4,4] + (L21**2) * (L12**2) * S_orthorombic [5,5] )

            self.EE3 = ( L31**4 + 2 * (L31**2) * (L32**2) * S_orthorombic [0,1] +  2 * (L31**2) * (L33**2) * S_orthorombic [0,2] + L32**4 * S_orthorombic [1,1]\
                         + 2 * (L32**2) * (L33**2) * S_orthorombic [1,2] +  L33**4 * S_orthorombic [2,2] + (L32**2) * (L33**2) * S_orthorombic [3,3]\
                         + (L31**2) * (L33**2) * S_orthorombic [4,4] + (L31**2) * (L32**2) * S_orthorombic [5,5] )


            dC11     = self.C11 - K0 - 4/3 * mue0    # Thanks to Python 3, we do not need to explicitly state "float" in these 'dC's
            dC22     = self.C22 - K0 - 4/3 * mue0
            dC33     = self.C33 - K0 - 4/3 * mue0
            dC12     = self.C12 - K0 + 2/3 * mue0
            dC13     = self.C13 - K0 + 2/3 * mue0
            dC23     = self.C23 - K0 + 2/3 * mue0

            a        = dC11 + dC22 + dC33
            b        = dC12 + dC13 + dC23
            c        = dC11 * dC22 + dC11 * dC33 + dC22 * dC33
            d        = dC12**2 + dC13**2 + dC23**2
            e        = dC12 * dC13 + dC12 * dC23 + dC13 * dC23 - dC11 * dC23 - dC22 * dC13 - dC33 * dC12
            deltak   = dC11 * dC22 * dC33 + 2 * dC12 * dC13 * dC23 - dC11 * dC23**2 - dC22 * dC13**2 - dC33 * dC12**2

            Knue     = 1./9 * ( self.C11 + self.C22 + self.C33 + 2 * ( self.C12 + self.C13 + self.C23 ))

            beta     = -3 * (K0 + 2 * mue0) / (5 * mue0 * (3 * K0 + 4 * mue0))
            eta      = -3./(3 * K0 + 4 * mue0)
            gamma    = 1./9 * (eta - 3 * beta)

            denom = (3 * (1 - a * beta - 9 * gamma * (Knue -K0) + beta * (beta + 2 * gamma)\
            * (c - d) - 2 * e * beta * gamma - 1./3 * eta * beta**2 * deltak))
            coeff  = 1 / 15    # Thanks to Python 3, we do not need to explicitly state "float"

            fa = (a - b + beta * (2 * d - 2 * c - e) + 3 * gamma * (d - c + e) + eta * beta * deltak)
            fb = (1 - a * beta - 9 * gamma * (Knue - K0) + beta * (beta + 2 * gamma) \
            * (c - d) - 2 * e * beta * gamma - 1./3 * eta * beta**2 * deltak)
            fc = ((self.C44 - mue0) / (1 - 2 * beta * (self.C44 - mue0)))
            fd = ((self.C55 - mue0) / (1 - 2 * beta * (self.C55 - mue0)))
            fe = ((self.C66 - mue0) / (1 - 2 * beta * (self.C66 - mue0)))

            self.local_1 = coeff *((fa / fb)  + 3 * (fc + fd + fe))
            self.local_2 = (9 * (Knue - K0) + 2 * beta * ( d - c + e) + 3 * beta**2 * deltak)/denom

            self.f1 = self.local_1
            self.f2 = self.local_2

            # -------------------------------  Upper & Lower Bound  ------------------------------- #
    #        G_low  = [C44, C66, C2 / (6* Knue)]
    #        G1 = max(G_low)
    #
    #        G_up = [C44, C66, (1/6. * M)]
    #        G2 = min(G_up)
    #        G3 = 0xBI6B00B5

            # Bulk Modulus


            # 1/ Voigt
            Kvoigt = 1/9. * (self.C11 + self.C22 + self.C33 + 2*self.C12 + 2*self.C13 + 2*self.C23)

            # 2/ Reuss
            kr1 = (self.C11 * self.C22 * self.C33 + 2*(self.C12*self.C13*self.C23) - self.C11*self.C23**2 - self.C22*self.C13**2 - self.C33*self.C12**2)
            kr2 = (self.C11*self.C22 + self.C11*self.C33 + self.C22*self.C33 - self.C12**2 - self.C13**2 - self.C23**2)
            kr3 = 2*(self.C12*self.C13 + self.C12*self.C23 + self.C13*self.C23 - self.C12*self.C33 - self.C13*self.C22 - self.C23*self.C11)
            Kreuss = kr1 * (kr2 + kr3)**-1

            # 3 / HS-
            # 4 / HS+

            # Shear Modulus

            # 1/ Voigt
            Gvoigt = 1/15. * (self.C11 + self.C22 + self.C33 - self.C12 -self.C13 -self.C23 + 3*self.C44 +3*self.C55 + 3*self.C66)

            # 2/ Reuss
            gr1 = self.C11*self.C22 + self.C11*self.C33 + self.C22*self.C33 + self.C11*self.C23 + self.C22*self.C13 +\
                  self.C33*self.C12 - self.C12**2 - self.C13**2 - self.C23**2 - self.C12*self.C13 - self.C12*self.C23 - self.C13*self.C23
            gr2 = (self.C11*self.C22*self.C33 + 2*(self.C12*self.C13*self.C23) - self.C11*self.C23**2 - self.C22*self.C13**2 - self.C33*self.C12**2)**-1
            gr3 = 3/4. *(self.C44**-1 + self.C55**-1 + self.C66**-1)

            Greuss = 15/4. * (gr1 * gr2 + gr3)**-1

            # 3 / HS-
            # 4 / HS+

            self.voigt_bulk   = Kvoigt
            self.reuss_bulk   = Kreuss
        #        self.HS_min_bulk  = KHSn
        #        self.HS_post_bulk = KHSp

            self.voigt_shear   = Gvoigt
            self.reuss_shear   = Greuss
        #        self.HS_min_shear  = KHSn # TODO: wrong equation in literature??
        #        self.HS_post_shear = KHSp
        else:
            print('warning' + '<br>')
            print()
            exit()

def getElasticDict(self):
    mydict = {'C11' : self.C11*100,    # Thanks to Python 3, we do not need to explicitly state "float"
              'C12' : self.C12*100,
              'C13' : self.C13*100,
              'C22' : self.C22*100,
              'C23' : self.C23*100,
              'C33' : self.C33*100,
              'C44' : self.C44*100,
              'C55' : self.C55*100,
              'C66' : self.C66*100}
    return mydict

def checkCond(self):
        # Check condition

        m_orthorombic = Matrix (([self.C11,self.C12,self.C13,       0,       0,       0],
                                 [self.C12,self.C22,self.C23,       0,       0,       0],
                                 [self.C13,self.C23,self.C33,       0,       0,       0],
                                 [       0,       0,       0,self.C44,       0,       0],
                                 [       0,       0,       0,       0,self.C55,       0],
                                 [       0,       0,       0,       0,       0,self.C66]))

        orthorombic_det         = det(m_orthorombic)

        orthorombic_condition_1 = self.C11

        orthorombic_condition_2 = det(Matrix(([self.C11,self.C12],
                                              [self.C12,self.C22])))

        orthorombic_condition_3 = det(Matrix(([self.C11,self.C12,self.C13],
                                              [self.C12,self.C22,self.C23],
                                              [self.C13,self.C23,self.C33])))

        orthorombic_condition_4 = det(Matrix(([self.C11,self.C12,self.C13,       0],
                                              [self.C12,self.C22,self.C23,       0],
                                              [self.C13,self.C23,self.C33,       0],
                                              [       0,       0,       0,self.C44])))

        orthorombic_condition_5= det(Matrix(([self.C11,self.C12,self.C13,       0,       0],
                                             [self.C12,self.C22,self.C23,       0,       0],
                                             [self.C13,self.C23,self.C33,       0,       0],
                                             [       0,       0,       0,self.C44,       0],
                                             [       0,       0,       0,       0,self.C55])))

        orthorombic_condition_6 = det(m_orthorombic)

        if orthorombic_condition_1 > 0:
            pass
        else:
#                u.inputError('Please enter right value: C11 > 0')    # It seems to me that this print statement should not end with a whitespace like it did before!
            return False, 'Please enter such a value that C11 > 0 (orthorombic phase).'
        if orthorombic_condition_2 > 0:
            pass
        else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12), (C12, C22) )')
            return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12), (C12, C22)] (orthorombic phase).'
        if orthorombic_condition_3 > 0:
            pass
        else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13), (C12, C22, C23), (C13, C23, C33) )')
            return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13), (C12, C22, C23), (C13, C23, C33)] (orthorombic phase).'
        if orthorombic_condition_4 > 0:
            pass
        else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13, 0), (C12, C22, C23, 0), (C13, C23, C33, 0), (0, 0, 0, C44) )')
            return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13, 0), (C12, C22, C23, 0), (C13, C23, C33, 0), (0, 0, 0, C44)] (orthorombic phase).'
        if orthorombic_condition_5 > 0:
            pass
        else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13, 0, 0), (C12, C22, C23, 0, 0), (C13, C23, C33, 0, 0), (0, 0, 0, C44, 0), (0, 0, 0, 0, C55) )')
            return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13, 0, 0), (C12, C22, C23, 0, 0), (C13, C23, C33, 0, 0), (0, 0, 0, C44, 0), (0, 0, 0, 0, C55)] (orthorombic phase).'
        if orthorombic_condition_6 > 0:
            pass
        else:
#                u.inputError('Please enter right value: det (matrix A) > 0!\n A = ( (C11, C12, C13, 0, 0, 0), (C12, C22, C23, 0, 0, 0), (C13, C23, C33, 0, 0, 0), (0, 0, 0, C44, 0, 0), (0, 0, 0, 0, C55, 0), (0, 0, 0, 0, 0, C66) )')
            return False, 'Please enter such a value that det (matrix A) > 0!  A = [(C11, C12, C13, 0, 0, 0), (C12, C22, C23, 0, 0, 0), (C13, C23, C33, 0, 0, 0), (0, 0, 0, C44, 0, 0), (0, 0, 0, 0, C55, 0), (0, 0, 0, 0, 0, C66)] (orthorombic phase).'

        return True, ''
