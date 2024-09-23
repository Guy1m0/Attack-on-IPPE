from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.ABEnc import ABEnc
from charm.toolbox.msp import MSP


import random

class KRK_ABE(ABEnc):

    def __init__(self, n, assump_size, group_obj, math):
        ABEnc.__init__(self)
        self.group = group_obj
        self.n = n
        self.assump_size = assump_size
        self.math = math

    def setup(self):
        k = self.assump_size
        group = self.group

        gamma = group.random(ZR)

