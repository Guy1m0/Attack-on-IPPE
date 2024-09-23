from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.ABEnc import ABEnc
from charm.toolbox.msp import MSP


import random,time

class DIPE_ABE(ABEnc):

    def __init__(self, n, assump_size, group_obj, math):
        ABEnc.__init__(self)
        self.group = group_obj
        self.n = n
        self.assump_size = assump_size
        self.math = math

    def setup(self):
        k = self.assump_size
        group = self.group

        g = group.random(G1)


        pp = {'g':g}

        return pp
    
    def auth_setup(self, pp):

        pks = {}
        sks = {}
        n = self.n

        for i in range(n):
            pk, sk = self._gen_pk_sk(pp)
            pks[str(i+1)] = pk
            sks[str(i+1)] = sk

        return pks, sks            


    def keygen(self, pp, pks, sks, GID, vec_x):
        n = self.n
        l = n - 1

        vec_x_str = self._vec2str(vec_x)
        D0 = self.group.hash(str(GID)+vec_x_str, G1)
        K = {}

        for i in range(n):
            sk = sks[str(i+1)]
            g_ala = sk['g^ala']
            alaphas = sk['alaphas']

            D1 = g_ala * D0 ** alaphas[0]
            Ki = []

            for j in range(1, l): #2 to l
                tmp = -vec_x[j]/vec_x[0] * alaphas[1]
                k_tmp = D0 ** tmp * D0 ** alaphas[j]
                Ki.append(k_tmp)

            K[str(i+1)] = Ki

        return K


    def _gen_pk_sk(self, pp):
        k = self.assump_size
        group = self.group
        n = self.n
        l = n-1

        g = pp['g']
        ala = group.random(ZR)
        alaphas = []
        g_alaphas = []

        for i in range(l + 1):
            tmp = group.random(ZR)
            alaphas.append(tmp)
            g_alaphas.append(g ** tmp)

        g_ala = g**ala

        pk = {'g^alapha': g_alaphas, 'Z':pair(g,g_ala)}
        sk = {'g^ala': g_ala, 'alaphas': alaphas}

        return pk, sk


    def _vec2str(self, vec):
        vec_str = ""
        for v in vec:
            vec_str += str(v)

        return vec_str