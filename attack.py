from abenc_ph_mj18 import PH_ABE, mat_math, Inner_Product
from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.ABEnc import ABEnc
from nizk import NIZK
import matplotlib.pyplot as plt
import numpy as np


class RogueKeyAtt():
    def __init__(self, n, assump_size):
        self.n = n
        self.assump_size = assump_size

    def pks_update(self, ad_index, pks):
        n = self.n
        k = self.assump_size

        pk_ad = [[1] * k] * (k+1)

        for i in range(n-1):
            if i + 1 == ad_index:
                continue
            g_XTA = pks[str(i+1)]['g_1^{X^T A}']

            tmp = []
            for j in range(len(g_XTA)):
                tmp.append([x / y for x,y in zip(pk_ad[j], g_XTA[j])])
            pk_ad = tmp

        pk = pks[str(ad_index)]
        pk['g_1^{X^T A}'] = pk_ad
        pks[str(ad_index)] = pk
        return pks
    
    # def keygen_update(self, pp, sks, vec_v, mus, h, ad_index, K):
    #     k = self.assump_size
    #     g2 = pp['g_2']

    #     sk = sks[str(ad_index+1)]
    #     v = vec_v[ad_index]
    #     mu = mus[str(ad_index+1)]

    #     X = sk['X']
    #     tau = sk['tau']

    #     tmp = self.math.mul_matrices(X, h)
    #     #print (len(mu))
    #     exponent = [x - v * y for x,y in zip(tau, tmp)]
    #     exponent = [x + y for x,y in zip(exponent, mu)]
    
    #     K_i = []
    #     #print (len(exponent))
    #     for j in range(k+1):
    #         K_i.append(g2 ** exponent[j])

    #     K[str(ad_index+1)] = K_i
    #     return K

    
    def gen_omega(self, K, C):
        n = self.n
        k = self.assump_size

        H = K['g_2^h']
        C_i_s = C['C_i']

        C_ = [1] * (k + 1)
        for i in range(n-1):
            C_ = [x / y for x,y in zip(C_, C_i_s[i])]

        return [pair(x,y) for x,y in zip(C_, H)]

        