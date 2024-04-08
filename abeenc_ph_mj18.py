from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.ABEnc import ABEnc
from charm.toolbox.msp import MSP

import random

class PH_ABE(ABEnc):

    def __init__(self, n, assump_size, group_obj):
        ABEnc.__init__(self)
        self._group = group_obj
        self._n = n
        self._assump_size = assump_size

    
    def setup(self):
        k = self._assump_size
        group = self._group

        # Generate Pub Param
        A = []
        for i in range(k):
            A.append(group.random(ZR))

        '''
        a_1  0   0   0   0
        0   a_2  0   0   0
        0    0  a_3  0   0
        0    0   0  a_4  0
        0    0   0   0  a_5
        1    1   1   1   1
        '''
        # a special vector to make A^T.dot(a) = 0
        # a = [1/a_1 1/a_2 1/a_3 1/a_4 1/a_5 -1]
        a = []
        for i in range(k):
            a.append(1/A[i])
        a.append(-1)

        U = []
        for j1 in range(k + 1): #col
            y = []
            for j2 in range(k + 1): #row
                y.append(group.random(ZR))
            U.append(y)

        g = group.random(G1)
        h = group.random(G2)
        e_gh = pair(g,h)

        g_A = []
        for j in range(k):
            g_A.append(g**A[j])

        g_UTA = []
        for j1 in range(k+1): #col
            y = []
            for j2 in range(k): #row
                prod = (A[j2] * U[j2][j1] + U[k][j1])
                y.append(g ** prod)
            g_UTA.append(y)    

        pp = {'g_1': g, 'g_2': h, 'e_g1_g2':e_gh, 'g_1^A': g_A, 'g_1^{U^T A}': g_UTA}
        msk = {'A': A, 'U': U, 'a': a}
        return pp, msk
    
    def auth_setup(self, pp):
        pks = {}
        sks = {}
        group = self._group
        k = self._assump_size

        g = pp['g_1']
        h = pp['g_2']
        e_gh = pp['e_g1_g2']
        g_A = pp['g_1^A']

        for i in range(self._n):
            X = [] # random samples <- F_p^{(k+1) * (k+1)}
            for j1 in range(k + 1): #col
                y = []
                for j2 in range(k + 1): #row
                    y.append(group.random(ZR))
                X.append(y)

            tau = []
            for i in range(k+1):
                tau.append(group.random(ZR))

            sigma = group.random(ZR)

            g_XTA = []
            for j1 in range(k+1):
                y = []
                for j2 in range(k):
                    y.append(g_A[j2]**X[j2][j1] * g ** X[k][j1])
                g_XTA.append(y)

            e_gAh = []
            for i in range(k):
                e_gAh.append(pair(g_A[i],h))
            e_gAh.append(e_gh)

            e_gh_tauTA = []
            for i in range(k):
                e_gh_tauTA.append(e_gAh[i]**tau[i] * e_gh**tau[k])

            y = h**sigma

            pk = {'g_1^{X^T A}': g_XTA, 'e(g_1,g_2)^{tau^T A}': e_gh_tauTA, 'y': y}
            sk = {'X': X, 'tau': tau, 'sigma': sigma}

            pks[str(i+1)] = pk
            sks[str(i+1)] = sk
        return pks, sks

    
    def _gen_mu_i(self, i, n, k, H):
        group = self._group
        mu = [0] * (k+1)
        for j in range(1, i):
            mu = [x + y for x,y in zip(mu, H[str(j)+str(i)])]

        if i == n:
            return mu

        for j in range(i+1,n+1):
            mu = [x - y for x,y in zip(mu, H[str(i)+str(j)])]

        return mu

    def _gen_mus(self, n,assump_size):
        group = self._group
        mus = {}
        k = assump_size

        H = self._rand_oracle(n,k)
        h = []

        for i in range(n):
            mus[str(i+1)] = self._gen_mu_i(i+1,n,k,H)
        
        for i in range(k+1):
            h.append(group.random(ZR))
        
        return mus, h, H
    
    # n refers to number of AA here
    # Generate a random oracl that only return constant value for each AA_i
    # @update: Use group.hash(), maps to Z_p^{k+1}
    def _rand_oracle(self, n, k):
        group = self._group

        H = {}
        for i in range(n):
            for j in range(i, n):
                H[str(i+1)+str(j+1)] = self._gen_Z_p(k+1,group)
        
        return H
    

class Inner_Product:
    def __init__(self, group):
        self._group = group

    def gen_x_v(self, n, assump_size, authorized = []):
        k = assump_size
        group = self._group

        if not authorized:
            size = random.randint(1,n-2)
            authorized = random.sample(range(n-1), size)

        #print (authorized)
        vec_x = [0] * n
        vec_v = [0] * n

        for index in authorized:
            if index == authorized:
                print ('incorrect index to be ', n)
                continue

            vec_x[index] = group.random(ZR)
            vec_v[index] = 1

        vec_x[-1] = -sum(vec_x[:-1])
        vec_v[-1] = 1

        return vec_x, vec_v

def main():
    ### Sum up
    n = 10
    assump_size = 3
    group = PairingGroup('MNT224')

    ph_abe = PH_ABE(group)

    pp, msk = ph_abe.sys_setup(assump_size)
    vec_x, vec_v = ph_abe.gen_x_v(n, assump_size)

    pks, sks = ph_abe.auth_setup(n, pp, assump_size)

    print ('Authorized list: ', vec_v)
    K = ph_abe.gen_keys(pp, n, sks, vec_v, assump_size)

    M = group.random(GT)
    print ('M:', M)
    C,vec_s = ph_abe.encrypt(pp, pks, vec_x, M, assump_size)

    M_ = ph_abe.decrypt(K, C, vec_v, n, pp)

    print ('M_:', M_)
