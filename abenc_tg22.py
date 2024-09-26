from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.ABEnc import ABEnc
from charm.toolbox.msp import MSP


import random,time

class DIPE_ABE(ABEnc):

    def __init__(self, n, group_obj):
        ABEnc.__init__(self)
        self.group = group_obj
        self.n = n

    def setup(self):
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


    def _gen_pk_sk(self, pp):
        group = self.group
        n = self.n
        l = n-1

        g = pp['g']
        ala_hat = group.random(ZR)
        ala_0 = group.random(ZR)

        alaphas = []
        g_alaphas = []

        for i in range(l):
            tmp = group.random(ZR)
            alaphas.append(tmp)
            g_alaphas.append(g ** tmp)

        g_ala_hat = g**ala_hat

        pk = {'g^ala_0':g ** ala_0,'g^alaphas': g_alaphas, 'Z':pair(g,g_ala_hat)}
        sk = {'g^ala_hat': g_ala_hat, 'ala_0':ala_0,'alaphas': alaphas}

        return pk, sk

    def encrypt(self, pp, pks, vec_y, M):
        group = self.group
        n = self.n
        l = n-1

        s = group.random(ZR)

        tmp_E0 = 1
        E1 = 1
        for i in range(n):
            pk = pks[str(i+1)]
            #g_alaphas = pk['g^alaphas']
            g_ala_0 = pk['g^ala_0']

            tmp_E0 = tmp_E0 * pk['Z']
            E1 = E1 * g_ala_0

        E0 = M * tmp_E0**s

        for j in range(l):
            tmp_E1 = 1
            for i in range(n):
                pk = pks[str(i+1)]
                g_alaphas = pk['g^alaphas']

                tmp_E1 = tmp_E1 * g_alaphas[j]

            E1 = E1 * tmp_E1 ** vec_y[j]

        E1 = E1 ** s

        E2 = pp['g']**s
        C = {'E0':E0, 'E1':E1, 'E2':E2}
        
        return C



    def keygen(self, pp, sks, GID, vec_x):
        n = self.n
        l = n - 1

        vec_v_str = self._vec2str(vec_x)
        D0 = self.group.hash(str(GID)+vec_v_str, G1)

        K = {}
        D1 = {}

        for i in range(n):
            sk = sks[str(i+1)]

            g_ala_hat = sk['g^ala_hat']
            ala_0 = sk['ala_0']
            alaphas = sk['alaphas']

            D1i = g_ala_hat * D0 ** ala_0

            Ki = []
            for j in range(1, l): #2 to l
                tmp = -vec_x[j]/vec_x[0] * alaphas[0]
                k_tmp = D0 ** tmp * D0 ** alaphas[j]
                Ki.append(k_tmp)

            K[str(i+1)] = Ki
            D1[str(i+1)] = D1i

        return D0, D1, K


    def _vec2str(self, vec):
        vec_str = ""
        for v in vec:
            vec_str += str(v)

        return vec_str
    

    def decrypt(self, D0, D1, K, C, vec_y):
        group = self.group
        n = self.n
        l = n-1

        pi_D = 1
        for i in range(n):
            
            pi_D = pi_D * D1[str(i+1)]

        tmp = pi_D
        for j in range(1,l):
            pi_K = 1
            for i in range(n):
                pi_K = pi_K * K[str(i+1)][j-1]

            pi_K = pi_K ** vec_y[j]

            tmp = tmp*pi_K

        upper = pair(tmp, C['E2'])
        lower = pair(C['E1'], D0)

        return C['E0'] * lower / upper

class Inner_Product_TG22:
    def __init__(self, group):
        self._group = group

    #x1 not equals 0
    def gen_x_y(self, n, authorized = []):
        group = self._group
        l = n - 1

        if not authorized:
            size = random.randint(1,n-2)
            authorized = random.sample(range(n-2), size)

        #print (authorized)
        vec_x = []
        for i in range(l):
            vec_x.append(group.random(ZR))

        vec_v = [0] * l 
        tmp = 0
        for index in authorized:
            if index == authorized:
                print ('incorrect index to be ', n)
                continue

            #vec_x[index] = group.random(ZR)
            vec_v[index] = 1
            tmp += vec_x[index]

        vec_x[-1] = -tmp
        vec_v[-1] = 1

        return vec_x, vec_v