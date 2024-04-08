from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.ABEnc import ABEnc
from charm.toolbox.msp import MSP

import random

class PH_ABE(ABEnc):

    def __init__(self, n, assump_size, group_obj, math):
        ABEnc.__init__(self)
        self._group = group_obj
        self._n = n
        self._assump_size = assump_size
        self._math = math

    
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
        n = self._n

        g = pp['g_1']
        h = pp['g_2']
        e_gh = pp['e_g1_g2']
        g_A = pp['g_1^A']

        for i in range(n):
            pk, sk = self._gen_pk_sk(pp)
            pks[str(i+1)] = pk
            sks[str(i+1)] = sk

        return pks, sks
    
    def _gen_pk_sk(self, pp):
        k = self._assump_size
        group = self._group

        g = pp['g_1']
        h = pp['g_2']
        e_gh = pp['e_g1_g2']
        g_A = pp['g_1^A']

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
        return pk, sk
    
    def keygen(self, pp, sks, vec_v):
        k = self._assump_size
        g2 = pp['g_2']
        n = self._n

        mus,h, _ = self._gen_mus()
        # mu_list = []
        # for i in range(n):
        #     mu_list.append(mus[str(i+1)])
        # print (vrf_mus(mu_list, k))

        K = {}

        for i in range(n):
            #pk = pks[str(i+1)]
            sk = sks[str(i+1)]
            v = vec_v[i]
            mu = mus[str(i+1)]


            X = sk['X']
            tau = sk['tau']

            tmp = self._math.mul_matrices(X, h)
            #print (len(mu))
            exponent = [x - v * y for x,y in zip(tau, tmp)]
            exponent = [x + y for x,y in zip(exponent, mu)]
        
            K_i = []
            #print (len(exponent))
            for j in range(k+1):
                K_i.append(g2 ** exponent[j])

            K[str(i+1)] = K_i
            #print (len(h))
            K['g_2^h'] = [g2 ** ele for ele in h]

        return K
    
    def encrypt(self, pp, pks, vec_x, M):
        k = self._assump_size
        group = self._group
        n = len(pks)

        g = pp['g_1']
        h = pp['g_2']
        e_gh = pp['e_g1_g2']
        g_A = pp['g_1^A']
        g_UTA = pp['g_1^{U^T A}']

        vec_s = []
        for i in range(k):
            vec_s.append(group.random(ZR))

        #print (vec_x)

        C_0 = []
        for i in range(k):
            C_0.append(g_A[i] ** vec_s[i])
        C_0.append(g ** sum(vec_s))
        #C_0 = power_mul(g_A, vec_s, type='vector')

        C_i = []
        for i in range(n):
            g_XTA = pks[str(i+1)]['g_1^{X^T A}']
            g_UTAs = self._math.power_mul(g_UTA, vec_s, type='m_vector')
            
            if vec_x[i] != 0:
                tmp = self._math.power_mul(g_UTAs, vec_x[i], type = 'v_scalar')
            else:
                tmp = [1] * (k+1)
                
            tmp2 = self._math.power_mul(g_XTA, vec_s, type = 'm_vector')
            #print (len(tmp), len(tmp2))
            C_i.append([x * y for x,y in zip(tmp, tmp2)])

        C_ = M
        #print (M)
        for i in range(n):
            e_gh_tauTA = pks[str(i+1)]['e(g_1,g_2)^{tau^T A}']
            #print (e_gh_tauTA)
            C_ = C_ * self._math.power_mul(e_gh_tauTA, vec_s, type = 'v_vector')

        C = {'C_0': C_0, 'C_i': C_i, 'C':C_}
        return C, vec_s
    
    def decrypt(self, K, C, vec_v, pp):
        H = K['g_2^h']
        group = self._group

        g1 = pp['g_1']
        g2 = pp['g_2']

        C_i_s = C['C_i']
        C_0 = C['C_0']
        ciphertext = C['C']

        K_ = K['1']
        C_ = C_i_s[-1]
        #print (len(C_))
        for i in range(1, self._n):
            K_ = [x * y for x,y in zip(K_, K[str(i+1)])]

            if vec_v[i-1]:
                #print (len(C_i_s[i-1]))
                C_ = [x * y for x,y in zip(C_, C_i_s[i-1])]

        #print (len(H), len(C_),len(C_0), len(K_))
        e_gh_C_0 = [pair(x, y) for x,y in zip(C_0, K_)]
        e_gh_C_i = [pair(x, y) for x,y in zip(C_, H)]
        #print ('cancel out?')
        #print (e_gh_C_i)

        rst = e_gh_C_0[0] * e_gh_C_i[0]
        #print (H, C_)
        for i in range(1, len(e_gh_C_0)):
            rst = rst * e_gh_C_0[i] * e_gh_C_i[i]
        #print (rst)
        #return rst
        return ciphertext/rst
    
    # n refers to number of AA here
    # Generate a random oracl that only return constant value for each AA_i
    # @update: Use group.hash(), maps to Z_p^{k+1}
    def _rand_oracle(self, n, k):
        group = self._group

        H = {}
        for i in range(n):
            for j in range(i, n):
                H[str(i+1)+str(j+1)] = self._gen_Z_p(k+1)
        
        return H
    
    def _gen_Z_p(self, k):
        group = self._group
        vec = []
        if k == 1:
            return group.random(ZR)
        
        for i in range(k):
            vec.append(group.random(ZR))
        
        return vec
    
    def _gen_mu_i(self, i, n, H):
        group = self._group
        k = self._assump_size

        mu = [0] * (k+1)
        for j in range(1, i):
            mu = [x + y for x,y in zip(mu, H[str(j)+str(i)])]

        if i == n:
            return mu

        for j in range(i+1,n+1):
            mu = [x - y for x,y in zip(mu, H[str(i)+str(j)])]

        return mu

    def _gen_mus(self):
        group = self._group
        n = self._n
        k = self._assump_size
        mus = {}

        H = self._rand_oracle(n,k)
        h = []

        for i in range(n):
            mus[str(i+1)] = self._gen_mu_i(i+1,n,H)
        
        for i in range(k+1):
            h.append(group.random(ZR))
        
        return mus, h, H
    

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
    

class mat_math():
    def prod(self, iterable, type = 'scalar'):
        if type == 'scalar':
            result = 1
            for item in iterable:
                result *= item

            return result
        
        if type == 'vector':
            result = [self.prod(values) for values in zip(*iterable)]
            return result

    def transpose(self, matrix):
        # If it's a vector, convert to a 1xN matrix
        if not isinstance(matrix[0], list):
            matrix = [matrix]
        
        # Transpose the matrix
        transposed = [[row[i] for row in matrix] for i in range(len(matrix[0]))]
        
        return transposed

    def mul_matrices(self, A, B):
        # Check if B is a vector
        if not isinstance(B[0], list):
            result = [sum(A[i][j] * B[j] for j in range(len(B))) for i in range(len(A))]
        else:  # B is a matrix
            result = [[sum(A[i][j] * B[j][k] for j in range(len(B))) for k in range(len(B[0]))] for i in range(len(A))]

        return result

    def add_matrices(self, A, B):
        # Check if they are vectors (1D) or matrices (2D)
        if isinstance(A[0], list) and isinstance(B[0], list):  # Both are matrices
            if len(A) != len(B) or len(A[0]) != len(B[0]):
                raise ValueError("Matrices dimensions do not match")
            return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
        elif not isinstance(A[0], list) and not isinstance(B[0], list):  # Both are vectors
            if len(A) != len(B):
                raise ValueError("Vectors lengths do not match")
            return [A[i] + B[i] for i in range(len(A))]
        else:
            raise ValueError("Incompatible dimensions: One is vector, the other is matrix")

    def subtract_matrices(self, A, B):
        # Check if they are vectors (1D) or matrices (2D)
        if isinstance(A[0], list) and isinstance(B[0], list):  # Both are matrices
            if len(A) != len(B) or len(A[0]) != len(B[0]):
                raise ValueError("Matrices dimensions do not match")
            return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
        elif not isinstance(A[0], list) and not isinstance(B[0], list):  # Both are vectors
            if len(A) != len(B):
                raise ValueError("Vectors lengths do not match")
            return [A[i] - B[i] for i in range(len(A))]
        else:
            raise ValueError("Incompatible dimensions: One is vector, the other is matrix")

    def power_mul(self, A, B, type = 'matrix'):
        if type == 'v_scalar':
            result=[]
            for elem in A:
                # print ('new elem:')
                # print (elem)
                result.append(elem ** B)
            return result

        if type == 'm_scalar':
            return [[elem ** B for elem in row] for row in A]
        
        if type == 'v_vector':
            if len(A) != len(B):
                raise ValueError("Number of elements in A must be equal to number of elements in B")

            return self.prod([a ** b for a,b in zip(A,B)])
        
        if type == 'm_vector':
            # print (A[0])
            # print (B)
            if len(A[0]) != len(B):
                raise ValueError("Number of columns in A must be equal to number of elements in B")
            
            result = [1] * len(A)

            for i in range(len(A)):
                for j in range(len(B)):
                    result[i] *= A[i][j] ** B[j]
            return result

        result = [[1 for _ in range(len(B[0]))] for _ in range(len(A))]
        for i in range(len(A)):
            for j in range(len(B[0])):
                for k in range(len(B)):
                    result[i][j] *= A[i][k] ** B[k][j]

        return result

if __name__ == "__main__":
    ### Sum up
    # print ("Start?")
    n = 10
    assump_size = 3
    group = PairingGroup('MNT224')
    math_lib = mat_math()

    ph_abe = PH_ABE(n, assump_size, group, math_lib)

    pp, msk = ph_abe.setup()
    attributes = Inner_Product(group)
    vec_x, vec_v = attributes.gen_x_v(n, assump_size)

    pks, sks = ph_abe.auth_setup(pp)
    #print (sks)

    print ('Authorized list: ', vec_v)
    K = ph_abe.keygen(pp, sks, vec_v)

    M = group.random(GT)
    print ('M:', M)
    C,vec_s = ph_abe.encrypt(pp, pks, vec_x, M)

    M_ = ph_abe.decrypt(K, C, vec_v, pp)

    print ('M_:', M_)
