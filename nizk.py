from charm.toolbox.pairinggroup import ZR, pair

class NIZK():
    def __init__(self, group):
        self.group = group

    def commit(self, rps, h):
        a = self.group.random(ZR)
        R = rps.A ** a
        #print (a, R)
        hash_input = str(R) + h
        c = self.group.hash(hash_input, ZR)
        #print ("add 1:", c * self.rps.s)
        #print ("add 2:", a)
        u = a + c * rps.s
        return R, u
    
    def commit_gA(self, pp, s, h):
        g_A = pp['g_1^A']
        a = self.group.random(ZR)
        #a = 1
        R = []
        k = len(g_A)
        for i in range(k):
            R.append(g_A[i] ** a)
        R.append(pp['g_1']**a)

        hash_input = str(R) + h
        c = self.group.hash(hash_input, ZR)
        u = []
        # test
        # c = 1
        for i in range(k + 1):
            y = []
            for j in range(k + 1):
                y.append(a + c * s[i][j])
            u.append(y)

        return R, u
    
    def commitPK(self, pp, pk, sk):
        e_gh = pp['e_g1_g2']
        g_A = pp['g_1^A']
        k = len(g_A)
        a = self.group.random(ZR)
        #a = 0

        #A = [pair(x, pp['g_2']) for x in g_A].append(e_gh)
        R = [] #A
        A = []
        for i in range(k):
            tmp = pair(g_A[i],pp['g_2'])
            A.append(tmp)
            R.append(tmp**a)
        A.append(e_gh)
        R.append(e_gh ** a)

        B = pk['e(g_1,g_2)^{tau^T A}']
        s = sk['tau']

        h_x = self.group.hash([A,B], ZR)
        #h_tau = self.group.hash(sk['tau'])
        #h_sig = self.group.hash(sk['sigma'])
        #h = h_x * h_tau * h_sig

        hash_input = str(R) + str(h_x)
        c = self.group.hash(hash_input, ZR)
        #c = 1
        #print ('c:', c)
        u = []
        for i in range(k + 1):
            u.append(a + c * s[i])
        #u = [a + c * x for x in s]

        return R, u
    
    def verify(self, rps, R, u, h):
        hash_input = str(R) + h
        c = self.group.hash(hash_input, ZR)

        return rps.A ** u == R * (rps.B ** c)
    
    def verify_gA(self, pp, B, R, u, h):
        g_A = pp['g_1^A'] #A
        g = pp['g_1']
        hash_input = str(R) + h
        c = self.group.hash(hash_input, ZR)

        #g_A to u (scaled XT)
        g_uTA = []
        RBc = []
        k = len(u) - 1
        for j1 in range(k + 1):
            for j2 in range(k):
                eq1 = g_A[j2]**u[j2][j1] * g ** u[k][j1]

                eq2 = (B[j1][j2] ** c) *  R[j2] * R[-1]
                if eq1 != eq2:# * R[j1]:
                    print ("False", j1, j2)
                    return False
                # else:
                #     print ("True ", j1, j2)

        return True


    def verifyPK(self, pp, pk, R, u):
        e_gh = pp['e_g1_g2']
        g_A = pp['g_1^A']
        k = len(g_A)

        e_gAh = [] #A
        for i in range(k):
            e_gAh.append(pair(g_A[i],pp['g_2']))
        e_gAh.append(e_gh)

        A = e_gAh
        B = pk['e(g_1,g_2)^{tau^T A}']

        h_x = self.group.hash([A,B], ZR)

        hash_input = str(R) + str(h_x)
        c = self.group.hash(hash_input, ZR)
        #c = 1
        #print ('c:',c)

        for i in range(k):
            eq1 = e_gAh[i] ** u[i] * e_gh ** u[k]
            eq2 = B[i] ** c * R[i] * R[k]
            if eq1 != eq2:
                print ("False", i)
                return False
            
        return True

class s_pair:
    def __init__(self, s, A, group, B = None):
        self.s = s
        if B :
             self.B = B
        else :
            self.B = A ** s
        self.A = A
        self.group = group

    def getA(self):
        return self.A
