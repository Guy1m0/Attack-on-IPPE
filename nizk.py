from charm.toolbox.pairinggroup import ZR

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

        hash_input = str(R[-1]) + h
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
       
    
    def verify(self, rps, R, u, h):
        hash_input = str(R) + h
        c = self.group.hash(hash_input, ZR)

        return rps.A ** u == R * (rps.B ** c)
    
    def verify_gA(self, pp, B, R, u, h):
        g_A = pp['g_1^A']
        g = pp['g_1']
        hash_input = str(R[-1]) + h
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
