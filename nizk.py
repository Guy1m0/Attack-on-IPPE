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
      
    def commit_pk(self, pp, pk, sk):
        g_A = pp['g_1^A']
        e_gh = pp['e_g1_g2']
        k = len(g_A)

        a1 = self.group.random(ZR)
        a2 = self.group.random(ZR)
        a3 = self.group.random(ZR)
        R1, R2, u1,u2, A2, hash_inputs = ([] for i in range(6))
        
        for i in range(k):
            R1.append(g_A[i] ** a1)
            tmp = pair(g_A[i], pp['g_2'])
            A2.append(tmp)
            R2.append(tmp**a2)
        R1.append(pp['g_1']**a1)
        A2.append(e_gh)
        R2.append(e_gh ** a2)

        A1 = g_A
        B1 = pk['g_1^{X^T A}']
        s1 = sk['X']
        B2 = pk['e(g_1,g_2)^{tau^T A}']
        s2 = sk['tau']
        A3 = pp['g_2']
        B3 = pk['y']
        R3 = A3 ** a3
        s3 = sk['sigma']

        h_x = self.group.hash([A1,B1], ZR)
        h_tau = self.group.hash([A2,B2], ZR)
        h_sig = self.group.hash([A3,B3], ZR)
        h = h_x * h_tau * h_sig

        hash_inputs.append(str(h) + str(h_x))
        hash_inputs.append(str(h) + str(h_tau))
        hash_inputs.append(str(h) + str(h_sig))

        c1 = self.group.hash(hash_inputs[0], ZR)
        c2 = self.group.hash(hash_inputs[1], ZR)
        c3 = self.group.hash(hash_inputs[2], ZR)

        for i in range(k + 1):
            y1 = []
            for j in range(k + 1):
                y1.append(a1 + c1 * s1[i][j])
            u1.append(y1)
            u2.append(a2+c2 * s2[i])
        u3 = a3 + c3 * s3

        pi1 = {'R': R1, 'u': u1}
        pi2 = {'R': R2, 'u': u2}
        pi3 = {'R': R3, 'u': u3}
        pis = {'X': pi1, 'tau': pi2, 'sig': pi3}

        rps1 = {'A': A1, 'B': B1}
        rps2 = {'A': A2, 'B': B2}
        rps3 = {'A': A3, 'B': B3}
        s_pairs = {'X': rps1, 'tau': rps2, 'sig': rps3}

        return s_pairs, pis

    def verify_pk(self, pp, s_pairs, pis):
        g_A = pp['g_1^A']
        e_gh = pp['e_g1_g2']
        k = len(g_A)
        hash_inputs = []

        h_x = self.group.hash([s_pairs['X']['A'], s_pairs['X']['B']], ZR)
        h_tau = self.group.hash([s_pairs['tau']['A'], s_pairs['tau']['B']], ZR)
        h_sig = self.group.hash([s_pairs['sig']['A'], s_pairs['sig']['B']], ZR)
        h = h_x * h_tau * h_sig

        hash_inputs.append(str(h) + str(h_x))
        hash_inputs.append(str(h) + str(h_tau))
        hash_inputs.append(str(h) + str(h_sig))

        c1 = self.group.hash(hash_inputs[0], ZR)
        c2 = self.group.hash(hash_inputs[1], ZR)
        c3 = self.group.hash(hash_inputs[2], ZR)

        for j1 in range(k + 1):
            for j2 in range(k):
                eq1 = g_A[j2]**pis['X']['u'][j2][j1] * pp['g_1'] ** pis['X']['u'][k][j1]
                eq2 = (s_pairs['X']['B'][j1][j2] ** c1) *  pis['X']['R'][j2] * pis['X']['R'][-1]
                if eq1 != eq2:# * R[j1]:
                    print ("False for X in index: ", j1, j2)
                    return False
                
        for i in range(k):
            eq1 = s_pairs['tau']['A'][i] ** pis['tau']['u'][i] * e_gh ** pis['tau']['u'][k]
            eq2 = s_pairs['tau']['B'][i] ** c2 * pis['tau']['R'][i] * pis['tau']['R'][k]
            if eq1 != eq2:
                print ("False for tau in index:", i)
                return False
            
        return s_pairs['sig']['A'] ** pis['sig']['u'] == pis['sig']['R'] * (s_pairs['sig']['B'] ** c3)

    
    def verify(self, rps, R, u, h):
        hash_input = str(R) + h
        c = self.group.hash(hash_input, ZR)

        return rps.A ** u == R * (rps.B ** c)
