from abenc_ph_mj18 import PH_ABE, mat_math, Inner_Product
from charm.toolbox.pairinggroup import GT
from charm.toolbox.ABEnc import ABEnc
from nizk import NIZK
from attack import RogueKeyAtt
from collections import defaultdict 

import time
import matplotlib.pyplot as plt
import numpy as np
from charm.toolbox.pairinggroup import ZR

class Benchmark():
    def __init__(self, group):
        self.group = group

    def benchmark_scheme(self, assump_size, list_n, datasets):
        k = assump_size
        group = self.group
        # result = []
        math_lib = mat_math()
        if not datasets[str(k)]:
            datasets[str(k)] = defaultdict(str)

        datasets[str(k)]['seq'] = []
        datasets[str(k)]['seq_n'] = list_n

        for n in list_n:
            print ("----------------", n, "-----------------------")
            data = {}
            data['total']  = 0
            data['n'] = n

            ph_abe = PH_ABE(n, assump_size, group, math_lib)
            attributes = Inner_Product(group)   
            nizk = NIZK(group)
            att = RogueKeyAtt(n, assump_size)

            # sys setup
            start_time = time.time()
            pp, msk = ph_abe.setup()
            elapsed_time = time.time() - start_time
            data['sys'] = elapsed_time
            data['total'] += elapsed_time

            # Auth setup
            start_time = time.time()
            pks, sks = ph_abe.auth_setup(pp)
            elapsed_time = time.time() - start_time
            data['auth'] = elapsed_time
            data['total'] += elapsed_time

            # Key Gen
            start_time = time.time()
            GID = group.random(ZR)
            vec_x, vec_v = attributes.gen_x_v(n, assump_size)
            K = ph_abe.keygen(pp, pks, sks, GID, vec_v)
            elapsed_time = time.time() - start_time
            data['keygen'] = elapsed_time
            data['total'] += elapsed_time

            print ("Each AA's cost: ", data['auth']/n + elapsed_time/n)

            # Worst Case
            start_time = time.time()
            for i in range(n):
                s_pairs, pis = nizk.prove_pk(pp, pks[str(i+1)], sks[str(i+1)])
            elapsed_time = time.time() - start_time
            data['prove'] = elapsed_time
            #print ("Each cost of proving: ", elapsed_time/n)

            start_time = time.time()
            for i in range(n):
                check = nizk.verify_pk(pp, s_pairs, pis)
                if not check:
                    print ("Check not pass for AA:", str(i+1))
                    break
            tmp = time.time()
            print ("Each AA's extra cost: ", elapsed_time/n + tmp - start_time)
            elapsed_time = tmp - start_time
            data['verify'] = elapsed_time
            #print ("Verification time: ", elapsed_time)


            #print ("user key_gen: ", elapsed_time)

            # AD's setup
            start_time = time.time()
            ad = vec_v.index(0) + 1
            ad_vec_v = [0] * (n-1) + [1]
            pks = att.pks_update(ad, pks)
            elapsed_time = time.time() - start_time
            data['ad_setup'] = elapsed_time
            #print ("Adv's setup: ", elapsed_time)

            # AD KeyGen
            start_time = time.time()
            GID_ad = group.random(ZR)
            K_ = ph_abe.keygen(pp, pks, sks, GID_ad, ad_vec_v)
            elapsed_time = time.time() - start_time
            data['ad_keygen'] = elapsed_time
            #print ("Adv's keyGen: ", elapsed_time)

            # Encryption
            start_time = time.time()
            M = group.random(GT)
            #print ('M:', M)
            C, vec_s = ph_abe.encrypt(pp, pks, vec_x, M)
            elapsed_time = time.time() - start_time
            #print (elapsed_time)
            data['encrypt'] = elapsed_time
            data['total'] += elapsed_time

            # Decryption
            start_time = time.time()
            M_ = ph_abe.decrypt(K, C, vec_v, pp)
            if M_ != M:
                print ('Error in decrypt M (usr): ', M_)
            elapsed_time = time.time() - start_time
            data['decrypt'] = elapsed_time
            data['total'] += elapsed_time

            # AD decrypt
            tmp = ph_abe.decrypt(K_, C, ad_vec_v, pp)

            start_time = time.time()
            omega = att.gen_omega(K_,C)
            M_ =  tmp * math_lib.prod(omega)
            if M_ != M:
                print ('Error in decrypt M (adv): ', M_)
            elapsed_time = time.time() - start_time
            data['ad_cancel_out'] = elapsed_time
            print ("Adv's extra cost: ", data['ad_setup'] + data['ad_cancel_out'])

            print ("Total cost: ", data['total'])
            print ("Total NIZK cost: ", data['prove'] + data['verify'])
            datasets[str(k)]['n']=data
            datasets[str(k)]['seq'].append(data)

        return datasets
            


class Plot():
    def __init__(self):
        return
    
    def plot_total(xs, ys, k_values, sp_ks, sp_xs, sp_ys):
        # Plotting:
        plt.figure(figsize=(10,6))

        # List of linestyles for variation. You can extend this if needed.
        linestyles = ['-', '-.']
        sp_linestyles = ['--',':']

        # Plotting the regular data points for each k_value
        for k, y, style in zip(k_values, ys, linestyles):
            plt.plot(xs, y, label=f'k={k}', linestyle=style)

        # Plotting the special cases
        for k, y, style in zip(sp_ks, sp_ys, sp_linestyles):
            plt.plot(sp_xs, y, linestyle=style, marker='o', label=f'Special k={k}')

        vlines_x = [45, 85]  # example x-coordinates for the vertical lines

        for vx in vlines_x:
            plt.axvline(x=vx, color='gray', linestyle='--', alpha=0.7)  # draw vertical line
            
            if vx in xs:  # If vx is a value in xs, get its intersection with the regular plots
                idx = xs.index(vx)
                for y in ys:
                    rounded_y = round(y[idx], 2)
                    plt.scatter(vx, rounded_y, color='red', zorder=5)
                    plt.annotate(f'({vx}, {rounded_y})', (vx, rounded_y), textcoords="offset points", xytext=(0,10), ha='center')
            
            if vx in sp_xs:  # If vx is a value in sp_x, get its intersection with the special plots
                idx_sp = sp_xs.index(vx)
                for y in sp_ys:
                    rounded_y_sp = round(y[idx_sp], 2)
                    plt.scatter(vx, rounded_y_sp, color='blue', zorder=5)
                    plt.annotate(f'({vx}, {rounded_y_sp})', (vx, rounded_y_sp), textcoords="offset points", xytext=(0,10), ha='center')

        # Setting labels, title, legend, etc.
        plt.xlabel('Size of attribute authority')
        plt.ylabel('Elapsed Time (seconds)')
        #plt.title('Performance based on k-LIN Assumption Size')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    def plot_extra(xs, ys, extra_costs, k_values):
        plt.figure(figsize=(10,6))
        linestyles = ['-', '-.']
        sp_linestyles = ['--', ':']

        # Ensure extra_costs has the same structure as ys
        assert len(ys) == len(extra_costs), "Each 'ys' must have a corresponding 'extra_costs' list."

        # Plotting the regular data points for each k_value and their adjusted times
        for k, y, extra, style in zip(k_values, ys, extra_costs, linestyles):
            adjusted_y = [original + extra for original, extra in zip(y, extra)]
            plt.plot(xs, y, label=f'Original k={k}', linestyle=style)
            plt.plot(xs, adjusted_y, label=f'Adjusted k={k}', linestyle=style, alpha=0.7)

            # Filling the area between the original and adjusted lines
            plt.fill_between(xs, y, adjusted_y, alpha=0.2, color='gray')

            # Annotating the percentage increase at specific points
            for xi, original, adj in zip(xs, y, adjusted_y):
                percentage_increase = ((adj - original) / original) * 100
                if percentage_increase > 0:
                    plt.annotate(f'{percentage_increase:.2f}%', (xi, original + (adj - original)/2), textcoords="offset points", xytext=(0,5), ha='center')

        plt.xlabel('Size of attribute authority')
        plt.ylabel('Total Time Cost (seconds)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

