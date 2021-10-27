import numpy as np 
from .Models_and_solvers import com_v
import time


class cal():
    def recovery_time(dynamics ,ss, error_rate=0.0001, time = 1200000):
        p_steady = float(ss)
        p = abs(dynamics - p_steady)
        error = p_steady*error_rate
        t = np.linspace(0, time, time*2)

        s = 0
        n = 0
        m = time*2
        for i in range(m):
            if i> 1/4000*m:
                if p[i] < error and n == 0:
                    s = t[i]
                    n +=1
        return s
    
class performance:
    def __init__(self, func):
        self.func = func
        
    def efficiency(self, fc=1, coop=1, Kr1=com_v.Kr1,Kr2= com_v.Kr2,Kt1=com_v.Kt1, Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1,c_syn_scale=1, K=com_v.K, time=1200000):
        try:
            P = self.func(coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
        except:
            P = self.func(coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale)
        try:
            p3 = P[2][-1]
            m1 = P[3][-1]
        except:
            p3 = P[2]
            m1 = P[3]
            
        return p3/m1
        
    def dimer_ratio(self, fc=1, coop=1, Kr1=com_v.Kr1,Kr2= com_v.Kr2,Kt1=com_v.Kt1, Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1,c_syn_scale=1, K=com_v.K):
        try:
            P = self.func(coop=coop, Kr1=fc*Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
        except:
            P = self.func(coop=coop, Kr1=fc*Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale)
        try:
            p3 = P[2][-1]
            p1 = P[0][-1] 
            p2 = P[1][-1]
        except:
            p3 = P[2]
            p1 = P[0]
            p2 = P[1]

        return p3/(p1+p2+p3)*100
    
    def controllability(self, fc=0.5, coop=1, Kr1=com_v.Kr1,Kr2= com_v.Kr2,Kt1=com_v.Kt1, Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1,c_syn_scale=1, K=com_v.K):
        try:
            P = self.func(coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
            P_d = self.func(coop=coop, Kr1=Kr1*fc, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)      
        except:
            P = self.func(coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale)
            P_d = self.func(coop=coop, Kr1=Kr1*fc, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale)      
        try:
            p3 = P[2][-1]
            p3_d = P_d[2][-1] 
        except:
            p3 = P[2]
            p3_d = P_d[2]

        return p3_d/p3
    
    def recovery_time(self, shot=5, translation=True, coop=1, Kr1=com_v.Kr1,Kr2= com_v.Kr2,Kt1=com_v.Kt1, Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1,c_syn_scale=1, K=com_v.K, time=1200000):
        try:
            p1, p2, p3, m1, m2, TF1, TF2, mt1, mt2 = self.func(others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
            init = [mt1[-1], mt2[-1], TF1[-1], TF2[-1], m1[-1], m2[-1], shot*p1[-1], p2[-1], p3[-1]]
            P = self.func(init=init, others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
        except:
            p1, p2, p3, m1, m2, TF1, TF2, mt1, mt2 = self.func(others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale)
            init = [mt1[-1], mt2[-1], TF1[-1], TF2[-1], m1[-1], m2[-1], shot*p1[-1], p2[-1], p3[-1]]
            P = self.func(init=init, others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale)
        recover_time = cal.recovery_time(P[0], p1[-1]*1.3, error_rate=0.001, time=time)
        
        return recover_time
    
    def turn_on(self, shot=5, translation=True, coop=1, Kr1=com_v.Kr1,Kr2= com_v.Kr2,Kt1=com_v.Kt1, Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1,c_syn_scale=1, K=com_v.K, time=1200000):
        
        try:
            P = self.func(others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
        except:
            P = self.func(others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale)
                                                 #50
        recover_time = cal.recovery_time(P[2], 200, error_rate=0.001, time=time) 
        
        return recover_time
    
    def turn_off(self,shot=5, translation=True,coop=1, Kr1=com_v.Kr1,Kr2= com_v.Kr2,Kt1=com_v.Kt1, Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1,c_syn_scale=1, K=com_v.K, time=1200000):
        try:
            p1, p2, p3, m1, m2, TF1, TF2, mt1, mt2 = self.func(others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
        except:
            p1, p2, p3, m1, m2, TF1, TF2, mt1, mt2 = self.func(others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale)
        
        init = [mt1[-1], mt2[-1], TF1[-1], TF2[-1], m1[-1], m2[-1], p1[-1], p2[-1], p3[-1]]
        
        if translation==True:
            try:
                P = self.func(init=init, others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=0, Kt2=0 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
            except:
                P = self.func(init=init, others=True, time=time, coop=coop, Kr1=Kr1, Kr2=Kr2, Kt1=0, Kt2=0 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale)
        else:
            try:
                P = self.func(init=init, others=True, time=time, coop=coop, Kr1=0, Kr2=0, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
            except:
                P = self.func(init=init, others=True, time=time, coop=coop, Kr1=0, Kr2=0, Kt1=Kt1, Kt2=Kt2 ,kon=kon, syn_scale=syn_scale,c_syn_scale=c_syn_scale, K=K)
        recover_time = cal.recovery_time(P[2], p3[-1]*0.5, error_rate=0.001, time=time)
        return recover_time

def evaluation_ss(open_func, closed_func, n=2, fc=1, l='_label'):
    start = time.time()
    n = n

    o_wo = []
    o_with = []
    c_with = []
    c_wo = []

    x = np.logspace(-6, -3, n) #Disscociation constant 
    y = np.logspace(np.log10(0.036), np.log10(36.3037), n) #Protein synthesis rate
    #z = np.linspace(1.1, 50, n)
    z = np.linspace(0.5, 22, 10) #cooperative stability
    z = z[::-1]
    z = 25/z
    w = np.logspace(0, 3, n)  #TF-DNA Kd
    w = w[::-1]


    X, Y = np.meshgrid(x, y)
    Z, W = np.meshgrid(z, w)

    first_layer = [np.ravel(X), np.ravel(Y)]

    m = 0
    for i in range(len(first_layer[0])):
        coop = 0
        K = 0
        
        kon = first_layer[0][i]
        syn_scale = first_layer[1][i]
        e = open_func(fc=fc, kon=kon, syn_scale=syn_scale)
        config = (kon, syn_scale, coop, K)
        o_wo.append((e, config))
        m += 1

        for j in range(len(z)):
            K = 0
            coop = z[j]
            e = open_func(fc=fc, coop=coop, kon=kon, syn_scale=syn_scale)
            config = (kon, syn_scale, coop, K)
            o_with.append((e, config))
            m += 1

            for s in range(len(w)):
                K = w[s]
                e = closed_func(fc=fc, coop=coop, kon=kon, syn_scale=syn_scale, K=K)
                config = (kon, syn_scale, coop, K)
                c_with.append((e, config))
                m += 1

        coop = 0
        for j in range(len(w)):
            K = w[j]
            e = closed_func(fc=fc, kon=kon, syn_scale=syn_scale, K=K)
            config = (kon, syn_scale, coop, K)
            c_wo.append((e, config))
            m += 1

    end = time.time()
    print(end-start)
    print(m)
    return o_wo, o_with, c_wo, c_with

def evaluation_d(open_func, closed_func, n=2, shot=1, translation=True, l='_label'):
    start = time.time()
    n = n

    o_wo = []
    o_with = []
    c_with = []
    c_wo = []


    x = np.logspace(-6, -3, n) #Disscociation constant 
    y = np.logspace(np.log10(0.036), np.log10(36.3037), n) #Protein synthesis rate
    z = np.linspace(0.5, 22, n) #Cooperative stability
    z = z[::-1]
    z = 25/z
    w = np.logspace(0, 3, n)  #TF-DNA Kd
    w = w[::-1]


    X, Y = np.meshgrid(x, y)
    Z, W = np.meshgrid(z, w)

    first_layer = [np.ravel(X), np.ravel(Y)]
    second_layer = [np.ravel(Z), np.ravel(W)]

    m = 0
    for i in range(len(first_layer[0])):
        coop = 0
        K = 0
        
        kon = first_layer[0][i]
        syn_scale = first_layer[1][i]
        e = open_func(shot=shot, translation=translation, kon=kon, syn_scale=syn_scale)
        config = (kon, syn_scale, coop, K)
        o_wo.append((e, config))
        m += 1

        for j in range(len(z)):
            K = 0
            coop = z[j]
            e = open_func(shot=shot, translation=translation, coop=coop, kon=kon, syn_scale=syn_scale)
            config = (kon, syn_scale, coop, K)
            o_with.append((e, config))
            m += 1

            for s in range(len(w)):
                K = w[s]
                e = closed_func(shot=shot, translation=translation, coop=coop, kon=kon, syn_scale=syn_scale, K=K)
                config = (kon, syn_scale, coop, K)
                c_with.append((e, config))
                m += 1

        coop = 0
        for j in range(len(w)):
            K = w[j]
            e = closed_func(shot=shot, translation=translation, kon=kon, syn_scale=syn_scale, K=K)
            config = (kon, syn_scale, coop, K)
            c_wo.append((e, config))
            m += 1

    end = time.time()
    print(end-start)
    print(m)

    return o_wo, o_with, c_wo, c_with

def save_data(filename, data):
    f = open(filename+'.txt', 'w')
    for i in data:
        f.write(str(i)+'\n')
    f.close()
    
def read_data(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    LS= []
    for i in lines:
        ls = i.strip().split(', ')
        ls = [float(x.strip('()')) for x in ls]
        ls = [ls[0], ls[1:]]
        LS.append(ls)

    f.close()
    return LS