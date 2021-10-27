from sympy import Symbol, solve, simplify, nsolve, symbols
from collections import Iterable
from scipy.integrate import odeint 
import numpy as np

class heterodimer():

    def __init__(self):
        #Set up the parameters and variables
        lam1 = Symbol('lam1')
        lam2 = Symbol('lam2')
        lam3 = Symbol('lam3')
        K = Symbol('K')
        c1 = Symbol('c1')
        c2 = Symbol('c2')
        p1 = Symbol('p1')
        p2 = Symbol('p2')
        p3 = Symbol('p3')
        #Assume association/dissociation rate >> degradation and production rate
        #at steady state
        eq1 = c1 - lam1*p1 - lam3*p3
        eq2 = c2 - lam2*p2 - lam3*p3
        eq3 = p3 - p1*p2*K
        #Solving the equations of the heterodimer formation system at steady-state
        ans1 = solve(eq3, p3) #p3
        ans2 = solve([p3-ans1[0], eq1],[p1, p3]) 
        ans3 = solve([p3-ans1[0], eq2],[p2, p3]) 
        ans4 = solve([p2-ans3[p2], p1-ans2[p1]],[p1, p2]) #ans4[1](positive) ans4[1][0]=p1, ans4[1][1]=p2
        ans5 = solve([p3-ans3[p3], p1-ans4[1][0]],[p1, p3]) #ans5[0][1] p3

        P3 = ans5[0][1]
        P3_equal = P3.subs(c1, c2)
        P3_norm = P3/P3_equal
        P1 = ans4[1][0]
        P1_norm = simplify(P1/P3_equal)
        P2 = ans4[1][1]
        P2_norm = simplify(P2/P3_equal)

        self.P3 = P3
        self.P1 = P1
        self.P2 = P2
        self.P3_equal = P3_equal
        self.P3_norm = P3_norm
        self.P1_norm = P1_norm
        self.P2_norm = P2_norm

def numerical_solver_trimer(c1=30, c2=30, c3=30, lam1=0.027, lam2=0.027, lam3=0.027, lam12=0.027, lam123=0.027,
              K1=0.05, K2=0.05, ini_guess=(1, 1, 1)):
    
    p1 = Symbol('p1')
    p2 = Symbol('p2')
    p3 = Symbol('p3')
    
    synthesis_ls = [c1, c2, c3]
    
    a = 0
    for i in range(3):
        if isinstance(synthesis_ls[i], Iterable):
            a = synthesis_ls[i] #save the variable
            c = i+1             #which is the variable
    
    if not isinstance(a, Iterable):  # if none is Iterable
        a = [c1]
        c = 1
    
    x = [] #p1
    y = [] #p2
    z = [] #p3
    p12 = []
    p123 = []
    x_, y_, z_ = ini_guess 
    
    for i in range(len(a)): 
        if c==1: 
            c1 = a[i]
        elif c==2:
            c2 = a[i]
        elif c==3:
            c3 = a[i]
        
        f1 = c1 - lam1*p1 - lam12*p1*p2*K1 - lam123*p1*p2*p3*K1*K2
        f2 = c2 - lam2*p2 - lam12*p1*p2*K1 - lam123*p1*p2*p3*K1*K2
        f3 = c3 - lam3*p3 - lam123*p1*p2*p3*K1*K2
        
        try:
            sol = nsolve((f1, f2, f3), (p1, p2, p3), (x_, y_, z_)) #update the ini_guess to approach the solution
            x_, y_, z_ = sol

        except:
            sol = nsolve((f1, f2, f3), (p1, p2, p3), (1, 1, 1))    #If error, go back to (1, 1, 1)
            x_, y_, z_ = sol

        p12_ = x_*y_*K1
        p123_ = x_*y_*z_*K1*K2

        x.append(x_)
        y.append(y_)
        z.append(z_)
        p12.append(p12_)
        p123.append(p123_)

    return x, y, z, p12, p123

time_scale='s' #second
bacteria=False 

#Common parameter value
class com_v():
    
    if time_scale == 's':
        sc = 1
    else:
        sc = 3600
    
    if bacteria == False:
        #Trancription
        Krt1, Krt2, Kr1, Kr2 = 1e-5*sc, 1e-5*sc, 3e-5*sc, 3e-5*sc
        #Translation 
        Ktt1, Ktt2, Kt1, Kt2 = 1e-2*sc, 1e-2*sc, 1e-2*sc, 1e-2*sc
        #mRNA degradation 
        Lmt1, Lmt2, Lm1, Lm2 = 2e-5*sc, 2e-5*sc, 2e-5*sc, 2e-5*sc
        #Protein degradation 
        Ltf1, Ltf2, L1, L2, L3 = 7e-5*sc, 7e-5*sc, 7e-6*sc, 7e-6*sc, 7e-6*sc
        #binding
        kon, koff = 1e-4, 1e-3
        #Km, Hill coeffi, binding 
        Ktf1, Ktf2, K, h = 70, 70, 50, 2 
    else:
        Krt1, Krt2, Kr1, Kr2 = 1.6e-2*sc, 1.6e-2*sc, 1.6e-2*sc, 1.6e-2*sc
        Ktt1, Ktt2, Kt1, Kt2 = 1.6e-2*sc, 1.6e-2*sc, 1.6e-2*sc, 1.6e-2*sc
        Lmt1, Lmt2, Lm1, Lm2 = 3e-3*sc, 3e-3*sc, 3e-3*sc, 3e-3*sc
        Ltf1, Ltf2, L1, L2, L3 = 5.7e-4*sc, 5.7e-4*sc, 5.7e-4*sc, 5.7e-4*sc, 5.7e-4*sc
        kon, koff = 1e-4, 1e-3
        Ktf1, Ktf2, K, h = 20, 20, 50, 2

class open_circuit():
    
    def steady_state(others=False,coop=1, Kr1=com_v.Kr1, Kr2=com_v.Kr2, Kt1=com_v.Kt1, Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1,c_syn_scale=1,K=com_v.K):
        #Assign values to model parameters
        Krt1, Krt2, Kr1, Kr2 = com_v.Krt1, com_v.Krt2, Kr1*c_syn_scale, Kr2*c_syn_scale
        Ktt1, Ktt2, Kt1, Kt2 = com_v.Ktt1, com_v.Ktt2, syn_scale*Kt1, syn_scale*Kt2
        Lmt1, Lmt2, Lm1, Lm2 = com_v.Lmt1, com_v.Lmt2, com_v.Lm1, com_v.Lm2
        Ltf1, Ltf2, L1, L2, L3 = com_v.Ltf1, com_v.Ltf2, coop*com_v.L1, coop*com_v.L2, com_v.L3
        kon, koff = kon, com_v.koff 
        Ktf1, Ktf2, K, h = com_v.Ktf1, com_v.Ktf2, K, com_v.h
        
        #Chemical species in the system
        mt1, mt2, TF1, TF2, m1, m2, p1, p2, p3 = symbols(('mt1', 'mt2', 'TF1', 'TF2', 'm1', 'm2',
                                                         'p1', 'p2', 'p3'))
    
        #Differential Equations for each chemical species
        eq1 = Krt1 - Lmt1*mt1  
        eq2 = Krt2 - Lmt2*mt2 
        eq3 = Ktt1*mt1 - Ltf1*TF1
        eq4 = Ktt2*mt2 - Ltf2*TF2

        eq5 = Kr1*(TF1**h/(Ktf1**h + TF1**h)) - Lm1*m1
        eq6 = Kr2*(TF2**h/(Ktf2**h + TF2**h)) - Lm2*m2
        eq7 = Kt1*m1 + koff*p3 - kon*p1*p2 - L1*p1
        eq8 = Kt2*m2 + koff*p3 - kon*p1*p2 - L2*p2
        eq9 = kon*p1*p2 - koff*p3 - L3*p3

        #Step1: Solve the TF1 and TF2 from TF circuit
        TF_circuit = solve([eq1, eq2, eq3, eq4], [mt1, mt2, TF1, TF2])
        #Step2: Solve the m1 and m2 by replacing TF1 and TF2 from step1
        mRNA_sol = solve([eq5.subs({TF1:float(TF_circuit[TF1])}), eq6.subs({TF2:float(TF_circuit[TF2])})], [m1, m2])
        #Step3: Solve the p3
        p3_sol = solve(eq9, p3) 
        #Step4: Solve p1 and p2 by replacing p3 from step3
        intermed_1 = solve([p3 - p3_sol[0], eq7], [p1, p3])
        intermed_2 = solve([p3 - p3_sol[0], eq8], [p2, p3])
        #Step5: Solve p1, p2 from step4 by replacing m1, m2 from step2
        intermed_3 = solve([m2 - mRNA_sol[m2], p2 - intermed_2[p2]], [m2, p2]) 
        intermed_4 = solve([m1 - mRNA_sol[m1], p1 - intermed_1[p1]], [m1, p1]) 
        #Step6: Solve p1 from step5 by replacing p2 from step5
        p1p2_sol = solve([p2 - intermed_3[p2], p1 - intermed_4[p1]], [p1, p2]) #p1p2_sol[1][1] p1正號解
        p3_ss = p3_sol[0].subs({p1:p1p2_sol[1][0], p2:p1p2_sol[1][1]})

        mt1 = TF_circuit[mt1] 
        mt2 = TF_circuit[mt2] 
        TF1 = TF_circuit[TF1]
        TF2 = TF_circuit[TF2]
        m1 = intermed_4[m1]
        m2 = intermed_3[m2]

        p1 = p1p2_sol[1][0]
        p2 = p1p2_sol[1][1]
        p3 = p3_ss

        if others == False:
            return p1, p2, p3
        else:
            return p1, p2, p3, m1, m2, TF1, TF2, mt1, mt2
        
    def dynamics(time=1200000, init=[0,0,0,0,0,0,0,0,0], others=False, coop=1, Kr1=com_v.Kr1, Kr2= com_v.Kr2,Kt1=com_v.Kt1,Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1, c_syn_scale=1):
        
        def open_(z, t, Krt1, Krt2, Ktt1, Ktt2, Kr1, 
                      Kr2, Kt2, Ktf1, Ktf2, kon, koff, 
                      Lmt1, Lmt2, Ltf1, Ltf2,
                      Lm1, Lm2, L1, L2, L3, h):

            #Molecules in transcription factor circuit
            mt1 = z[0]
            mt2 = z[1]
            TF1 = z[2]
            TF2 = z[3]

            #Molecules in Protein Complex formation circuit 
            m1 = z[4]
            m2 = z[5]
            p1 = z[6]
            p2 = z[7]
            p3 = z[8]

            #Differential equations 
            dmt1dt = Krt1 - Lmt1*mt1 
            dmt2dt = Krt2 - Lmt2*mt2 
            dTF1dt = Ktt1*mt1 - Ltf1*TF1
            dTF2dt = Ktt2*mt2 - Ltf2*TF2
            dm1dt = Kr1*(TF1**h/(Ktf1**h + TF1**h)) - Lm1*m1
            dm2dt = Kr2*(TF2**h/(Ktf2**h + TF2**h)) - Lm2*m2
            dp1dt = Kt1*m1 + koff*p3 - kon*p1*p2 - L1*p1
            dp2dt = Kt2*m2 + koff*p3 - kon*p1*p2 - L2*p2
            dp3dt = kon*p1*p2 - koff*p3 - L3*p3

            return [dmt1dt, dmt2dt, dTF1dt, dTF2dt, dm1dt, dm2dt, dp1dt, dp2dt, dp3dt]
        #Assign values to model parameters
        Krt1, Krt2, Kr1, Kr2 = com_v.Krt1, com_v.Krt2, Kr1*c_syn_scale, Kr2*c_syn_scale
        Ktt1, Ktt2, Kt1, Kt2 = com_v.Ktt1, com_v.Ktt2, syn_scale*Kt1, syn_scale*Kt2
        Lmt1, Lmt2, Lm1, Lm2 = com_v.Lmt1, com_v.Lmt2, com_v.Lm1, com_v.Lm2
        Ltf1, Ltf2, L1, L2, L3 = com_v.Ltf1, com_v.Ltf2, coop*com_v.L1, coop*com_v.L2, com_v.L3
        kon, koff = kon, com_v.koff 
        Ktf1, Ktf2, h = com_v.Ktf1, com_v.Ktf2, com_v.h
        
        z0 = init
        
        time_period = time
        time_points = time_period*2
        time_plot = 1
        end = int(time_points*time_plot)

        #time period for simulation
        t = np.linspace(0, time_period, time_points)
        
        #numerically solve the system
        z = odeint(open_, z0, t, args=(Krt1, Krt2, Ktt1, Ktt2, Kr1,
                      Kr2, Kt2, Ktf1, Ktf2, kon, koff, 
                      Lmt1, Lmt2, Ltf1, Ltf2,
                      Lm1, Lm2, L1, L2, L3, h))

        #obtain the concentration of each molecule all the time points
        mt1 = z[:end, 0]
        mt2 = z[:end, 1]
        TF1 = z[:end, 2]
        TF2 = z[:end, 3]
        m1 = z[:end, 4]
        m2 = z[:end, 5]
        p1 = z[:end, 6]
        p2 = z[:end, 7]
        p3 = z[:end, 8]
        
        if others==False:
            return p1, p2, p3, m1, m2
        else:
            return p1, p2, p3, m1, m2, TF1, TF2, mt1, mt2

class closed_circuit():
    
    def steady_state(return_others=False, coop=1, Kr1=com_v.Kr1, Kr2= com_v.Kr2,Kt1=com_v.Kt1, Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1,c_syn_scale=1, K=com_v.K):
        #Assign values to model parameters
        Krt1, Krt2, Kr1, Kr2 = com_v.Krt1, com_v.Krt2, Kr1*c_syn_scale, Kr2*c_syn_scale
        Ktt1, Ktt2, Kt1, Kt2 = com_v.Ktt1, com_v.Ktt2, syn_scale*Kt1, syn_scale*Kt2
        Lmt1, Lmt2, Lm1, Lm2 = com_v.Lmt1, com_v.Lmt2, com_v.Lm1, com_v.Lm2
        Ltf1, Ltf2, L1, L2, L3 = com_v.Ltf1, com_v.Ltf2, coop*com_v.L1, coop*com_v.L2, com_v.L3
        kon, koff = kon, com_v.koff 
        Ktf1, Ktf2, K, h = com_v.Ktf1, com_v.Ktf2, K, com_v.h
        
        '''Major players in the system'''
        mt1, mt2, TF1, TF2, m1, m2, p1, p2, p3 = symbols(('mt1', 'mt2', 'TF1', 'TF2', 'm1', 'm2',
                                                         'p1', 'p2', 'p3'))
        
        '''Differential Equations for each major player'''
        eq1 = Krt1 - Lmt1*mt1  
        eq2 = Krt2 - Lmt2*mt2 
        eq3 = Ktt1*mt1 - Ltf1*TF1
        eq4 = Ktt2*mt2 - Ltf2*TF2

        eq5 = Kr1*(TF1**h/(Ktf1**h + TF1**h))*(K/(K + p1)) - Lm1*m1
        eq6 = Kr2*(TF2**h/(Ktf2**h + TF2**h))*(K/(K + p2)) - Lm2*m2
        eq7 = Kt1*m1 + koff*p3 - kon*p1*p2 - L1*p1
        eq8 = Kt2*m2 + koff*p3 - kon*p1*p2 - L2*p2
        eq9 = kon*p1*p2 - koff*p3 - L3*p3

        '''Step1: Solve the TF1 and TF2 from TF circuit'''
        TF_circuit = solve([eq1, eq2, eq3, eq4], [mt1, mt2, TF1, TF2])
        '''Step2: Solve the m1 and m2 by replacing TF1 and TF2 from step1'''
        mRNA_sol = solve([eq5.subs({TF1:float(TF_circuit[TF1])}), eq6.subs({TF2:float(TF_circuit[TF2])})], [m1, m2])
        '''Step3: Solve the p3'''
        p3_sol = solve(eq9, p3) 
        '''Step4: Solve p1 and p2 by replacing p3 from step3'''
        intermed_1 = solve([p3 - p3_sol[0], eq7], [p1, p3])
        intermed_2 = solve([p3 - p3_sol[0], eq8], [p2, p3])
        '''Step5: Solve p1, p2 from step4 by replacing m1, m2 from step2'''
        intermed_3 = solve([m2 - mRNA_sol[m2], p2 - intermed_2[p2]], [m2, p2]) #intermed_3[1][1] p2正號解
        intermed_4 = solve([m1 - mRNA_sol[m1], p1 - intermed_1[p1]], [m1, p1]) #intermed_3[1][1] p1正號解
        '''Step6: Solve p1 from step5 by replacing p2 from step5'''
        p1p2_sol = solve([p2 - intermed_3[1][1], p1 - intermed_4[1][1]], [p1, p2])
        '''Step7: Judge if the solution is negative'''
        try:
            if p1p2_sol[0][0] > 0:
                p1_s = p1p2_sol[0][0]
                p2_s = p1p2_sol[0][1]
            else:
                p1_s = p1p2_sol[1][0]
                p2_s = p1p2_sol[1][1]
        except:
            p1_s = p1p2_sol[1][0]
            p2_s = p1p2_sol[1][1]
        p3_ss = p3_sol[0].subs({p1:p1_s, p2:p2_s})
        '''Step8: Solve m1 and m2 by replacing p1 and p2 from step6'''
        m1_sol = intermed_4[1][0].subs({p2:p2_s})
        m2_sol = intermed_3[1][0].subs({p1:p1_s})

        mt1 = TF_circuit[mt1]
        mt2 = TF_circuit[mt2]
        TF1 = TF_circuit[TF1]
        TF2 = TF_circuit[TF2]
        p3_s = p3_ss


        if return_others == False:
            return p1_s, p2_s, p3_s
        else:
            return p1_s, p2_s, p3_s, m1_sol, m2_sol, TF1, TF2, mt1, mt2

    def dynamics(time=1200000,init=[0,0,0,0,0,0,0,0,0], others=False, coop=1, Kr1=com_v.Kr1,Kr2= com_v.Kr2,Kt1=com_v.Kt1, Kt2=com_v.Kt2 ,kon=com_v.kon, syn_scale=1,c_syn_scale=1, K=com_v.K):
        '''Return p1, p2, p3, m1, m2, TF1, TF2, mt1, mt2'''
        
        def close_(z, t, Krt1, Krt2, Ktt1, Ktt2,
                      Kr1, Kr2, Kt1, Kt2, Ktf1, Ktf2, kon, koff, 
                      Lmt1, Lmt2, Ltf1, Ltf2,
                      Lm1, Lm2, L1, L2, L3, h, K):
            #Molecules in transcription factor circuit
            mt1 = z[0]
            mt2 = z[1]
            TF1 = z[2]
            TF2 = z[3]

            #Molecules in Protein Complex formation circuit 
            m1 = z[4]
            m2 = z[5]
            p1 = z[6]
            p2 = z[7]
            p3 = z[8]

            #Differential equations 
            dmt1dt = Krt1 - Lmt1*mt1 
            dmt2dt = Krt2 - Lmt2*mt1 
            dTF1dt = Ktt1*mt1 - Ltf1*TF1
            dTF2dt = Ktt2*mt2 - Ltf2*TF2
            dm1dt = Kr1*(TF1**h/(Ktf1**h + TF1**h))*(K/(K + p1)) - Lm1*m1
            dm2dt = Kr2*(TF2**h/(Ktf2**h + TF2**h))*(K/(K + p2)) - Lm2*m2
            dp1dt = Kt1*m1 + koff*p3 - kon*p1*p2 - L1*p1
            dp2dt = Kt2*m2 + koff*p3 - kon*p1*p2 - L2*p2
            dp3dt = kon*p1*p2 - koff*p3 - L3*p3

            return [dmt1dt, dmt2dt, dTF1dt, dTF2dt, dm1dt, dm2dt, dp1dt, dp2dt, dp3dt]

        #Assign values to model parameters    
        Krt1, Krt2, Kr1, Kr2 = com_v.Krt1, com_v.Krt2, Kr1*c_syn_scale, Kr2*c_syn_scale
        Ktt1, Ktt2, Kt1, Kt2 = com_v.Ktt1, com_v.Ktt2, syn_scale*Kt1, syn_scale*Kt2
        Lmt1, Lmt2, Lm1, Lm2 = com_v.Lmt1, com_v.Lmt2, com_v.Lm1, com_v.Lm2
        Ltf1, Ltf2, L1, L2, L3 = com_v.Ltf1, com_v.Ltf2, coop*com_v.L1, coop*com_v.L2, com_v.L3
        kon, koff = kon, com_v.koff 
        Ktf1, Ktf2, K, h = com_v.Ktf1, com_v.Ktf2, K, com_v.h
        
        z0 = init
        
        time_period = time
        time_points = time_period*2
        time_plot = 1
        end = int(time_points*time_plot)

        #time period for simulation
        t = np.linspace(0, time_period, time_points)

        #numerically solve the system
        z = odeint(close_, z0, t, args=(Krt1, Krt2, Ktt1, Ktt2,
                          Kr1, Kr2, Kt1, Kt2, Ktf1, Ktf2, kon, koff, 
                          Lmt1, Lmt2, Ltf1, Ltf2,
                          Lm1, Lm2, L1, L2, L3, h, K))

        #obtain the concentration of each molecule all the time points
        mt1 = z[:end, 0]
        mt2 = z[:end, 1]
        TF1 = z[:end, 2]
        TF2 = z[:end, 3]
        m1 = z[:end, 4]
        m2 = z[:end, 5]
        p1 = z[:end, 6]
        p2 = z[:end, 7]
        p3 = z[:end, 8]
        
        if others==False:
            return p1, p2, p3, m1, m2
        else:
            return p1, p2, p3, m1, m2, TF1, TF2, mt1, mt2

