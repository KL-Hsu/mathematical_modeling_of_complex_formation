from sympy import Polygon
import numpy as np 
from itertools import combinations as comb

def owo_1000(ls1, ls2):
    n = 0
    m = 0
    t = 0
    for i in ls1:
        key = i[0]
        ls = ls2[n:n+10]
        for j in ls:
            if key > j[0] and key!=0 and j[0]!=0:
                m+=1
            elif key==0 or j[0]==0:
                t+=1
        n+=10

    return m*10, t*10
    
def owo_10000(ls1, ls2):
    n = 0
    m = 0
    t = 0
    for i in ls1:
        key = i[0]
        ls = ls2[n:n+100]
        for j in ls:
            if key > j[0] and key!=0 and j[0]!=0:
                m+=1
            elif key==0 or j[0]==0:
                t+=1
        n+=100
    return m, t
    
def owith_cwo(ls1, ls2):
    n = 0
    m = 0
    t = 0
    for i in range(100):
        l1 = ls1[n:n+10]
        l2 = ls2[n:n+10]
        for j in l1:
            for k in l2:
                if j[0] == 0 or k[0] == 0:
                    t+=1 
                else:
                    if j[0]>k[0]:
                        m+=1
        n+=10
        
    return m, t

def owith_cwith(ls1, ls2):
    n = 0
    m = 0
    t = 0
    for i in ls1:
        key = i[0]
        ls = ls2[n:n+10]
        for j in ls:
            if key > j[0] and key!=0 and j[0]!=0:
                m+=1
            elif key ==0 or j[0] == 0:
                t+=1
        n +=10
    return m, t
    
def cwo_cwith(ls1, ls2):
    
    def count_100(n):
        if n ==0:
            return 0
        elif n%100 == 0:
            return 100
        else:
            return 0

    def each_10(n):
        return n%100//10
    
    n = 0
    b = 0
    k = 0
    m = 0
    t = 0
    for i in range(len(ls1)):
        key = ls1[i][0]
        b += count_100(n)
        p = each_10(n)
        k = 0
        for j in range(10):
            index = b + k + p
            e1, e2 = ls1[i][0], ls2[index][0]
            if e1 ==0 or e2==0:
                t += 1
            elif e1 > e2:
                m += 1                
            k+=10
            n +=1
    return m, t

def convert(ls):
    L = []
    for i in ls:
        L.append((i[0], 30000-i[1]))
    return L

def convert_com(compare):
    for i in compare:
        ori = compare[i]
        compare[i] = 10000-ori

def multi_comparison(l, result_dict):
    lst = ['owo', 'owith', 'cwo', 'cwith']
    score = {}
    compare = {}
    for i in lst:
        score[i] = 0
    
    owo = result_dict['o_wo_'+l]
    owith = result_dict['o_with_'+l]
    cwo = result_dict['c_wo_'+l]
    cwith = result_dict['c_with_'+l]
    
    s1, t = owo_1000(owo, owith)
    s2 = 10000 - s1 - t
    tot = s1 + s2
    s1 *= 10000/tot
    s2 *= 10000/tot
    s1 = int(s1)
    s2 = int(s2)
    if s1 + s2 == 10000:
        pass
    elif s1 > s2:
        s1+=1
    elif s1 < s2:
        s2+=1
    compare['owoowith'] = s1
    compare['owithowo'] = s2
    score['owo'] += s1
    score['owith'] += s2
    
    
    s1, t  = owo_1000(owo, cwo)
    s2 = 10000 - s1 - t
    tot = s1 + s2
    s1 *= 10000/tot
    s2 *= 10000/tot
    s1 = int(s1)
    s2 = int(s2)
    if s1 + s2 == 10000:
        pass
    elif s1 > s2:
        s1+=1
    elif s1 < s2:
        s2+=1      
    compare['owocwo'] = s1
    compare['cwoowo'] = s2
    score['owo'] += s1
    score['cwo'] += s2
    
    s1, t = owo_10000(owo, cwith)
    s2 = 10000 - s1 - t
    tot = s1 + s2
    s1 *= 10000/tot
    s2 *= 10000/tot
    s1 = int(s1)
    s2 = int(s2)
    if s1 + s2 == 10000:
        pass
    elif s1 > s2:
        s1+=1
    elif s1 < s2:
        s2+=1
    compare['owocwith'] = s1
    compare['cwithowo'] = s2
    score['owo'] += s1
    score['cwith'] += s2   
    
    s1, t = owith_cwo(owith, cwo)
    s2 = 10000 - s1 - t
    tot = s1 + s2
    s1 *= 10000/tot
    s2 *= 10000/tot
    s1 = int(s1)
    s2 = int(s2)
    if s1 + s2 == 10000:
        pass
    elif s1 > s2:
        s1+=1
    elif s1 < s2:
        s2+=1
    compare['owithcwo'] = s1
    compare['cwoowith'] = s2
    score['owith'] += s1
    score['cwo'] += s2   

    s1, t = owith_cwith(owith, cwith)
    s2 = 10000 - s1 - t
    tot = s1 + s2
    s1 *= 10000/tot
    s2 *= 10000/tot
    s1 = int(s1)
    s2 = int(s2)
    if s1 + s2 == 10000:
        pass
    elif s1 > s2:
        s1+=1
    elif s1 < s2:
        s2+=1
    compare['owithcwith'] = s1
    compare['cwithowith'] = s2
    score['owith'] += s1
    score['cwith'] += s2   
    
    s1, t = cwo_cwith(cwo, cwith)
    s2 = 10000 - s1 - t
    tot = s1 + s2
    s1 *= 10000/tot
    s2 *= 10000/tot
    s1 = int(s1)
    s2 = int(s2)
    if s1 + s2 == 10000:
        pass
    elif s1 > s2:
        s1+=1
    elif s1 < s2:
        s2+=1
    compare['cwocwith'] = s1
    compare['cwithcwo'] = s2
    score['cwo'] += s1
    score['cwith'] += s2 
    
    LS = sorted(score.items(), key=lambda x:x[1], reverse=True)
    
    if l.find('turn')!=-1 or l.find('recov')!= -1 or l.find('downward')!= -1:
        LS = convert(LS)
        convert_com(compare)
        LS.reverse()

    ls = []
    for i in range(len(LS)-1):
        key1 = LS[i][0] + LS[i+1][0]
        key2 = LS[i+1][0] + LS[i][0]
        try:
            k = compare[key1]
        except:
            k = compare[key2]
        if k<=5000:
            ls.append(10000-k)
        else:
            ls.append(k)

    return LS, compare

def obtain_N(l, result_dict, ranking_ls):
    LS, compare_x = multi_comparison(l, result_dict)
    P = 0
    for j in range(len(ranking_ls)):
        ls = ranking_ls[j]
        z = comb(ls, 2)

        p = 1
        for i in z:
            key = i[0]+i[1]
            proba = compare_x[key]/10000
            p *= proba
        P += p
    return P, compare_x

def P_per_rank(compare_x, N, ranking_ls):
    LS_ = []
    for j in range(len(ranking_ls)):
        ls = ranking_ls[j]
        z = comb(ls, 2)

        p = 1
        for i in z:
            key = i[0]+i[1]
            proba = compare_x[key]/10000
            p *= proba
        LS_.append((ranking_ls[j], p/N))
    return LS_


def obtain_ranking_prob_per_circuit(l, circuit, variables, normalize=True):
    P_ = []
    for j in range(4):
        p = 0
        if normalize == True:
            for i in variables.P_rank_N[l]:
                if i[0][j] == circuit:
                    p+=i[1]
            P_.append(p)
        else:
            for i in variables.P_rank[l]:
                if i[0][j] == circuit:
                    p+=i[1]
            P_.append(p)
    return P_

def gen_coordinate(P_):
    y_ls = []
    y = 0
    for i in P_:
        y+=i
        y_ls.append(y)

    x_ls = np.linspace(0, 3, 4)
    Lst = list(zip(x_ls, y_ls))
    return Lst

def cal_area(Lst):
    a = Polygon((0,0),(1,0),(2,0),(3,0),
                Lst[3], Lst[2], Lst[1], Lst[0]).area
    s = float(a.evalf())
    return s

def rescale(s):
    N = 3*(s-0.5)/2.5
    return N

def scores_for_radar(l, variables,normalize=True, rescaling=True):
    circuit_name = ['owo', 'owith', 'cwo','cwith']
    Ls = []
    for j in circuit_name:
        P_ = obtain_ranking_prob_per_circuit(l,j, variables, normalize=normalize)
        Lst = gen_coordinate(P_)
        if rescaling == True:
            score = rescale(cal_area(Lst))
        else:
            score = cal_area(Lst)
        Ls.append(score)
    return Ls 


class variables():
    suffix_n= ['efficiency_n', 'recover_n', 'downward_n', 'ratio_n',
                'turn_off_n','ratio_2_n','upward_n','turn_on_200_n']
    circuits = ['o_wo', 'o_with', 'c_wo','c_with']

    def __init__(self, result_dict, ranking_ls):
        P_rank = {}
        P_rank_N = {} 
        for i in self.suffix_n:
            N, compare_x = obtain_N(i, result_dict, ranking_ls)
            P_rank_N[i] = P_per_rank(compare_x, N, ranking_ls)
            P_rank[i] = P_per_rank(compare_x, 1, ranking_ls)
        
        self.P_rank = P_rank
        self.P_rank_N = P_rank_N 