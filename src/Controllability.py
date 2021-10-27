from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from numpy import sqrt, arange
import numpy as np
import pylab as pl
from pylab import meshgrid, imshow, contour, clabel, title, show
from .Models_and_solvers import numerical_solver_trimer
from .Robustness import color_obj

def draw_dimer_without(lam1, lam2, lam3, K, c2, fold, normalize):
    c1 = np.linspace(0.000001*c2, fold*c2, 101)
    
    if normalize == True:
        p1 = -(lam1*lam2 + sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2)))*(K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(2*K*c2*lam1*(lam1*lam2 - sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2))))
        p2 = (lam1*lam2 + sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2)))*(K*c1*lam3 - K*c2*lam3 + lam1*lam2 - sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(2*K*c2*lam2*(lam1*lam2 - sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2))))
        p3 = (lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))*(K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/((-lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))*(K*c1*lam3 - K*c2*lam3 + lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
        
    else:
        p1 = (K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(2*K*lam1*lam3)
        p2 = (-K*c1*lam3 + K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(2*K*lam2*lam3)
        p3 = c2*(K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(lam3*(K*c1*lam3 - K*c2*lam3 + lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
        
    dimer_balance = c2*(-lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(lam3*(lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
        
    return p1, p2, p3, dimer_balance

def draw_dimer_response(lam1=0.027, lam2=0.027, lam3=0.027, K=0.05, c2=98, fold=2, lam12_equal=True, normalize=True, compare=True):
    
    if lam12_equal==True:
        lam2 = lam1      #lam1
    
    linewidth = 7
    
    font2 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 25}
    
    
    figure, ax = plt.subplots(figsize=(9, 6))
    
    if normalize == True:

        c1 = np.linspace(0.000001*c2, fold*c2, 101)
        x = np.linspace(0.000001*1, fold*100, 101)
 
        p3 = (lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))*(K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/((-lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))*(K*c1*lam3 - K*c2*lam3 + lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
    
        dimer_balance = c2*(-lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(lam3*(lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
            
        
        if compare == True:
            p1_control, p2_control, p3_control, ref_dimer = draw_dimer_without(lam1=0.027, lam2=0.027, lam3=0.027, K=0.05, c2=98,fold=fold, normalize=normalize)
            plt.plot(x, p3_control, linewidth=linewidth, color='k',label='Dimer w/o cooperative stability')
            parameter_without ='\nwithout: \u03BB1 = \u03BB3 = %.3f, K=%.3f, c2=%.2f'%(lam3, K, c2)
            print('Ref dimer amount: ', ref_dimer)

        plt.plot(x, p3, color=(0.1, 0.6, 0.1),linewidth=linewidth, label='Dimer with cooperative stability', zorder=2)   
       
        plt.xlabel('Relative $p_{1}$ synthesis rate (%)', font2, labelpad=10)
        plt.ylabel('Relative abundance', font2, labelpad=10)
        plt.ylim(0, 1.5)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    else:
        c1 = np.linspace(0.000001*c2, fold*c2, 101)

        p3 = c2*(K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(lam3*(K*c1*lam3 - K*c2*lam3 + lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
        
        dimer_balance = c2*(-lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(lam3*(lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
        
        if compare == True:
            p1_control, p2_control, p3_control, ref_dimer = draw_dimer_without(lam1=0.027, lam2=0.027, lam3=0.027, K=0.05, c2=98,fold=fold, normalize=normalize)

            plt.plot(c1, p3_control,'--', linewidth=3 ,label='dimer ref')
            parameter_without ='\nwithout: \u03BB1 = \u03BB3 = %.3f'%lam3
            print('Ref dimer amount: ', ref_dimer)
            

        plt.plot(c1, p3, color=(0.1, 0.5, 0.1),  linewidth=3 ,label='dimer')
        
        plt.xlabel('c1 % of c2', font2)
        plt.ylabel('Protein concentration (nM)', font2)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1000))
        plt.ylim(0, 3800)
        
    plt.legend(bbox_to_anchor=[1.3, 0.8])
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
    plt.tick_params(width=5, length=6, labelsize=20)
    plt.title('Heterodimer response', font2)
    ax1=plt.gca()
    
    ax_linewidth = 5
    
    ax1.spines['bottom'].set_linewidth(ax_linewidth)
    ax1.spines['left'].set_linewidth(ax_linewidth)
    ax1.spines['right'].set_linewidth(ax_linewidth)
    ax1.spines['top'].set_linewidth(ax_linewidth)
    
    plt.show()
    
    print('Control')
    print('Downward response (%):\t',(1 - (p3_control[50]-p3_control[25])/p3_control[50])*100)
    print('Upward response (%):\t',100+(p3_control[75]-p3_control[50])/p3_control[50]*100)

    print('With cooperative stability')
    print('Downward response (%):\t', (1 - (p3[50]-p3[25])/p3[50])*100)
    print('Upward response (%):\t', 100+(p3[75]-p3[50])/p3[50]*100)

class scenario:
    one = (0.13, 0.13, 0.13, 0.017, 0.017)       #all monomer
    two = (0.13, 0.13, 0.13, 0.13, 0.017)        #all intermediate(only second step)
    three = (0.13, 0.13, 0.13, 0.0735, 0.017)    #hierarchy (arithmetic) 
    four = (0.13, 0.13, 0.017, 0.017, 0.017)     #only first step
    five = (0.13, 0.13, 0.047, 0.047, 0.017)     #hierarchy (geometric)
    stable = (0.017, 0.017, 0.017, 0.017, 0.017) #all stable
    unstable = (0.13, 0.13, 0.13, 0.13, 0.13)    #all unstable
    test = (0.0017,0.0017,0.0017,0.0017,0.0017)
    test1 = (0.006,0.006,0.006,0.006,0.006)
    
def decide_tri(c_low, c_high, xmiddle, who):
    ls = ['c1', 'c2', 'c3']
    index = ls.index(who)
    globals()['_'+ls[index]] = c_low
    globals()[ls[index]+'_'] = c_high
    
    ls.remove(ls[index])    
        
    for i in ls:
        globals()['_'+i] = xmiddle
        globals()[i+'_'] = xmiddle
    
    return _c1, c1_, _c2, c2_, _c3, c3_

def p123_draw_trimer(scenario, code, xmiddle=98, who='c1',linewidth=3, color=color_obj.rbg[0]):
    lam1, lam2, lam3, lam12, lam123 = scenario
    
    _arr = np.linspace(0.000000001, xmiddle, 101)
    arr_ = np.linspace(xmiddle, 2*xmiddle, 101)
    
    c_low = _arr[::-1]
    c_high = arr_
    
    _c1, c1_, _c2, c2_, _c3, c3_ = decide_tri(c_low, c_high, xmiddle, who)
    
    _x, _y, _z, _p12, _p123 = numerical_solver_trimer(c1=_c1,c2=_c2,c3=_c3, lam1=lam1, lam2=lam2, lam3=lam3, lam12=lam12, lam123=lam123, ini_guess=(1, 1, 1))
    x_, y_, z_, p12_, p123_ = numerical_solver_trimer(c1=c1_,c2=c2_,c3=c3_,lam1=lam1, lam2=lam2, lam3=lam3, lam12=lam12, lam123=lam123, ini_guess=(1, 1, 1))

    plt.plot(np.array(list(_arr)+list(arr_))/xmiddle*100, np.array(_p123[::-1]+p123_)/p123_[0], linewidth=linewidth,label=code, color=color)

    print('(lam1, lam2, lam3, lam12, lam123)=', scenario)
    print('Downward response (%):\t', (1 - (_p123[0]-_p123[50])/_p123[0])*100)
    print('Upward response (%):\t',(p123_[50])/p123_[0]*100)

    
def draw_trimer_response(scenario_ls, xmiddle=98, who='c1'):
    
    linewidth=7
    ax_linewidth = 5
    figure, ax = plt.subplots(figsize=(9, 6))
    color_ls= ['k', color_obj.rbg[0]]
    
    #try:
    for i in range(len(scenario_ls)):
        p123_draw_trimer(scenario_ls[i], xmiddle=xmiddle, who=who, code='scenario '+str(i+1), linewidth=linewidth, color=color_ls[i])

    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
    plt.tick_params(width=ax_linewidth, labelsize=20,length=6)
    #plt.grid()
    
    font2 = {'family':'Times New Roman',
    'weight' : 'bold',
    'size'   : 25}
    
    ax1=plt.gca();#obtain the axis
    
    ax1.spines['bottom'].set_linewidth(ax_linewidth)
    ax1.spines['left'].set_linewidth(ax_linewidth)
    ax1.spines['right'].set_linewidth(ax_linewidth)
    ax1.spines['top'].set_linewidth(ax_linewidth)

    plt.title('Heteromeric trimer response', font2)
    plt.ylabel('Normalized Abundance', font2)
    plt.xlabel('Relative $p_%s$ synthesis rate (%%)'%who[-1], font2)
    plt.ylim(0, 1.5)
    plt.legend(bbox_to_anchor=(1.3, 1))
    plt.show()

class controllability:
    def downward_response(lam1, lam3, K, c2):
        lam2 = lam1 
        return (lam1*lam2 + sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2)))*(0.5*K*c2*lam3 + lam1*lam2 - sqrt(0.25*K**2*c2**2*lam3**2 + 3.0*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/((lam1*lam2 - sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2)))*(-0.5*K*c2*lam3 + lam1*lam2 + sqrt(0.25*K**2*c2**2*lam3**2 + 3.0*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
    def upward_response(lam1, lam3, K,  c2):
        lam2 = lam1 
        return -(lam1*lam2 + sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2)))*(0.5*K*c2*lam3 - lam1*lam2 + sqrt(0.25*K**2*c2**2*lam3**2 + 5.0*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/((lam1*lam2 - sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2)))*(0.5*K*c2*lam3 + lam1*lam2 + sqrt(0.25*K**2*c2**2*lam3**2 + 5.0*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
        
def controllability_heatmap(lam1=0.027, lam3=0.027,c_min=1,c_max=1000,k_min=0.001,k_max=1,resolution1 = 500, resolution2=500, up=True):
    
    font2 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 20}
        
    k_plot_log = np.linspace(np.log10(k_min), np.log10(k_max), resolution1)
    k_plot_linear = 10**k_plot_log

    c_plot_log = np.linspace(np.log10(c_min), np.log10(c_max), resolution2)
    c_plot_linear = 10**c_plot_log
    
    #generate x and y array within range
    x = c_plot_linear
    y = k_plot_linear
    #generate X and Y matrix for ploting in 3D or contour plot
    c_range, k_range = meshgrid(x, y)
    #evaluation of the function
    if up == True:
        Z = controllability.upward_response(c2=c_range, K=k_range, lam1=lam1, lam3=lam3)
        response_text = 'Upward '
    else:
        Z = controllability.downward_response(c2=c_range, K=k_range, lam1=lam1, lam3=lam3)
        response_text = 'Downward '
    Z = Z*100 

    mi = Z.min()
    ma = Z.max()
    diff = ma - mi

    #set the contour line stpes based on the range of the data
    if diff<=10:
        con_steps = 2
    else:
        con_steps = 10
    #norm = Normalize(vmin=0., vmax=1., clip=False)
    figure, ax = plt.subplots(figsize=(9, 7))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.tick_params(width=5, labelsize=15)
    ax.axis()
    
    ax1=plt.gca();#obtain the axis
    
    ax_linewidth = 3
    
    ax1.spines['bottom'].set_linewidth(ax_linewidth)
    ax1.spines['left'].set_linewidth(ax_linewidth)
    ax1.spines['right'].set_linewidth(ax_linewidth)
    ax1.spines['top'].set_linewidth(ax_linewidth)
    
    #drawing the function
    im = imshow(Z,cmap='coolwarm', aspect='auto',extent=[np.log10(c_min), np.log10(c_max), np.log10(k_min), np.log10(k_max)], origin='lower')

    #adding the Contour lines with labels
    cset = contour(Z,arange(0, ma, con_steps),extent=[np.log10(c_min), np.log10(c_max), np.log10(k_min), np.log10(k_max)],linewidths=2, cmap='binary_r')
    clabel(cset, colors='k',inline=True, fmt='%d', fontsize=20)
    #adding the colorbar on the right

    cb = plt.colorbar(im)
    cb.ax.tick_params(labelsize=25)
    title(response_text+'response (%%), $\u03BB_{monomer}/\u03BB_{dimer}=%d$\n'%(lam1/lam3), font2)
    pl.xlabel('$log_{10}$(synthesis rate) $(nM/hr)$', font2)
    pl.ylabel('$log_{10}K$ $(nM⁻¹)$', font2)

    show()

