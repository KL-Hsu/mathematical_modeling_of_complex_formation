from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from numpy import sqrt, arange
import numpy as np
import pylab as pl
from pylab import meshgrid, imshow, contour, clabel, title, show


# Coler code
class colors():
    def __init__(self):
        self.acedemia = ['k',(0.8, 0.8, 0.65), (0.2, 0.55, 0.75)]
        self.rbg = [(0.8, 0.2, 0.2),(0.3, 0.4, 0.7), (0.1, 0.6, 0.1)]
        self.obg = [(0.95, 0.5, 0.2),(0.3, 0.4, 0.7), (0.1, 0.6, 0.1)]

    def show(self, color_ls, name):
        plt.figure(figsize=(6, 1.5))
    
        plt.plot([1], [1], 'o', color=color_ls[0],markersize=70)
        plt.plot([2], [1], 'o', color=color_ls[1], markersize=70)
        plt.plot([3], [1], 'o', color=color_ls[2], markersize=70)
        
        plt.title(name)
        plt.axis('off')
        plt.xlim(0, 4)
        plt.ylim(0.75, 1.25)
        plt.show()
        
color_obj = colors()

def cal_tolerance_score(lam1, lam2, lam3, K, c2):
    upper_limit = c2*lam1/lam3 + c2 - lam1*lam2/(K*lam3) - lam2/K
    lower_limit = (K*c2*lam3 + lam1*lam2 + lam1*lam3)/(K*(lam2 + lam3))

    if lower_limit >= upper_limit:
        return 0 
    else:
        return np.log2(upper_limit/lower_limit)
    
def sqrt(x):
    return x**0.5


def draw_dimer_explore(lam1=0.027, lam2=0.027, lam3=0.027, K=0.05, c2=98, fold=2, lam12_equal=True, normalize=True):
    
    if lam12_equal==True:
        lam2 = lam1     
    
    linewidth = 7
    
    font2 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 25}
    
    figure, ax = plt.subplots(figsize=(9, 6))
    
    if normalize == True:

        c1 = np.linspace(0.000001*c2, fold*c2, 101)
        x = np.linspace(0.000001*1, 100*fold, 101)

        p1 = -(lam1*lam2 + sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2)))*(K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(2*K*c2*lam1*(lam1*lam2 - sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2))))
        p2 = (lam1*lam2 + sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2)))*(K*c1*lam3 - K*c2*lam3 + lam1*lam2 - sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(2*K*c2*lam2*(lam1*lam2 - sqrt(lam1*lam2*(4*K*c2*lam3 + lam1*lam2))))
        p3 = (lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))*(K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/((-lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))*(K*c1*lam3 - K*c2*lam3 + lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
    
        dimer_balance = c2*(-lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(lam3*(lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
            
        plt.plot(x, p3, color=color_obj.obg[2],linewidth=linewidth, label='Dimer', zorder=2)   
        plt.plot(x, p1, color=color_obj.obg[1],linewidth=linewidth, label='p1 monomer')
        plt.plot(x, p2, color=color_obj.obg[0],linewidth=linewidth, label='p2 monomer')
       
        plt.xlabel('Relative $p_{1}$ synthesis rate (%)', font2, labelpad=10)
        plt.ylabel('Relative abundance', font2, labelpad=10)
        plt.ylim(0, 1.5)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    else:
        c1 = np.linspace(0.000001*c2, fold*c2, 101)
        x = np.linspace(0.000001*1, 100*fold, 101)
        
        p1 = (K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(2*K*lam1*lam3)
        p2 = (-K*c1*lam3 + K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(2*K*lam2*lam3)
        p3 = c2*(K*c1*lam3 - K*c2*lam3 - lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(lam3*(K*c1*lam3 - K*c2*lam3 + lam1*lam2 + sqrt(K**2*c1**2*lam3**2 - 2*K**2*c1*c2*lam3**2 + K**2*c2**2*lam3**2 + 2*K*c1*lam1*lam2*lam3 + 2*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))
        
        dimer_balance = c2*(-lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2))/(lam3*(lam1*lam2 + sqrt(4*K*c2*lam1*lam2*lam3 + lam1**2*lam2**2)))

        plt.plot(x, p3, color=color_obj.obg[2],  linewidth=linewidth ,label='Dimer', zorder=2)  
        plt.plot(x, p1, color=color_obj.obg[1], linewidth=linewidth ,label='p1 monomer')
        plt.plot(x, p2, color=color_obj.obg[0],  linewidth=linewidth ,label='p2 monomer')
      
        plt.xlabel('Relative $p_{1}$ synthesis rate (%)', font2)
        plt.ylabel('Protein concentration (nM)', font2)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1000))
        plt.ylim(0, 3800)
   
    plt.legend(bbox_to_anchor=[1.3, 0.8])
    plt.tick_params(width=5, length=6, labelsize=20)    
    plt.title('Heterodimer', font2)
    
    ax1=plt.gca();#obtain the axis
    
    ax_linewidth = 5
    
    ax1.spines['bottom'].set_linewidth(ax_linewidth)
    ax1.spines['left'].set_linewidth(ax_linewidth)
    ax1.spines['right'].set_linewidth(ax_linewidth)
    ax1.spines['top'].set_linewidth(ax_linewidth)
    
    plt.show()
    
    lower_limit = (K*c2*lam3 + lam1*lam2 + lam1*lam3)/(K*(lam2 + lam3))
    upper_limit = c2*lam1/lam3 + c2 - lam1*lam2/(K*lam3) - lam2/K
    tolerance_score = cal_tolerance_score(lam1, lam2, lam3, K, c2)
    
    print('Normalized dimer (nM):\t', dimer_balance)
    print('Upper limit:', upper_limit/c2)
    print('Lower limit:', lower_limit/c2)
    print('Tolerance score:', tolerance_score)

def tolerance_score(c2=30, lam1=0.027, lam3=0.027, K=0.02):
    lam2 = lam1
    lower = ((K*c2*lam3 + lam1*lam2 + lam1*lam3)/(K*(lam2 + lam3)))/c2
    upper = (c2*lam1/lam3 + c2 - lam1*lam2/(K*lam3) - lam2/K)/c2
    
    #label those beyond the limit
    lower[lower>1] = 2**1000
    upper[upper<1] = 2**1000
    
    log2_lower = np.log2(lower)
    log2_upper = np.log2(upper)
    
    #Transfer label into zero
    log2_lower[log2_lower==1000] = 0
    log2_upper[log2_upper==1000] = 0
    
    #If beyond the limits, return 0
    if np.all(log2_lower==0) or np.all(log2_upper==0):
        R = log2_lower*log2_upper
    else:
        R = log2_upper - log2_lower

    return R

def tolerance_score_heatmap(lam1=0.027, lam3=0.027,c_min=1,c_max=1000,k_min=0.001,k_max=1,resolution1 = 500, resolution2=500):
    
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
    Z = tolerance_score(c2=c_range, K=k_range, lam1=lam1, lam3=lam3)

    lst = []
    for i in Z:
        for ii in i:
            if ii !=0:
                lst.append(ii)

    mi = min(lst)
    ma = max(lst)
    diff = ma - mi

    #set the contour line stpes based on the range of the data
    con_steps = 2
  
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
    pl.clim(0, 7)
    #adding the Contour lines with labels
    cset = contour(Z,arange(0, 7, con_steps),extent=[np.log10(c_min), np.log10(c_max), np.log10(k_min), np.log10(k_max)],linewidths=2, cmap='binary_r')
    
    clabel(cset, colors='k',inline=True, fmt='%d', fontsize=20)
    #adding the colorbar on the right

    cb = plt.colorbar(im)
    cb.ax.tick_params(labelsize=25)
    title('Heatmap of Tolerance Scores, $\u03BB_{monomer}/\u03BB_{dimer}=%d$\n'%(lam1/lam3), font2)
    pl.xlabel('$log_{10}$(synthesis rate) $(nM/hr)$', font2)
    pl.ylabel('$log_{10}K$ $(nM⁻¹)$', font2)

    show()

def tolerance_score_heatmap3D(lam_min=0.027, lam_max=0.7, k_min=0.001, k_max=1, c2=10, offset=0, resolution=100):
    font1 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 25}
    
    font2 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 18}


    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')

    k_plot_log = np.linspace(np.log10(k_min), np.log10(k_max), resolution)
    k_plot_linear = 10**k_plot_log


    #generate x and y array within range
    x = np.linspace(lam_min, lam_max, resolution)
    y = np.logspace(np.log10(k_min), np.log10(k_max), resolution)

    #generate X and Y matrix for ploting in 3D or contour plot
    lam_range, k_range = np.meshgrid(x, y)
    #evaluation of the function
    Z = tolerance_score(lam1=lam_range, K=k_range, c2=c2)
    lam_ratio = lam_range/lam_min
    
    plt.title('Tolerance Score\n', font1)

    ax.plot_surface(lam_ratio, np.log10(k_range), Z, linewidth=0,cmap=plt.get_cmap('jet'),zorder=5,
                   rcount=100, ccount=100)    

    plt.xlabel('$\u03BB_{mono}/\u03BB_{dimer}$',font2)
    plt.ylabel('$log_{10}(K) (nM⁻¹)$', font2)
    ax.set_zlabel('Tolerance Score',font2)

    ax.set_zlim(0, 10)
    ax.view_init(30, -120)

    plt.show()