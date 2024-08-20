import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.ticker as mticker
from matplotlib import gridspec
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable


f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.1e'%x))
fmt = mticker.FuncFormatter(g)

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

plt.rcParams.update({'backend' : 'Qt5Agg'})
plt.rcParams.update({'text.usetex' : True})

plt.rcParams.update({'font.size' : 11.0})
plt.rcParams.update({'axes.titlesize' : 14.0})  # Font size of title
plt.rcParams.update({'axes.titlepad'  : 10.0})
plt.rcParams.update({'axes.labelsize' : 14.0})  # Axes label sizes
plt.rcParams.update({'axes.labelpad'  : 10.0})
plt.rcParams.update({'xtick.labelsize'  : 14.0})
plt.rcParams.update({'ytick.labelsize'  : 14.0})
plt.rcParams.update({'xtick.labelsize'  : 10.0})
plt.rcParams.update({'ytick.labelsize'  : 10.0})

plt.rcParams.update({'axes.spines.left'  : True})
plt.rcParams.update({'axes.spines.right'  : True})
plt.rcParams.update({'axes.spines.top'  : True})
plt.rcParams.update({'axes.spines.bottom'  : True})
plt.rcParams.update({'savefig.format'     : 'pdf'})
plt.rcParams.update({'savefig.bbox'       : 'tight'})
plt.rcParams.update({'savefig.pad_inches' : 0.1})
plt.rcParams.update({'pdf.compression' : 6})
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#672766', '#a6228c', '#ed7625', '#1b988c', '#75bc43'])
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#b974c6', '#510c5e', '#f96a35', '#96ca4f', '#f45523'])
#plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#001219', '#005f73', '#0a9396', '#94d2bd', '#e9d8a6', '#ee9b00', '#ca6702', '#bb3e03', '#ae2012', '#9b2226'])


tableau1 = [(50, 162, 81), (172,217,141), (255,127,15), (255,185,119), (60,183,204), (152,217,228), (184,90,13), (255,217,74), (57,115,124), (134,180,169), (130,133,59), (204,201,77)]
for i in range(len(tableau1)):
    r, g, b = tableau1[i]    
    tableau1[i] = (r / 255., g / 255., b / 255.)

tableau2 = [(44,105,176), (181,200,226), (240,39,32), (255,182,176), (172,97,60), (233,195,155), (107,163,214), (181,223,253), (172,135,99), (221,201,180), (189,10,54), (244,115,122)]
for i in range(len(tableau2)):
    r, g, b = tableau2[i]    
    tableau2[i] = (r / 255., g / 255., b / 255.)

linestyle = ['solid', 'dashed', 'dotted']
tabl = [tableau2[0],tableau2[2],tableau1[0]]
tabl2 = [tableau2[1],tableau2[3],tableau1[1]]


def simple_imshow(simulation, extent_floats, title='Sim', cmp='PRGn'):
    fig, ax0 = plt.subplots(1, 1, figsize=(6,5))
    im0 = plt.imshow(simulation, aspect='auto', interpolation='none', extent=extent_floats, cmap=cmp, origin='lower')
    clb = plt.colorbar(im0, ax = ax0)
    ax0.set_title(title)
    return ax0


def simulation_imshow(simulation, extent_floats, cmp='PRGn'):
    fig, ax0 = plt.subplots(1, 1, figsize=(15,4))
    im0 = plt.imshow(simulation, aspect='auto', interpolation='none', extent=extent_floats, cmap=cmp, origin='lower')
    return ax0


def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter


    
def plot_dmdz(ms, zs, func, count=10, title='No Title'):
    fig, ax = plt.subplots(2, 1, figsize=(12, 10))
    for mi, mm in enumerate(ms):
        if mm in ms[:count]:
            ax[0].plot(zs, func[:,mi], linewidth=1, label=('$m=${}'.format(fmt(mm)) if mm in ms[:count//count] else None), color='k')
        elif mm in ms[-count:]:
            ax[0].plot(zs, func[:,mi], linewidth=1, label=('$m=${}'.format(fmt(mm)) if mm in ms[-count//count:] else None), color='g')
    for zi, zz in enumerate(zs):
        if zz in zs[:count]:
            ax[1].plot(ms, func[zi,:], linewidth=1, label=('$z=${}'.format(fmt(zz)) if zz in zs[:count//count] else None), color='r')
        elif zz in zs[-count:]:
            ax[1].plot(ms, func[zi,:], linewidth=1, label=('$z=${}'.format(fmt(zz)) if zz in zs[-count//count:] else None), color='b')

    ax[0].set_xlabel('z')
    ax[1].set_xlabel('m')
    for axx in ax:
        axx.set_yscale('log')
        axx.set_xscale('log')
        axx.set_ylabel(title)
        axx.legend(); axx.grid()
    return ax

def plot_ucosth(ms, zs, angs, ucosth, prob, title, count=10):
    fig, ax = plt.subplots(2, 1, figsize=(12, 10))
    nMs, nZs = len(ms), len(zs)
    for mi, mm in enumerate(ms):
        if mm in ms[:count]:
            lab = lambda zi: (r'$m=${}'.format(fmt(mm))+f', $z=%5.2f$'%(zs[zi]) if mi==0 else None)
            zi = 0
            ax[0].plot(angs[:,zi,mi], ucosth[:,zi,mi], ms=3, marker='o', label=lab(zi), color='k')
            zi = nZs-1
            ax[0].plot(angs[:,zi,mi], ucosth[:,zi,mi], ms=3, marker='o', label=lab(zi), color='r')
        elif mm in ms[-count:]:
            lab = lambda zi: (r'$m=${}'.format(fmt(mm))+f', $z=%5.2f$'%(zs[zi]) if mi==len(ms)-1 else None)
            zi = 0
            ax[0].plot(angs[:,zi,mi], ucosth[:,zi,mi], ms=3, marker='o', label=lab(zi), color='g')
            zi = nZs-1
            ax[0].plot(angs[:,zi,mi], ucosth[:,zi,mi], ms=3, marker='o', label=lab(zi), color='b')

    for mi, mm in enumerate(ms):
        if mm in ms[:count]:
            lab = lambda zi: (r'$m=${}'.format(fmt(mm))+f', $z=%5.2f$'%(zs[zi]) if mi==0 else None)
            zi = 0
            ax[1].plot(angs[:,zi,mi], ucosth[:,zi,mi]*prob[zi,mi], ms=3, marker='*', label=lab(zi), color='k')
            zi = nZs-1
            ax[1].plot(angs[:,zi,mi], ucosth[:,zi,mi]*prob[zi,mi], ms=3, marker='*', label=lab(zi), color='r')
        elif mm in ms[-count:]:
            lab = lambda zi: (r'$m=${}'.format(fmt(mm))+f', $z=%5.2f$'%(zs[zi]) if mi==len(ms)-1 else None)
            zi = 0
            ax[1].plot(angs[:,zi,mi], ucosth[:,zi,mi]*prob[zi,mi], ms=3, marker='*', label=lab(zi), color='g')
            zi = nZs-1
            ax[1].plot(angs[:,zi,mi], ucosth[:,zi,mi]*prob[zi,mi], ms=3, marker='*', label=lab(zi), color='b')

    for axx in ax:
        axx.set_yscale('log'); axx.set_xscale('log')
        axx.set_ylabel(r'$u(\cos(\theta))$'); axx.set_xlabel(r'$\theta$')
        axx.legend(); axx.grid()
    return ax

#Label line with line2D label data
def labelLine(line,x,label=None,align=True,**kwargs):

    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    #Find corresponding y co-ordinate and angle of the line
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1])*(x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = line.get_label()

    if align:
        #Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy,dx))

        #Transform to screen co-ordinates
        pt = np.array([x,y]).reshape((1,2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)),pt)[0]

    else:
        trans_angle = 0

    #Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    ax.text(x,y,label,rotation=trans_angle,**kwargs)

def labelLines(lines,align=True,xvals=None,**kwargs):

    ax = lines[0].axes
    labLines = []
    labels = []

    #Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin,xmax = ax.get_xlim()
        xvals = np.linspace(xmin,xmax,len(labLines)+2)[1:-1]

    for line,x,label in zip(labLines,xvals,labels):
        labelLine(line,x,label,align,**kwargs)

def annot_max(xmax, ymax, lab, col, xcap, mind, ax):
    text = lab
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec=col, lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60",color=col)
    kw = dict(xycoords='data',textcoords="data",
              arrowprops=arrowprops, bbox=bbox_props, ha="left", va="top")
    ax.annotate(text, xy=(xmax, ymax), xytext=(xmax-0.01, ymax-5e6), **kw)

class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))
    

def plot_slices(t, areal):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    tabl = cycle((tableau1[::]))
    line = cycle(('solid', 'dashed', 'dashdot', 'dotted'))
    i0 = i
    for slice in areal:
        plt.plot(slice, label=r'$\tau = \,$'+str((i0-i)//t+1), color=next(tabl), linestyle=next(line))
        i0 = i0 + t
#    labelLines(plt.gca().get_lines(), xvals=(2*len(slice)/3, len(slice)), align=False)
    plt.grid(alpha=0.8, linestyle='dashed', linewidth=0.5)
    plt.legend(); plt.xlabel( r'$\phi_0^{-1} \sqrt{V_0} \; r$'); plt.ylabel(r'$\bar{\phi}$')
    ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 2))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
    plt.savefig('./bubble_evolution.pdf')
    return

   

def animate_field(bubble, step, filename='Animation'):
    nT, nN = np.shape(bubble)
    dx2plot = np.sqrt(4*nu)*dx
    xlist = np.arange(nN)
    tlist = np.arange(0, nT, step)
    philist = np.linspace(1.4*np.pi/2., 7.4, 100)
    
    fig, ax = plt.subplots(1, 2, figsize = (9, 4.5), dpi=(1920/16))
    camera = Camera(fig)
    for tt in tlist:
        slice = ax[0].plot(xlist*dx2plot, bubble[tt], color='b', linewidth=2, ls='-')
        #ax[1].errorbar(np.mean(bubble[tt]), V(np.mean(bubble[tt]), 2), xerr=np.std(bubble[tt]), ecolor='r', elinewidth=2)
        
        ax[1].plot(philist, V(philist, 2), color='k', linewidth=2, ls='-')
        xerr1 = np.linspace(np.mean(bubble[tt]), np.mean(bubble[tt])-np.std(bubble[tt]), 100)
        xerr2 = np.linspace(np.mean(bubble[tt]), np.mean(bubble[tt])+np.std(bubble[tt]), 100)
        ax[1].plot(xerr1, V(xerr1, 2), color='r', linewidth=2)
        ax[1].plot(xerr2, V(xerr2, 2), color='r', linewidth=2)
        ax[1].plot(np.mean(bubble[tt]), V(np.mean(bubble[tt]), 2), color='b', marker='o')

        for aa in ax:
            aa.grid(ls=':', color='darkgray', alpha=0.5)
        ax[0].legend(slice, [f't = {tt}'], loc=1)
        ax[0].yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
        ax[0].yaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
        ax[1].xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
        ax[1].xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
        ax[1].get_yaxis().set_ticks([])
        ax[0].get_xaxis().set_ticks([])

        ax[0].set_xlabel(r'$r$')
        ax[0].set_ylabel(r'$\phi(r)$')
        ax[1].set_xlabel(r'$\phi(r)$')
        ax[1].set_ylabel(r'$V(\phi)$')
        camera.snap()
    animation = camera.animate(interval = 0.05);
    animation.save('./'+filename+'.gif', writer = 'imagemagick')
    return
