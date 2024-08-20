from bubble_correlations import *

np_load_old = np.load
np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, **k)

### Params
nLat  = 8192
nTime = 512

phi0 = 2.*np.pi/3.6
temp = 0.
lamb = 1.5

#### CONSTANTS
nu     = 2e-3
lenLat = 1400./(2.*nu)**0.5

phieq  = np.pi
alph   = 8.

dx    = lenLat/nLat
dk    = 2.*np.pi/lenLat
knyq  = nLat//2+1

dt         = dx/alph
dtout      = dt*alph
light_cone = dtout//dx

# Lattice
lattice     = np.arange(nLat)
xlist       = lattice*dx
klist       = np.roll((lattice - nLat//2)*dk, nLat//2)
inv_phases  = np.exp(1j*np.outer(xlist, klist))
dit_phases  = nLat**-1. * np.exp(-1j*np.outer(xlist, klist))

#### SPECTRA
# Free field (constant mass term) field modes \phi_k
m2   = lambda la: 4.*nu*(-1.+la**2.)
norm = lambda ph0: 1./ ph0 / np.sqrt(2. * lenLat)
w2   = lambda la: klist**2. + m2(la)
free_eigenbasis  = lambda la, ph0: np.asarray([norm(ph0)/(w2(la)[k]**0.25) if kk!=0. else 0. for k,kk in enumerate(klist)])
free_pspec       = lambda la, ph0: np.abs(free_eigenbasis(la, ph0))**2.
free_fluct_stdev = lambda la, ph0: np.sqrt(np.sum(free_pspec(la, ph0)))

thermal_eigenbasis  = lambda la, ph0, te: free_eigenbasis(la, ph0) * np.sqrt(2./(np.exp(w2(la)**0.5/te)-1.))
thermal_pspec       = lambda la, ph0, te: np.abs(thermal_eigenbasis(la, ph0, te))**2.
thermal_fluct_stdev = lambda la, ph0, te: np.sqrt(np.sum(thermal_pspec(la, ph0, te)))

pspec       = lambda la, ph0, te: thermal_pspec(la, ph0, te) if te!=0 else free_pspec(la, ph0)
fluct_stdev = lambda la, ph0, te: thermal_fluct_stdev(la, ph0, te) if te!=0 else free_fluct_stdev(la, ph0)


### POTENTIAL
V    = lambda x, la:  ( -np.cos(x) + 0.5 * la**2. * np.sin(   x)**2. + 1. ) * 4. * nu
Vinv = lambda x, la: -( -np.cos(x) + 0.5 * la**2. * np.sin(   x)**2. + 1. ) * 4. * nu
dV   = lambda x, la:  (  np.sin(x) + 0.5 * la**2. * np.sin(2.*x)          ) * 4. * nu

### Paths to files
root_dir = '/gpfs/dpirvu/bubble_correlations/'
batch_params = lambda nL, la, ph, te: 'x'+str(nL)+'_phi0'+str('%.4f'%ph)+'_lambda'+str('%.4f'%la)+'_T'+str('%.4f'%te) 
sim_location = lambda nL, la, ph, te, sim: root_dir + batch_params(nL,la,ph,te) + '_sim'+str(sim)+'_fields.dat'

corrs_prefix = lambda sim, multi, window: '_sim'+str(sim)+'_multi'+str(multi)+'_filter{:.4f}'.format(window)
positive_data_file = lambda sim, multi, window: root_dir+'positive_targets'+corrs_prefix(sim, multi, window)+'.npy'
negative_data_file = lambda sim, multi, window: root_dir+'negative_targets'+corrs_prefix(sim, multi, window)+'.npy'

clean_sim_location = lambda nL, la, ph, te, sim: root_dir + 'clean_' + batch_params(nL,la,ph,te) + '_sim'+str(sim)+'_fields'
bubble_at_rest = lambda nL, la, ph, te, sim: root_dir + 'rest_bubble_' + batch_params(nL,la,ph,te) + '_sim'+str(sim)+'_fields'

velocities_bubbles_file = lambda nL, la, ph, te: root_dir + 'velocitiesCOM_' + batch_params(nL,la,ph,te)
average_bubble_file = lambda nL, la, ph, te: root_dir + 'average_bubble_' + batch_params(nL,la,ph,te)

triage_pref = lambda minS, maxS, nTM: '_minSim'+str(minS)+'_maxSim'+str(maxS)+'_up_to_nTMax'+str(nTM)
path_nodecay_sims  = lambda nL, la, ph, te, minS, maxS, nTM: root_dir+'sims_no_decay' + triage_pref(minS,maxS,nTM) + batch_params(nL,la,ph,te)
path_decayed_sims  = lambda nL, la, ph, te, minS, maxS, nTM: root_dir+'sims_good_decay' + triage_pref(minS,maxS,nTM) + batch_params(nL,la,ph,te)

#pickle_corr_location = lambda type: root_dir + 'th_correlator_type' + type
#thrcorr_file = lambda type, threshold: root_dir + 'thcorr_type' + str(type) + '_threshold' + str(threshold)+'.npy'

data_file_1D = lambda corrtype, xstep, deltat, multi, filter: root_dir+'corr1d_type'+str(corrtype)+'_xstep'+str(xstep)+'_deltat'+str(deltat)+'_multiplier'+str(multi)+'_filter{:.4f}'.format(filter)

data_file_2D = lambda corrtype, xstep, tstep, multi, filter: root_dir+'corr2d_type'+str(corrtype)+'_xstep'+str(xstep)+'_tstep'+str(tstep)+'_multiplier'+str(multi)+'_filter{:.4f}'.format(filter)


path_decaytimes = lambda la, ph, te: './data/tdecaylists_lamb'+str('%.4f'%la)+'_phi0'+str('%.4f'%ph)+'_temp'+str('%.4f'%te)
pspec_path = lambda la, ph, te: './data/powspec_lamb'+str('%.4f'%la)+'_phi0'+str('%.4f'%ph)+'_temp'+str('%.4f'%te)
en_path = lambda la, ph, te: './data/energy_lamb'+str('%.4f'%la)+'_phi0'+str('%.4f'%ph)+'_temp'+str('%.4f'%te)

instanton_location = '/home/dpirvu/inst/instantons/dev/bubcorr_instanton_sim.dat'

path_inst = lambda nL, la, ph, te: './data/instanton_profile'+batch_params(nL, la, ph, te)
path_indcrit = lambda nL, la, ph, te: './data/critical_timeslice'+batch_params(nL, la, ph, te)
path_encrit = lambda nL, la, ph, te: './data/critical_energy'+batch_params(nL, la, ph, te)

def import_all_data(minS, maxS, multiplier, tcheck, filter):
    targets_pos, targets_neg = [], []

    for sim in range(minS, maxS):
        simdatp = np.load(positive_data_file(sim, multiplier, filter))
        targets_pos.append(simdatp)
        
        simdatn = np.load(negative_data_file(sim, multiplier, filter))
        targets_neg.append(simdatn)

    return targets_pos, targets_neg

def import_all_data2(minS, maxS, multiplier, tcheck, filter):
    targets_pos, targets_neg = [], []

    for sim in range(minS, maxS):
        simdatp = np.load(positive_data_file(sim, multiplier, filter))
        for tgs in simdatp:
            targets_pos.append(tgs)
        
        simdatn = np.load(negative_data_file(sim, multiplier, filter))
        for tgs in simdatn:
            targets_neg.append(tgs)

    return targets_pos, targets_neg

titles = [r'$\phi(x)$', r'$\partial_t \phi(x)$', r'$|\nabla \phi(x)|^2$', r'$V(\phi(x))$']

list_title_type = [r'$\xi^{++}_{bb}(t,r)$', r'$\xi^{--}_{bb}(t,r)$', r'$\xi^{overall}_{bb}(t,r)$', r'$\xi^{+-}_{bb}(t,r)$', r'$\xi^{-+}_{bb}(t,r)$']


# Standard order or columns in .dat files is:
# field, momentum, gradient
normal = np.asarray([phieq, 0., 0.])

tcheck = int(m2(lamb)**-0.5/dtout)

right_Vmax  = sco.minimize_scalar(Vinv, args=lamb, bounds=(np.pi, 2*np.pi), method='bounded')
left_Vmax   = sco.minimize_scalar(Vinv, args=lamb, bounds=(0    ,   np.pi), method='bounded')

sigmafld    = fluct_stdev(lamb, phi0, temp)
crit_thresh = right_Vmax.x+2.*sigmafld
crit_rad    = 20


### DATA ANALYSIS Params
minSim = 0
maxSim = 1000

list_multiplier = [0., 0.5, 1.]
collision_thr = right_Vmax.x + 7.*sigmafld
list_threshold = [right_Vmax.x + mult*sigmafld for mult in list_multiplier]

filter_size = 5

tMinCutoff = 30
tMaxCutoff = nTime
tFutureCheck = 3
lightc = light_cone
