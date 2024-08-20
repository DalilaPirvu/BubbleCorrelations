from bubble_correlations import *
from experiment import *

for sim in range(minSim, maxSim):
    print('Simulation:', sim)
    path_sim = sim_location(nLat, lamb, phi0, temp, sim)
    real = extract_data(nLat, path_sim)
    nC, nT, nN = np.shape(real)

    smreal = np.asarray([gaussian_filter1d(real[0,tt], sigma=filter_size, mode='wrap') for tt in range(nT)])
    anti_smreal = 2.*phieq - smreal

    for mind, multiplier in enumerate(list_multiplier):
        threshold = right_Vmax.x + multiplier * sigmafld

        targets_pos = identify_bubble_sites(smreal, tMinCutoff, tMaxCutoff, tFutureCheck, light_cone, crit_thresh, threshold, collision_thr)
        targets_neg = identify_bubble_sites(anti_smreal, tMinCutoff, tMaxCutoff, tFutureCheck, light_cone, crit_thresh, threshold, collision_thr)

        np.save(positive_data_file(sim, multiplier, filter_size), targets_pos)
        np.save(negative_data_file(sim, multiplier, filter_size), targets_neg)

print('Done.')
