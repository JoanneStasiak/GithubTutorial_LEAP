import os
import pandas as pd
from mne.stats import permutation_cluster_test
import numpy as np

#####################
#########################################################################################################

trev=pd.read_csv("/Volumes/labshare/Joanne/ASAP/AE_physio_data/TREV/processed/may25/allSubs_trevCountdown17s_interpolated_noOutliers.csv")
timepoints = np.sort(trev['time'].unique())

condition_arrays = {}
for cond, group in trev.groupby('condition'):
    trial_matrices = []
    for (sub, run, trial), trial_df in group.groupby(['sub', 'run', 'trial']):
        trial_series = trial_df.sort_values('time')['z_contractility'].values
        if len(trial_series) == len(timepoints):  # Ensure consistency
            trial_matrices.append(trial_series)
    if trial_matrices:
        condition_arrays[cond] = np.array(trial_matrices)
min_trials = min([arr.shape[0] for arr in condition_arrays.values()])
aligned_conditions = {k: v[:min_trials] for k, v in condition_arrays.items()}
aligned_data = [aligned_conditions[k] for k in sorted(aligned_conditions.keys())]
aligned_shapes = {k: v.shape for k, v in aligned_conditions.items()}
aligned_shapes

T_obs, clusters, cluster_p_values, H0 = permutation_cluster_test(
    aligned_data,
    n_permutations=1000,
    tail=0,  # two-tailed
    threshold=None,
    out_type='mask'
)
results = {
    'T_obs': T_obs,
    'clusters': clusters,
    'p_values': cluster_p_values,
    'timepoints': timepoints
}
results['significant_clusters'] = np.where(cluster_p_values < 0.05)[0]
results['n_significant'] = len(results['significant_clusters'])
results['significant_clusters']

#########################################################################################################
trev=pd.read_csv("/Volumes/labshare/Joanne/ASAP/AE_physio_data/TREV/processed/may25/allSubs_trevCountdown17s_interpolated_noOutliers.csv")
timepoints = np.sort(trev['time'].unique())

shock_arrays = {}
for shockz, group in trev.groupby('shock'):
    trial_matrices = []
    for (sub, run, trial), trial_df in group.groupby(['sub', 'run', 'trial']):
        trial_series = trial_df.sort_values('time')['z_contractility'].values
        if len(trial_series) == len(timepoints):  # Ensure consistency
            trial_matrices.append(trial_series)
    if trial_matrices:
        shock_arrays[shockz] = np.array(trial_matrices)
min_trials = min([arr.shape[0] for arr in shock_arrays.values()])
aligned_conditions = {k: v[:min_trials] for k, v in shock_arrays.items()}
aligned_data = [aligned_conditions[k] for k in sorted(aligned_conditions.keys())]
aligned_shapes = {k: v.shape for k, v in aligned_conditions.items()}
aligned_shapes

T_obs, clusters, cluster_p_values, H0 = permutation_cluster_test(
    aligned_data,
    n_permutations=1000,
    tail=0,  # two-tailed
    threshold=None,
    out_type='mask'
)
results = {
    'T_obs': T_obs,
    'clusters': clusters,
    'p_values': cluster_p_values,
    'timepoints': timepoints
}
results['significant_clusters'] = np.where(cluster_p_values < 0.05)[0]
results['n_significant'] = len(results['significant_clusters'])
results['significant_clusters']

#################### plot timeseries of F statistics 
import matplotlib.pyplot as plt
T_obs = results['T_obs']
clusters = results['clusters']
p_values = results['p_values']
times = results['timepoints']

plt.figure(figsize=(9, 5))
plt.plot(times, T_obs, label='Observed F-values', color='black')
for i_c, (cluster_mask,) in enumerate(clusters):  # each cluster is a tuple (mask,)
    if p_values[i_c] < 0.05:
        time_indices = np.where(cluster_mask)[0]
        if len(time_indices) > 0:
            start_time = times[time_indices[0]]
            end_time = times[time_indices[-1]]
            plt.axvspan(start_time, end_time, color='red', alpha=0.3,
                        label='Significant cluster' if i_c == results['significant_clusters'][0] else None)

plt.title('Cluster-based Permutation Test (Shock Conditions)')
plt.xlabel('Time (s)')
plt.ylabel('F-statistic')
plt.legend()
plt.tight_layout()
plt.show()


plt.figure(figsize=(9, 5))
plt.plot(times, T_obs, label='Observed F-values', color='black')

for i_c, cluster in enumerate(clusters):
    if p_values[i_c] < 0.05:
        start_idx = cluster[0].start
        end_idx = cluster[0].stop
        start_time = times[start_idx]
        end_time = times[end_idx - 1]  # inclusive end
        plt.axvspan(start_time, end_time, color='#cbb5e2', alpha=0.3,
                    label='Significant cluster' if i_c == 0 else None)

plt.title('Cluster-based Permutation Test (Shock Conditions)')
plt.xlabel('Time (s)')
plt.ylabel('F-statistic')
plt.legend()
plt.tight_layout()
plt.show()
## from time 8.5 - 9.7s
################################# plot timeseries of shock conditions
mean_timecourses = trev.groupby(['shock', 'time'])['z_contractility'].mean().unstack(level=0)
plt.figure(figsize=(9, 5))
for shock_type in mean_timecourses.columns:
    plt.plot(mean_timecourses.index, mean_timecourses[shock_type], label=shock_type)

plt.title("Average z_contractility Over Time by Shock Condition")
plt.xlabel("Time (s)")
plt.ylabel("z_contractility")
plt.tight_layout()
plt.show()

############################
########################control#################################################################################
trev=pd.read_csv("/Volumes/labshare/Joanne/ASAP/AE_physio_data/TREV/processed/may25/allSubs_trevCountdown17s_interpolated_noOutliers.csv")
timepoints = np.sort(trev['time'].unique())

control_arrays = {}
for controlz, group in trev.groupby('control'):
    trial_matrices = []
    for (sub, run, trial), trial_df in group.groupby(['sub', 'run', 'trial']):
        trial_series = trial_df.sort_values('time')['z_contractility'].values
        if len(trial_series) == len(timepoints):  # Ensure consistency
            trial_matrices.append(trial_series)
    if trial_matrices:
        control_arrays[controlz] = np.array(trial_matrices)
min_trials = min([arr.shape[0] for arr in control_arrays.values()])
aligned_conditions = {k: v[:min_trials] for k, v in control_arrays.items()}
aligned_data = [aligned_conditions[k] for k in sorted(aligned_conditions.keys())]
aligned_shapes = {k: v.shape for k, v in aligned_conditions.items()}
aligned_shapes

T_obs, clusters, cluster_p_values, H0 = permutation_cluster_test(
    aligned_data,
    n_permutations=1000,
    tail=0,  # two-tailed
    threshold=None,
    out_type='mask'
)
results = {
    'T_obs': T_obs,
    'clusters': clusters,
    'p_values': cluster_p_values,
    'timepoints': timepoints
}
cluster_p_values
################################# plot timeseries of control conditions
mean_timecourses = trev.groupby(['control', 'time'])['z_contractility'].mean().unstack(level=0)
plt.figure(figsize=(9, 5))
for control_type in mean_timecourses.columns:
    plt.plot(mean_timecourses.index, mean_timecourses[control_type], label=control_type)

plt.xlabel("Time (s)")
plt.ylabel("z_contractility")
plt.tight_layout()
plt.show()

#########################################
##### Make new df with averaged contractility during the significant cluster:

trev=pd.read_csv("/Volumes/labshare/Joanne/ASAP/AE_physio_data/TREV/processed/may25/allSubs_trevCountdown17s_interpolated_noOutliers.csv")

trev_cluster = trev[(trev['time'] >= 8.5) & (trev['time'] <= 9.7)]
# Group by sub, run, trial and compute mean z_contractility
cluster_avg_df = trev_cluster.groupby(['sub', 'run', 'trial','condition','shock', 'control'], as_index=False)['z_contractility'].mean()

# Rename for clarity
cluster_avg_df.rename(columns={'z_contractility': 'z_contractility_cluster'}, inplace=True)
cluster_avg_df.to_csv("/Volumes/labshare/Joanne/ASAP/AE_physio_data/TREV/processed/may25/cluster_avg_z_contractility.csv", index=False)

#### adding secret text down here to show version control!