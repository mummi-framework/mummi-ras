#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------

import logging
import numpy as np
import glob
import os
import random
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import distance
from itertools import combinations
import MDAnalysis as mda
from MDAnalysis.analysis import align
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

LOGGER = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
def load_validation_data(tpath, glb_rank=-1, keys_to_load = []):

    # need to function to decide whether to read a key or not
    read_key = lambda k: k in keys_to_load if len(keys_to_load) > 0 else True

    # first, figure out of all files need to be loaded or just for one rank
    if glb_rank != -1:
        fnames = [os.path.join(tpath, f'validation_data_r{glb_rank}.npz')]
    else:
        fnames = sorted(glob.glob(os.path.join(tpath, 'validation_data_*.npz')))

    # read all files
    LOGGER.info(f'Loading {len(fnames)} files from ({tpath})')
    vdata = {k: [] for k in keys_to_load}
    for f in fnames:
        with np.load(f, allow_pickle=True) as npz:
            for k in np.intersect1d(npz.files, keys_to_load):
                vdata[k].append(npz[k])

    vdata = {k: np.concatenate(v) for k,v in vdata.items()}
    return vdata

def plot_labeled_scatters(pos, labels, colors, ax=None):

    if np.unique(labels).shape[0] > len(colors):
        LOGGER.info('please add more colors')

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(12, 12))

    for _l, _c in zip(np.unique(labels), colors):
        _ids = labels == _l
        _lcnt = np.count_nonzero(_ids)
        ax.scatter(pos[_ids, 0], pos[_ids, 1], c=_c, alpha=0.5, label='{:02d}={}'.format(_l, _lcnt))

    ax.legend()
    return ax

def plot_labeled_scatters3(pos, labels, colors, ax=None):

    if np.unique(labels).shape[0] > len(colors):
        LOGGER.info('please add more colors')

    if ax is None:
        fig= plt.subplots(1, 1, figsize=(12, 12))
        ax = plt.axes(projection ="3d")

    for _l, _c in zip(np.unique(labels), colors):
        _ids = labels == _l
        _lcnt = np.count_nonzero(_ids)
        ax.scatter3D(pos[_ids, 0], pos[_ids, 1], pos[_ids, 2], c=_c, alpha=0.25, label='{:02d}={}'.format(_l, _lcnt))

    ax.legend()
    return ax
# ------------------------------------------------------------------------------
# PLOT PCA
# ------------------------------------------------------------------------------
def plot_pca(encodings, label1, filename, filename_xz, filename_xy, filename_yz, title, show):
    colors1 = ['k',
               '#a6cee3', 'orange', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
               '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928',
               '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462',
               '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f',
               'orange', 'r']
    colors2 = ['k', '#a6cee3', '#b2df8a', '#fb9a99']
    colors = ['k', '#a6cee3', '#b2df8a']

    from sklearn.decomposition import PCA

    pca = PCA(n_components=3)
    x = pca.fit_transform(encodings)

    fig, axs = plt.subplots(1, 1, figsize=(16, 8))

    axs = plot_labeled_scatters3(x, labels1, colors1)


    plt.xlabel('PC 1')
    plt.ylabel('PC 2')
    axs.set_zlabel('PC 3')

    if len(title) > 0:
        fig.suptitle(title)

    plt.savefig(filename, bbox_inches='tight')
    axs.view_init(elev=0, azim=-90)
    plt.savefig(filename_xz, bbox_inches='tight')
    axs.view_init(elev=90, azim=-90)
    plt.savefig(filename_xy, bbox_inches='tight')
    axs.view_init(elev=0, azim=0)
    plt.savefig(filename_yz, bbox_inches='tight')

    if show:
        plt.show()
    plt.close()
# ------------------------------------------------------------------------------
def plot_project_pca(encodings, encodings_new, timeframe_ids, labels1, labels_new, uniq_structure_ids, filename, filename_xz, filename_xy, filename_yz, title, show):
    colors1 = ['cyan',
               'hotpink', 'orange', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
               '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928',
               '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462',
               '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f',
               'orange', 'r']
    colors2 = ['k', '#a6cee3', '#b2df8a', '#fb9a99']
    colors = ['k', '#a6cee3', '#b2df8a']
    cmaps = ['Blues', 'Greens', 'Greys', 'Oranges', 'Purples', 'Reds']
    # cmaps = ['Greens', 'Greys', 'Oranges', 'Purples', 'Reds']

    from sklearn.decomposition import PCA

    pca = PCA(n_components=3)

    x = pca.fit_transform(encodings)
    print(pca.explained_variance_ratio_)
    # exit(1)
    x_new = pca.transform(encodings_new)

    LOGGER.info(f'x {x.shape} x_new {x_new.shape}')
    # x = np.append(x, x_new, axis=0)
    LOGGER.info(f'x {x.shape} x_new {x_new.shape}')
    LOGGER.info(f'labels1 {labels1.shape} labels_new {labels_new.shape}')
    # labels1 = np.append(labels1, labels_new, axis=0)
    LOGGER.info(f'labels1 {labels1.shape} labels_new {labels_new.shape}')
    valid = 0
    for id in range(108):
        fig, axs = plt.subplots(1, 1, figsize=(16, 8), constrained_layout=True)
        # axs.grid(axis='x', color='0.95')
        # axs = plot_labeled_scatters3(x, labels1, colors1)
        axs = plot_labeled_scatters(x, labels1, colors1)
        print(np.where(x_new[:,1]< -20))
        print(len(np.where(x_new[:,1]< -20)[0]))
        for sid in range(6*id,6*(id+1)):
            a = [i for i,item in enumerate(timeframe_ids) if uniq_structure_ids[sid] in str(item)]
            structure_id = timeframe_ids[a][0].split("_")[-2].split(".")[-2]
            print(structure_id)
            # traj = axs.scatter3D(x_new[a, 0], x_new[a, 1], x_new[a, 2], c=np.arange(len(a)), alpha=1,
            traj = axs.scatter(x_new[a, 0], x_new[a, 1], c=np.arange(len(a)), alpha=1,
            cmap=cmaps[sid-(6*id)], zorder=2.5)
            axins = inset_axes(axs,
                        width="3%",
                        height="50%",
                        loc='right',
                        borderpad=0,
                        bbox_to_anchor=(1.05*((sid - (6*id)+1)/10), 0., 1, 1),
                        bbox_transform=axs.transAxes,
                    )
            cba = fig.colorbar(traj, cax=axins)
            cba.set_label('{:03d}={}'.format(int(structure_id), len(a)))
            axs.set_zorder(1)
        axs.set_xlabel('PC 1')
        axs.set_ylabel('PC 2')
        # axs.set_zlabel('PC 3')
        if len(title) > 0:
            fig.suptitle(title)

        res = filename.find(".png")
        filename_new = filename[:res-1] + str(id) + filename[res :]
        plt.savefig(filename_new, bbox_inches='tight')

        # axs.view_init(elev=0, azim=-90)
        # plt.savefig(filename_xz, bbox_inches='tight')
        # axs.view_init(elev=90, azim=-90)
        # plt.savefig(filename_xy, bbox_inches='tight')
        # axs.view_init(elev=0, azim=0)
        # plt.savefig(filename_yz, bbox_inches='tight')

        if show:
            plt.show()
        plt.close()
# ------------------------------------------------------------------------------
# PLOT TSNE
# ------------------------------------------------------------------------------
def plot_tsne(encodings, labels1, filename, title, show):
    colors1 = ['k',
               '#a6cee3', 'orange', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
               '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928',
               '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462',
               '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f',
               'orange', 'r']
    colors2 = ['k', '#a6cee3', '#b2df8a', '#fb9a99']
    colors = ['k', '#a6cee3', '#b2df8a']

    from sklearn.manifold import TSNE

    tsne = TSNE(n_components=2, random_state=0)
    x = tsne.fit_transform(encodings)

    fig, axs = plt.subplots(1, 1, figsize=(16, 8))

    plot_labeled_scatters(x, labels1, colors1)

    plt.xlabel('tsne dim 1')
    plt.ylabel('tsne dim 2')
    if len(title) > 0:
        fig.suptitle(title)

    plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()
# --------------------------------------------------------------------------
def plot_crd_dist(data, title, filename, nbins=30, ref=None, show=True):

    assert isinstance(data, np.ndarray)

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(9, 5))
    arr = plt.hist(data, density=False, bins=nbins)
    for i in range(nbins):
        plt.text(arr[1][i],arr[0][i],"%.f" % (arr[0][i]), fontsize=8)

    if ref !=None:
        plt.axvline(x = ref, color = 'r')
    axs.set_ylabel('Count')
    axs.set_xlabel('CRD distance (Ångström) ({} structures)'.format(data.shape[0]))
    axs.set_title(title)
    plt.xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
    plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()

# --------------------------------------------------------------------------
def plot_positions(positions_real, positions_gen, filename, title, show):
    fig = plt.figure(figsize=plt.figaspect(0.5))

    #---- First subplot
    ax1 = fig.add_subplot(1, 2, 2, projection='3d')
    l = 100
    ax1.scatter3D(positions_real[:, 0], positions_real[:, 1], positions_real[:, 2], color="navy")
    ax1.set_xlabel('X axis')
    ax1.set_ylabel('Y axis')
    ax1.set_zlabel('Z axis')
    ax1.title.set_text('Real')
    

    #---- Second subplot
    ax2 = fig.add_subplot(1, 2, 1, projection='3d')
    ax2.scatter3D(positions_gen[:, 0], positions_gen[:, 1], positions_gen[:, 2], color="navy")
    ax2.set_xlabel('X axis')
    ax2.set_ylabel('Y axis')
    ax2.set_zlabel('Z axis')
    ax2.title.set_text('Generated')

    max_x = np.amax([np.amax(positions_real[:, 0]), np.amax(positions_gen[:, 0])])
    min_x = np.amin([np.amin(positions_real[:, 0]), np.amin(positions_gen[:, 0])])

    max_y = np.amax([np.amax(positions_real[:, 1]), np.amax(positions_gen[:, 1])])
    min_y = np.amin([np.amin(positions_real[:, 1]), np.amin(positions_gen[:, 1])])

    max_z = np.amax([np.amax(positions_real[:, 2]), np.amax(positions_gen[:, 2])])
    min_z = np.amin([np.amin(positions_real[:, 2]), np.amin(positions_gen[:, 2])])

    ax1.set_xlim([min_x, max_x])
    ax1.set_ylim([min_y, max_y])
    ax1.set_zlim([min_z, max_z])

    ax2.set_xlim([min_x, max_x])
    ax2.set_ylim([min_y, max_y])
    ax2.set_zlim([min_z, max_z])

    fig.suptitle(title)
    plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()

# --------------------------------------------------------------------------
def prepare_hierarchical_distances(clusters_idxs, all_clusters_idxs, mean_mt, std_mt, data, dist_matrices_rep):

    data_orig = data
    # create a filter for removing repeated distances
    cluster_range =np.arange(0,120+1,3)
    filter_matrix = np.ones((120,120), dtype=np.float32)
    for cli in range(len(cluster_range) - 1):
        comb=list(combinations(np.arange(cluster_range[cli],cluster_range[cli+1]),2))
        for ci in range(len(comb)):
            # set negative the value if we have to remove this element
            filter_matrix[comb[ci][0]][comb[ci][1]] = -1
    filter_matrix = filter_matrix[np.triu_indices(120,1)]
    filter_ids = np.where(filter_matrix==1)


    # 120 representative beads
    mean_mt_rep = mean_mt[clusters_idxs, :][:, clusters_idxs]
    std_mt_rep = std_mt[clusters_idxs, :][:, clusters_idxs]
    mean_ut_rep = mean_mt_rep[np.triu_indices(120,1)]
    std_ut_rep = std_mt_rep[np.triu_indices(120,1)]
    data_rep, dmean_rep, dstd_rep = standardize(dist_matrices_rep, 'distance matrices',
                                            mean_ut_rep, std_ut_rep)

    # 40 UCG clusters
    for clid in range(len(all_clusters_idxs)):
        beads_ids = all_clusters_idxs[clid]
        n_beads_cluster = len(beads_ids)
        mean_mt_ucg = mean_mt[beads_ids, :][:, beads_ids]
        std_mt_ucg = std_mt[beads_ids, :][:, beads_ids]
        mean_ut_ucg = mean_mt_ucg[np.triu_indices(n_beads_cluster,1)]
        std_ut_ucg = std_mt_ucg[np.triu_indices(n_beads_cluster,1)]

        dist_matrices_ucg = data_orig[:, beads_ids, :][:, :, beads_ids]
        dist_matrices_ucg_new = np.empty((dist_matrices_ucg.shape[0],) + mean_ut_ucg.shape, dtype=np.float32)
        for ig in range(dist_matrices_ucg.shape[0]):
            dist_matrices_ucg_group = dist_matrices_ucg[ig,:,:]
            dist_matrices_ucg_new[ig] = dist_matrices_ucg_group[np.triu_indices(n_beads_cluster,1)] 
        data_ucg, dmean_ucg, dstd_ucg = standardize(dist_matrices_ucg_new, 'distance matrices',
                                            mean_ut_ucg, std_ut_ucg)
        if clid == 0:
            data = data_ucg
        else:
            data = np.concatenate((data, data_ucg), axis=1)

    # Add distance matrix of representative beads
    # find and remove distances which already exist in clusters distance matrices using a filter matrix
    data_rep_new = np.empty((data_rep.shape[0], int(data_rep.shape[1] - 120)), dtype=np.float32)
    for ri in range(data_rep.shape[0]):
        data_rep_new[ri] = data_rep[ri][filter_ids]
    data = np.concatenate((data, data_rep_new), axis=1)
    return data

# ------------------------------------------------------------------------------
def produce_hierarchical_distances(INPATH, positions):
    clusters_path = os.path.join(INPATH, 'rep_indices_per_cluster.npz')
    all_clusters_path = os.path.join(INPATH, 'all_indices_per_cluster.npz')
    std_fname = os.path.join(INPATH, 'CG_data_summary_pos_dis.npz')
    # std_fname = os.path.join(INPATH, 'CG-data_summary_rasraf_new.npz')

    clusters_file = np.load(clusters_path, allow_pickle=True)
    all_clusters_file = np.load(all_clusters_path, allow_pickle=True)
    clusters_idxs = clusters_file["indices_per_cluster"]
    all_clusters_idxs = all_clusters_file["indices_per_cluster"]
    clusters_idxs = clusters_idxs.flatten()

    outshp = np.shape(positions)[1:]
    n_structures = np.shape(positions)[0]
    dist_ms = np.empty((n_structures,outshp[0],outshp[0]), dtype=np.float32)
    dist_ms_rep = np.empty((n_structures,7140), dtype=np.float32)
    for tidx in range(n_structures):
        patch_pos = positions[tidx]
        # LOGGER.info(f' patch_pos {patch_pos.shape}')
        patch_pos_rep = patch_pos[clusters_idxs]
        distance_m_rep = distance.cdist(patch_pos_rep,patch_pos_rep)
        distance_m = distance.cdist(patch_pos,patch_pos)

        dist_m_u_rep = distance_m_rep[np.triu_indices(120,1)]
        dist_ms_rep[tidx] = dist_m_u_rep
        dist_ms[tidx] = distance_m
        # print(tidx)

    data = dist_ms
    
    with np.load(std_fname) as npz:
        mean_mt = npz['mean_dm']
        std_mt = npz['std_dm']
        data = prepare_hierarchical_distances(clusters_idxs, all_clusters_idxs, mean_mt, std_mt, data, dist_ms_rep)
    # adjust shapes and data types
    _do_cast = lambda d,t: d.astype(t) if d.dtype != t else d
    data = _do_cast(data, np.float32)

    return data
# ------------------------------------------------------------------------------
def standardize(data, label, dmean=None, dstd=None):

    LOGGER.info(f'Standardizing {data.shape} {label}')
    data -= dmean
    data /= dstd
    return data, dmean, dstd

# --------------------------------------------------------------------------      
def destandardize(data, label, dmean=None, dstd=None):

    LOGGER.info(f'De-standardizing {data.shape} {label}')
    data *= dstd
    data += dmean
    return data

# --------------------------------------------------------------------------   
def standardize_zero_std(data, label, ind, dmean=None, dstd=None):

    LOGGER.info(f'Standardizing {data.shape} {label}')
    data -= dmean

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if j == ind:
                data[i,j,0:2] = data[i,j,0:2] / dstd[j,0:2]
            else:
                data[i,j,:] = data[i,j,:] / dstd[j,:]
    return data, dmean, dstd
# --------------------------------------------------------------------------  
def rmsd(x, y):
    assert x.shape == y.shape
    d = np.square(x - y).sum(axis=-1)   # compute the squared distance
    return np.sqrt(d.mean(axis=-1))     # compute the mean over all points
# ------------------------------------------------------------------------------
# Class to analyze the protein latent space and select new coordinates in ls
# ------------------------------------------------------------------------------
class CGSampler():

    def __init__(self):
        self.random = True

    def analyze_ls(self, mdu):
        LOGGER.info(f'Analyze latent space')
        
    def draw_samples(self, nsamples, model_trainpath):
        # load validation data
        _keys = ['predicted_positions', 'x_test_encoded', 'labels', 'patch_ids']
        validation_data = load_validation_data(model_trainpath, keys_to_load=_keys)
        
        encodings = validation_data['x_test_encoded']
        predicted_positions = validation_data['predicted_positions']
        labels = validation_data['labels']

        n = nsamples/2
        class_a_ids = np.where(labels==0)[0]
        class_b_ids = np.where(labels==1)[0]

        # find the mean of each classs
        mean_encod_a = np.mean(encodings[class_a_ids], axis=0)
        mean_encod_b = np.mean(encodings[class_b_ids], axis=0)

        # Find the distance between the mean of the two classes
        dist = mean_encod_a - mean_encod_b
        LOGGER.info(f'class_a_ids {class_a_ids.shape} class_b_ids: {class_b_ids.shape} ')

        # For random encodings in each class produce new encodings
        rdm_ids_a = random.sample(range(0, class_a_ids.shape[0]), int(n))
        rdm_ids_b = random.sample(range(0, class_b_ids.shape[0]), int(n))
        encodings_a = encodings[class_a_ids[rdm_ids_a]] - (0.5 * dist)
        encodings_b = encodings[class_b_ids[rdm_ids_b]] + (0.5 * dist)

        new_encodings = np.concatenate((encodings_a, encodings_b), axis=0)

        LOGGER.info(f'Draw {nsamples} samples: {new_encodings.shape} ')

        new_labels = np.full((new_encodings.shape[0]), 2, dtype=int)

        return new_encodings, encodings, new_labels, labels

    # --------------------------------------------------------------------------
    def calculate_crd_distance_of_already_produced(self, model_trainpath, analysis_path, datamean, datastd, gen_path):
        # load generated data
        files = glob.glob(os.path.join(gen_path, 'structure*.gro'))
        n = len(files)

        gen_pos = np.empty((n,751,3), dtype=np.float32)
        for fidx in range(n):
            f = files[fidx]
            f = f[:-3] + "npz"
            with np.load(f, allow_pickle=True) as npz:
                patch_pos = npz['data']
                gen_pos[fidx] = patch_pos
        # gen_pos = np.concatenate(gen_pos,axis=0)


        # load validation data
        _keys = ['positions', 'labels', 'patch_ids']
        validation_data = load_validation_data(model_trainpath, keys_to_load=_keys)

        positions = validation_data['positions']
        labels = validation_data['labels']

        RESIDXS_CRD = np.array([628,630,634,636,659,663,665,667])

        # real_positions =p ositions
        class_a_ids = np.where(labels==0)[0]
        class_b_ids = np.where(labels==1)[0]
        do_destandardization = True
        if do_destandardization:
            # Apply destandardization
            real_positions = destandardize(positions, 'real_val_positions',
                                    datamean, datastd)

        pos_crd = real_positions[:,RESIDXS_CRD,2]
        LOGGER.info(f'pos_crd {pos_crd.shape}')
        mean_crd = np.mean(pos_crd, axis=1)
        CYF_bead_sc = real_positions[:,413,2]
        LOGGER.info(f'CYF_bead_sc {CYF_bead_sc.shape} mean_crd: {mean_crd.shape} ')
        # crd_dist = np.linalg.norm(mean_crd - CYF_bead_sc, axis=1)
        crd_dist = np.absolute(mean_crd - CYF_bead_sc)

        LOGGER.info(f'crd_dist {crd_dist.shape} mean_crd: {mean_crd.shape} ')


        name_crd_distance = os.path.join(os.path.join(analysis_path, 'distribution_crd_dis.png'))
        plot_crd_dist(crd_dist, "Distibution of CRD distance in real structures",
                name_crd_distance, nbins=30, show=False)


        LOGGER.info(f'gen_pos {gen_pos.shape} ')
        pos_crd_g = gen_pos[:,RESIDXS_CRD,2]
        LOGGER.info(f'pos_crd {pos_crd_g.shape}')
        mean_crd_g = np.mean(pos_crd_g, axis=1)
        CYF_bead_sc_g = gen_pos[:,413,2]
        LOGGER.info(f'CYF_bead_sc {CYF_bead_sc_g.shape} mean_crd: {mean_crd_g.shape} ')
        # crd_dist_g = np.linalg.norm(mean_crd_g - CYF_bead_sc_g, axis=1)
        crd_dist_g = np.absolute(mean_crd_g - CYF_bead_sc_g)
        LOGGER.info(f'crd_dist {crd_dist_g.shape} mean_crd: {mean_crd_g.shape} ')


        name_crd_distance_g = os.path.join(os.path.join(analysis_path, 'distribution_crd_dis_generated.png'))
        plot_crd_dist(crd_dist_g, "Distibution of CRD distance in generated structures",
                name_crd_distance_g, nbins=30, show=False)

# --------------------------------------------------------------------------
    def calculate_crd_distance(self, model_trainpath, analysis_path, gen_pos, datamean, datastd, ):

        # load real validation data
        _keys = ['positions', 'labels', 'patch_ids']
        validation_data = load_validation_data(model_trainpath, keys_to_load=_keys)

        positions = validation_data['positions']
        labels = validation_data['labels']

        RESIDXS_CRD = np.array([628,630,634,636,659,663,665,667])

        # real_positions =p ositions
        class_a_ids = np.where(labels==0)[0]
        class_b_ids = np.where(labels==1)[0]
        do_destandardization = True
        if do_destandardization:
            # Apply destandardization
            real_positions = destandardize(positions, 'real_val_positions',
                                    datamean, datastd)

        pos_crd = real_positions[:,RESIDXS_CRD,2]
        LOGGER.info(f'pos_crd {pos_crd.shape}')
        mean_crd = np.mean(pos_crd, axis=1)
        CYF_bead_sc = real_positions[:,413,2]
        LOGGER.info(f'CYF_bead_sc {CYF_bead_sc.shape} mean_crd: {mean_crd.shape} ')
        # crd_dist = np.linalg.norm(mean_crd - CYF_bead_sc, axis=1)
        crd_dist = np.absolute(mean_crd - CYF_bead_sc)

        LOGGER.info(f'crd_dist {crd_dist.shape} mean_crd: {mean_crd.shape} ')


        name_crd_distance = os.path.join(os.path.join(analysis_path, 'distribution_crd_dis.png'))
        plot_crd_dist(crd_dist, "Distibution of CRD distance in real structures",
                name_crd_distance, nbins=30, show=False)

        # generated structures
        LOGGER.info(f'gen_pos {gen_pos.shape} ')
        pos_crd_g = gen_pos[:,RESIDXS_CRD,2]
        LOGGER.info(f'pos_crd {pos_crd_g.shape}')
        mean_crd_g = np.mean(pos_crd_g, axis=1)
        CYF_bead_sc_g = gen_pos[:,413,2]
        LOGGER.info(f'CYF_bead_sc {CYF_bead_sc_g.shape} mean_crd: {mean_crd_g.shape} ')
        # crd_dist_g = np.linalg.norm(mean_crd_g - CYF_bead_sc_g, axis=1)
        crd_dist_g = np.absolute(mean_crd_g - CYF_bead_sc_g)
        LOGGER.info(f'crd_dist {crd_dist_g.shape} mean_crd: {mean_crd_g.shape} ')


        name_crd_distance_g = os.path.join(os.path.join(analysis_path, 'distribution_crd_dis_generated.png'))
        plot_crd_dist(crd_dist_g, "Distibution of CRD distance in generated structures",
                name_crd_distance_g, nbins=30, show=False)

        # Return the ids of middle patches
        middle_patch_ids = np.where((crd_dist_g>17) & (crd_dist_g<25))[0]
        return middle_patch_ids

    # --------------------------------------------------------------------------
    def find_closest_real_patches_to_generated(self, model, gen_path, MODELPATH_LS, OUTPATH):
        
        datamean = model.datamean
        datastd = model.datastd 
        REAL_MIDDLE_PATH = os.path.join(OUTPATH, 'real_middle_val')

        bb_ids_filename = os.path.join(MODELPATH_LS, 'bb_ids.npz')
        bb_ids_file = np.load(bb_ids_filename, allow_pickle=True)
        bb_ids = bb_ids_file["bb_ids"]

        # load generated data
        files = glob.glob(os.path.join(gen_path, 'structure*.gro'))
        RESIDXS_CRD = np.array([628,630,634,636,659,663,665,667])
        n = len(files)

        gen_pos = np.empty((n,751,3), dtype=np.float32)
        for fidx in range(n):
            f = files[fidx]
            f = f[:-3] + "npz"
            with np.load(f, allow_pickle=True) as npz:
                patch_pos = npz['data']
                gen_pos[fidx] = patch_pos

        LOGGER.info(f'gen_pos {gen_pos.shape} ')
        pos_crd_g = gen_pos[:,RESIDXS_CRD,2]
        # LOGGER.info(f'pos_crd {pos_crd_g.shape}')
        mean_crd_g = np.mean(pos_crd_g, axis=1)
        CYF_bead_sc_g = gen_pos[:,413,2]
        LOGGER.info(f'CYF_bead_sc {CYF_bead_sc_g.shape} mean_crd: {mean_crd_g.shape} ')
        crd_dist_g = np.absolute(mean_crd_g - CYF_bead_sc_g)
        LOGGER.info(f'crd_dist {crd_dist_g.shape} mean_crd: {mean_crd_g.shape} ')

        middle_patch_ids = np.where((crd_dist_g>17) & (crd_dist_g<25))[0]
        LOGGER.info(f'crd_dist {crd_dist_g[middle_patch_ids].shape} {crd_dist_g[middle_patch_ids]} ')
        middle_gen_pos = gen_pos[middle_patch_ids]
        files = np.array(files)
        middle_files_gen = files[middle_patch_ids]
        crd_distance_middle_gen = crd_dist_g[middle_patch_ids] 
        # LOGGER.info(f'files {files[middle_patch_ids]}')
        LOGGER.info(f'middle_gen_pos {middle_gen_pos.shape}')

        # Prepare the imput data for the autoencoder
        distances_gen = produce_hierarchical_distances(MODELPATH_LS, middle_gen_pos)
        # LOGGER.info(f'distances_gen {distances_gen.shape}') 
        do_standardization = True
        if do_standardization:
            # Apply destandardization
            # index 413 has z =0
            middle_gen_pos_g, dmean, dstd = standardize_zero_std(middle_gen_pos, 'gen_middle_positions',413,
                                    datamean, datastd) 
        n_patches = middle_gen_pos_g.shape[0]
        n_beads = middle_gen_pos_g.shape[1]
        flatten_data_gen = np.empty((n_patches, n_beads*3), dtype=np.float32)
        for pidx in range(n_patches):
            patch_pos_cg = middle_gen_pos_g[pidx]
            flatten_data_gen[pidx] = patch_pos_cg.flatten()

        data_gen = np.concatenate((flatten_data_gen, distances_gen), axis=1)

        encodings_gen = model.predict_ae(data_gen)
        LOGGER.info(f'encodings_gen {encodings_gen.shape}') 

 
        # load real validation data in the middle
        _keys = ['positions', 'labels', 'patchids']
        validation_data = load_validation_data(REAL_MIDDLE_PATH, keys_to_load=_keys) 
        positions = validation_data['positions']
        labels = validation_data['labels']
        patch_ids = validation_data['patchids'] 

        
        positions = np.reshape(positions, (positions.shape[0], 751, 3))
        LOGGER.info(f'Real positions {positions.shape}')

        # Prepare the imput data for the autoencoder 
        do_destandardization = True
        if do_destandardization:
            # Apply destandardization
            positions_raw = destandardize(positions, 'real_positions',
                                    datamean, datastd) 
        distances = produce_hierarchical_distances(MODELPATH_LS, positions_raw)
        do_standardization = True
        if do_standardization:
            # Apply destandardization
            # index 413 has z =0
            positions_d, dmean, dstd = standardize_zero_std(positions_raw, 'real_positions',413,
                                    datamean, datastd) 

        n_patches = positions_d.shape[0]
        n_beads = positions_d.shape[1]
        flatten_data = np.empty((n_patches, n_beads*3), dtype=np.float32)
        for pidx in range(n_patches):
            patch_pos_cg = positions_d[pidx]
            flatten_data[pidx] = patch_pos_cg.flatten()
        data = np.concatenate((flatten_data, distances), axis=1)
        encodings_real = model.predict_ae(data)
        LOGGER.info(f'encodings_real {encodings_real.shape}')


        CLOSES_PATCHES_TPATH = os.path.join(gen_path, 'closest_real_patches')
        if not os.path.isdir(CLOSES_PATCHES_TPATH):
            os.makedirs(CLOSES_PATCHES_TPATH, exist_ok=True)
        # Find closest points of real patches in between ensembles
        n_patches = encodings_gen.shape[0]
        closest_real_patches = np.empty((n_patches, 128), dtype=np.float32)
        closest_positions = np.empty((n_patches, 751, 3), dtype=np.float32)
        closest_real_patches_ids = np.empty((n_patches), dtype='object')
        for pidx in range(n_patches):
            encoding_g = encodings_gen[pidx]
            encodings_real = np.asarray(encodings_real)

            dist_2 = np.sum((encodings_real - encoding_g)**2, axis=1)
            closest_index = np.argmin(dist_2)

            closest_real_patches[pidx] = encodings_real[closest_index]
            closest_real_patches_ids[pidx] = patch_ids[closest_index]
            closest_positions[pidx] = positions_raw[closest_index]
        LOGGER.info(f'closest_real_patches {closest_real_patches.shape} encodings_gen {encodings_gen.shape}')

        LOGGER.info(f'middle_gen_pos= {middle_gen_pos.shape} closest_positions = {closest_positions.shape}')
     
        for pidx in range(n_patches):
            gen_struc = middle_gen_pos[pidx]
            real_struc = closest_positions[pidx]
            gromacs_ref_file = os.path.join(MODELPATH_LS, 'RAS-RAF_state_ma_0.gro')
            u_real = mda.Universe(gromacs_ref_file)
            u_real.atoms.positions = real_struc

            u_gen = mda.Universe(gromacs_ref_file)
            u_gen.atoms.positions = gen_struc
            
            align.alignto(u_gen, u_real)

            gen_struc = u_gen.atoms.positions
            real_struc = u_real.atoms.positions
            rmsd_v = rmsd(gen_struc, real_struc)
            rmsd_bb = rmsd(gen_struc[bb_ids], real_struc[bb_ids])
            filename = os.path.join(os.path.join(CLOSES_PATCHES_TPATH, f'positions_{pidx}.png'))
            title = f"{pidx}: Generated {middle_files_gen[pidx]} - Real {closest_real_patches_ids[pidx]} (RMSD: {rmsd_v:.2f}; RMSD_BB: {rmsd_bb:.2f}) (CRD dis {crd_distance_middle_gen[pidx]:.2f})"
            plot_positions(real_struc, gen_struc, filename, title, show=False)
        
        encodings = np.append(encodings_real, encodings_gen, axis=0)
        new_labels = np.full((encodings_gen.shape[0]), 3, dtype=int)
        labels = np.append(labels, new_labels)
        name_pca_validation = os.path.join(os.path.join(CLOSES_PATCHES_TPATH, 'pca_space_middle.png'))
        name_pca_validation_xz = os.path.join(os.path.join(CLOSES_PATCHES_TPATH, 'pca_space_middle_xz.png'))
        name_pca_validation_xy = os.path.join(os.path.join(CLOSES_PATCHES_TPATH, 'pca_space_middle_xy.png'))
        name_pca_validation_yz = os.path.join(os.path.join(CLOSES_PATCHES_TPATH, 'pca_space_middle_yz.png'))
        plot_pca(encodings, labels, name_pca_validation, name_pca_validation_xz, name_pca_validation_xy, name_pca_validation_yz, [], False)
        
        name_tsne_validation = os.path.join(os.path.join(CLOSES_PATCHES_TPATH, 'tsne_space_middle.png'))
        plot_tsne(encodings, labels, name_tsne_validation, [], False)

        filename_files = os.path.join(os.path.join(CLOSES_PATCHES_TPATH, f'files_generated_real.npz'))
        np.savez_compressed(filename_files, generated_files=middle_files_gen, real_files=closest_real_patches_ids )

