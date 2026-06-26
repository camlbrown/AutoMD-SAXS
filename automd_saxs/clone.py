"""CLoNe clustering (Clustering of Local Neighborhoods).

Faithful port of the legacy ``slurms/post_processing/python_scripts/clone.py``
into the package. The numerical algorithm is unchanged -- it still uses numpy,
scipy, and scikit-learn -- so the scientific behaviour is preserved exactly. The
only difference is that those heavy dependencies are imported lazily inside
:meth:`CLoNe.fit`, so this module imports on a bare interpreter and raises a
clear :class:`~automd_saxs.command_runner.MissingDependencyError` only if you
actually try to cluster without them.

Algorithm (unchanged): kNN distance matrix -> Gaussian-kernel local densities ->
candidate cluster cores -> real cluster centres -> density-ordered assignation ->
border-point merging via Bhattacharyya distance -> optional outlier/noise removal.
"""

import operator

from .command_runner import MissingDependencyError


def _require_numerics():
    try:
        import numpy as np
        from sklearn.neighbors import NearestNeighbors
        from scipy.stats import gaussian_kde
    except ImportError as exc:  # pragma: no cover - exercised only without deps
        raise MissingDependencyError(
            "CLoNe clustering requires numpy, scipy, and scikit-learn. Install "
            "them (they are in automdsaxs.yml) or use the parsing/summary helpers "
            "in automd_saxs.structural instead. Original import error: {0}".format(exc)
        )
    return np, NearestNeighbors, gaussian_kde


class CLoNe(object):
    """CLoNe clustering. See module docstring; ``pdc`` is the main parameter."""

    def __init__(self, pdc=4, n_resize=4, filt=0.1, verbose=False):
        self.pdc = pdc
        self.n_resize = n_resize
        self.filt = filt
        self.verbose = verbose

    def fit(self, vectors):
        np, NearestNeighbors, gaussian_kde = _require_numerics()

        # 1. Build neighbour distance matrix
        nb_el = len(vectors)
        effective_nb_el = int(nb_el / self.n_resize)

        if self.verbose:
            print("> Computing kNN...")

        vectors = vectors.astype(np.float32)
        nbrs = NearestNeighbors(n_neighbors=effective_nb_el, algorithm="auto",
                                n_jobs=1, metric="l2").fit(vectors)
        knn_dist, knn_idx = nbrs.kneighbors(vectors)
        knn_dist = knn_dist.astype(np.float32)
        knn_idx = knn_idx.astype(np.int32)

        # Get Gaussian sigma for local density based on neighbour matrix
        position = int(round(nb_el ** 2 * self.pdc * 0.01))
        if position >= np.prod(knn_dist.shape):
            position = -1
        dc = np.sort(knn_dist.flatten())[position]

        # 3. Compute local densities
        inv_dc = 1 / dc
        if self.verbose:
            print("> Computing local densities...")

        ordered_knn_dens = np.exp(-1 * np.square(knn_dist * inv_dc))
        summed_okd = ordered_knn_dens.cumsum(axis=1) - 1
        rho = summed_okd[:, -1].copy()
        summed_okd /= int(round(self.pdc * 0.01 * nb_el))
        core_card = np.unravel_index(
            np.argmax(summed_okd > ordered_knn_dens, axis=1), summed_okd.shape)[1]

        # 4. Get real cluster centers
        if self.verbose:
            print("> Identifying real centers...")

        nneigh = [0] * nb_el
        centers = []
        labels = [-1] * nb_el
        nb_clust = 0
        for x in range(0, nb_el):
            cur_idx_list = knn_idx[x]
            is_set = False
            for y in range(1, effective_nb_el):
                if rho[cur_idx_list[y]] > rho[x]:
                    is_set = True
                    nneigh[x] = cur_idx_list[y]
                    if y > core_card[x] and core_card[x]:
                        centers.append(x)
                        nb_clust += 1
                    break
            if not is_set:
                centers.append(x)
                nb_clust += 1

        if nb_clust == 1:
            print("Unique cluster found with set parameters.")
            self.labels_ = labels
            self.labels_all = labels
            self.centers = centers
            self.core_card = core_card
            self.rho = rho
            return

        to_sort = []
        for x in range(nb_clust):
            to_sort.append([rho[centers[x]], centers[x]])
        to_sort = sorted(to_sort, key=operator.itemgetter(0))
        centers = [el[1] for el in to_sort]
        for x in range(nb_clust):
            labels[centers[x]] = x
        idx_rho_sorted = np.argsort(-np.array(rho))

        # 4. Cluster assignation
        if self.verbose:
            print("> Assigning points to clusters...")

        labels = np.array(labels)
        for x in range(nb_el):
            if labels[idx_rho_sorted[x]] == -1:
                labels[idx_rho_sorted[x]] = labels[nneigh[idx_rho_sorted[x]]]

        # 5. Merging
        outliers = np.empty(0, dtype=int)
        core_points = np.empty(0, dtype=int)

        if self.verbose:
            print("> Merging centers...")

        core_dict = {}
        for cl in centers:
            lab = labels[cl]
            core_dict[lab] = []
            for j in knn_idx[cl, :core_card[cl]]:
                core_dict[lab].append(j)

        border_points = {}
        for c1_idx in range(1, len(centers) - 1):
            c1 = centers[c1_idx]
            lab1 = labels[c1]
            mask1 = np.where(labels == lab1)[0]
            for p1 in mask1:
                for neigh in range(effective_nb_el):
                    if knn_dist[p1, neigh] > dc:
                        break
                    p2 = knn_idx[p1, neigh]
                    lab2 = labels[p2]
                    if lab2 == -1:
                        continue
                    if lab1 != lab2:
                        if rho[c1] > rho[centers[lab2]]:
                            key = tuple((lab2, lab1))
                        else:
                            key = tuple((lab1, lab2))
                        if key not in border_points.keys():
                            border_points[key] = []
                        border_points[key].append(p1)
                        border_points[key].append(p2)

        # Merge based on Bhattacharyya distance > 0.65
        for k in sorted(border_points.keys(), reverse=True):
            k0, k1 = k
            if labels[centers[k0]] == labels[centers[k1]]:
                continue
            border_points[k] = np.unique(border_points[k])
            mask0 = np.where(labels == labels[centers[k0]])[0]
            np.where(labels == labels[centers[k1]])
            try:
                kde_bord = gaussian_kde(core_card[border_points[k]])
                kde_core0 = gaussian_kde(core_card[core_dict[k0]])
                kde_core1 = gaussian_kde(core_card[core_dict[k1]])
            except Exception:
                continue
            mb = min(min(core_card[border_points[k]]), min(core_card[core_dict[k0]]),
                     min(core_card[core_dict[k1]]))
            mc = max(max(core_card[border_points[k]]), max(core_card[core_dict[k0]]),
                     max(core_card[core_dict[k1]]))

            support = np.linspace(mb, mc, 100)
            c2_0 = kde_core0.evaluate(support)
            c2_1 = kde_core1.evaluate(support)
            b2 = kde_bord.evaluate(support)

            c2_1 /= np.sum(c2_1)
            c2_0 /= np.sum(c2_0)
            b2 /= np.sum(b2)
            bc_0 = np.sum(np.sqrt(np.multiply(b2, c2_0)))
            bc_1 = np.sum(np.sqrt(np.multiply(b2, c2_1)))
            bc_01 = 0.5 * (bc_0 + bc_1)

            if bc_01 > 0.65:
                labels[mask0] = labels[centers[k1]]

        # Clean centers
        cur_lab = 0
        centers = []
        for ulab in np.unique(labels):
            if ulab == -1:
                continue
            mask_l = np.where(labels == ulab)[0]
            centers.append(mask_l[np.argmax(rho[mask_l])])
            labels[mask_l] = cur_lab
            cur_lab += 1

        for cl in centers:
            lab = labels[cl]
            for j in knn_idx[cl, :core_card[cl]]:
                core_points = np.append(core_points, j)

        # Removing outliers
        labels_all = labels.copy()
        outliers = np.array([], dtype=int)
        outlier_clusters = []
        if self.filt:
            if self.verbose:
                print("> Removing outliers...")

            for cl in centers:
                lab = labels[cl]
                mask = np.where(labels == lab)[0]
                cl_v = vectors[mask]
                core_dc = knn_dist[cl, core_card[cl]]
                nbrs = NearestNeighbors(n_neighbors=len(cl_v), algorithm="auto",
                                        n_jobs=1, metric="l2").fit(cl_v)
                cl_knn_dist, cl_knn_idx = nbrs.kneighbors(cl_v)
                cl_ordered_knn_dens = np.exp(-1 * np.square(cl_knn_dist / core_dc))
                cl_summed_okd = cl_ordered_knn_dens.cumsum(axis=1) - 1
                cl_rho = cl_summed_okd[:, -1]
                rho[mask] = cl_rho
                rho_center = rho[cl]
                rho_thresh = self.filt * rho_center
                noise_mask = np.where(cl_rho < rho_thresh)[0]
                outliers = np.append(outliers, mask[noise_mask])
                labels[mask[noise_mask]] = -1

            if len(outliers) > 2:
                for cl in centers:
                    lab = labels[cl]
                    for j in knn_idx[cl, :core_card[cl]]:
                        core_points = np.append(core_points, j)

                kde_noise = gaussian_kde(rho[outliers])
                kde_core = gaussian_kde(rho[core_points])

                mb = np.amin(rho)
                mc = np.amax(rho)
                support = np.linspace(mb, mc, 1000)

                log_pN = np.log(len(outliers) / len(vectors))
                log_pC = np.log(len(core_points) / len(vectors))

                for cl in centers:
                    lab = labels[cl]
                    mask_core = knn_idx[cl, :core_card[cl]]
                    rho_core = [rho[i] for i in mask_core if labels[i] == lab]
                    core_size = len(rho_core)
                    d_core = kde_core.evaluate(rho_core)
                    d_noise = kde_noise.evaluate(rho_core)

                    if len(d_noise[d_noise > 1e-16]) < core_size:
                        log_pN_pNX = np.NINF
                    else:
                        log_pN_pNX = np.sum(np.log(d_noise[d_noise > 1e-16])) + log_pN

                    if not len(d_core[d_core > 1e-16]):
                        log_pC_pCX = np.NINF
                    else:
                        log_pC_pCX = np.sum(np.log(d_core[d_core > 1e-16])) + log_pC

                    if log_pN_pNX > log_pC_pCX:
                        outlier_clusters.append(np.where(labels == lab)[0])
                        labels[labels == lab] = -1

                init_c_idx = 0
                updated_centers = []
                for c in centers:
                    lab = labels[c]
                    if lab > -1:
                        labels_all[labels_all == lab] = init_c_idx
                        labels[labels == lab] = init_c_idx
                        init_c_idx += 1
                        updated_centers.append(c)
                centers = updated_centers

        self.rho = rho
        self.core_card = core_card
        self.centers = centers
        self.labels_ = labels
        self.labels_all = labels_all
