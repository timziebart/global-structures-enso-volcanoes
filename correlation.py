
import numpy as np
import functools as ft

# always flush print output
print = ft.partial(print, flush=True)

def corr_coeff(A,B = None):
    # makes sure they are arrays and have the same shape
    A = np.asarray(A)
    if B is None:
        B = np.asarray(A)
    else:
        B = np.asarray(B)
    assert A.shape == B.shape
    assert A.ndim <= 2

    # shape of the reduced arrays
    # flat_shape = list(A.shape)
    # flat_shape[axis] = 1

    # mean of input along axis 'axis' & subtract from input arrays themeselves
    A_mA = A - A.mean(0)
    B_mB = B - B.mean(0)

    # calculate covariance
    numerator = np.dot(A_mA.T,B_mB)
    # axis = 0
    # numerator = np.tensordot(A_mA,B_mB, axes=[(axis,), (axis,)])

    # Sum of squares across rows
    ssA = (A_mA**2).sum(0)
    ssB = (B_mB**2).sum(0)

    # calcualte sqrt of variances
    denominator = np.sqrt(np.outer(ssA, ssB))

    # Finally get corr coeff
    return numerator / denominator

def thresholding_matrix(C, percentage, save_histo=False):
    '''
    In this function we will set a threshold
    in the correlation matrix in order to keep
    only the strongest links. 
    '''

    C = np.nan_to_num(C)
    C[np.diag_indices_from(C)] = 0.0
    il = np.tril_indices(C.shape[0], k=-1)

    #         flat = C[il]
    #         hist, bins = np.histogram(flat,5000)
    # ##         print len(hist), hist
    # ##         print len(bins), bins
    #         h = cumulative_distribution(1.0*hist)
    #
    #         v = (h>=1.0-percentage).nonzero()
    #         thr = bins[v[0][0]]
    flat = np.sort(C[il])
    del il
    assert flat.ndim == 1, "something went wrong with flattening"
    thr_index = int((1 - percentage) * flat.shape[0])
    thr = flat[thr_index]
    del flat
    print("at", thr, end=" ... ")

    # create weighted adjacency matrix
    # mask = C>thr
    # Cbin = np.zeros_like(C)
    # Cbin[mask] = C[mask]

    # create unweighted adjacency matrix
    Cbin = (C > thr).astype(np.int8)


    return Cbin