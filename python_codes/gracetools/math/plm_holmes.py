import numpy as np


def plm_holmes(x, lmax, astype=np.float64):
    """

    Args:
        x: scalar or array of value to determine the plm (between -1 and 1), usually cos of the latitude
        lmax: maximum degree to evaluate
        astype: output type variable

    Returns:
        plm: legender polynomials
        dplm: first derivateve on x

    References:
         S. A. Holmes and W. E. Featherstone, 2002, https://doi.org/10.1007/s00190-002-0216-2

    """
    x = np.atleast_1d(x).flatten().astype(astype)
    jm = len(x)
    lmax = np.int64(lmax)
    scalef = 1.0e-280

    #f1 is a
    #f2 is b
    f1 = np.zeros(((lmax + 1) * (lmax + 2) // 2), dtype=astype)
    f2 = np.zeros(((lmax + 1) * (lmax + 2) // 2), dtype=astype)
    p = np.zeros(((lmax + 1) * (lmax + 2) // 2, jm), dtype=astype)
    plm = np.zeros((lmax + 1, lmax + 1, jm), dtype=astype)
    dplm = np.zeros((lmax + 1, lmax + 1, jm), dtype=astype)

    k = 2 #eq 12
    for l in range(2, lmax + 1):
        k += 1
        f1[k] = np.sqrt(2.0 * l - 1.0) * np.sqrt(2.0 * l + 1.0) / np.float128(l)
        f2[k] = np.float128(l - 1.0) * np.sqrt(2.0 * l + 1.0) / (np.sqrt(2.0 * l - 3.0) * np.float128(l))
        for m in range(1, l - 1):
            k += 1
            f1[k] = np.sqrt(2.0 * l + 1.0) * np.sqrt(2.0 * l - 1.0) / (np.sqrt(l + m) * np.sqrt(l - m))
            f2[k] = np.sqrt(2.0 * l + 1.0) * np.sqrt(l - m - 1.0) * np.sqrt(l + m - 1.0) / \
                    (np.sqrt(2.0 * l - 3.0) * np.sqrt(l + m) * np.sqrt(l - m))
        k += 2

    u = np.sqrt(1.0 - x ** 2)
    u[u == 0] = np.finfo(u.dtype).eps

    p[0, :] = 1.0
    p[1, :] = np.sqrt(3.0) * x
    k = 1
    for l in range(2, lmax + 1):
        k += l
        p[k, :] = f1[k] * x * p[k - l, :] - f2[k] * p[k - 2 * l + 1, :]

    pmm = np.sqrt(2.0) * scalef
    rescalem = 1.0 / scalef
    kstart = 0

    for m in range(1, lmax):
        rescalem = rescalem * u
        kstart += m + 1
        pmm = pmm * np.sqrt(2 * m + 1) / np.sqrt(2 * m)
        p[kstart, :] = pmm
        k = kstart + m + 1
        p[k, :] = x * np.sqrt(2 * m + 3) * pmm
        for l in range(m + 2, lmax + 1):
            k += l
            p[k, :] = x * f1[k] * p[k - l, :] - f2[k] * p[k - 2 * l + 1, :]
            p[k - 2 * l + 1, :] = p[k - 2 * l + 1, :] * rescalem
        p[k, :] = p[k, :] * rescalem
        p[k - lmax, :] = p[k - lmax, :] * rescalem

    rescalem = rescalem * u
    kstart += m + 2
    p[kstart, :] = pmm * np.sqrt(2 * lmax + 1) / np.sqrt(2 * lmax) * rescalem
    for m in range(lmax + 1):
        for l in range(m, lmax + 1):
            lm = (l * (l + 1)) // 2 + m
            plm[l, m, :] = p[lm, :]
            if l == m:
                dplm[l, m, :] = np.float128(m) * (x / u) * plm[l, m, :]
            else:
                flm = np.sqrt(((l ** 2.0 - m ** 2.0) * (2.0 * l + 1.0)) / (2.0 * l - 1.0))
                dplm[l, m, :] = (1.0 / u) * (l * x * plm[l, m, :] - flm * plm[l - 1, m, :])

    return plm, dplm
