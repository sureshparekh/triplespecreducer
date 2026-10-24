"""
The MPFIT algorithm (Markwardt 2009; More 1978, MINPACK-1 LMDIF), in Python.

Spextool fits most of its line and peak positions with MPFIT, and MPFIT
works in the precision of the residuals it is given. Where Spextool builds
its starting values with FLOAT() -- MC_FINDLINES1D, the profile peaks of
MC_FINDPEAKS -- the whole fit runs in single precision: float machine
epsilon, finite-difference steps of sqrt(eps)|p|, float sums. That fit
stops somewhere a double-precision solver does not, and on a weak line the
difference decides whether the line is kept. So the algorithm is followed
step by step, in the same dtype IDL would use.

Only the parts Spextool reaches are here: no parameter limits, ties or
external state; analytic derivatives when the model supplies them
(MPFITPEAK), forward differences otherwise (MPFITEXPR).
"""

import numpy as np


class _Mach:
    def __init__(self, dtype):
        fi = np.finfo(dtype)
        self.D = dtype
        self.machep = dtype(fi.eps)
        self.minnum = dtype(fi.tiny)
        self.maxnum = dtype(fi.max)
        self.rdwarf = dtype(np.sqrt(self.minnum * dtype(1.5)) * dtype(10))
        # FLOAT literals (0.1, 0.001, 0.0001) keep their single-precision
        # value in a double fit, as in IDL: L(v).
        self.rgiant = dtype(np.sqrt(self.maxnum) * dtype(np.float32(0.1)))


def _L(D, v):
    """An IDL FLOAT literal in precision D: its single-precision value."""
    return D(np.float32(v))


def _sq(x):
    """x^2 as IDL computes it, x*x (numpy's scalar ** 2 is not always that)."""
    return x * x


def _total(a, D):
    """IDL TOTAL below the thread-pool threshold: a running sum in D."""
    a = np.asarray(a, dtype=D).ravel()
    return D(np.cumsum(a, dtype=D)[-1]) if a.size else D(0)


def _enorm(vec, M):
    D = M.D
    vec = np.asarray(vec, dtype=D)
    agiant = D(M.rgiant / D(vec.size))
    adwarf = D(M.rdwarf * D(vec.size))
    mx = D(max(abs(vec.max()), abs(vec.min())))
    if mx == 0:
        return D(0)
    if mx > agiant or mx < adwarf:
        return D(mx * np.sqrt(_total(_sq((vec / mx)), D)))
    return D(np.sqrt(_total(_sq(vec), D)))


def _qrfac(a, M):
    """MPFIT_QRFAC with pivoting. `a` is [m, n] and is overwritten."""
    D = M.D
    m, n = a.shape
    acnorm = np.array([_enorm(a[:, j], M) for j in range(n)], dtype=D)
    rdiag = acnorm.copy()
    wa = rdiag.copy()
    ipvt = np.arange(n)
    for j in range(min(m, n)):
        rmax = rdiag[j:].max()
        kmax = int(np.where(rdiag[j:] == rmax)[0][0]) + j
        if kmax != j:
            ipvt[j], ipvt[kmax] = ipvt[kmax], ipvt[j]
            rdiag[kmax] = rdiag[j]
            wa[kmax] = wa[j]
        lj = ipvt[j]
        ajj = a[j:, lj].copy()
        ajnorm = _enorm(ajj, M)
        if ajnorm == 0:
            rdiag[j] = -ajnorm
            continue
        if a[j, lj] < 0:
            ajnorm = -ajnorm
        ajj = (ajj / ajnorm).astype(D)
        ajj[0] = ajj[0] + D(1)
        a[j:, lj] = ajj
        for k in range(j + 1, n):
            lk = ipvt[k]
            ajk = a[j:, lk].copy()
            if a[j, lj] != 0:
                a[j:, lk] = (ajk - (ajj * _total(ajk * ajj, D)).astype(D) / a[j, lj]).astype(D)
            if rdiag[k] != 0:
                temp = D(a[j, lk] / rdiag[k])
                rdiag[k] = D(rdiag[k] * np.sqrt(max(D(1.) - temp * temp, D(0))))
                temp = D(rdiag[k] / wa[k])
                if 0.05 * float(temp) * float(temp) <= float(M.machep):
                    rdiag[k] = _enorm(a[j + 1:, lk], M)
                    wa[k] = rdiag[k]
        rdiag[j] = -ajnorm
    return ipvt, rdiag, acnorm


def _qrsolv(r, ipvt, diag, qtb, sdiag, M):
    D = M.D
    n = r.shape[1]
    for j in range(n):
        r[j:n, j] = r[j, j:n]
    x = np.array([r[i, i] for i in range(n)], dtype=D)
    wa = qtb.astype(D).copy()
    half, quart = D(0.5), D(0.25)
    for j in range(n):
        l = ipvt[j]
        if diag[l] != 0:
            sdiag[j:] = 0
            sdiag[j] = diag[l]
            qtbpj = D(0)
            for k in range(j, n):
                if sdiag[k] == 0:
                    continue
                if abs(r[k, k]) < abs(sdiag[k]):
                    cotan = D(r[k, k] / sdiag[k])
                    sine = D(half / np.sqrt(quart + quart * cotan * cotan))
                    cosine = D(sine * cotan)
                else:
                    tang = D(sdiag[k] / r[k, k])
                    cosine = D(half / np.sqrt(quart + quart * tang * tang))
                    sine = D(cosine * tang)
                r[k, k] = D(cosine * r[k, k] + sine * sdiag[k])
                temp = D(cosine * wa[k] + sine * qtbpj)
                qtbpj = D(-sine * wa[k] + cosine * qtbpj)
                wa[k] = temp
                if n > k + 1:
                    temp = (cosine * r[k + 1:n, k] + sine * sdiag[k + 1:n]).astype(D)
                    sdiag[k + 1:n] = (-sine * r[k + 1:n, k] + cosine * sdiag[k + 1:n]).astype(D)
                    r[k + 1:n, k] = temp
        sdiag[j] = r[j, j]
        r[j, j] = x[j]
    nsing = n
    wh = np.where(sdiag == 0)[0]
    if wh.size:
        nsing = wh[0]
        wa[nsing:] = 0
    if nsing >= 1:
        wa[nsing - 1] = D(wa[nsing - 1] / sdiag[nsing - 1])
        for j in range(nsing - 2, -1, -1):
            s = _total(r[j + 1:nsing, j] * wa[j + 1:nsing], D)
            wa[j] = D((wa[j] - s) / sdiag[j])
    x = np.zeros(n, dtype=D)
    x[ipvt] = wa
    return x


def _lmpar(r, ipvt, diag, qtb, delta, sdiag, par, M):
    D = M.D
    n = r.shape[1]
    dg = np.array([r[i, i] for i in range(n)], dtype=D)
    nsing = n
    wa1 = qtb.astype(D).copy()
    rthresh = D(np.abs(dg).max() * M.machep)
    wh = np.where(np.abs(dg) < rthresh)[0]
    if wh.size:
        nsing = wh[0]
        wa1[wh[0]:] = 0
    if nsing >= 1:
        for j in range(nsing - 1, -1, -1):
            wa1[j] = D(wa1[j] / r[j, j])
            if j - 1 >= 0:
                wa1[0:j] = (wa1[0:j] - r[0:j, j] * wa1[j]).astype(D)
    x = np.zeros(n, dtype=D)
    x[ipvt] = wa1
    it = 0
    wa2 = (diag * x).astype(D)
    dxnorm = _enorm(wa2, M)
    fp = D(dxnorm - delta)
    if fp <= _L(D, 0.1) * delta:
        return D(0), x
    parl = D(0)
    if nsing >= n:
        wa1 = (diag[ipvt] * wa2[ipvt] / dxnorm).astype(D)
        wa1[0] = D(wa1[0] / r[0, 0])
        for j in range(1, n):
            s = _total(r[0:j, j] * wa1[0:j], D)
            wa1[j] = D((wa1[j] - s) / r[j, j])
        temp = _enorm(wa1, M)
        parl = D(((fp / delta) / temp) / temp)
    for j in range(n):
        s = _total(r[0:j + 1, j] * qtb[0:j + 1], D)
        wa1[j] = D(s / diag[ipvt[j]])
    gnorm = _enorm(wa1, M)
    paru = D(gnorm / delta)
    if paru == 0:
        paru = D(M.minnum / min(delta, _L(D, 0.1)))
    par = D(max(par, parl))
    par = D(min(par, paru))
    if par == 0:
        par = D(gnorm / dxnorm)
    while True:
        it += 1
        if par == 0:
            par = D(max(M.minnum, paru * _L(D, 0.001)))
        temp = D(np.sqrt(par))
        wa1 = (temp * diag).astype(D)
        x = _qrsolv(r, ipvt, wa1, qtb, sdiag, M)
        wa2 = (diag * x).astype(D)
        dxnorm = _enorm(wa2, M)
        temp = fp
        fp = D(dxnorm - delta)
        if (abs(float(fp)) <= 0.1 * float(delta)) or \
                ((parl == 0) and (fp <= temp) and (temp < 0)) or it == 10:
            return par, x
        wa1 = (diag[ipvt] * wa2[ipvt] / dxnorm).astype(D)
        for j in range(n - 1):
            wa1[j] = D(wa1[j] / sdiag[j])
            wa1[j + 1:n] = (wa1[j + 1:n] - r[j + 1:n, j] * wa1[j]).astype(D)
        wa1[n - 1] = D(wa1[n - 1] / sdiag[n - 1])
        temp = _enorm(wa1, M)
        parc = D(((fp / delta) / temp) / temp)
        if fp > 0:
            parl = D(max(parl, par))
        if fp < 0:
            paru = D(min(paru, par))
        par = D(max(parl, par + parc))


def mpfit(fcn, p0, dtype=np.float64, jac=None, ftol=1e-10, xtol=1e-10,
          gtol=1e-10, maxiter=200, factor=100., trace=None, dbg=None):
    """
    MPFIT. `fcn(p)` returns the residual vector; `jac(p)`, when given,
    returns d(residual)/dp as [m, n] (MPFIT negates the model derivative).
    Returns (p, status).
    """
    D = dtype
    M = _Mach(D)
    x = np.array(p0, dtype=D)
    n = x.size
    fvec = np.asarray(fcn(x), dtype=D)
    m = fvec.size
    fnorm = _enorm(fvec, M)
    par = D(0)
    it = 1
    qtf = np.zeros(n, dtype=D)
    info = 0
    eps = D(np.sqrt(max(M.machep, M.machep)))
    diag = None
    delta = D(0)
    xnorm = D(0)
    fnorm1 = D(-1)
    while True:
        if trace is not None:                    # what ITERPROC is shown
            trace.append((it, x.copy(), D(fnorm * fnorm)))
        # FDJAC2
        if jac is not None:
            fjac = np.asarray(jac(x), dtype=D).copy()
        else:
            fjac = np.zeros((m, n), dtype=D)
            h = (eps * np.abs(x)).astype(D)
            h[h == 0] = eps
            for j in range(n):
                xp = x.copy()
                xp[j] = D(xp[j] + h[j])
                fp = np.asarray(fcn(xp), dtype=D)
                fjac[:, j] = ((fp - fvec) / h[j]).astype(D)
        if dbg is not None and it == dbg.get('it'):
            dbg['fjac0'] = fjac.copy()
        ipvt, wa1, wa2 = _qrfac(fjac, M)
        if dbg is not None and it == dbg.get('it'):
            dbg.update(ipvt=ipvt.copy(), wa1=wa1.copy(), wa2=wa2.copy(), fjacqr=fjac.copy())
        if it == 1:
            diag = wa2.copy()
            diag[diag == 0] = D(1)
            wa3 = (diag * x).astype(D)
            xnorm = _enorm(wa3, M)
            delta = D(D(factor) * xnorm)
            if delta == 0:
                delta = D(factor)
        wa4 = fvec.copy()
        for j in range(n):
            lj = ipvt[j]
            temp3 = fjac[j, lj]
            if temp3 != 0:
                fj = fjac[j:, lj]
                wj = wa4[j:]
                wa4[j:] = (wj - (fj * _total(fj * wj, D)).astype(D) / temp3).astype(D)
            fjac[j, lj] = wa1[j]
            qtf[j] = wa4[j]
        fjac = fjac[:n, :n][:, ipvt].copy()
        if not np.all(np.isfinite(fjac)):
            return x, -16
        gnorm = D(0)
        if fnorm != 0:
            for j in range(n):
                l = ipvt[j]
                if wa2[l] != 0:
                    s = D(_total(fjac[0:j + 1, j] * qtf[0:j + 1], D) / fnorm)
                    gnorm = D(max(gnorm, abs(s / wa2[l])))
        if gnorm <= gtol:
            info = 4
        if info != 0:
            return x, info
        diag = np.maximum(diag, wa2).astype(D)
        if dbg is not None and it == dbg.get('it'):
            dbg.update(qtf=qtf.copy(), gnorm=gnorm, diag=diag.copy(), delta=delta)
        while True:
            sdiag = np.zeros(n, dtype=D)
            par, wa1 = _lmpar(fjac.copy(), ipvt, diag, qtf, delta, sdiag, par, M)
            if dbg is not None and it == dbg.get('it'):
                dbg.setdefault('lm', []).append((par, wa1.copy(), delta))
            wa1 = -wa1
            alpha = D(1)
            wa2 = (x + wa1).astype(D)
            wa3 = (diag * wa1).astype(D)
            pnorm = _enorm(wa3, M)
            if it == 1:
                delta = D(min(delta, pnorm))
            wa4 = np.asarray(fcn(wa2), dtype=D)
            fnorm1 = _enorm(wa4, M)
            actred = D(-1)
            if 0.1 * float(fnorm1) < float(fnorm):
                actred = D(-((fnorm1 / fnorm) * (fnorm1 / fnorm)) + D(1.))   # x^2 is x*x; a numpy scalar's x ** 2 is not
            wa3 = np.zeros(n, dtype=D)
            for j in range(n):
                wa3[j] = 0
                wa3[0:j + 1] = (wa3[0:j + 1] + fjac[0:j + 1, j] * wa1[ipvt[j]]).astype(D)
            temp1 = D(_enorm(alpha * wa3, M) / fnorm)
            temp2 = D((np.sqrt(alpha * par) * pnorm) / fnorm)
            half = D(0.5)
            prered = D(temp1 * temp1 + (temp2 * temp2) / half)
            dirder = D(-(temp1 * temp1 + temp2 * temp2))
            ratio = D(0)
            tenth = _L(D, 0.1)
            if prered != 0:
                ratio = D(actred / prered)
            if dbg is not None and it == dbg.get('it'):
                dbg.setdefault('steps', []).append(dict(fnorm=fnorm, pnorm=pnorm, fnorm1=fnorm1, actred=actred,
                                                        prered=prered, dirder=dirder, ratio=ratio,
                                                        temp1=temp1, temp2=temp2))
            if ratio <= 0.25:
                if actred >= 0:
                    temp = half
                else:
                    temp = D(half * dirder / (dirder + half * actred))
                if (0.1 * float(fnorm1) >= float(fnorm)) or (temp < 0.1):
                    temp = tenth
                delta = D(temp * min(delta, D(pnorm / tenth)))
                par = D(par / temp)
            else:
                if par == 0 or ratio >= D(0.75):
                    delta = D(pnorm / half)
                    par = D(half * par)
            if ratio >= _L(D, 0.0001):
                x = wa2
                wa2 = (diag * x).astype(D)
                fvec = wa4
                xnorm = _enorm(wa2, M)
                fnorm = fnorm1
                it += 1
            if abs(float(actred)) <= ftol and float(prered) <= ftol and 0.5 * float(ratio) <= 1:
                info = 1
            if float(delta) <= xtol * float(xnorm):
                info = 2
            if abs(float(actred)) <= ftol and float(prered) <= ftol and \
                    0.5 * float(ratio) <= 1 and info == 2:
                info = 3
            if info != 0:
                return x, info
            if it >= maxiter:
                info = 5
            if abs(actred) <= M.machep and prered <= M.machep and 0.5 * ratio <= 1:
                info = 6
            if delta <= M.machep * xnorm:
                info = 7
            if gnorm <= M.machep:
                info = 8
            if info != 0:
                return x, info
            if ratio < _L(D, 0.0001):
                continue
            if not (np.all(np.isfinite(wa1)) and np.all(np.isfinite(wa2))
                    and np.all(np.isfinite(x)) and np.isfinite(ratio)):
                return x, -16
            break


# ------------------------------------------------------------- callers

def _is_double(*arrays):
    return any(np.asarray(a).dtype.kind == "f" and np.asarray(a).dtype.itemsize == 8
               for a in arrays)


def gaussian_f(x, p, D):
    """astrolib GAUSSIAN in dtype D, with its underflow guard."""
    x = np.asarray(x, dtype=D)
    p = np.asarray(p, dtype=D)
    z = ((x - p[1]) / p[2]).astype(D)
    zz = (z * z).astype(D)
    from triplespec_reducer.idl import exp, alog      # IDL's float EXP/ALOG
    minexp = alog(D(np.finfo(D).tiny))
    g = np.zeros(x.size, dtype=D)
    w = zz < -2 * minexp
    g[w] = exp((-zz[w] / D(2)).astype(D))
    return (p[0] * g).astype(D)


def fit_gauss_const(x, y, start, trace=None):
    """
    MPFITEXPR('p[0] + gaussian(x,p[1:3])', x, y, err, FLOAT(start)) as
    MC_FINDLINES1D calls it: no errors, so plain residuals, in float.
    """
    D = np.float64 if _is_double(x, y) else np.float32
    x = np.asarray(x, dtype=D)
    y = np.asarray(y, dtype=D)

    def fcn(p):
        return (y - (p[0] + gaussian_f(x, p[1:4], D))).astype(D)

    p, status = mpfit(fcn, np.asarray(start, dtype=np.float32), dtype=D, trace=trace)
    return p.astype(float), status


def _exp(a):
    from triplespec_reducer.idl import exp
    return exp(a)


def fit_peak(x, y, estimates, nterms=4, lorentzian=False, weights=None):
    """
    MPFITPEAK with ESTIMATES: MPFITFUN with the peak's analytic derivatives.

    The precision follows the data: double when x or y is double, as when
    Spextool wraps them in DOUBLE(), otherwise float.
    """
    D = np.float64 if _is_double(x, y, estimates) else np.float32
    x = np.asarray(x, dtype=D)
    y = np.asarray(y, dtype=D)
    w = np.ones(x.size, dtype=D) if weights is None else np.asarray(weights, dtype=D)
    wc = np.sqrt(np.abs(w)).astype(D)
    est = np.asarray(estimates, dtype=float)
    p0 = np.zeros(max(nterms, est.size))
    p0[:est.size] = est
    p0 = p0[:nterms]
    pdtype = D
    smax = D(26.) if D == np.float64 else D(13.)

    def model(p):
        p = p.astype(D)
        wid = D(max(abs(p[2]), D(1e-20)))
        u = ((x - p[1]) / wid).astype(D)
        f = np.full(x.size, p[3] if p.size >= 4 else D(0), dtype=D)
        if p.size >= 5:
            f = (f + p[4] * x).astype(D)
        if lorentzian:
            denom = (D(1) / (_sq(u) + D(1))).astype(D)
            pder = np.zeros((x.size, p.size), dtype=D)
            pder[:, 0] = denom
            if p[2] != 0:
                pder[:, 1] = 2 * p[0] * u * denom * denom / p[2]
            pder[:, 2] = pder[:, 1] * u
            if p.size >= 4:
                pder[:, 3] = 1.
            if p.size >= 5:
                pder[:, 4] = x
            return (f + p[0] * denom).astype(D), pder
        mask = (np.abs(u) < smax).astype(D)
        uz = (mask * _exp((-D(0.5) * _sq(u) * mask).astype(D))).astype(D)
        pder = np.zeros((x.size, p.size), dtype=D)
        pder[:, 0] = uz
        if p[2] != 0:
            pder[:, 1] = p[0] * uz * u / p[2]
        pder[:, 2] = pder[:, 1] * u
        if p.size >= 4:
            pder[:, 3] = 1.
        if p.size >= 5:
            pder[:, 4] = x
        return (f + p[0] * uz).astype(D), pder

    def fcn(p):
        return ((y - model(p)[0]) * wc).astype(D)

    def jac(p):
        return (-(model(p)[1] * wc[:, None])).astype(D)

    p, status = mpfit(fcn, p0.astype(pdtype), dtype=D, jac=jac)
    p = p.astype(float)
    p[2] = abs(p[2])
    return model(p.astype(D))[0], p, status
