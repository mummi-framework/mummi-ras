"""Function sample lipids on a grid base on provided lipid densities. Written by Tomas Oppelstrup <oppelstrup2@llnl.gov>."""

import numpy as np


def pdcglob(nx,ny,nchan,C):
    # 40, 40, 8/6 etc, C consentration array

    ng2 = nx*ny

    conc = np.copy(C).reshape(ng2,nchan)
    csum = np.sum(conc,1);
    for i in range(ng2):
        conc[i,:] = conc[i,:] / csum[i]
    csum = np.sum(conc,0)

    crem = np.sum(conc,0)
    gsum = np.zeros(nchan)
    #print "crem = ",crem
    #print "gsum = ",gsum
    gamma = np.zeros(nchan)
    gamma[np.where(crem > 0)] = 1
    #print "gamma = ",gamma

    rnd = np.random.rand(ng2)
    p = np.random.permutation(ng2)

    grid = np.zeros([ng2,nchan])

    for ii in range(ng2):
        i = p[ii]
        c = conc[i,:] * gamma
        c0 = c
        c0[np.where(c0 < 0)] = 0
        s = np.cumsum(c0)

        r = rnd[i] * s[-1]
        j = 0
        while s[j] < r:
            j = j+1
        grid[i,j] = grid[i,j] + 1

        if c[j] <= 0:
            raise Exception("Picked zero concentration channel")
        #if gsum[j] > np.ceil(csum[j])+1e-9:
        #    raise("Oversampling")

        gsum[j] += 1
        crem -= conc[i,:]
        idx = np.where(crem > 0)
        gamma[idx] = (csum[idx] - gsum[idx]) / crem[idx];
        gamma[np.where(crem <= 0)] = 0

    if False:
        idx1 = np.asarray(np.where(gsum > np.ceil(csum)+1e-9))
        idx2 = np.asarray(np.where(gsum < np.floor(csum)-1e-9))
        print("idx1 = ",idx1," len idx1 = ",np.size(idx1))
        print("idx2 = ",idx2," len idx2 = ",len(idx2))
        if len(idx1) > 0 or len(idx2) > 0:
            print("gsum = ",gsum)
            print("csum = ",csum)
            print("crem = ",crem)
            print("idx1 = ",idx1)
            print("idx2 = ",idx2)
            raise Exception("Over- or under-sampling!");

    # check if ny,nx is correct vs nx,ny
    return grid.reshape([ny,nx,nchan])

def reinterp(nxin,nyin,nchan,C,nxout,nyout):
    # We assume the boundary is included in the grid, so that
    # e.g. points 0 and (nx-1) are at left and right edge,
    # respectively.

    # print "%% ",nxin,nyin,nchan,nxout,nyout

    # We want to sample points in the interior of the patch,
    # So let's devide it into nxout x nyout cells, and assume
    # the average concentration for each cell is defined by
    # bilinear interpolation from the input grid to the center
    # of the output cell. The cells are then at e.g.:
    #   (0.5, 1.5, 2.5, ...)*dxout
    dxout = (nxin-1)/float(nxout)
    dyout = (nyin-1)/float(nyout)
    # print "%% ",dxout,dyout

    Cout = np.zeros([nyout,nxout,nchan])
    for iy in range(nyout):
        for ix in range(nxout):

            x = (ix + 0.5)*dxout
            y = (iy + 0.5)*dyout
            xptr = int(x)
            yptr = int(y)
            xfrac = x - xptr
            yfrac = y - yptr

            for j in range(nchan):
                Cout[iy,ix,j] = \
                    (1-xfrac)*(1-yfrac)*C[yptr  ,xptr  ,j] + \
                       xfrac *(1-yfrac)*C[yptr  ,xptr+1,j] + \
                    (1-xfrac)*   yfrac *C[yptr+1,xptr  ,j] + \
                       xfrac *   yfrac *C[yptr+1,xptr+1,j]

    return Cout


def main():
    """Test sampler."""

    nx = 40
    ny = 40
    nchan = 8

    C = np.random.rand(ny,nx,nchan)
    csum = np.sum(C,2);
    for ix in range(nx):
        for iy in range(ny):
            C[iy,ix,:] /= csum[iy,ix]

    print("C = ",C)

    iterlist = [10 , 100, 1000, 10000, 100000]
    #iterlist = [10000 , 100000]
    for niter in iterlist:
        G0 = np.zeros([ny,nx,nchan])
        for iter in range(niter):
            G = pdcglob(nx,ny,nchan,C)
            G0 = G0 + G
        G0 /= niter
        #print G0 - C
        print("niter = %8d, dev = %12.4e, relerr = %12.4e" % \
            (niter,np.max(np.abs(G0 - C)),np.max(np.abs(G0-C)/C)))


if __name__ == "__main__":
    main()
