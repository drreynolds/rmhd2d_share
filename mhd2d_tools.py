# matplotlib-based plotting utilities for RMHD2d tests
# Daniel R. Reynolds, reynolds@smu.edu


##########
def load_cons():
    """Returns times, energy, mass arrays from the conservation.txt file"""
    import shlex
    import numpy as np
    # load the conservation.txt file
    # (columns are:  time, total energy, total mass)
    f = open("conservation.txt")
    times=[];  tenergy=[];  totmass=[]; 
    for line in f:
        txt = shlex.split(line)
        times.append(float(txt[0]));
        tenergy.append(float(txt[1]));
        totmass.append(float(txt[2]));
    f.close()
    return [np.array(times), np.array(tenergy), np.array(totmass)]

##########
def load_energies():
    """Returns kenergy, menergy, ienergy from the energyHistory.txt file"""
    import shlex
    import numpy as np
    # load the energyHistory.txt file
    # (columns are:  time, kinetic energy, magnetic energy, internal energy)
    f = open("energyHistory.txt")
    kenergy=[];  menergy=[];  ienergy=[]; 
    for line in f:
        txt = shlex.split(line)
        kenergy.append(float(txt[1]));
        menergy.append(float(txt[2]));
        ienergy.append(float(txt[3]));
    f.close()
    return [np.array(kenergy), np.array(menergy), np.array(ienergy)]

##########
def load_props():
    """Returns mu, eta, kappa, gamma from prop.inp"""
    import shlex
    # load the physical simulation properties
    f = open("prop.inp")
    mu=0.0;  eta=0.0;  kappa=0.0;  gamma=5.0/3.0;
    for line in f:
        txt = shlex.split(line)
        if ("mu" in txt):
            mu = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("eta" in txt):
            eta = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("kappa" in txt):
            kappa = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("gamma" in txt):
            gamma = float(txt[len(txt)-1].replace('d','e').replace(',',''))
    f.close()
    return [mu, eta, kappa, gamma]    

##########
def load_mhd():
    """Returns xl,xr,yl,yr,zl,zr,ndump from mhd.inp"""
    import shlex
    # load the grid bounds
    f = open("mhd.inp")
    xl=0.0;  xr=1.0;  yl=0.0;  yr=1.0;  zl=0.0;  zr=1.0;  ndump=1;
    for line in f:
        txt = shlex.split(line)
        if ("xl" in txt):
            xl = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("xr" in txt):
            xr = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("yl" in txt):
            yl = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("yr" in txt):
            yr = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("zl" in txt):
            zl = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("zr" in txt):
            zr = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("ndump" in txt):
            ndump = int(txt[len(txt)-1].replace(',',''))
    f.close()
    return [xl, xr, yl, yr, zl, zr, ndump]

##########
def load_mesh():
    """Returns nx,ny,nz,xprocs,yprocs,zprocs,xbc,ybc,zbc from mesh.inp"""
    import shlex
    # load the mesh information
    f = open("mesh.inp")
    nx=1; ny=1; nz=1; xbc=0; ybc=0; zbc=0;
    xprocs=1; yprocs=1; zprocs=1; 
    for line in f:
        txt = shlex.split(line)
        if ("nx" in txt):
            nx = int(txt[len(txt)-1].replace(',',''))
        elif ("ny" in txt):
            ny = int(txt[len(txt)-1].replace(',',''))
        elif ("nz" in txt):
            nz = int(txt[len(txt)-1].replace(',',''))
        elif ("xprocs" in txt):
            xprocs = int(txt[len(txt)-1].replace(',',''))
        elif ("yprocs" in txt):
            yprocs = int(txt[len(txt)-1].replace(',',''))
        elif ("zprocs" in txt):
            zprocs = int(txt[len(txt)-1].replace(',',''))
        elif ("xbc" in txt):
            xbc = int(txt[len(txt)-1].replace(',',''))
        elif ("ybc" in txt):
            ybc = int(txt[len(txt)-1].replace(',',''))
        elif ("zbc" in txt):
            zbc = int(txt[len(txt)-1].replace(',',''))
    f.close()
    return [nx, ny, nz, xprocs, yprocs, zprocs, xbc, ybc, zbc];

##########
def load_RunHistory():
    """Returns Newton, Krylov, Fnorm, divB, SolTime from Run.history"""
    import shlex
    import numpy as np
    # load the Run history information
    f = open("Run.history")
    Newton=[];  Krylov=[];  Fnorm=[];  divB=[];  SolTime=0.0;
    for line in f:
        txt = shlex.split(line)
        if ("Num. Newton Iterations" in txt):
            Newton.append(int(txt[len(txt)-1]))
        elif ("Number Krylov Iterations" in txt):
            Krylov.append(int(txt[len(txt)-1]))
        elif ("Final scaled norm of f" in txt):
            Fnorm.append(float(txt[len(txt)-1]).replace('d','e'))
        elif ("divB" in txt):
            divB = float(txt[len(txt)-4].replace('d','e'))
        elif ("Total Solution Time" in txt):
            SolTime = float(txt[len(txt)-1].replace('d','e'))
    f.close()
    return [np.array(Newton), np.array(Krylov), np.array(Fnorm), 
            np.array(divB), np.array(SolTime)];
    
##########
def load_vals(tdump):
    """Returns x,y,rho,u,v,w,bx,by,bz,p,dB,j,te from a given data dump"""
    import shlex
    import numpy as np
    sdump = repr(tdump).zfill(6)
    outfile = 'output.001.' + sdump
    nx,ny,nz,xprocs,yprocs,zprocs,xbc,ybc,zbc = load_mesh()
    f = open(outfile)
    x   = np.zeros(nx*ny,dtype=float); y  = np.zeros(nx*ny,dtype=float); 
    rho = np.zeros(nx*ny,dtype=float); u  = np.zeros(nx*ny,dtype=float); 
    v   = np.zeros(nx*ny,dtype=float); w  = np.zeros(nx*ny,dtype=float); 
    bx  = np.zeros(nx*ny,dtype=float); by = np.zeros(nx*ny,dtype=float); 
    bz  = np.zeros(nx*ny,dtype=float); p  = np.zeros(nx*ny,dtype=float); 
    dB  = np.zeros(nx*ny,dtype=float); j  = np.zeros(nx*ny,dtype=float); 
    te  = np.zeros(nx*ny,dtype=float); 
    idx = 0;
    for line in f:
        txt = shlex.split(line)
        if (("#" not in txt) and (txt != [])):
            x[idx]   = float(txt[0])
            y[idx]   = float(txt[1])
            rho[idx] = float(txt[2])
            u[idx]   = float(txt[3])
            v[idx]   = float(txt[4])
            w[idx]   = float(txt[5])
            bx[idx]  = float(txt[6])
            by[idx]  = float(txt[7])
            bz[idx]  = float(txt[8])
            p[idx]   = float(txt[9])
            dB[idx]  = float(txt[10])
            j[idx]   = float(txt[11])
            te[idx]  = float(txt[12])
            idx += 1
    f.close()
    # reshape 1D arrays to 2D
    xarr = np.array(x[0:nx]);
    yarr = np.zeros(ny);
    idx = 0;
    for k in range(ny):
        yarr[idx] = y[k*nx];
        idx += 1;
    rho2D = np.reshape(rho, (nx,ny), order='F')
    u2D   = np.reshape(u,   (nx,ny), order='F')
    v2D   = np.reshape(v,   (nx,ny), order='F')
    w2D   = np.reshape(w,   (nx,ny), order='F')
    bx2D  = np.reshape(bx,  (nx,ny), order='F')
    by2D  = np.reshape(by,  (nx,ny), order='F')
    bz2D  = np.reshape(bz,  (nx,ny), order='F')
    p2D   = np.reshape(p,   (nx,ny), order='F')
    dB2D  = np.reshape(dB,  (nx,ny), order='F')
    j2D   = np.reshape(j,   (nx,ny), order='F')
    te2D  = np.reshape(te,  (nx,ny), order='F')
    return [xarr,yarr,rho2D,u2D,v2D,w2D,bx2D,by2D,bz2D,p2D,dB2D,j2D,te2D]
##########

