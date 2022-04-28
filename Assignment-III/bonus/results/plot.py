import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

for lbl in ['2048', '4096', '4096_avg', '4096_custom_dist']:
    fdata = np.genfromtxt(f'results_runtime_{lbl}')

    ts = fdata.T[2:]

    # mean
    t = ts.mean(axis=0)

    nprocs = fdata.T[0].astype(int)
    M     = fdata.T[1].astype(int)

    plt.plot(nprocs, t, '+-')

    # --- linear regression ---
    y = t * np.sqrt(nprocs) / (2*M*M)
    x = M /np.sqrt(nprocs)
    
    # linregress for less than 64 processes
    mask = nprocs < 64
    res = linregress(x[mask], y[mask])
    t_comm, t_flop = res.intercept, res.slope
    print(f'{lbl}: t_comm={t_comm}, t_flop={t_flop} (R2={res.rvalue})')

    n = np.linspace(16, 49)
    t = 2*M[0]**3 / n * t_flop + 2*M[0]**2 / np.sqrt(n) * t_comm

    plt.plot(n, t, 'r--')
    
    # linregress for more than 64 processes
    mask = np.logical_and(nprocs >= 64, nprocs != 256)
    res = linregress(x[mask], y[mask])
    t_comm, t_flop = res.intercept, res.slope
    print(f'{lbl}: t_comm={t_comm}, t_flop={t_flop} (R2={res.rvalue})')

    n = np.linspace(64, 225)
    t = 2*M[0]**3 / n * t_flop + 2*M[0]**2 / np.sqrt(n) * t_comm

    plt.plot(n, t, 'r--')
    
    # -- -- 



    plt.xlabel('number of processes')
    plt.ylabel('runtime (s)')
    plt.title(f'runtime of Fox algorithm on Dardel ($M\\approx {M[0]}$)')
    plt.xticks(nprocs)
    plt.tight_layout()
    plt.savefig(f'runtime_{lbl}.svg')
    plt.cla()

    