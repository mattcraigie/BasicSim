import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.fft import rfftn, irfftn
import matplotlib as mpl
mpl.rcParams['animation.ffmpeg_path'] = r'C:\\Program Files (x86)\\FFMpeg\\bin\\ffmpeg.exe'


def apply_boundary(p):
    p = p % n_grid_cells
    return p


def assign_densities(p):
    rho = np.array(np.zeros((n_grid_cells, n_grid_cells, n_grid_cells), dtype=float))
    cell_centers = np.floor(p)
    cell_centers = cell_centers.astype(int)
    d = p - cell_centers
    t = 1 - d

    # assigning mass using CIC
    # might be able to vectorise somehow
    for i in range(n_particles):
        l, m, n = cell_centers[:, i]
        mp = mass_particle

        lp1 = (l + 1) % n_grid_cells
        mp1 = (m + 1) % n_grid_cells
        np1 = (n + 1) % n_grid_cells

        rho[l, m, n] += mp * t[0, i] * t[1, i] * t[2, i]
        rho[l, mp1, n] += mp * t[0, i] * d[1, i] * t[2, i]
        rho[l, m, np1] += mp * t[0, i] * t[1, i] * d[2, i]
        rho[l, mp1, np1] += mp * t[0, i] * d[1, i] * d[2, i]
        rho[lp1, m, n] += mp * d[0, i] * t[1, i] * t[2, i]
        rho[lp1, mp1, n] += mp * d[0, i] * d[1, i] * t[2, i]
        rho[lp1, m, np1] += mp * d[0, i] * t[1, i] * d[2, i]
        rho[lp1, mp1, np1] += mp * d[0, i] * d[1, i] * d[2, i]

    return rho


def apply_force_pm(p, v):
    rho = assign_densities(p)
    rhobar = (np.sum(rho) / n_grid_cells**3)
    delta = (rho - rhobar) / rhobar
    delta_ft = rfftn(delta, axes=(0, 1, 2)) #/ n_grid_cells**(2/3)

    phi_ft = delta_ft * greens
    phi = irfftn(phi_ft, axes=(0, 1, 2))

    g_field_slow = get_acceleration_field(phi)

    g_particle = assign_accelerations(p, g_field_slow)

    g_norm = np.linalg.norm(g_particle, axis=0)
    cond = g_norm > 50
    rhats = g_particle[:, cond] / g_norm[cond]
    g_particle[:, cond] = (50 - g_norm[cond]) * 0.1 * rhats + 50 * rhats


    if plot_powerspectrum:

        # fast
        k_bin, P = power_spectrum_fast(rho, k_norm, n_bins=50, delta=delta)
        r, xi = correlation_function_fast(k_bin, P)

        # slow
        # r, xi = correlation_function_slow(p, n_bins=30)
        # k, P = power_spectrum_slow(r, xi)
        # r = r * n_grid_cells


        for rect, h in zip(powspec, P * k_bin):
            rect.set_height(h)
        if np.min(P*k_bin) < 0:
            ax2.set_ylim([np.min(P* k_bin) * 1.2, np.max(P * k_bin) * 1.2])
        else:
            ax2.set_ylim([0, np.max(P * k_bin) * 1.2])

        for rect, h in zip(corrfunc, xi):
            rect.set_height(h)
        if np.min(xi) < 0:
            ax3.set_ylim([np.min(xi) * 1.2, np.max(xi) * 1.2])
        else:
            ax3.set_ylim([0, np.max(xi) * 1.2])

    v += g_particle * da

    return v


def assign_accelerations(p, g_field):
    g = np.array(np.zeros((3, n_particles), dtype=float))
    cell_centers = np.floor(p)
    cell_centers = cell_centers.astype(int)
    d = p - cell_centers
    t = 1 - d

    # assigning acceleration using CIC
    # might be able to vectorise this somehow
    for i in range(n_particles):
        l, m, n = cell_centers[:, i]

        lp1 = (l + 1) % n_grid_cells
        mp1 = (m + 1) % n_grid_cells
        np1 = (n + 1) % n_grid_cells

        g[:, i] = g_field[:, l, m, n] * t[0, i] * t[1, i] * t[2, i] \
                  + g_field[:, l, mp1, n] * t[0, i] * d[1, i] * t[2, i] \
                  + g_field[:, l, m, np1] * t[0, i] * t[1, i] * d[2, i] \
                  + g_field[:, l, mp1, np1] * t[0, i] * d[1, i] * d[2, i] \
                  + g_field[:, lp1, m, n] * d[0, i] * t[1, i] * t[2, i] \
                  + g_field[:, lp1, mp1, n] * d[0, i] * d[1, i] * t[2, i] \
                  + g_field[:, lp1, m, np1] * d[0, i] * t[1, i] * d[2, i] \
                  + g_field[:, lp1, mp1, np1] * d[0, i] * d[1, i] * d[2, i]

    return g


def get_acceleration_field(phi):
    g = np.zeros((3, n_grid_cells, n_grid_cells, n_grid_cells), dtype=float)

    for i in range(n_grid_cells):
        ip1 = (i + 1) % n_grid_cells
        im1 = (i - 1) % n_grid_cells
        g[0, i, :, :] = -(phi[ip1, :, :] - phi[im1, :, :]) / 2

    for j in range(n_grid_cells):
        jp1 = (j + 1) % n_grid_cells
        jm1 = (j - 1) % n_grid_cells
        g[1, :, j, :] = -(phi[:, jp1, :] - phi[:, jm1, :]) / 2

    for k in range(n_grid_cells):
        kp1 = (k + 1) % n_grid_cells
        km1 = (k - 1) % n_grid_cells
        g[2, :, :, k] = -(phi[:, :, kp1] - phi[:, :, km1]) / 2

    return g


def greens_function(k_vec):
    greens = -3 * Omega_0 / (8 * a * np.linalg.norm(np.sin(k_vec/2), axis=0))
    greens[0, 0, 0] = 0
    return greens


### ~~~ Fast Functions ~~~ ###


def power_spectrum_fast(rho, k_norm, n_bins):
    bin_edges = np.linspace(np.min(k_norm), np.max(k_norm), n_bins + 1)
    bin_width = bin_edges[1] - bin_edges[0]
    bin_mids = (bin_edges + bin_width)[:-1]

    G_r = (rho - np.mean(rho)) / np.mean(rho)

    F_rvec = G_r #- 1 / alpha * R_r
    F_kvec = rfftn(G_r, axes=(0, 1, 2))

    # F_kvec = delta_ft #- alpha * delta_ft_random
    P_kvec = np.real(F_kvec * np.conj(F_kvec))
    P = np.zeros(n_bins)

    for i in range(n_bins):
        condition = np.where(np.abs(k_norm - bin_mids[i]) < bin_width/2)
        P[i] = np.sum(P_kvec[condition]) / np.sum(condition)

    return bin_mids, P


def correlation_function_fast(k, P):

    r = n_grid_cells * k / (2 * np.pi)
    xi = np.zeros(len(r))

    for i in range(len(r)):
        integrand = P * np.sin(k*r[i]) * k
        xi[i] = 1 / (2 * np.pi**2 * r[i]) * np.trapz(integrand, k)

    return r, xi


### ~~~ Slow Functions ~~~ ###

def get_random_catalogue(alpha):
    rand_pos = n_grid_cells * np.random.random((3, alpha * n_particles))
    return rand_pos

def correlation_function_slow(p, n_bins):
    n_random = n_particles * alpha
    bin_mids, DD = separation_counts(p, p, n_bins)
    xi = (n_random / n_particles)**2 * DD / RR - 1
    xi[np.isnan(xi)] = 0
    xi[np.isinf(xi)] = 0

    return bin_mids, xi

def power_spectrum_slow(bin_mids, xi):
    half = int(len(bin_mids)/2 + 1)
    k = 2 * np.pi * bin_mids[:half]
    P = np.fft.rfft(xi)
    return k, np.real(P)

def separation_counts(p1, p2, n_bins):
    r_vec = p1[:, None, :] - p2[:, :, None]  # find N x N x Nd matrix of particle separations
    r = np.linalg.norm(r_vec, axis=0)
    r = np.ravel(r)
    r = r[np.nonzero(r)]

    bin_edges = np.linspace(0, n_grid_cells * np.sqrt(3), n_bins + 1)
    bin_width = bin_edges[1] - bin_edges[0]
    bin_mids = (bin_edges + bin_width)[:-1]
    counts, _ = np.histogram(r, bins=bin_edges)

    return bin_mids, counts


def fieldshow(field, n):
    plt.close()
    plt.imshow(np.sum(field, axis=0))
    plt.savefig(n)

def histshow(field, n):
    plt.close()
    plt.hist(field.flatten())
    plt.savefig(n)



def update(i):
    global position, velocity, frame, a

    if frame % int(n_frames / 10) == 0:
        print(str(int(frame / int(n_frames / 10))) + '%')


    position += velocity * da
    position = apply_boundary(position)
    velocity = apply_force_pm(position, velocity)

    points_sim.set_data(position[0, :], position[1, :])

    a += da
    frame += 1

    return points_sim,


if __name__ == '__main__':
    print('Do not run this file. Import using >> from simulation import run_sim_pm')






