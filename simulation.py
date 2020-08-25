import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.fft import rfftn, irfftn
import matplotlib as mpl
mpl.rcParams['animation.ffmpeg_path'] = r'C:\\Program Files (x86)\\FFMpeg\\bin\\ffmpeg.exe'


def apply_boundary(p):
    """Defines the boundary conditions - we want a periodic boundary"""
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

    if frame % 10 == 0:
        print(frame)

    position += velocity * da
    position = apply_boundary(position)
    velocity = apply_force_pm(position, velocity)

    points_sim.set_data(position[0, :], position[1, :])

    a += da
    frame += 1

    return points_sim,


def run_sim_pm(filename, n_particles=1000, n_grid=30, n_frames=1000, show_ps=False, ps_bins=50):
    """
    Creates and runs a simulation using a particle-mesh N-body implementation. An mp4 file of the simulation is saved at
    the location specified by filename. The particle-mesh method scales with O(N + G*log[G]), and operations are done in
    place so large number of particles can be computed in reasonable times. Setting show_ps=True will slow computation
    but is also calculated with O(G*log[G]).

    :param filename: name of the file to save the simulation video
    :param n_particles: the number of particles to use in the simulation.
    :param n_grid: number of grid cells along a single dimension, total grid cells will be n_grid**3
    :param n_frames: number of frames to run the simulation for
    :param show_ps: set to true to plot the power spectrum (k * P[k]) and correlation function (xi[r]) alongside the
    simulation
    :param ps_bins: the number of bins to use for the power spectrum and correlation function if show_ps=True
    :return: returns None
    """

    np.random.seed(42)

    # Sim setup
    n_particles = 15 ** 3
    n_frames = 10
    frame_duration = 100

    n_grid_cells = 30

    L_box = 1
    H_0 = 1
    G = 1
    Omega_m0 = 1
    mass_particle = 1
    Omega_0 = 1

    r_0 = L_box / n_grid_cells
    t_0 = 1 / H_0
    rho_0 = 3 * H_0 ** 2 * Omega_m0 / (8 * np.pi * G)
    v_0 = r_0 / t_0
    phi_0 = v_0 ** 2

    # Set initial positions and velocities
    position = n_grid_cells * np.random.random((3, n_particles))

    x = np.linspace(0, 1, 15)
    position = np.array([[i, j, k] for i in x for j in x for k in x]).transpose() * n_grid_cells
    velocity = 0 * (0.5 - np.random.random((3, 15 ** 3)))

    # position = np.array([[25., 35.], [0., 0.], [0., 0.]])
    # velocity = np.array([[0., 0.], [10., 10.], [0., 0.]])

    # frame counter
    frame = 0

    #
    a = 0.1
    da = 1e-2

    plot_powerspectrum = show_ps

    x = np.arange(0, n_grid_cells)
    pos_grid = np.array(np.meshgrid(x, x, x))
    k_grid = 2 * np.pi * pos_grid / n_grid_cells
    k_norm = np.linalg.norm(k_grid, axis=0)

    half = int(n_grid_cells / 2 + 1)
    k_grid = k_grid[:, :, :, :half]
    k_norm = k_norm[:, :, :half]

    greens = greens_function(k_grid)

    if plot_powerspectrum:
        alpha = 100
        pos_rand = get_random_catalogue(alpha)

        # For fast
        R_r = assign_densities(pos_rand)

        # For slow:
        # _, RR = separation_counts(pos_rand, pos_rand, n_bins=30)

        plt.ion()

        fig = plt.figure(figsize=(12, 8), constrained_layout=True)
        gs = fig.add_gridspec(2, 3)
        ax1 = fig.add_subplot(gs[0:, 0:2])
        ax1.set_title('Main Simulation')

        ax2 = fig.add_subplot(gs[0, 2])
        ax2.set_title('Power Spectrum [$k\\ P(k)$]')

        ax3 = fig.add_subplot(gs[1, 2])
        ax3.set_title('Correlation Function [$\\xi(r)$]')

        ax1.set_xlim(0, n_grid_cells)
        ax1.set_ylim(0, n_grid_cells)

        ax2.set_xlim(np.min(k_norm), np.max(k_norm))
        ax2.set_ylim(0, 0.5)

        ax3.set_xlim(0, n_grid_cells * np.sqrt(3))
        ax3.set_ylim(-1.2, 1.2)

        points_sim, = ax1.plot([], [], 'o', markersize=1)

        n_bins = 60
        bin_edges = np.linspace(np.min(k_norm), np.max(k_norm), n_bins + 1)
        bin_width = bin_edges[1] - bin_edges[0]
        bin_mids = (bin_edges + bin_width)[:-1]
        powspec = ax2.bar(bin_mids, np.zeros(len(bin_mids)), align='center', width=9 / n_bins, color='white',
                          edgecolor='royalblue')

        n_bins = 60
        bin_edges = np.linspace(0, n_grid_cells * np.sqrt(3), n_bins + 1)
        bin_width = bin_edges[1] - bin_edges[0]
        bin_mids = (bin_edges + bin_width)[:-1]
        corrfunc = ax3.bar(bin_mids, np.zeros(len(bin_mids)), align='center',
                           width=(n_grid_cells * np.sqrt(3)) / n_bins, color='white', edgecolor='royalblue')

        ani = FuncAnimation(fig, update, frames=n_frames)
        f = r"sim_pm.mp4"
        writervideo = FFMpegWriter(fps=60)
        ani.save(f, writer=writervideo)

    if not plot_powerspectrum:

        # Set the axes on which the points will be shown
        plt.ion()  # Set interactive mode on
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlim(0, n_grid_cells)
        ax.set_ylim(0, n_grid_cells)

        # Create command which will plot the positions of the particles
        points_sim, = plt.plot([], [], 'o', markersize=1)

        ani = FuncAnimation(fig, update, frames=n_frames)  # interval=frame_duration)
        writervideo = FFMpegWriter(fps=60)

        if filename[:-3] != '.mp4':
            f = filename + '.mp4'
        else:
            f = filename

        ani.save(filename, writer=writervideo)

    return None


if __name__ == '__main__':
    print('Do not run this file. Import using >> from simulation import run_sim_pm')






