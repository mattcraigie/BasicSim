import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.fft import rfftn, irfftn


class SimPM():
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


    def __init__(self, filename, n_particles=100, n_grid=30, n_frames=1000, show_ps=False, ps_bins=50, gravity_factor=1):

        self.filename = filename
        self.n_particles = n_particles
        self.n_grid = n_grid
        self.n_frames = n_frames
        self.show_ps = show_ps
        self.ps_bins = ps_bins


        self.L_box = 1
        self.H_0 = 1
        self.G = 1
        self.Omega_m0 = 0.3
        self.Omega_k0 = 0
        self.Omega_L0 = 0.7
        self.mass_particle = 1
        self.Omega_0 = 1

        self.r_0 = self.L_box / self.n_grid
        self.t_0 = 1 / self.H_0
        rho_0 = 3 * self.H_0 ** 2 * self.Omega_m0 / (8 * np.pi * self.G)
        self.v_0 = self.r_0 / self.t_0
        phi_0 = self.v_0 ** 2

        # Set initial positions and velocities
        self.position = n_grid * np.random.random((3, n_particles))

        # x = np.linspace(0, 1, 7)[:-1]
        # self.n_particles = 6 ** 3
        # self.position = np.array([[i, j, k] for i in x for j in x for k in x]).transpose() * (n_grid -1)
        self.velocity = 1 * (0.5 - np.random.random((3, self.n_particles)))

        # frame counter
        self.frame = 0

        # time steps
        self.a = 0.4
        self.da = 1e-2

        x = np.arange(0, n_grid)
        self.pos_grid = np.array(np.meshgrid(x, x, x))
        self.k_grid = 2 * np.pi * self.pos_grid / n_grid
        self.k_norm = np.linalg.norm(self.k_grid, axis=0)

        half = int(n_grid / 2 + 1)
        self.half_k_grid = self.k_grid[:, :, :, :half]
        self.half_k_norm = self.k_norm[:, :, :half]

        if show_ps:
            # slowfunc setup
            self.alpha = 10
            self.DD = None

            plt.ion()

            self.fig = plt.figure(figsize=(12, 8), constrained_layout=True)
            gs = self.fig.add_gridspec(2, 3)

            self.ax1 = self.fig.add_subplot(gs[0:, 0:2])
            self.ax1.set_title('Main Simulation [a={}]'.format(self.a))
            self.ax1.set_xlabel('Normalised Comoving $x$')
            self.ax1.set_ylabel('Normalised Comoving $y$')
            self.ax1.set_xlim(0, 1)
            self.ax1.set_ylim(0, 1)

            self.ax2 = self.fig.add_subplot(gs[0, 2])
            self.ax2.set_title('Power Spectrum [$k\\ P(k)$]')
            self.ax2.set_xlim(np.min(self.k_norm), np.max(self.k_norm))
            self.ax2.set_ylim(0, 0.5)

            self.ax3 = self.fig.add_subplot(gs[1, 2])
            self.ax3.set_title('Correlation Function [$\\xi(r)$]')
            self.ax3.set_xlim(0, n_grid * np.sqrt(3))
            self.ax3.set_ylim(-1.2, 1.2)

            self.points_sim, = self.ax1.plot([], [], 'o', markersize=1)

            bin_edges = np.linspace(np.min(self.k_norm), np.max(self.k_norm), self.ps_bins + 1)
            bin_width = bin_edges[1] - bin_edges[0]
            bin_mids = (bin_edges + bin_width)[:-1]
            self.powspec = self.ax2.bar(bin_mids, np.zeros(len(bin_mids)), align='center',
                                   width= np.max(bin_edges) / self.ps_bins, color='white', edgecolor='royalblue')

            bin_edges = np.linspace(0, self.n_grid * np.sqrt(3), self.ps_bins + 1)
            bin_width = bin_edges[1] - bin_edges[0]
            bin_mids = (bin_edges + bin_width)[:-1]
            self.corrfunc = self.ax3.bar(bin_mids, np.zeros(len(bin_mids)), align='center',
                               width=(self.n_grid * np.sqrt(3)) / self.ps_bins, color='white', edgecolor='royalblue')

            print('Running simulation. \nProgress:')

            ani = FuncAnimation(self.fig, self.update, frames=self.n_frames)
            f = r"sim_pm.mp4"
            writervideo = FFMpegWriter(fps=60)

            if filename[:-4] != '.mp4':
                f = filename + '.mp4'
            else:
                f = filename

            ani.save(filename, writer=writervideo)

        if not show_ps:

            # Set the axes on which the points will be shown
            plt.ion()  # Set interactive mode on
            fig, self.ax1 = plt.subplots(figsize=(8, 8))
            self.ax1.set_title('Main Simulation [a={}]'.format(self.a))
            self.ax1.set_xlabel('Normalised Comoving $x$')
            self.ax1.set_ylabel('Normalised Comoving $y$')
            self.ax1.set_xlim(0, 1)
            self.ax1.set_ylim(0, 1)

            # Create command which will plot the positions of the particles
            self.points_sim, = plt.plot([], [], 'o', markersize=1)

            ani = FuncAnimation(fig, self.update, frames=n_frames)  # interval=frame_duration)
            writervideo = FFMpegWriter(fps=60)

            if filename[:-3] != '.mp4':
                f = filename + '.mp4'
            else:
                f = filename

            ani.save(filename, writer=writervideo)


    def apply_boundary(self, p):
        p = p % self.n_grid
        return p

    def assign_densities(self, p):
        rho = np.array(np.zeros((self.n_grid, self.n_grid, self.n_grid), dtype=float))

        cell_centers = np.floor(p)
        cell_centers = cell_centers.astype(int)
        d = p - cell_centers
        t = 1 - d


        for i in range(self.n_particles):
            l, m, n = cell_centers[:, i]
            mp = self.mass_particle

            lp1 = (l + 1) % (self.n_grid)
            mp1 = (m + 1) % (self.n_grid)
            np1 = (n + 1) % (self.n_grid)

            rho[l, m, n] += mp * t[0, i] * t[1, i] * t[2, i]
            rho[l, mp1, n] += mp * t[0, i] * d[1, i] * t[2, i]
            rho[l, m, np1] += mp * t[0, i] * t[1, i] * d[2, i]
            rho[l, mp1, np1] += mp * t[0, i] * d[1, i] * d[2, i]
            rho[lp1, m, n] += mp * d[0, i] * t[1, i] * t[2, i]
            rho[lp1, mp1, n] += mp * d[0, i] * d[1, i] * t[2, i]
            rho[lp1, m, np1] += mp * d[0, i] * t[1, i] * d[2, i]
            rho[lp1, mp1, np1] += mp * d[0, i] * d[1, i] * d[2, i]

        return rho

    def f(self, a):
        return (a**-1 * self.Omega_m0 * self.Omega_k0 * a + self.Omega_L0 * a**3) ** (-0.5)

    def apply_force_pm(self, p, v):
        rho = self.assign_densities(p)
        rhobar = (np.sum(rho) / self.n_grid ** 3)
        delta = (rho - rhobar) / rhobar
        delta_ft = rfftn(delta, axes=(0, 1, 2))  # / n_grid_cells**(2/3)

        phi_ft = delta_ft * self.greens_function()
        phi = irfftn(phi_ft, axes=(0, 1, 2))

        g_field = self.get_acceleration_field(phi)

        g_particle = self.assign_accelerations(p, g_field)

        g_norm = np.linalg.norm(g_particle, axis=0)
        cond = g_norm > 50
        rhats = g_particle[:, cond] / g_norm[cond]
        g_particle[:, cond] = (50 - g_norm[cond]) * 0.1 * rhats + 50 * rhats

        if self.show_ps:

            # fast
            k_bin, P = self.power_spectrum_fast(rho, self.k_norm, self.ps_bins)
            r, xi = self.correlation_function_fast(k_bin, P)

            P = P * k_bin

            for rect, h in zip(self.powspec, P):
                rect.set_height(h)
            if np.min(P) < 0:
                self.ax2.set_ylim([np.min(P) * 1.2, np.max(P) * 1.2])
            else:
                self.ax2.set_ylim([0, np.max(P) * 1.2])

            for rect, h in zip(self.corrfunc, xi):
                rect.set_height(h)
            if np.min(xi) < 0:
                self.ax3.set_ylim([np.min(xi) * 1.2, np.max(xi) * 1.2])
            else:
                self.ax3.set_ylim([0, np.max(xi) * 1.2])

        v += self.f(self.a) * g_particle * self.da

        return v

    def assign_accelerations(self, p, g_field):
        g = np.array(np.zeros((3, self.n_particles), dtype=float))
        cell_centers = np.floor(p)
        cell_centers = cell_centers.astype(int)
        d = p - cell_centers
        t = 1 - d

        # assigning acceleration using CIC
        # might be able to vectorise this somehow
        for i in range(self.n_particles):
            l, m, n = cell_centers[:, i]

            lp1 = (l + 1) % self.n_grid
            mp1 = (m + 1) % self.n_grid
            np1 = (n + 1) % self.n_grid

            g[:, i] = g_field[:, l, m, n] * t[0, i] * t[1, i] * t[2, i] \
                      + g_field[:, l, mp1, n] * t[0, i] * d[1, i] * t[2, i] \
                      + g_field[:, l, m, np1] * t[0, i] * t[1, i] * d[2, i] \
                      + g_field[:, l, mp1, np1] * t[0, i] * d[1, i] * d[2, i] \
                      + g_field[:, lp1, m, n] * d[0, i] * t[1, i] * t[2, i] \
                      + g_field[:, lp1, mp1, n] * d[0, i] * d[1, i] * t[2, i] \
                      + g_field[:, lp1, m, np1] * d[0, i] * t[1, i] * d[2, i] \
                      + g_field[:, lp1, mp1, np1] * d[0, i] * d[1, i] * d[2, i]

        return g

    def get_acceleration_field(self, phi):
        g = np.zeros((3, self.n_grid, self.n_grid, self.n_grid), dtype=float)

        for i in range(self.n_grid):
            ip1 = (i + 1) % self.n_grid
            im1 = (i - 1) % self.n_grid
            g[0, i, :, :] = -(phi[ip1, :, :] - phi[im1, :, :]) / 2

        for j in range(self.n_grid):
            jp1 = (j + 1) % self.n_grid
            jm1 = (j - 1) % self.n_grid
            g[1, :, j, :] = -(phi[:, jp1, :] - phi[:, jm1, :]) / 2

        for k in range(self.n_grid):
            kp1 = (k + 1) % self.n_grid
            km1 = (k - 1) % self.n_grid
            g[2, :, :, k] = -(phi[:, :, kp1] - phi[:, :, km1]) / 2

        return g

    def greens_function(self):
        # greens = - G_factor / (self.a * np.linalg.norm(np.sin(self.half_k_grid / 2), axis=0))
        greens = -3 * self.Omega_0 / (8 * self.a * np.linalg.norm(np.sin(self.half_k_grid / 2), axis=0))
        greens[0, 0, 0] = 0
        return greens

    ### ~~~ Fast Functions ~~~ ###

    def power_spectrum_fast(self, rho, k_norm, n_bins):
        bin_edges = np.linspace(np.min(k_norm), np.max(k_norm), n_bins + 1)
        bin_width = bin_edges[1] - bin_edges[0]
        bin_mids = (bin_edges + bin_width)[:-1]

        G_r = (rho - np.mean(rho)) / np.mean(rho)

        F_rvec = G_r
        F_kvec = rfftn(F_rvec, axes=(0, 1, 2))

        P_kvec = np.real(F_kvec * np.conj(F_kvec))
        P = np.zeros(n_bins)

        things = []

        for i in range(n_bins):
            for i in range(n_bins):
                condition = np.where(np.abs(k_norm - bin_mids[i]) < bin_width / 2)
            print(condition)

            if np.sum(condition) != 0:
                P[i] = np.sum(P_kvec[condition]) / np.sum(condition)
            else:
                P[i] = 0

        return bin_mids, P

    def correlation_function_fast(self, k, P):

        r = self.n_grid * k / (2 * np.pi)
        xi = np.zeros(len(r))

        for i in range(len(r)):
            integrand = P * np.sin(k * r[i]) * k
            xi[i] = 1 / (2 * np.pi ** 2 * r[i]) * np.trapz(integrand, k)

        return r, xi

    ### ~~~ Slow Functions ~~~ ###

    def get_random_catalogue(self, alpha):
        rand_pos = self.n_grid * np.random.random((3, alpha * self.n_particles))
        return rand_pos

    def correlation_function_slow(self, p, n_bins):
        n_random = self.n_grid * self.alpha
        bin_mids, DD = self.separation_counts(p, p, n_bins)
        xi = (n_random / self.n_particles) ** 2 * DD / self.RR - 1
        xi[np.isnan(xi)] = 0
        xi[np.isinf(xi)] = 0

        return bin_mids, xi

    def power_spectrum_slow(self, bin_mids, xi):
        half = int(len(bin_mids) / 2 + 1)
        k = 2 * np.pi * bin_mids[:half]
        P = np.fft.rfft(xi)
        return k, np.real(P)

    def separation_counts(self, p1, p2, n_bins):
        r_vec = p1[:, None, :] - p2[:, :, None]  # find N x N x Nd matrix of particle separations
        r = np.linalg.norm(r_vec, axis=0)
        r = np.ravel(r)
        r = r[np.nonzero(r)]

        bin_edges = np.linspace(0, self.n_grid * np.sqrt(3), n_bins + 1)
        bin_width = bin_edges[1] - bin_edges[0]
        bin_mids = (bin_edges + bin_width)[:-1]
        counts, _ = np.histogram(r, bins=bin_edges)

        return bin_mids, counts


    def update(self, i):
        try:
            if self.frame % int(self.n_frames / 10) == 0:
                print(str(int(self.frame / int(self.n_frames / 10)) * 10) + '%')
        except ZeroDivisionError:
            pass


        self.position += self.velocity * self.da * self.f(self.a) / self.a**2
        self.position = self.apply_boundary(self.position)
        self.velocity = self.apply_force_pm(self.position, self.velocity)

        self.points_sim.set_data(self.position[0, :] / self.n_grid, self.position[1, :] / self.n_grid)

        self.a += self.da
        self.frame += 1

        self.ax1.set_title('Main Simulation [a={:.2f}]'.format(self.a))

        return self.points_sim,

