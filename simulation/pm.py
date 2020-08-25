def run_sim_pm(filename, n_particles=100, n_grid=30, n_frames=1000, show_ps=False, ps_bins=50):
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
