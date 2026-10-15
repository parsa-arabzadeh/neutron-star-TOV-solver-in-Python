import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Physical Constants (SI Units)
# =============================================================================
G = 6.67430e-11          # Gravitational constant [m^3 kg^-1 s^-2]
c = 299792458.0          # Speed of light [m s^-1]
M_SUN = 1.98847e30       # Solar mass [kg]

# =============================================================================
# Equation of State (Polytropic EOS: P = K * rho^Gamma)
# =============================================================================
K = 2.0e-2               # Polytropic constant [SI units for Gamma=2]
GAMMA = 2.0              # Adiabatic index


def density_from_pressure(P):
    """Rest-mass density rho [kg/m^3] from pressure P [Pa]."""
    if P <= 0.0:
        return 0.0
    return (P / K) ** (1.0 / GAMMA)


def energy_density(P):
    """Total energy density epsilon [J/m^3] from pressure P [Pa]."""
    if P <= 0.0:
        return 0.0
    rho = density_from_pressure(P)
    return rho * c**2 + P / (GAMMA - 1.0)


def sound_speed_squared(P):
    """
    Squared speed of sound (c_s / c)^2 = dP / d(epsilon).
    Used to check the relativistic causality condition (c_s <= c)
    and required for the tidal perturbation equation.
    """
    if P <= 0.0:
        return 0.0
    rho = density_from_pressure(P)
    # d(epsilon)/dP = (c^2 / (GAMMA * K)) * P^(1/GAMMA - 1) + 1/(GAMMA - 1)
    deps_dP = (c**2 * rho) / (GAMMA * P) + 1.0 / (GAMMA - 1.0)
    return 1.0 / deps_dP


# =============================================================================
# Coupled Relativistic Structure & Tidal Perturbation Equations (TOV + Love)
# =============================================================================
def tov_tidal_equations(r, state):
    """
    Coupled ODE system for:
    state[0] = m     : Enclosed gravitational mass [kg]
    state[1] = P     : Pressure [Pa]
    state[2] = m_b   : Enclosed baryonic (rest) mass [kg]
    state[3] = phi   : Metric potential Phi (where g_tt = -e^(2*Phi)) [dimensionless]
    state[4] = y_tid : Relativistic logarithmic derivative of quadrupole
                       metric perturbation H(r), i.e., y = r * H'(r) / H(r)
    """
    m, P, m_b, phi, y_tid = state

    if P <= 0.0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]

    rho = density_from_pressure(P)
    eps = energy_density(P)
    cs2 = sound_speed_squared(P)  # (c_s / c)^2

    # Relativistic compactness parameter at radius r: u = 2Gm / (r c^2)
    u = 2.0 * G * m / (r * c**2)
    if u >= 1.0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]

    e_lambda = 1.0 / (1.0 - u)  # Metric component g_rr = e^(2*Lambda)

    # 1. Gravitational Mass Equation: dm/dr
    dm_dr = 4.0 * np.pi * r**2 * eps / c**2

    # 2. Metric Potential Equation: dPhi/dr
    dphi_dr = (G * (m + 4.0 * np.pi * r**3 * P / c**2)) / (c**2 * r**2 * (1.0 - u))

    # 3. Hydrostatic Equilibrium (TOV) Equation: dP/dr
    dP_dr = -(eps + P) * dphi_dr

    # 4. Baryonic Mass Equation: dm_b/dr
    dm_b_dr = 4.0 * np.pi * r**2 * rho * np.sqrt(e_lambda)

    # 5. Tidal Perturbation Equation: dy/dr (Hinderer 2008 / Damour & Nagar 2009)
    # Describes quadrupolar (l=2) static tidal field perturbations for GW signals
    term_F = 1.0 - (4.0 * np.pi * G * r**2 / c**4) * (eps - P)
    F_r = e_lambda * term_F

    term_Q1 = (4.0 * np.pi * G * r**2 / c**4) * (
        5.0 * eps + 9.0 * P + (eps + P) / cs2
    )
    Q_r = e_lambda * (term_Q1 - 6.0) - 4.0 * (r * dphi_dr) ** 2

    dy_dr = -(y_tid**2 + y_tid * F_r + Q_r) / r

    return [dm_dr, dP_dr, dm_b_dr, dphi_dr, dy_dr]


# =============================================================================
# Surface Event Detection
# =============================================================================
def surface_event(r, state):
    """Stop integration when pressure drops to a negligible surface threshold."""
    return state[1] - 1.0e15  # Cut-off at P = 10^15 Pa to prevent stiffness at P -> 0


surface_event.terminal = True
surface_event.direction = -1


# =============================================================================
# Tidal Love Number (k2) and Deformability (Lambda) Calculator
# =============================================================================
def compute_tidal_parameters(C, y_R):
    """
    Computes the second-order quadrupolar Love number k2 and dimensionless
    tidal deformability Lambda from compactness C = G*M/(R*c^2) and surface y(R).
    Reference: Hinderer, T., ApJ 677, 1216 (2008).
    """
    if C <= 0.0 or C >= 0.5:
        return 0.0, 0.0

    num = (8.0 / 5.0) * (C**5) * ((1.0 - 2.0 * C) ** 2) * (
        2.0 + 2.0 * C * (y_R - 1.0) - y_R
    )

    term1 = 2.0 * C * (
        6.0 - 3.0 * y_R + 3.0 * C * (5.0 * y_R - 8.0)
    )
    term2 = 4.0 * (C**3) * (
        13.0 - 11.0 * y_R + C * (3.0 * y_R - 2.0) + 2.0 * (C**2) * (1.0 + y_R)
    )
    term3 = 3.0 * ((1.0 - 2.0 * C) ** 2) * (
        2.0 - y_R + 2.0 * C * (y_R - 1.0)
    ) * np.log(1.0 - 2.0 * C)

    den = term1 + term2 + term3

    k2 = num / den
    # Dimensionless tidal deformability: Lambda = (2/3) * k2 * C^(-5)
    Lambda = (2.0 / 3.0) * k2 / (C**5)

    return k2, Lambda


# =============================================================================
# Single Star Integration
# =============================================================================
def solve_star(central_pressure):
    """Integrates the coupled TOV + Tidal equations for a given central pressure."""
    r_start = 1.0        # Start at r = 1 m to avoid coordinate singularity at r = 0
    r_end = 100000.0     # Max integration domain: 100 km

    rho_c = density_from_pressure(central_pressure)
    eps_c = energy_density(central_pressure)

    # Small-radius Taylor expansion for initial conditions near r -> 0
    m0 = (4.0 / 3.0) * np.pi * (r_start**3) * (eps_c / c**2)
    P0 = central_pressure
    mb0 = (4.0 / 3.0) * np.pi * (r_start**3) * rho_c
    phi0 = 0.0           # Arbitrary integration constant; shifted at surface later
    y0 = 2.0             # Regularity condition at center for l=2 perturbation

    initial_state = [m0, P0, mb0, phi0, y0]

    sol = solve_ivp(
        tov_tidal_equations,
        (r_start, r_end),
        initial_state,
        method="RK45",
        events=surface_event,
        max_step=50.0,
        rtol=1e-8,
        atol=[1e10, 1e20, 1e10, 1e-10, 1e-8],
    )

    return sol


# =============================================================================
# Stellar Sequence Simulation & Post-Processing
# =============================================================================
if __name__ == "__main__":
    central_pressures = np.logspace(33.0, 36.2, 60)

    results = {
        "Pc": [],
        "rho_c": [],
        "R_km": [],
        "M_sun": [],
        "Mb_sun": [],
        "C": [],
        "z_surf": [],
        "k2": [],
        "Lambda": [],
        "cs_c_max": [],
    }

    sample_star_sol = None
    target_mass_for_profile = 1.4  # Save internal profile of a ~1.4 M_sun star

    for Pc in central_pressures:
        sol = solve_star(Pc)

        if len(sol.t_events[0]) == 0:
            continue

        R = sol.t[-1]
        M = sol.y[0, -1]
        Mb = sol.y[2, -1]
        y_R = sol.y[4, -1]

        # Dimensionless compactness C = GM / (Rc^2)
        C = G * M / (R * c**2)
        if C >= 0.49:
            continue

        # Surface gravitational redshift: z = (1 - 2C)^(-1/2) - 1
        z_surf = (1.0 - 2.0 * C) ** (-0.5) - 1.0

        # Tidal Love number k2 and dimensionless tidal deformability Lambda
        k2, Lambda = compute_tidal_parameters(C, y_R)

        # Central speed of sound (c_s / c)
        cs_c = np.sqrt(sound_speed_squared(Pc))

        results["Pc"].append(Pc)
        results["rho_c"].append(density_from_pressure(Pc))
        results["R_km"].append(R / 1000.0)
        results["M_sun"].append(M / M_SUN)
        results["Mb_sun"].append(Mb / M_SUN)
        results["C"].append(C)
        results["z_surf"].append(z_surf)
        results["k2"].append(k2)
        results["Lambda"].append(Lambda)
        results["cs_c_max"].append(cs_c)

        # Pick a representative star near 1.4 M_sun for internal profile plotting
        if sample_star_sol is None and (M / M_SUN) >= target_mass_for_profile:
            sample_star_sol = sol

    # Convert lists to NumPy arrays
    for key in results:
        results[key] = np.array(results[key])

    # Identify Maximum Mass (TOV Limit) and stability boundary (dM/d(rho_c) > 0)
    idx_max = np.argmax(results["M_sun"])
    M_max = results["M_sun"][idx_max]
    R_max = results["R_km"][idx_max]

    # Interpolate canonical 1.4 M_sun neutron star properties (on the stable branch)
    stable_slice = slice(0, idx_max + 1)
    M_stable = results["M_sun"][stable_slice]
    R_14 = np.interp(1.4, M_stable, results["R_km"][stable_slice])
    Lambda_14 = np.interp(1.4, M_stable, results["Lambda"][stable_slice])
    k2_14 = np.interp(1.4, M_stable, results["k2"][stable_slice])
    z_14 = np.interp(1.4, M_stable, results["z_surf"][stable_slice])

    # =========================================================================
    # Console Summary Output
    # =========================================================================
    print("=" * 60)
    print("   RELATIVISTIC NEUTRON STAR & GRAVITATIONAL WAVE SOLVER")
    print("=" * 60)
    print(f"Computed Models         : {len(results['M_sun'])} stars")
    print(f"Maximum Stable Mass     : {M_max:.3f} M_sun (at R = {R_max:.2f} km)")
    print("-" * 60)
    print("Canonical 1.4 M_sun Neutron Star Observables (GW170817 ref):")
    print(f"  Radius (R_1.4)        : {R_14:.2f} km")
    print(f"  Love Number (k2)      : {k2_14:.4f}")
    print(f"  Tidal Deformability   : Lambda_1.4 = {Lambda_14:.1f}")
    print(f"  Surface Redshift (z)  : {z_14:.3f}")
    print("=" * 60)

    # =========================================================================
    # Multi-Panel Publication Quality Plots
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    plt.subplots_adjust(hspace=0.3, wspace=0.25)

    # 1. Mass-Radius Relation (Stable vs Unstable Branches)
    ax1 = axes[0, 0]
    ax1.plot(
        results["R_km"][: idx_max + 1],
        results["M_sun"][: idx_max + 1],
        color="navy",
        lw=2.2,
        label="Stable Branch ($dM/d\\rho_c > 0$)",
    )
    ax1.plot(
        results["R_km"][idx_max:],
        results["M_sun"][idx_max:],
        color="crimson",
        ls="--",
        lw=1.8,
        label="Unstable Branch ($dM/d\\rho_c < 0$)",
    )
    ax1.scatter(
        [R_max],
        [M_max],
        color="black",
        zorder=5,
        label=f"$M_{{\\rm max}} = {M_max:.2f}\\,M_\\odot$",
    )
    ax1.axhline(1.4, color="gray", ls=":", alpha=0.7, label="$1.4\\,M_\\odot$ Reference")
    ax1.set_xlabel("Radius $R$ (km)", fontsize=11)
    ax1.set_ylabel("Gravitational Mass $M$ ($M_\\odot$)", fontsize=11)
    ax1.set_title("Mass–Radius Relation & Stability", fontsize=12, fontweight="bold")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(True, alpha=0.3)

    # 2. Dimensionless Tidal Deformability (Lambda vs Mass - GW Observable)
    ax2 = axes[0, 1]
    ax2.plot(
        results["M_sun"][: idx_max + 1],
        results["Lambda"][: idx_max + 1],
        color="darkorange",
        lw=2.2,
        label="$\\Lambda(M)$ (Stable Branch)",
    )
    ax2.axvline(1.4, color="gray", ls=":", alpha=0.7)
    ax2.axhline(
        800,
        color="seagreen",
        ls="-.",
        lw=1.5,
        label="LIGO/Virgo GW170817 Limit ($\\Lambda_{1.4} \\lesssim 800$)",
    )
    ax2.scatter([1.4], [Lambda_14], color="red", zorder=5, label=f"$\\Lambda_{{1.4}} = {Lambda_14:.0f}$")
    ax2.set_yscale("log")
    ax2.set_ylim(1, 1e5)
    ax2.set_xlabel("Gravitational Mass $M$ ($M_\\odot$)", fontsize=11)
    ax2.set_ylabel("Tidal Deformability $\\Lambda$ (log scale)", fontsize=11)
    ax2.set_title("GW Tidal Deformability $\\Lambda$ vs. Mass", fontsize=12, fontweight="bold")
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(True, which="both", alpha=0.3)

    # 3. Gravitational vs Baryonic Mass & Redshift
    ax3 = axes[1, 0]
    binding_energy = (results["Mb_sun"][: idx_max + 1] - results["M_sun"][: idx_max + 1]) / results["M_sun"][: idx_max + 1]
    ax3.plot(
        results["C"][: idx_max + 1],
        results["z_surf"][: idx_max + 1],
        color="purple",
        lw=2.2,
        label="Surface Redshift $z_{\\rm surf}$",
    )
    ax3.plot(
        results["C"][: idx_max + 1],
        binding_energy,
        color="teal",
        lw=2.2,
        ls="--",
        label="Fractional Binding Energy $(M_b - M)/M$",
    )
    ax3.plot(
        results["C"][: idx_max + 1],
        results["k2"][: idx_max + 1],
        color="brown",
        lw=2.0,
        ls=":",
        label="Tidal Love Number $k_2$",
    )
    ax3.set_xlabel("Compactness $\\mathcal{C} = GM / (Rc^2)$", fontsize=11)
    ax3.set_ylabel("Dimensionless Observables", fontsize=11)
    ax3.set_title("Relativistic Observables vs. Compactness", fontsize=12, fontweight="bold")
    ax3.legend(loc="upper left", fontsize=9)
    ax3.grid(True, alpha=0.3)

    # 4. Internal Structure Profile for a Representative Star
    ax4 = axes[1, 1]
    if sample_star_sol is not None:
        r_prof = sample_star_sol.t / 1000.0
        m_prof = sample_star_sol.y[0] / sample_star_sol.y[0, -1]
        P_prof = sample_star_sol.y[1] / sample_star_sol.y[1, 0]
        M_sample = sample_star_sol.y[0, -1] / M_SUN

        ax4.plot(r_prof, P_prof, color="royalblue", lw=2.2, label="Normalized Pressure $P(r)/P_c$")
        ax4.plot(r_prof, m_prof, color="firebrick", lw=2.2, label="Normalized Mass $m(r)/M$")
        ax4.set_xlabel("Radial Coordinate $r$ (km)", fontsize=11)
        ax4.set_ylabel("Normalized Profile", fontsize=11)
        ax4.set_title(
            f"Internal Structure Profile ($M = {M_sample:.2f}\\,M_\\odot$)",
            fontsize=12,
            fontweight="bold",
        )
        ax4.legend(loc="center right", fontsize=9)
        ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("tov_gw_analysis.png", dpi=300)
    plt.show()
