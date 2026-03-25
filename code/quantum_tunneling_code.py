#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quantum Tunneling: Microscopic Marvel and Transdimensional Revolution
Complete Code Collection for the Book

Author: Liu Ming (pen name)
Contact: huanxinmengmeng@126.com

This file contains all teaching demonstration code and science fiction speculations.
Real physics code can be used for teaching; science fiction code is marked with warnings
and is intended for thought experiments only.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e, e

# ============================================================================
# 1. Real Physics Modules
# ============================================================================

def transmission_rectangle(V0, E, d):
    """
    Rectangular barrier tunneling probability (real physics).

    Parameters:
        V0 : float
            Barrier height (eV).
        E : float
            Particle energy (eV).
        d : float
            Barrier width (nm).

    Returns:
        float
            Transmission coefficient.
    """
    V0_J = V0 * e
    E_J = E * e
    d_m = d * 1e-9
    if E >= V0:
        k = np.sqrt(2 * m_e * (E_J - V0_J)) / hbar
        T = 1 / (1 + (V0_J**2 * np.sin(k * d_m)**2) / (4 * E_J * (E_J - V0_J)))
    else:
        kappa = np.sqrt(2 * m_e * (V0_J - E_J)) / hbar
        T = 1 / (1 + (V0_J**2 * np.sinh(kappa * d_m)**2) / (4 * E_J * (V0_J - E_J)))
    return T

def stm_current(distance, phi=4.5, V=0.1):
    """
    STM tunneling current simulation (real physics).

    Parameters:
        distance : float
            Tip-sample distance (nm).
        phi : float, optional
            Average barrier height (eV). Default is 4.5.
        V : float, optional
            Bias voltage (V). Default is 0.1.

    Returns:
        float
            Relative tunneling current (arbitrary units).
    """
    kappa = np.sqrt(2 * m_e * phi * e) / hbar
    kappa_nm = kappa * 1e9
    return V * np.exp(-2 * kappa_nm * distance)

def tfet_transfer(Vg, SS=45, Vth=0.2, Ioff=1e-12, Ion=1e-3):
    """
    TFET transfer characteristic simulation (based on literature).

    Parameters:
        Vg : array_like
            Gate voltage (V).
        SS : float, optional
            Subthreshold swing (mV/dec). Default is 45.
        Vth : float, optional
            Threshold voltage (V). Default is 0.2.
        Ioff : float, optional
            Off-state current (A/µm). Default is 1e-12.
        Ion : float, optional
            On-state current (A/µm). Default is 1e-3.

    Returns:
        ndarray
            Drain current Id (A/µm).
    """
    Vg = np.asarray(Vg)
    subthreshold = Ioff * 10**((Vg - Vth) * 1000 / SS)
    above_threshold = np.zeros_like(Vg)
    mask = Vg >= Vth
    above_threshold[mask] = Ion * (Vg[mask] - Vth)**1.5
    Id = np.where(Vg < Vth, subthreshold, above_threshold)
    return Id

def nonhermitian_lz(alpha, Omega, gamma):
    """
    Non-Hermitian Landau-Zener tunneling probability (PRL 132, 156802).

    Parameters:
        alpha : float
            Sweep rate.
        Omega : float
            Coupling strength.
        gamma : float
            Dissipation rate.

    Returns:
        float
            Tunneling probability P (clipped to [0, 1]).
    """
    gamma_c = 0.46 * Omega
    if gamma < gamma_c:
        P = 1 - np.exp(-np.pi * Omega**2 / (2 * alpha))
    else:
        delta = Omega**2 - gamma**2
        if delta < 0:
            delta = 0.0
        P = 1 - 0.5 * (1 + np.cos(2 * np.pi * np.sqrt(delta) / alpha))
    return np.clip(P, 0, 1)

def topological_tunneling(delta, xi, theta=0):
    """
    Topological insulator surface state tunneling (JAP 137, 124301).

    Parameters:
        delta : float
            Domain wall width (nm).
        xi : float
            Topological coherence length (nm).
        theta : float, optional
            Incidence angle (rad). Default is 0.

    Returns:
        float
            Tunneling probability T.
    """
    return 1 / (1 + np.sinh(delta / xi)**2 * np.cos(theta)**2)

def electron_transfer_rate(distance, beta=1.2):
    """
    Electron transfer rate vs distance (simplified model).

    Parameters:
        distance : float
            Donor-acceptor distance (Å).
        beta : float, optional
            Decay coefficient (Å⁻¹). Default is 1.2.

    Returns:
        float
            Relative electron transfer rate.
    """
    return np.exp(-beta * distance)

def exciton_transfer_simple(n_sites=7):
    """
    Simple rate equation model for exciton transfer (teaching only, no quantum coherence).

    Parameters:
        n_sites : int, optional
            Number of pigment molecules. Default is 7.

    Returns:
        ndarray
            Probability history array with shape (time_steps, n_sites).
    """
    prob = np.zeros(n_sites)
    prob[0] = 1.0
    transfer_rate = 0.3
    n_steps = 20
    prob_history = [prob.copy()]
    for _ in range(n_steps):
        new_prob = prob.copy()
        for i in range(n_sites):
            if prob[i] > 0:
                if i > 0:
                    new_prob[i-1] += prob[i] * transfer_rate
                if i < n_sites-1:
                    new_prob[i+1] += prob[i] * transfer_rate
                new_prob[i] -= prob[i] * 2 * transfer_rate
        prob = new_prob
        prob_history.append(prob.copy())
    return np.array(prob_history)

def transmission_wkb_quad(V_func, E, x1, x2):
    """
    WKB tunneling probability for arbitrary barrier using scipy.integrate.quad.

    Parameters:
        V_func : callable
            Potential function V(x) in eV, x in meters.
        E : float
            Particle energy (eV).
        x1, x2 : float
            Classical turning points (m).

    Returns:
        float
            Transmission coefficient.
    """
    try:
        from scipy.integrate import quad
    except ImportError:
        raise ImportError("scipy is required for this function. Install with: pip install scipy")
    E_J = E * e
    def integrand(x):
        Vx = V_func(x) * e
        if Vx <= E_J:
            return 0
        kappa = np.sqrt(2 * m_e * (Vx - E_J)) / hbar
        return kappa
    integral, _ = quad(integrand, x1, x2)
    return np.exp(-2 * integral)

def transmission_wkb_trapezoid(V_func, E, x1, x2, n=1000):
    """
    WKB tunneling probability for arbitrary barrier using trapezoidal rule (no scipy).

    Parameters:
        V_func : callable
            Potential function V(x) in eV, x in meters.
        E : float
            Particle energy (eV).
        x1, x2 : float
            Classical turning points (m).
        n : int, optional
            Number of integration points. Default is 1000.

    Returns:
        float
            Transmission coefficient.
    """
    E_J = E * e
    x = np.linspace(x1, x2, n)
    dx = (x2 - x1) / (n - 1)
    integral = 0
    for i, xi in enumerate(x):
        Vx = V_func(xi) * e
        if Vx <= E_J:
            kappa = 0
        else:
            kappa = np.sqrt(2 * m_e * (Vx - E_J)) / hbar
        weight = 0.5 if i == 0 or i == n-1 else 1.0
        integral += weight * kappa * dx
    return np.exp(-2 * integral)


# ============================================================================
# 2. Plotting Functions
# ============================================================================

def plot_transmission_rectangle():
    """Plot rectangular barrier tunneling probability."""
    V0, d = 3.0, 1.0
    energies = np.linspace(0.1, 5.0, 500)
    transmissions = [transmission_rectangle(V0, E, d) for E in energies]
    plt.figure(figsize=(8,5))
    plt.semilogy(energies, transmissions, 'b-', lw=2)
    plt.axvline(x=V0, color='r', linestyle='--', label=f'Barrier Height V0={V0}eV')
    plt.xlabel('Particle Energy (eV)')
    plt.ylabel('Tunneling Probability')
    plt.title('Rectangular Barrier Tunneling Probability')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()

def plot_stm_current():
    """Plot STM current vs tip-sample distance."""
    distances = np.linspace(0.3, 1.0, 100)
    currents = stm_current(distances)
    if np.any(currents <= 0):
        print("Warning: STM current contains non-positive values; adjusting for log plot.")
        currents = np.maximum(currents, 1e-20)
    plt.figure(figsize=(8,5))
    plt.semilogy(distances, currents, 'b-', lw=2)
    plt.xlabel('Tip-Sample Distance (nm)')
    plt.ylabel('Tunneling Current (arb. units)')
    plt.title('STM Current vs Distance')
    plt.grid(alpha=0.3)
    plt.show()

def plot_tfet_transfer():
    """Plot TFET transfer characteristics."""
    Vg = np.linspace(-0.2, 0.8, 200)
    Id_ss45 = tfet_transfer(Vg, SS=45)
    Id_ss60 = tfet_transfer(Vg, SS=60)
    plt.figure(figsize=(8,5))
    plt.semilogy(Vg, Id_ss45 * 1e6, 'r-', lw=2, label='SS=45 mV/dec (TFET)')
    plt.semilogy(Vg, Id_ss60 * 1e6, 'b--', lw=2, label='SS=60 mV/dec (MOSFET limit)')
    plt.xlabel('Gate Voltage Vg (V)')
    plt.ylabel('Drain Current Id (µA/µm)')
    plt.title('TFET Transfer Characteristics')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()

def plot_nonhermitian_lz():
    """Plot non-Hermitian Landau-Zener tunneling probability."""
    alpha, Omega = 0.5, 1.0
    gammas = np.linspace(0, 1.2, 500)
    probs = [nonhermitian_lz(alpha, Omega, g) for g in gammas]
    plt.figure(figsize=(8,5))
    plt.plot(gammas, probs, 'b-', lw=2)
    plt.axvline(x=0.46, color='r', linestyle='--', label=r'$\gamma_c = 0.46$ (experimental)')
    plt.xlabel(r'Dissipation rate $\gamma$')
    plt.ylabel('Tunneling Probability $P$')
    plt.title('Non-Hermitian Landau-Zener Tunneling (PRL 132, 156802)')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()

def plot_topological_tunneling():
    """Plot topological insulator surface state tunneling probability."""
    xi = 30.0
    delta_vals = np.linspace(0.1, 100, 500)
    T_vals = [topological_tunneling(d, xi) for d in delta_vals]
    plt.figure(figsize=(8,5))
    plt.semilogy(delta_vals, T_vals, 'g-', lw=2)
    plt.xlabel('Domain Wall Width δ (nm)')
    plt.ylabel('Tunneling Probability T')
    plt.title('Topological Insulator Surface State Tunneling (JAP 137, 124301)')
    plt.grid(alpha=0.3)
    plt.show()

def plot_electron_transfer_rate():
    """Plot electron transfer rate vs donor-acceptor distance."""
    distances = np.linspace(5, 20, 100)
    rates = electron_transfer_rate(distances)
    plt.figure(figsize=(8,5))
    plt.semilogy(distances, rates, 'b-', lw=2)
    plt.xlabel('Donor-Acceptor Distance (Å)')
    plt.ylabel('Relative Electron Transfer Rate')
    plt.title('Distance Dependence of Electron Transfer Rate')
    plt.grid(alpha=0.3)
    plt.show()

def plot_exciton_transfer():
    """Plot exciton transfer simplified model (no quantum coherence)."""
    prob_history = exciton_transfer_simple()
    plt.figure(figsize=(10,5))
    plt.imshow(prob_history.T, aspect='auto', cmap='hot', extent=[0, 20, 0, 7])
    plt.xlabel('Time Step')
    plt.ylabel('Pigment Molecule Index')
    plt.title('Exciton Transfer Simplified Model (No Quantum Coherence)')
    plt.colorbar(label='Probability')
    plt.show()

def demo_wkb():
    """Demonstrate WKB integration for a linear barrier."""
    def linear_barrier(x):
        if 0 <= x <= 1e-9:
            return 3.0 * (1 - x/1e-9)
        return 0
    prob_trap = transmission_wkb_trapezoid(linear_barrier, 1.5, 0, 1e-9, n=2000)
    print(f"Linear barrier tunneling probability (trapezoid): {prob_trap:.2e}")
    try:
        prob_quad = transmission_wkb_quad(linear_barrier, 1.5, 0, 1e-9)
        print(f"Linear barrier tunneling probability (quad): {prob_quad:.2e}")
    except ImportError:
        print("scipy not installed, skipping quad example.")


# ============================================================================
# 3. Science Fiction / Thought Experiment Functions (No Experimental Basis)
# ============================================================================

def fictional_barrier_collision(V0, E, d, spin='up'):
    """
    Fictional barrier collision probability (science fiction setting).
    No experimental basis; for thought experiments only.
    """
    print("WARNING: This function is a science fiction setting, no experimental basis.")
    spin_factor = 1.03 if spin == 'down' else 1.0
    prob = 0.03 * spin_factor
    return prob

def fictional_fusion_gain(T, catalyst='MoS2_WSe2'):
    """
    Fictional quantum catalytic fusion gain (science fiction setting).
    No experimental basis; for thought experiments only.
    """
    print("WARNING: This function is a science fiction setting, no experimental basis.")
    Q_base = 0.28 * T**2 / (T + 6.8)
    enhancement = 1.8 if catalyst == 'MoS2_WSe2' else 1.0
    Q = Q_base * enhancement / (1 + 0.1 * T**0.5)
    return min(Q, 3.2)

def fictional_interstellar_efficiency(distance_ly):
    """
    Fictional interstellar quantum communication efficiency (science fiction setting).
    No experimental basis; for thought experiments only.
    """
    print("WARNING: This function is a science fiction setting, no experimental basis.")
    return np.exp(-distance_ly / 10)

def fictional_lunar_storage_lifetime(radiation_level):
    """
    Fictional lunar data storage half-life (science fiction setting).
    No experimental basis; for thought experiments only.
    """
    print("WARNING: This function is a science fiction setting, no experimental basis.")
    damage_rate = 1e-6 * radiation_level
    protection = np.exp(-0.1 * 7)
    effective_damage = damage_rate * protection
    half_life = -np.log(0.5) / effective_damage
    return half_life

def fictional_wormhole_transfer():
    """
    Fictional wormhole information transfer (science fiction setting).
    No experimental basis; for thought experiments only.
    """
    print("WARNING: This function is a science fiction setting, no experimental basis.")
    return {"correlation": 0.95}

def fictional_neuro_interface(current, distance):
    """
    Fictional quantum neural interface risk score (science fiction setting).
    No experimental basis; for thought experiments only.
    """
    print("WARNING: This function is a science fiction setting, no experimental basis.")
    max_current = 0.005
    min_distance = 5.0
    risk = (current/max_current)**2 + (min_distance/distance)**2
    return min(risk, 1.0)

def fictional_compliance_index(capabilities):
    """
    Fictional national quantum compliance index (science fiction setting).
    No experimental basis; for thought experiments only.
    """
    print("WARNING: This function is a science fiction setting, no experimental basis.")
    weights = {'bio': 0.30, 'compute': 0.25, 'comm': 0.20, 'sensing': 0.15, 'energy': 0.10}
    score = sum(capabilities[field] * weights[field] for field in weights)
    return score * 100


# ============================================================================
# 4. Main Execution
# ============================================================================

if __name__ == "__main__":
    print("=== Quantum Tunneling Book Code Collection ===")
    print("1. Rectangular barrier tunneling probability")
    plot_transmission_rectangle()
    
    print("2. STM current vs distance")
    plot_stm_current()
    
    print("3. TFET transfer characteristics")
    plot_tfet_transfer()
    
    print("4. Non-Hermitian Landau-Zener tunneling")
    plot_nonhermitian_lz()
    
    print("5. Topological tunneling probability")
    plot_topological_tunneling()
    
    print("6. Electron transfer rate vs distance")
    plot_electron_transfer_rate()
    
    print("7. Exciton transfer simplified model")
    plot_exciton_transfer()
    
    print("8. WKB integration example")
    demo_wkb()
    
    print("\n=== Science Fiction Code Examples (Thought Experiments) ===")
    print("SF: Spin-up collision probability =", fictional_barrier_collision(5, 3, 0.8))
    print("SF: Spin-down collision probability =", fictional_barrier_collision(5, 3, 0.8, 'down'))
    print("SF: Fusion gain (T=150) =", fictional_fusion_gain(150))
    print("SF: Interstellar comm efficiency (4.24 ly) =", fictional_interstellar_efficiency(4.24))
    print("SF: Lunar data half-life =", fictional_lunar_storage_lifetime(2.5)/1e6, "million years")
    print("SF: Wormhole correlation =", fictional_wormhole_transfer()['correlation'])
    print("SF: Neural interface risk =", fictional_neuro_interface(0.002, 10))
    
    countries = {
        'Fictional Country A': {'bio': 0.85, 'compute': 0.92, 'comm': 0.88, 'sensing': 0.90, 'energy': 0.85},
        'Fictional Country B': {'bio': 0.88, 'compute': 0.95, 'comm': 0.82, 'sensing': 0.85, 'energy': 0.88}
    }
    print("\nFictional Compliance Index:")
    for country, caps in countries.items():
        print(f"{country}: {fictional_compliance_index(caps):.1f}")
    
    plt.show()