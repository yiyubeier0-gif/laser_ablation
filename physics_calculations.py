"""
physics_calculations.py

Core physics calculations for laser-material interaction
Includes: intensity, fluence, plasma density, Keldysh parameter, thermal diffusion, etc.
"""

import numpy as np
from scipy import constants
from typing import Tuple, Dict, Optional

# Physical constants
c = constants.c              # Speed of light [m/s]
h = constants.h              # Planck constant [J·s]
hbar = constants.hbar        # Reduced Planck constant [J·s]
e = constants.e              # Elementary charge [C]
m_e = constants.m_e          # Electron mass [kg]
epsilon_0 = constants.epsilon_0  # Vacuum permittivity [F/m]
k_B = constants.k            # Boltzmann constant [J/K]


# ═══════════════════════════════════════════════════════════════════════════
# INTENSITY AND FLUENCE CALCULATIONS
# ═══════════════════════════════════════════════════════════════════════════

def calculate_peak_intensity(pulse_energy: float, 
                            spot_diameter: float, 
                            pulse_duration: float,
                            beam_profile: str = 'gaussian') -> float:
    """
    Calculate peak intensity from pulse parameters
    
    Args:
        pulse_energy: Pulse energy [J]
        spot_diameter: Beam diameter (1/e² for Gaussian) [m]
        pulse_duration: Pulse duration (FWHM) [s]
        beam_profile: 'gaussian' or 'tophat'
    
    Returns:
        Peak intensity [W/cm²]
    """
    spot_radius = spot_diameter / 2  # m
    spot_area = np.pi * spot_radius**2  # m²
    
    if beam_profile == 'gaussian':
        # Gaussian: Peak intensity = 2 * Average intensity
        peak_power = pulse_energy / pulse_duration  # W
        peak_intensity = 2 * peak_power / spot_area  # W/m²
    elif beam_profile == 'tophat':
        # Top-hat: Uniform intensity
        peak_power = pulse_energy / pulse_duration
        peak_intensity = peak_power / spot_area
    else:
        raise ValueError(f"Unknown beam profile: {beam_profile}")
    
    # Convert to W/cm²
    return peak_intensity * 1e-4


def calculate_fluence(intensity: float, pulse_duration: float) -> float:
    """
    Calculate fluence from intensity and pulse duration
    
    Args:
        intensity: Peak intensity [W/cm²]
        pulse_duration: Pulse duration (FWHM) [s]
    
    Returns:
        Fluence [J/cm²]
    """
    return intensity * pulse_duration


def calculate_pulse_energy(intensity: float, 
                          spot_diameter: float, 
                          pulse_duration: float,
                          beam_profile: str = 'gaussian') -> float:
    """
    Calculate required pulse energy for given intensity
    
    Args:
        intensity: Desired peak intensity [W/cm²]
        spot_diameter: Beam diameter [m]
        pulse_duration: Pulse duration [s]
        beam_profile: 'gaussian' or 'tophat'
    
    Returns:
        Pulse energy [J]
    """
    spot_radius = spot_diameter / 2
    spot_area_cm2 = np.pi * (spot_radius * 100)**2  # cm²
    
    if beam_profile == 'gaussian':
        # For Gaussian: I_peak = 2·E / (π·w₀²·τ)
        # Therefore: E = I_peak·π·w₀²·τ / 2
        pulse_energy = intensity * np.pi * (spot_radius * 100)**2 * pulse_duration / 2
    else:
        pulse_energy = intensity * spot_area_cm2 * pulse_duration
    
    return pulse_energy


# ═══════════════════════════════════════════════════════════════════════════
# PLASMA PHYSICS
# ═══════════════════════════════════════════════════════════════════════════

def calculate_critical_density(wavelength: float) -> float:
    """
    Calculate critical plasma density
    
    Args:
        wavelength: Laser wavelength [m]
    
    Returns:
        Critical density [cm⁻³]
    """
    omega = 2 * np.pi * c / wavelength
    n_c = epsilon_0 * m_e * omega**2 / e**2
    return n_c * 1e-6  # Convert m⁻³ to cm⁻³


def calculate_keldysh_parameter(wavelength: float, 
                                intensity: float, 
                                ionization_potential: float) -> float:
    """
    Calculate Keldysh parameter
    
    γ = ω·√(2·m_e·Φ) / (e·E_field)
    
    Args:
        wavelength: Laser wavelength [m]
        intensity: Peak intensity [W/cm²]
        ionization_potential: Ionization potential [eV]
    
    Returns:
        Keldysh parameter (dimensionless)
    """
    omega = 2 * np.pi * c / wavelength
    E_field = np.sqrt(2 * intensity * 1e4 / (c * epsilon_0))  # Convert W/cm² to W/m²
    Phi = ionization_potential * e  # Convert eV to J
    
    gamma = omega * np.sqrt(2 * m_e * Phi) / (e * E_field)
    
    return gamma


def calculate_multiphoton_order(bandgap_eV: float, wavelength: float) -> int:
    """
    Calculate minimum number of photons required for ionization
    
    Args:
        bandgap_eV: Material bandgap [eV]
        wavelength: Laser wavelength [m]
    
    Returns:
        Multiphoton order (integer)
    """
    photon_energy_eV = (h * c / wavelength) / e
    n_photon = int(np.ceil(bandgap_eV / photon_energy_eV))
    return max(n_photon, 1)  # At least 1


def estimate_electron_density(intensity: float,
                              pulse_duration: float,
                              wavelength: float,
                              bandgap_eV: float,
                              mechanism: str = 'combined') -> float:
    """
    Estimate electron density evolution (simplified model)
    
    Args:
        intensity: Peak intensity [W/cm²]
        pulse_duration: Pulse duration [s]
        wavelength : Laser wavelength [m]
        bandgap_eV: Material bandgap [eV]
        mechanism: 'mpi' (multiphoton), 'avalanche', or 'combined'
    
    Returns:
        Estimated peak electron density [cm⁻³]
    """
    n_photon = calculate_multiphoton_order(bandgap_eV, wavelength)
    n_critical = calculate_critical_density(wavelength)
    
    if mechanism == 'mpi':
        # Multiphoton ionization rate (simplified)
        # W_MPI ∝ I^n (very approximate)
        sigma_n = 1e-50 * n_photon**(-3)  # Rough estimate [cm^(2n) s^-1]
        W_MPI = sigma_n * (intensity)**n_photon  # cm⁻³ s⁻¹
        n_e = W_MPI * pulse_duration
        
    elif mechanism == 'avalanche':
        # Avalanche ionization (exponential growth)
        # dn_e/dt = α·n_e, solution: n_e(t) = n_0·exp(α·t)
        E_field = np.sqrt(2 * intensity * 1e4 / (c * epsilon_0))
        nu_collision = 1e15  # Collision frequency [Hz] (rough estimate)
        alpha = (e * E_field) / (m_e * nu_collision * bandgap_eV * e)  # Growth rate
        
        n_seed = 1e10  # Background seed density [cm⁻³]
        n_e = n_seed * np.exp(alpha * pulse_duration)
        
    else:  # combined
        # Combined MPI + avalanche (sequential)
        # Phase 1: MPI seeds electrons
        sigma_n = 1e-50 * n_photon**(-3)
        W_MPI = sigma_n * (intensity)**n_photon
        t_seed = pulse_duration * 0.1  # First 10% of pulse
        n_seed = W_MPI * t_seed
        
        # Phase 2: Avalanche amplifies
        E_field = np.sqrt(2 * intensity * 1e4 / (c * epsilon_0))
        nu_collision = 1e15
        alpha = (e * E_field) / (m_e * nu_collision * bandgap_eV * e)
        t_avalanche = pulse_duration * 0.9
        n_e = n_seed * np.exp(alpha * t_avalanche)
    
    # Saturate at critical density
    return min(n_e, n_critical * 10)  # Allow some overshoot


# ═══════════════════════════════════════════════════════════════════════════
# THERMAL CALCULATIONS
# ═══════════════════════════════════════════════════════════════════════════

def calculate_thermal_diffusion_length(thermal_diffusivity: float, 
                                      time: float) -> float:
    """
    Calculate thermal diffusion length
    
    l_diff = √(D·t)
    
    Args:
        thermal_diffusivity: Thermal diffusivity [m²/s]
        time: Time scale [s]
    
    Returns:
        Diffusion length [m]
    """
    return np.sqrt(thermal_diffusivity * time)


def calculate_thermal_diffusion_time(thermal_diffusivity: float, 
                                    length_scale: float) -> float:
    """
    Calculate thermal diffusion time
    
    τ_thermal = l² / (4·D)
    
    Args:
        thermal_diffusivity: Thermal diffusivity [m²/s]
        length_scale: Characteristic length [m]
    
    Returns:
        Diffusion time [s]
    """
    return length_scale**2 / (4 * thermal_diffusivity)


def calculate_stress_confinement_time(sound_velocity: float, 
                                      length_scale: float) -> float:
    """
    Calculate stress confinement time
    
    τ_stress = l / c_sound
    
    Args:
        sound_velocity: Speed of sound in material [m/s]
        length_scale: Characteristic length [m]
    
    Returns:
        Stress confinement time [s]
    """
    return length_scale / sound_velocity


def estimate_temperature_rise(fluence: float,
                              absorption_coeff: float,
                              density: float,
                              specific_heat: float,
                              pulse_duration: float,
                              thermal_diffusivity: float) -> float:
    """
    Estimate peak temperature rise (simplified 1D model)
    
    Args:
        fluence: Laser fluence [J/cm²]
        absorption_coeff: Absorption coefficient [cm⁻¹]
        density: Material density [kg/m³]
        specific_heat: Specific heat [J/(kg·K)]
        pulse_duration: Pulse duration [s]
        thermal_diffusivity: Thermal diffusivity [m²/s]
    
    Returns:
        Temperature rise [K]
    """
    # Penetration depth
    l_abs = 1 / absorption_coeff * 1e-2  # Convert to m
    
    # Thermal diffusion length during pulse
    l_diff = calculate_thermal_diffusion_length(thermal_diffusivity, pulse_duration)
    
    # Effective heated volume depth
    l_eff = max(l_abs, l_diff)
    
    # Energy per unit volume
    energy_density = fluence * 1e4 / l_eff  # J/m³
    
    # Temperature rise (neglecting phase transitions)
    delta_T = energy_density / (density * specific_heat)
    
    return delta_T


def calculate_heat_accumulation(pulse_energy_absorbed: float,
                               rep_rate: float,
                               spot_size: float,
                               thermal_conductivity: float) -> float:
    """
    Estimate steady-state temperature rise from heat accumulation
    
    Args:
        pulse_energy_absorbed: Absorbed energy per pulse [J]
        rep_rate: Repetition rate [Hz]
        spot_size: Beam spot diameter [m]
        thermal_conductivity: Thermal conductivity [W/(m·K)]
    
    Returns:
        Steady-state temperature rise [K]
    """
    average_power = pulse_energy_absorbed * rep_rate  # W
    
    # Point source approximation
    # ΔT = P / (4π·κ·r)
    delta_T = average_power / (4 * np.pi * thermal_conductivity * (spot_size/2))
    
    return delta_T


# ═══════════════════════════════════════════════════════════════════════════
# ABLATION DEPTH MODELS
# ═══════════════════════════════════════════════════════════════════════════

def calculate_ablation_depth_single_pulse(fluence: float,
                                         threshold_fluence: float,
                                         absorption_coeff: float) -> float:
    """
    Calculate ablation depth per pulse (logarithmic model)
    
    d = (1/α) × ln(F / F_th)
    
    Args:
        fluence: Laser fluence [J/cm²]
        threshold_fluence: Ablation threshold fluence [J/cm²]
        absorption_coeff: Absorption coefficient [cm⁻¹]
    
    Returns:
        Ablation depth [μm]
    """
    if fluence <= threshold_fluence:
        return 0.0
    
    # Logarithmic regime
    if fluence < 5 * threshold_fluence:
        depth_cm = (1 / absorption_coeff) * np.log(fluence / threshold_fluence)
    else:
        # Saturation regime (plasma shielding)
        depth_max = (1 / absorption_coeff) * np.log(5)
        depth_cm = depth_max * (1 + 0.1 * np.log(fluence / (5 * threshold_fluence)))
    
    return depth_cm * 1e4  # Convert cm to μm


def calculate_ablation_depth_multipulse(fluence: float,
                                       threshold_fluence_single: float,
                                       absorption_coeff: float,
                                       num_pulses: int,
                                       incubation_parameter: float = 0.8) -> float:
    """
    Calculate total ablation depth with incubation effects
    
    Args:
        fluence: Fluence per pulse [J/cm²]
        threshold_fluence_single: Single-pulse threshold [J/cm²]
        absorption_coeff: Absorption coefficient [cm⁻¹]
        num_pulses: Number of pulses
        incubation_parameter: S parameter (0.7-0.95 typical)
    
    Returns:
        Total ablation depth [μm]
    """
    total_depth = 0.0
    
    for n in range(1, num_pulses + 1):
        # Incubation: F_th(N) = F_th(1) × N^(S-1)
        F_th_effective = threshold_fluence_single * n**(incubation_parameter - 1)
        
        # Ablation depth for this pulse
        depth = calculate_ablation_depth_single_pulse(
            fluence, F_th_effective, absorption_coeff
        )
        
        total_depth += depth
    
    return total_depth


# ═══════════════════════════════════════════════════════════════════════════
# NONLINEAR OPTICS
# ═══════════════════════════════════════════════════════════════════════════

def calculate_self_focusing_critical_power(wavelength: float,
                                          refractive_index: float,
                                          n2: float) -> float:
    """
    Calculate critical power for self-focusing
    
    P_cr = λ² / (2π·n₀·n₂)
    
    Args:
        wavelength: Laser wavelength [m]
        refractive_index: Linear refractive index
        n2: Nonlinear refractive index [m²/W]
    
    Returns:
        Critical power [W]
    """
    P_critical = wavelength**2 / (2 * np.pi * refractive_index * n2)
    return P_critical


def calculate_self_focusing_distance(beam_power: float,
                                    critical_power: float,
                                    rayleigh_range: float) -> Optional[float]:
    """
    Calculate self-focusing distance (Marburger formula)
    
    Args:
        beam_power: Beam power [W]
        critical_power: Critical power for self-focusing [W]
        rayleigh_range: Rayleigh range [m]
    
    Returns:
        Self-focusing distance [m], or None if no self-focusing
    """
    if beam_power <= critical_power:
        return None
    
    # Marburger formula
    z_sf = 0.367 * rayleigh_range / np.sqrt((beam_power / critical_power) - 0.852)
    
    return z_sf


def calculate_two_photon_absorption(intensity: float,
                                   beta: float) -> float:
    """
    Calculate two-photon absorption coefficient
    
    α_eff = α_linear + β·I
    
    Args:
        intensity: Laser intensity [W/cm²]
        beta: Two-photon absorption coefficient [cm/W]
    
    Returns:
        Effective absorption coefficient [cm⁻¹]
    """
    return beta * intensity


# ═══════════════════════════════════════════════════════════════════════════
# BEAM PROPAGATION
# ═══════════════════════════════════════════════════════════════════════════

def calculate_rayleigh_range(spot_size: float, wavelength: float) -> float:
    """
    Calculate Rayleigh range (depth of focus)
    
    z_R = π·w₀² / λ
    
    Args:
        spot_size: Beam waist diameter (1/e²) [m]
        wavelength: Wavelength [m]
    
    Returns:
        Rayleigh range [m]
    """
    w0 = spot_size / 2
    return np.pi * w0**2 / wavelength


def calculate_diffraction_limit(wavelength: float, numerical_aperture: float) -> float:
    """
    Calculate diffraction-limited spot size
    
    w₀ = λ / (π·NA)
    
    Args:
        wavelength: Wavelength [m]
        numerical_aperture: NA of focusing optic
    
    Returns:
        Minimum spot diameter [m]
    """
    return wavelength / (np.pi * numerical_aperture)


def calculate_numerical_aperture(focal_length: float, beam_diameter: float) -> float:
    """
    Calculate numerical aperture
    
    NA = sin(arctan(D/(2f))) ≈ D/(2f) for small angles
    
    Args:
        focal_length: Focal length of lens [m]
        beam_diameter: Beam diameter at lens [m]
    
    Returns:
        Numerical aperture
    """
    half_angle = np.arctan(beam_diameter / (2 * focal_length))
    return np.sin(half_angle)


def calculate_depth_of_focus(spot_size: float, wavelength: float) -> float:
    """
    Calculate depth of focus (DOF = 2·z_R)
    
    Args:
        spot_size: Beam waist diameter [m]
        wavelength: Wavelength [m]
    
    Returns:
        Depth of focus [m]
    """
    return 2 * calculate_rayleigh_range(spot_size, wavelength)


# ═══════════════════════════════════════════════════════════════════════════
# PULSE CHARACTERISTICS
# ═══════════════════════════════════════════════════════════════════════════

def calculate_pulse_overlap(scan_speed: float,
                           rep_rate: float,
                           spot_size: float) -> float:
    """
    Calculate pulse overlap during scanning
    
    Overlap = 1 - v / (f·d)
    
    Args:
        scan_speed: Scan speed [m/s]
        rep_rate: Repetition rate [Hz]
        spot_size: Spot diameter [m]
    
    Returns:
        Overlap fraction [0-1]
    """
    pulse_spacing = scan_speed / rep_rate
    overlap = 1 - (pulse_spacing / spot_size)
    return np.clip(overlap, 0, 1)


def calculate_effective_pulses(spot_size: float,
                              scan_speed: float,
                              rep_rate: float) -> float:
    """
    Calculate effective number of pulses hitting a point
    
    Args:
        spot_size: Spot diameter [m]
        scan_speed: Scan speed [m/s]
        rep_rate: Repetition rate [Hz]
    
    Returns:
        Effective number of pulses
    """
    dwell_time = spot_size / scan_speed
    n_pulses = dwell_time * rep_rate
    return n_pulses


def convert_fwhm_to_1e2(pulse_duration_fwhm: float,
                       pulse_shape: str = 'gaussian') -> float:
    """
    Convert FWHM pulse duration to 1/e² duration
    
    Args:
        pulse_duration_fwhm: FWHM pulse duration [s]
        pulse_shape: 'gaussian' or 'sech2'
    
    Returns:
        1/e² pulse duration [s]
    """
    if pulse_shape == 'gaussian':
        # For Gaussian: τ_1/e² = τ_FWHM / (2·√ln(2))
        return pulse_duration_fwhm / (2 * np.sqrt(np.log(2)))
    elif pulse_shape == 'sech2':
        # For sech²: τ_1/e² = τ_FWHM / 1.763
        return pulse_duration_fwhm / 1.763
    else:
        return pulse_duration_fwhm


# ═══════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

def photon_energy_eV(wavelength: float) -> float:
    """
    Calculate photon energy
    
    Args:
        wavelength: Wavelength [m]
    
    Returns:
        Photon energy [eV]
    """
    return (h * c / wavelength) / e


def wavelength_from_energy(energy_eV: float) -> float:
    """
    Calculate wavelength from photon energy
    
    Args:
        energy_eV: Photon energy [eV]
    
    Returns:
        Wavelength [m]
    """
    return h * c / (energy_eV * e)


def intensity_to_electric_field(intensity: float) -> float:
    """
    Convert intensity to electric field amplitude
    
    Args:
        intensity: Intensity [W/cm²]
    
    Returns:
        Electric field [V/m]
    """
    intensity_SI = intensity * 1e4  # Convert to W/m²
    E_field = np.sqrt(2 * intensity_SI / (c * epsilon_0))
    return E_field


def electric_field_to_intensity(E_field: float) -> float:
    """
    Convert electric field to intensity
    
    Args:
        E_field: Electric field amplitude [V/m]
    
    Returns:
        Intensity [W/cm²]
    """
    intensity_SI = 0.5 * c * epsilon_0 * E_field**2
    return intensity_SI * 1e-4  # Convert to W/cm²


# ═══════════════════════════════════════════════════════════════════════════
# COMPREHENSIVE PARAMETER CALCULATION
# ═══════════════════════════════════════════════════════════════════════════

def calculate_all_parameters(pulse_energy: float,
                            spot_diameter: float,
                            pulse_duration: float,
                            wavelength: float,
                            material_props: dict) -> Dict:
    """
    Calculate all relevant laser-material interaction parameters
    
    Args:
        pulse_energy: Pulse energy [J]
        spot_diameter: Spot diameter [m]
        pulse_duration: Pulse duration [s]
        wavelength: Wavelength [m]
        material_props: Dictionary of material properties
    
    Returns:
        Dictionary with all calculated parameters
    """
    results = {}
    
    # Basic laser parameters
    results['pulse_energy'] = pulse_energy
    results['spot_diameter'] = spot_diameter
    results['pulse_duration'] = pulse_duration
    results['wavelength'] = wavelength
    
    # Intensity and fluence
    results['peak_intensity'] = calculate_peak_intensity(
        pulse_energy, spot_diameter, pulse_duration
    )
    results['fluence'] = calculate_fluence(
        results['peak_intensity'], pulse_duration
    )
    
    # Photon properties
    results['photon_energy_eV'] = photon_energy_eV(wavelength)
    
    # Plasma parameters (if dielectric/semiconductor)
    if 'bandgap' in material_props and material_props['bandgap'] > 0:
        results['critical_density'] = calculate_critical_density(wavelength)
        results['keldysh_parameter'] = calculate_keldysh_parameter(
            wavelength, results['peak_intensity'], material_props['bandgap']
        )
        results['multiphoton_order'] = calculate_multiphoton_order(
            material_props['bandgap'], wavelength
        )
        results['estimated_electron_density'] = estimate_electron_density(
            results['peak_intensity'], pulse_duration, wavelength,
            material_props['bandgap']
        )
    
    # Thermal parameters
    if 'thermal_diffusivity' in material_props:
        D = material_props['thermal_diffusivity']
        feature_size = spot_diameter
        
        results['thermal_diffusion_time'] = calculate_thermal_diffusion_time(D, feature_size)
        results['thermal_diffusion_length'] = calculate_thermal_diffusion_length(D, pulse_duration)
        
        if 'sound_velocity' in material_props:
            results['stress_confinement_time'] = calculate_stress_confinement_time(
                material_props['sound_velocity'], feature_size
            )
        
        # Temperature rise estimate
        if all(k in material_props for k in ['absorption_coeff', 'density', 'specific_heat']):
            results['estimated_temp_rise'] = estimate_temperature_rise(
                results['fluence'],
                material_props['absorption_coeff'],
                material_props['density'],
                material_props['specific_heat'],
                pulse_duration,
                D
            )
    
    # Beam propagation
    results['rayleigh_range'] = calculate_rayleigh_range(spot_diameter, wavelength)
    results['depth_of_focus'] = calculate_depth_of_focus(spot_diameter, wavelength)
    
    # Nonlinear optics (if applicable)
    if 'n2' in material_props and material_props['n2'] is not None:
        results['self_focusing_critical_power'] = calculate_self_focusing_critical_power(
            wavelength, material_props['refractive_index'], material_props['n2']
        )
        beam_power = pulse_energy / pulse_duration
        if beam_power > results['self_focusing_critical_power']:
            results['self_focusing_distance'] = calculate_self_focusing_distance(
                beam_power, results['self_focusing_critical_power'],
                results['rayleigh_range']
            )
            results['self_focusing_warning'] = True
        else:
            results['self_focusing_warning'] = False
    
    return results


if __name__ == "__main__":
    # Test physics calculations
    print("="*70)
    print("PHYSICS CALCULATIONS TEST")
    print("="*70)
    
    # Test case: Femtosecond glass processing
    print("\nTest Case: Femtosecond Glass Processing")
    print("-"*70)
    
    pulse_energy = 100e-6  # 100 μJ
    spot_diameter = 20e-6  # 20 μm
    pulse_duration = 100e-15  # 100 fs
    wavelength = 800e-9  # 800 nm
    
    # Glass properties (simplified)
    material_props = {
        'bandgap': 9.0,  # eV
        'refractive_index': 1.45,
        'thermal_diffusivity': 8e-7,  # m²/s
        'absorption_coeff': 1e5,  # cm⁻¹ (effective, nonlinear)
        'density': 2200,  # kg/m³
        'specific_heat': 740,  # J/(kg·K)
        'sound_velocity': 5900,  # m/s
        'n2': 2.5e-20,  # m²/W
    }
    
    results = calculate_all_parameters(
        pulse_energy, spot_diameter, pulse_duration, wavelength, material_props
    )
    
    print(f"\nLaser Parameters:")
    print(f"  Pulse Energy:     {results['pulse_energy']*1e6:.1f} μJ")
    print(f"  Spot Diameter:    {results['spot_diameter']*1e6:.1f} μm")
    print(f"  Pulse Duration:   {results['pulse_duration']*1e15:.0f} fs")
    print(f"  Wavelength:       {results['wavelength']*1e9:.0f} nm")
    
    print(f"\nIntensity & Fluence:")
    print(f"  Peak Intensity:   {results['peak_intensity']:.2e} W/cm²")
    print(f"  Fluence:          {results['fluence']:.2f} J/cm²")
    
    print(f"\nPhoton Properties:")
    print(f"  Photon Energy:    {results['photon_energy_eV']:.2f} eV")
    print(f"  Multiphoton Order: {results['multiphoton_order']}")
    
    print(f"\nPlasma Parameters:")
    print(f"  Critical Density: {results['critical_density']:.2e} cm⁻³")
    print(f"  Keldysh Parameter: {results['keldysh_parameter']:.2f}")
    print(f"  Est. Electron Density: {results['estimated_electron_density']:.2e} cm⁻³")
    
    print(f"\nTemporal Hierarchy:")
    print(f"  τ_pulse:          {pulse_duration*1e15:.0f} fs")
    print(f"  τ_stress:         {results['stress_confinement_time']*1e9:.1f} ns")
    print(f"  τ_thermal:        {results['thermal_diffusion_time']*1e6:.1f} μs")
    print(f"  l_diff (thermal): {results['thermal_diffusion_length']*1e9:.1f} nm")
    
    print(f"\nTemperature Estimate:")
    print(f"  ΔT:               {results['estimated_temp_rise']:.0f} K")
    
    print(f"\nBeam Propagation:")
    print(f"  Rayleigh Range:   {results['rayleigh_range']*1e6:.1f} μm")
    print(f"  Depth of Focus:   {results['depth_of_focus']*1e6:.1f} μm")
    
    print(f"\nNonlinear Effects:")
    print(f"  P_critical (self-focus): {results['self_focusing_critical_power']*1e-3:.1f} kW")
    if results['self_focusing_warning']:
        print(f"  ⚠ Self-focusing at z = {results['self_focusing_distance']*1e3:.2f} mm")
    else:
        print(f"  ✓ No self-focusing")
    
    print("\n" + "="*70)
