"""
materials_database.py

Complete material properties database for laser-material interaction framework
"""

from dataclasses import dataclass, field
from typing import Dict, Optional

@dataclass
class MaterialProperties:
    """
    Comprehensive material property specification
    """
    # Identification
    name: str
    category: str  # 'metal', 'semiconductor', 'dielectric', 'polymer'
    
    # Electronic properties
    bandgap: float  # eV (0 for metals)
    refractive_index: float  # at 800 nm
    
    # Optical properties (wavelength-dependent)
    absorption_coefficient: Dict[str, float]  # {wavelength_str: alpha [cm^-1]}
    reflectivity: Dict[str, float] = field(default_factory=dict)  # {wavelength_str: R [0-1]}
    
    # Thermal properties
    thermal_diffusivity: float  # m²/s
    thermal_conductivity: float  # W/(m·K)
    specific_heat: float  # J/(kg·K)
    density: float  # kg/m³
    
    # Mechanical properties
    sound_velocity: float  # m/s
    youngs_modulus: Optional[float] = None  # Pa
    thermal_expansion_coeff: Optional[float] = None  # 1/K
    
    # Phase transition temperatures
    melting_temp: Optional[float] = None  # K
    boiling_temp: Optional[float] = None  # K
    decomposition_temp: Optional[float] = None  # K (for polymers)
    
    # Laser interaction parameters
    tau_e_phonon: float = 1e-12  # s (electron-phonon coupling time)
    incubation_parameter: float = 0.8  # S parameter (0.7-0.95 typical)
    
    # Nonlinear optics (if applicable)
    n2: Optional[float] = None  # m²/W (nonlinear refractive index)
    beta_TPA: Optional[float] = None  # m/W (two-photon absorption coefficient)
    multiphoton_cross_section: Optional[Dict[int, float]] = None  # {n_photon: sigma_n}
    
    # Notes
    notes: str = ""


# ═══════════════════════════════════════════════════════════════════════════
# MATERIAL DATABASE
# ═══════════════════════════════════════════════════════════════════════════

MATERIALS_DB = {}

# ───────────────────────────────────────────────────────────────────────────
# GLASS / DIELECTRICS
# ───────────────────────────────────────────────────────────────────────────

MATERIALS_DB['Glass'] = MaterialProperties(
    name='Fused Silica',
    category='dielectric',
    bandgap=9.0,
    refractive_index=1.45,
    absorption_coefficient={
        '193nm': 1e4,
        '248nm': 1e3,
        '355nm': 100.0,
        '532nm': 10.0,
        '800nm': 1.0,
        '1064nm': 0.1,
    },
    reflectivity={
        '800nm': 0.035,  # ~3.5% at normal incidence
        '1064nm': 0.035,
    },
    thermal_diffusivity=8e-7,
    thermal_conductivity=1.38,
    specific_heat=740,
    density=2200,
    sound_velocity=5900,
    youngs_modulus=73e9,
    thermal_expansion_coeff=0.55e-6,
    melting_temp=1983,
    tau_e_phonon=1e-11,
    incubation_parameter=0.85,
    n2=2.5e-20,
    notes="Requires femtosecond pulses for crack-free processing. Critical: τ < 10 ps"
)

MATERIALS_DB['BK7'] = MaterialProperties(
    name='BK7 Glass',
    category='dielectric',
    bandgap=4.0,
    refractive_index=1.51,
    absorption_coefficient={
        '355nm': 1e3,
        '532nm': 10.0,
        '800nm': 1.0,
        '1064nm': 0.1,
    },
    thermal_diffusivity=5.5e-7,
    thermal_conductivity=1.114,
    specific_heat=858,
    density=2510,
    sound_velocity=5100,
    melting_temp=830,
    tau_e_phonon=1e-11,
    incubation_parameter=0.80,
    n2=3.0e-20,
    notes="Lower quality than fused silica, lower melting point"
)

MATERIALS_DB['Sapphire'] = MaterialProperties(
    name='Sapphire (Al2O3)',
    category='dielectric',
    bandgap=9.9,
    refractive_index=1.76,
    absorption_coefficient={
        '355nm': 1e3,
        '532nm': 10.0,
        '800nm': 1.0,
        '1064nm': 0.1,
    },
    thermal_diffusivity=1.1e-5,
    thermal_conductivity=35,
    specific_heat=750,
    density=3980,
    sound_velocity =11000,
    youngs_modulus=400e9,
    thermal_expansion_coeff=5.0e-6,
    melting_temp=2323,
    tau_e_phonon=1e-11,
    incubation_parameter=0.90,
    n2=2.8e-20,
    notes="Very hard, high melting point. Requires very high intensity (>10^14 W/cm²)"
)

# ───────────────────────────────────────────────────────────────────────────
# POLYMERS
# ───────────────────────────────────────────────────────────────────────────

MATERIALS_DB['ABF_Polymer'] = MaterialProperties(
    name='ABF (Ajinomoto Build-up Film)',
    category='polymer',
    bandgap=4.5,
    refractive_index=1.55,
    absorption_coefficient={
        '193nm': 1e6,
        '248nm': 5e5,
        '266nm': 1e5,
        '355nm': 5e4,
        '532nm': 1e4,
        '800nm': 5e2,
        '1064nm': 1e2,
    },
    reflectivity={
        '355nm': 0.05,
        '532nm': 0.05,
    },
    thermal_diffusivity=1e-7,
    thermal_conductivity=0.2,
    specific_heat=1200,
    density=1100,
    sound_velocity=2000,
    decomposition_temp=623,  # ~350°C
    tau_e_phonon=1e-12,
    incubation_parameter=0.75,
    notes="Epoxy-based. UV optimal (355nm). Use N2 atmosphere to prevent carbonization."
)

MATERIALS_DB['Polyimide'] = MaterialProperties(
    name='Polyimide (Kapton)',
    category='polymer',
    bandgap=4.0,
    refractive_index=1.70,
    absorption_coefficient={
        '193nm': 1e6,
        '248nm': 5e5,
        '266nm': 1e5,
        '355nm': 8e4,
        '532nm': 1e4,
        '800nm': 1e3,
        '1064nm': 5e2,
    },
    thermal_diffusivity=1.5e-7,
    thermal_conductivity=0.12,
    specific_heat=1090,
    density=1420,
    sound_velocity=2200,
    decomposition_temp=823,  # ~550°C (higher than ABF)
    tau_e_phonon=1e-12,
    incubation_parameter=0.78,
    notes="Higher thermal stability than ABF. Aromatic structure."
)

MATERIALS_DB['PMMA'] = MaterialProperties(
    name='PMMA (Poly(methyl methacrylate))',
    category='polymer',
    bandgap=5.5,
    refractive_index=1.49,
    absorption_coefficient={
        '193nm': 5e5,
        '248nm': 1e5,
        '266nm': 5e4,
        '355nm': 1e4,
        '532nm': 1e3,
        '800nm': 1e2,
        '1064nm': 50,
    },
    thermal_diffusivity=1.2e-7,
    thermal_conductivity=0.19,
    specific_heat=1420,
    density=1180,
    sound_velocity=2690,
    decomposition_temp=433,  # ~160°C (low)
    tau_e_phonon=1e-12,
    incubation_parameter=0.70,
    notes="Low decomposition temperature. Photochemical ablation at UV."
)

# ───────────────────────────────────────────────────────────────────────────
# METALS
# ───────────────────────────────────────────────────────────────────────────

MATERIALS_DB['Copper'] = MaterialProperties(
    name='Copper',
    category='metal',
    bandgap=0.0,
    refractive_index=0.27,  # Complex refractive index (n for metals)
    absorption_coefficient={
        '193nm': 2e6,
        '248nm': 1.5e6,
        '355nm': 1.2e6,
        '532nm': 1e6,
        '800nm': 8e5,
        '1064nm': 7e5,
    },
    reflectivity={
        '355nm': 0.60,
        '532nm': 0.70,
        '800nm': 0.95,
        '1064nm': 0.98,
    },
    thermal_diffusivity=1.1e-4,
    thermal_conductivity=400,
    specific_heat=385,
    density=8960,
    sound_velocity=4700,
    youngs_modulus=130e9,
    thermal_expansion_coeff=16.5e-6,
    melting_temp=1358,
    boiling_temp=2835,
    tau_e_phonon=1e-12,
    incubation_parameter=0.90,
    notes="High reflectivity at NIR. Use green (532nm) or UV. τ < 1ps for melt-free."
)

MATERIALS_DB['Aluminum'] = MaterialProperties(
    name='Aluminum',
    category='metal',
    bandgap=0.0,
    refractive_index=1.44,
    absorption_coefficient={
        '193nm': 1.5e6,
        '248nm': 1e6,
        '355nm': 8e5,
        '532nm': 6e5,
        '800nm': 4e5,
        '1064nm': 3e5,
    },
    reflectivity={
        '355nm': 0.92,
        '532nm': 0.92,
        '800nm': 0.95,
        '1064nm': 0.97,
    },
    thermal_diffusivity=9.7e-5,
    thermal_conductivity=237,
    specific_heat=897,
    density=2700,
    sound_velocity=6420,
    youngs_modulus=70e9,
    thermal_expansion_coeff=23.1e-6,
    melting_temp=933,
    boiling_temp=2792,
    tau_e_phonon=1e-12,
    incubation_parameter=0.88,
    notes="Very high reflectivity. Rapid oxidation in air."
)

MATERIALS_DB['Stainless_Steel'] = MaterialProperties(
    name='Stainless Steel 316',
    category='metal',
    bandgap=0.0,
    refractive_index=2.50,
    absorption_coefficient={
        '355nm': 7e5,
        '532nm': 5e5,
        '800nm': 3e5,
        '1064nm': 2e5,
    },
    reflectivity={
        '355nm': 0.55,
        '532nm': 0.60,
        '800nm': 0.65,
        '1064nm': 0.70,
    },
    thermal_diffusivity=4.2e-6,
    thermal_conductivity=16.3,
    specific_heat=500,
    density=8000,
    sound_velocity=5790,
    youngs_modulus=193e9,
    thermal_expansion_coeff=17.3e-6,
    melting_temp=1673,
    boiling_temp=3023,
    tau_e_phonon=1.5e-12,
    incubation_parameter=0.85,
    notes="Lower thermal conductivity than Cu/Al. Better for laser processing."
)

MATERIALS_DB['Titanium'] = MaterialProperties(
    name='Titanium',
    category='metal',
    bandgap=0.0,
    refractive_index=2.90,
    absorption_coefficient={
        '355nm': 6e5,
        '532nm': 4e5,
        '800nm': 2.5e5,
        '1064nm': 2e5,
    },
    reflectivity={
        '355nm': 0.50,
        '532nm': 0.55,
        '800nm': 0.60,
        '1064nm': 0.65,
    },
    thermal_diffusivity=8.9e-6,
    thermal_conductivity=21.9,
    specific_heat=523,
    density=4506,
    sound_velocity=6070,
    youngs_modulus=116e9,
    thermal_expansion_coeff=8.6e-6,
    melting_temp=1941,
    boiling_temp=3560,
    tau_e_phonon=1.2e-12,
    incubation_parameter=0.87,
    notes="Biocompatible. Good for medical devices. Reactive (use inert atmosphere)."
)

MATERIALS_DB['Gold'] = MaterialProperties(
    name='Gold',
    category='metal',
    bandgap=0.0,
    refractive_index=0.47,
    absorption_coefficient={
        '355nm': 8e5,
        '532nm': 5e5,
        '800nm': 3e5,
        '1064nm': 2e5,
    },
    reflectivity={
        '355nm': 0.40,
        '532nm': 0.50,
        '800nm': 0.97,
        '1064nm': 0.98,
    },
    thermal_diffusivity=1.27e-4,
    thermal_conductivity=318,
    specific_heat=129,
    density=19300,
    sound_velocity=3240,
    youngs_modulus=78e9,
    thermal_expansion_coeff=14.2e-6,
    melting_temp=1337,
    boiling_temp=3129,
    tau_e_phonon=1e-12,
    incubation_parameter=0.92,
    notes="Very high reflectivity at NIR. Use UV/green. Inert (no oxidation)."
)

# ───────────────────────────────────────────────────────────────────────────
# SEMICONDUCTORS
# ───────────────────────────────────────────────────────────────────────────

MATERIALS_DB['Silicon'] = MaterialProperties(
    name='Silicon (crystalline)',
    category='semiconductor',
    bandgap=1.12,  # Indirect bandgap at 300K
    refractive_index=3.5,
    absorption_coefficient={
        '193nm': 1e6,
        '248nm': 5e5,
        '266nm': 1e6,
        '355nm': 1e6,
        '532nm': 1e5,
        '800nm': 1e3,  # Two-photon absorption
        '1064nm': 10,   # Very weak (free carrier absorption)
    },
    reflectivity={
        '532nm': 0.35,
        '800nm': 0.30,
        '1064nm': 0.30,
    },
    thermal_diffusivity=8.8e-5,
    thermal_conductivity=150,
    specific_heat=700,
    density=2329,
    sound_velocity=8400,
    youngs_modulus=170e9,
    thermal_expansion_coeff=2.6e-6,
    melting_temp=1687,
    tau_e_phonon=1e-12,
    incubation_parameter=0.85,
    beta_TPA=0.5e-11,  # m/W at 800 nm
    notes="Use λ<550nm for linear absorption. Amorphization risk with ns pulses."
)

MATERIALS_DB['GaAs'] = MaterialProperties(
    name='Gallium Arsenide',
    category='semiconductor',
    bandgap=1.42,  # Direct bandgap
    refractive_index=3.3,
    absorption_coefficient={
        '355nm': 1e6,
        '532nm': 5e5,
        '800nm': 1e4,
        '1064nm': 1e3,
    },
    reflectivity={
        '800nm': 0.32,
        '1064nm': 0.30,
    },
    thermal_diffusivity=2.3e-5,
    thermal_conductivity=55,
    specific_heat=330,
    density=5316,
    sound_velocity=5150,
    youngs_modulus=85e9,
    melting_temp=1511,
    tau_e_phonon=1e-12,
    incubation_parameter=0.83,
    notes="Direct bandgap. Brittle. Lower thermal conductivity than Si."
)

MATERIALS_DB['GaN'] = MaterialProperties(
    name='Gallium Nitride',
    category='semiconductor',
    bandgap=3.4,  # Wide bandgap
    refractive_index=2.3,
    absorption_coefficient={
        '193nm': 1e6,
        '248nm': 5e5,
        '266nm': 1e6,
        '355nm': 1e5,
        '532nm': 1e3,
        '800nm': 100,
        '1064nm': 10,
    },
    thermal_diffusivity=4.5e-5,
    thermal_conductivity=130,
    specific_heat=490,
    density=6150,
    sound_velocity=8000,
    youngs_modulus=295e9,
    melting_temp=2791,
    tau_e_phonon=1e-12,
    incubation_parameter=0.88,
    notes="Wide bandgap. Use UV for linear absorption. Very hard material."
)

MATERIALS_DB['SiC'] = MaterialProperties(
    name='Silicon Carbide',
    category='semiconductor',
    bandgap=3.26,  # 4H-SiC
    refractive_index=2.6,
    absorption_coefficient={
        '193nm': 1e6,
        '248nm': 5e5,
        '266nm': 1e6,
        '355nm': 1e5,
        '532nm': 1e3,
        '800nm': 100,
        '1064nm': 10,
    },
    thermal_diffusivity=1.1e-4,
    thermal_conductivity=370,
    specific_heat=750,
    density=3210,
    sound_velocity=11000,
    youngs_modulus=450e9,
    melting_temp=3003,
    tau_e_phonon=1e-12,
    incubation_parameter=0.90,
    notes="Very hard, high thermal conductivity. Requires high intensity."
)


def get_material(name: str) -> MaterialProperties:
    """
    Retrieve material properties by name
    
    Args:
        name: Material name (case-insensitive)
    
    Returns:
        MaterialProperties object
    
    Raises:
        KeyError: If material not found
    """
    # Normalize name
    for key in MATERIALS_DB.keys():
        if key.lower() == name.lower() or MATERIALS_DB[key].name.lower() == name.lower():
            return MATERIALS_DB[key]
    
    raise KeyError(f"Material '{name}' not found in database. Available: {list(MATERIALS_DB.keys())}")


def list_materials(category: Optional[str] = None) -> list:
    """
    List available materials, optionally filtered by category
    
    Args:
        category: Filter by category ('metal', 'semiconductor', 'dielectric', 'polymer')
    
    Returns:
        List of material names
    """
    if category is None:
        return list(MATERIALS_DB.keys())
    else:
        return [name for name, mat in MATERIALS_DB.items() if mat.category == category]


def get_absorption_coefficient(material: MaterialProperties, wavelength_nm: float) -> float:
    """
    Get absorption coefficient at specified wavelength (with interpolation)
    
    Args:
        material: MaterialProperties object
        wavelength_nm: Wavelength in nanometers
    
    Returns:
        Absorption coefficient in cm^-1
    """
    wavelength_str = f'{int(wavelength_nm)}nm'
    
    if wavelength_str in material.absorption_coefficient:
        return material.absorption_coefficient[wavelength_str]
    
    # Simple interpolation (log-log)
    wavelengths = sorted([int(w.replace('nm', '')) for w in material.absorption_coefficient.keys()])
    alphas = [material.absorption_coefficient[f'{w}nm'] for w in wavelengths]
    
    # Find nearest wavelengths
    import numpy as np
    
    if wavelength_nm <= wavelengths[0]:
        return alphas[0]
    if wavelength_nm >= wavelengths[-1]:
        return alphas[-1]
    
    # Log-log interpolation
    return float(np.exp(np.interp(
        np.log(wavelength_nm),
        np.log(wavelengths),
        np.log(alphas)
    )))


if __name__ == "__main__":
    # Test material database
    print("="*70)
    print("MATERIAL DATABASE TEST")
    print("="*70)
    
    print("\nAvailable materials by category:")
    for cat in ['metal', 'semiconductor', 'dielectric', 'polymer']:
        mats = list_materials(cat)
        print(f"\n{cat.upper()}: {len(mats)} materials")
        for mat in mats:
            print(f"  - {mat}")
    
    print("\n" + "="*70)
    print("SAMPLE MATERIAL: Glass")
    print("="*70)
    
    glass = get_material('Glass')
    print(f"Name: {glass.name}")
    print(f"Category: {glass.category}")
    print(f"Bandgap: {glass.bandgap} eV")
    print(f"Refractive index: {glass.refractive_index}")
    print(f"Thermal diffusivity: {glass.thermal_diffusivity:.2e} m²/s")
    print(f"τ_e-phonon: {glass.tau_e_phonon*1e12:.1f} ps")
    print(f"\nAbsorption coefficients:")
    for wl, alpha in glass.absorption_coefficient.items():
        print(f"  {wl}: {alpha:.1e} cm⁻¹")
    print(f"\nNotes: {glass.notes}")
    
    print("\n" + "="*70)
    print("INTERPOLATION TEST")
    print("="*70)
    
    test_wavelengths = [400, 600, 900]
    for wl in test_wavelengths:
        alpha = get_absorption_coefficient(glass, wl)
        print(f"α({wl} nm) = {alpha:.2e} cm⁻¹")
