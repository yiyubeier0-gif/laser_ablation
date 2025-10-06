"""
thresholds_database.py

Intensity threshold database for laser-material interaction framework
Defines I_modification, I_ablation, I_damage, I_catastrophic for each material
"""

from dataclasses import dataclass
from typing import Optional, Tuple, Dict

@dataclass
class ProcessingThresholds:
    """
    Material-specific intensity and temporal thresholds
    """
    # Material identification
    material_name: str
    
    # Intensity thresholds [W/cm²]
    I_modification: float      # Permanent change without removal
    I_ablation: float         # Material removal onset
    I_damage: float           # Damage onset (ablation WITH damage)
    I_catastrophic: float     # Catastrophic breakdown
    
    # Temporal constraints [s]
    tau_critical: float       # Above this: damage likely
    tau_optimal: Tuple[float, float]  # (min, max) optimal pulse duration range
    
    # Clean ablation windows by pulse duration regime [W/cm²]
    clean_window_fs: Optional[Tuple[float, float]] = None  # Femtosecond regime
    clean_window_ps: Optional[Tuple[float, float]] = None  # Picosecond regime
    clean_window_ns: Optional[Tuple[float, float]] = None  # Nanosecond regime
    
    # Additional constraints
    max_rep_rate: Optional[float] = None  # Hz (for heat accumulation concerns)
    recommended_wavelength: Optional[str] = None  # Optimal wavelength
    
    # Fluence thresholds [J/cm²] (calculated or measured)
    F_ablation_single_pulse: Optional[float] = None
    F_damage: Optional[float] = None
    
    # Notes
    notes: str = ""


# ═══════════════════════════════════════════════════════════════════════════
# THRESHOLDS DATABASE
# ═══════════════════════════════════════════════════════════════════════════

THRESHOLDS_DB: Dict[str, ProcessingThresholds] = {}

# ───────────────────────────────────────────────────────────────────────────
# GLASS / DIELECTRICS
# ───────────────────────────────────────────────────────────────────────────

THRESHOLDS_DB['Glass'] = ProcessingThresholds(
    material_name='Glass',
    
    # Intensity thresholds
    I_modification=1e12,      # Refractive index change (Type I/II modification)
    I_ablation=1e13,         # Plasma formation threshold
    I_damage=1e14,           # Cracking, shock damage onset
    I_catastrophic=5e14,     # Optical breakdown, explosive removal
    
    # Temporal constraints
    tau_critical=10e-12,     # HARD LIMIT: Above this → cracking
    tau_optimal=(50e-15, 500e-15),  # 50-500 fs optimal
    
    # Clean ablation windows
    clean_window_fs=(3e13, 1e15),   # VERY WIDE (>2 orders of magnitude)
    clean_window_ps=(1e13, 5e13),   # NARROW
    clean_window_ns=None,           # NO CLEAN WINDOW (cracking inevitable)
    
    # Additional constraints
    max_rep_rate=100e3,      # 100 kHz (higher → heat accumulation → cracking)
    recommended_wavelength='800nm',  # Ti:Sapphire standard
    
    # Fluence thresholds
    F_ablation_single_pulse=2.0,  # J/cm² (at 100 fs)
    F_damage=20.0,                # J/cm²
    
    notes="CRITICAL: τ < 10 ps MANDATORY for crack-free processing. Femtosecond strongly recommended."
)

THRESHOLDS_DB['BK7'] = ProcessingThresholds(
    material_name='BK7',
    I_modification=5e11,
    I_ablation=5e12,
    I_damage=5e13,
    I_catastrophic=2e14,
    tau_critical=10e-12,
    tau_optimal=(100e-15, 500e-15),
    clean_window_fs=(1e13, 5e13),
    clean_window_ps=(5e12, 1e13),
    clean_window_ns=None,
    max_rep_rate=100e3,
    recommended_wavelength='800nm',
    F_ablation_single_pulse=1.5,
    notes="Lower quality than fused silica. Lower thresholds due to lower melting point."
)

THRESHOLDS_DB['Sapphire'] = ProcessingThresholds(
    material_name='Sapphire',
    I_modification=5e12,
    I_ablation=5e13,
    I_damage=5e14,
    I_catastrophic=1e15,
    tau_critical=10e-12,
    tau_optimal=(50e-15, 300e-15),
    clean_window_fs=(1e14, 5e14),
    clean_window_ps=(5e13, 1e14),
    clean_window_ns=None,
    max_rep_rate=100e3,
    recommended_wavelength='800nm',
    F_ablation_single_pulse=5.0,
    notes="Very high thresholds due to high bandgap and thermal conductivity. Requires high pulse energy."
)

# ───────────────────────────────────────────────────────────────────────────
# POLYMERS
# ───────────────────────────────────────────────────────────────────────────

THRESHOLDS_DB['ABF_Polymer'] = ProcessingThresholds(
    material_name='ABF_Polymer',
    
    # Intensity thresholds
    I_modification=1e6,      # Color change, chemical modification
    I_ablation=5e6,         # Material removal onset
    I_damage=5e9,           # Carbonization, excessive thermal damage
    I_catastrophic=1e11,    # Complete destruction
    
    # Temporal constraints
    tau_critical=50e-9,     # Above this: significant thermal effects
    tau_optimal=(1e-9, 20e-9),  # 1-20 ns optimal for nanosecond systems
    
    # Clean ablation windows
    clean_window_fs=(1e8, 1e10),    # WIDE (femtosecond best quality)
    clean_window_ps=(5e7, 5e8),     # MODERATE
    clean_window_ns=(1e7, 5e7),     # NARROW (but acceptable for cost)
    
    # Additional constraints
    max_rep_rate=50e3,       # 50 kHz (heat accumulation at higher rates)
    recommended_wavelength='355nm',  # UV for photochemical ablation
    
    # Fluence thresholds
    F_ablation_single_pulse=0.1,  # J/cm² (at 10 ns, 355 nm)
    F_damage=10.0,                # J/cm²
    
    notes="Use UV (355nm) for best results. N2 atmosphere prevents carbonization. Nanosecond cost-effective."
)

THRESHOLDS_DB['Polyimide'] = ProcessingThresholds(
    material_name='Polyimide',
    I_modification=5e6,
    I_ablation=1e7,
    I_damage=1e10,
    I_catastrophic=5e11,
    tau_critical=50e-9,
    tau_optimal=(1e-9, 20e-9),
    clean_window_fs=(5e8, 5e10),
    clean_window_ps=(1e8, 1e9),
    clean_window_ns=(5e7, 1e8),
    max_rep_rate=50e3,
    recommended_wavelength='355nm',
    F_ablation_single_pulse=0.15,
    notes="Higher thermal stability than ABF. Aromatic structure more resistant to carbonization."
)

THRESHOLDS_DB['PMMA'] = ProcessingThresholds(
    material_name='PMMA',
    I_modification=5e5,
    I_ablation=1e6,
    I_damage=1e9,
    I_catastrophic=1e10,
    tau_critical=100e-9,
    tau_optimal=(5e-9, 50e-9),
    clean_window_fs=(5e7, 5e9),
    clean_window_ps=(1e7, 1e8),
    clean_window_ns=(5e6, 5e7),
    max_rep_rate=20e3,
    recommended_wavelength='248nm',  # Deep UV or 355nm
    F_ablation_single_pulse=0.05,
    notes="Low decomposition temperature. Photochemical ablation at UV optimal."
)

# ───────────────────────────────────────────────────────────────────────────
# METALS
# ───────────────────────────────────────────────────────────────────────────

THRESHOLDS_DB['Copper'] = ProcessingThresholds(
    material_name='Copper',
    
    # Intensity thresholds
    I_modification=1e7,      # Surface oxidation
    I_ablation=5e7,         # Melt expulsion (ns), plasma (fs)
    I_damage=5e13,          # Excessive plasma shielding, substrate damage
    I_catastrophic=1e15,    # Complete destruction
    
    # Temporal constraints
    tau_critical=100e-12,   # τ_e-phonon scale (melt-free boundary)
    tau_optimal=(100e-15, 10e-12),  # 100 fs - 10 ps
    
    # Clean ablation windows
    clean_window_fs=(1e12, 5e13),   # WIDE (plasma ablation, melt-free)
    clean_window_ps=(1e11, 1e12),   # MODERATE (some lattice heating)
    clean_window_ns=(5e7, 2e8),     # NARROW (thermal, HAZ present)
    
    # Additional constraints
    max_rep_rate=1e6,        # 1 MHz (metals have high thermal diffusivity)
    recommended_wavelength='532nm',  # Green (lower reflectivity than NIR)
    
    # Fluence thresholds
    F_ablation_single_pulse=0.5,  # J/cm² (at 10 ps, 532 nm)
    F_damage=50.0,                # J/cm²
    
    notes="Avoid NIR (high reflectivity). Use green (532nm) or UV (355nm). τ<1ps for melt-free."
)

THRESHOLDS_DB['Aluminum'] = ProcessingThresholds(
    material_name='Aluminum',
    I_modification=1e7,
    I_ablation=5e7,
    I_damage=5e13,
    I_catastrophic=1e15,
    tau_critical=100e-12,
    tau_optimal=(100e-15, 10e-12),
    clean_window_fs=(1e12, 5e13),
    clean_window_ps=(1e11, 1e12),
    clean_window_ns=(5e7, 2e8),
    max_rep_rate=1e6,
    recommended_wavelength='355nm',  # UV best (high reflectivity at all wavelengths)
    F_ablation_single_pulse=0.4,
    notes="Very high reflectivity. Rapid oxidation. Use UV and inert atmosphere."
)

THRESHOLDS_DB['Stainless_Steel'] = ProcessingThresholds(
    material_name='Stainless_Steel',
    I_modification=5e6,
    I_ablation=2e7,
    I_damage=1e13,
    I_catastrophic=5e14,
    tau_critical=200e-12,
    tau_optimal=(100e-15, 20e-12),
    clean_window_fs=(5e11, 1e13),
    clean_window_ps=(1e11, 5e11),
    clean_window_ns=(2e7, 1e8),
    max_rep_rate=1e6,
    recommended_wavelength='532nm',
    F_ablation_single_pulse=0.3,
    notes="Lower thermal conductivity than Cu/Al. Easier to process. Less reflective."
)

THRESHOLDS_DB['Titanium'] = ProcessingThresholds(
    material_name='Titanium',
    I_modification=5e6,
    I_ablation=2e7,
    I_damage=1e13,
    I_catastrophic=5e14,
    tau_critical=150e-12,
    tau_optimal=(100e-15, 15e-12),
    clean_window_fs=(5e11, 1e13),
    clean_window_ps=(1e11, 5e11),
    clean_window_ns=(2e7, 1e8),
    max_rep_rate=500e3,
    recommended_wavelength='532nm',
    F_ablation_single_pulse=0.35,
    notes="Biocompatible. Reactive (use inert atmosphere). Good for medical devices."
)

THRESHOLDS_DB['Gold'] = ProcessingThresholds(
    material_name='Gold',
    I_modification=1e7,
    I_ablation=5e7,
    I_damage=5e13,
    I_catastrophic=1e15,
    tau_critical=100e-12,
    tau_optimal=(100e-15, 10e-12),
    clean_window_fs=(1e12, 5e13),
    clean_window_ps=(1e11, 1e12),
    clean_window_ns=(5e7, 2e8),
    max_rep_rate=1e6,
    recommended_wavelength='532nm',  # Green (lower reflectivity)
    F_ablation_single_pulse=0.6,
    notes="Very high NIR reflectivity. Use green/UV. Inert (no oxidation). Biocompatible."
)

# ───────────────────────────────────────────────────────────────────────────
# SEMICONDUCTORS
# ───────────────────────────────────────────────────────────────────────────

THRESHOLDS_DB['Silicon'] = ProcessingThresholds(
    material_name='Silicon',
    
    # Intensity thresholds
    I_modification=1e6,      # Free carrier generation (transient)
    I_ablation=5e6,         # Surface melting (ns), plasma (fs)
    I_damage=5e12,          # Amorphization, thermal stress cracking
    I_catastrophic=1e14,    # Complete breakdown
    
    # Temporal constraints
    tau_critical=1e-9,      # Carrier recombination scale
    tau_optimal=(500e-15, 10e-12),  # 500 fs - 10 ps
    
    # Clean ablation windows
    clean_window_fs=(1e11, 5e12),   # WIDE (athermal)
    clean_window_ps=(1e10, 5e11),   # MODERATE
    clean_window_ns=(1e7, 5e7),     # NARROW (amorphization risk)
    
    # Additional constraints
    max_rep_rate=500e3,      # 500 kHz
    recommended_wavelength='532nm',  # Above bandgap (1.12 eV = 1107 nm)
    
    # Fluence thresholds
    F_ablation_single_pulse=0.2,  # J/cm² (at 10 ps, 532 nm)
    F_damage=10.0,                # J/cm²
    
    notes="Use λ<550nm for linear absorption. Amorphization risk with ns pulses. Check with Raman."
)

THRESHOLDS_DB['GaAs'] = ProcessingThresholds(
    material_name='GaAs',
    I_modification=5e5,
    I_ablation=5e6,
    I_damage=1e12,
    I_catastrophic=5e13,
    tau_critical=1e-9,
    tau_optimal=(500e-15, 10e-12),
    clean_window_fs=(5e10, 1e12),
    clean_window_ps=(1e10, 5e11),
    clean_window_ns=(5e6, 2e7),
    max_rep_rate=500e3,
    recommended_wavelength='800nm',  # Near bandgap (1.42 eV = 873 nm)
    F_ablation_single_pulse=0.15,
    notes="Direct bandgap. Brittle. Lower thermal conductivity than Si."
)

THRESHOLDS_DB['GaN'] = ProcessingThresholds(
    material_name='GaN',
    I_modification=1e7,
    I_ablation=5e7,
    I_damage=5e12,
    I_catastrophic=1e14,
    tau_critical=1e-9,
    tau_optimal=(500e-15, 10e-12),
    clean_window_fs=(1e11, 5e12),
    clean_window_ps=(5e10, 1e11),
    clean_window_ns=(5e7, 1e8),
    max_rep_rate=500e3,
    recommended_wavelength='355nm',  # Above wide bandgap (3.4 eV = 365 nm)
    F_ablation_single_pulse=0.5,
    notes="Wide bandgap. Use UV for linear absorption. Very hard material."
)

THRESHOLDS_DB['SiC'] = ProcessingThresholds(
    material_name='SiC',
    I_modification=1e7,
    I_ablation=1e8,
    I_damage=1e13,
    I_catastrophic=5e14,
    tau_critical=1e-9,
    tau_optimal=(500e-15, 10e-12),
    clean_window_fs=(5e11, 1e13),
    clean_window_ps=(1e11, 5e12),
    clean_window_ns=(1e8, 5e8),
    max_rep_rate=500e3,
    recommended_wavelength='355nm',
    F_ablation_single_pulse=1.0,
    notes="Very hard. High thermal conductivity. Requires high intensity."
)


def get_thresholds(material_name: str) -> ProcessingThresholds:
    """
    Retrieve processing thresholds by material name
    
    Args:
        material_name: Material name (case-insensitive)
    
    Returns:
        ProcessingThresholds object
    
    Raises:
        KeyError: If material not found
    """
    # Normalize name
    for key in THRESHOLDS_DB.keys():
        if key.lower() == material_name.lower():
            return THRESHOLDS_DB[key]
    
    raise KeyError(f"Thresholds for '{material_name}' not found. Available: {list(THRESHOLDS_DB.keys())}")


def get_clean_window(material_name: str, pulse_duration: float) -> Tuple[float, float]:
    """
    Get clean ablation window for specified material and pulse duration
    
    Args:
        material_name: Material name
        pulse_duration: Pulse duration in seconds
    
    Returns:
        Tuple (I_min, I_max) in W/cm², or (None, None) if no clean window exists
    """
    thresh = get_thresholds(material_name)
    
    # Classify pulse duration regime
    if pulse_duration < 1e-12:  # Femtosecond
        return thresh.clean_window_fs if thresh.clean_window_fs else (None, None)
    elif pulse_duration < 100e-12:  # Picosecond
        return thresh.clean_window_ps if thresh.clean_window_ps else (None, None)
    else:  # Nanosecond
        return thresh.clean_window_ns if thresh.clean_window_ns else (None, None)


def check_temporal_constraint(material_name: str, pulse_duration: float) -> dict:
    """
    Check if pulse duration satisfies material constraints
    
    Args:
        material_name: Material name
        pulse_duration: Pulse duration in seconds
    
    Returns:
        Dictionary with constraint status and warnings
    """
    thresh = get_thresholds(material_name)
    
    result = {
        'within_optimal': False,
        'within_critical': False,
        'warnings': []
    }
    
    # Check critical constraint
    if pulse_duration > thresh.tau_critical:
        result['warnings'].append(
            f"CRITICAL: τ = {pulse_duration*1e12:.1f} ps > τ_critical = {thresh.tau_critical*1e12:.1f} ps"
        )
        result['warnings'].append("Damage likely. Strongly recommend shorter pulses.")
        result['within_critical'] = False
    else:
        result['within_critical'] = True
    
    # Check optimal range
    tau_min, tau_max = thresh.tau_optimal
    if tau_min <= pulse_duration <= tau_max:
        result['within_optimal'] = True
    else:
        if pulse_duration < tau_min:
            result['warnings'].append(
                f"Below optimal range. τ = {pulse_duration*1e15:.0f} fs < {tau_min*1e15:.0f} fs"
            )
            result['warnings'].append("Very short pulse - ensure sufficient energy delivery.")
        else:
            result['warnings'].append(
                f"Above optimal range. τ = {pulse_duration*1e12:.1f} ps > {tau_max*1e12:.1f} ps"
            )
            result['warnings'].append("Quality may be reduced. Consider shorter pulses.")
    
    return result


if __name__ == "__main__":
    # Test thresholds database
    print("="*70)
    print("THRESHOLDS DATABASE TEST")
    print("="*70)
    
    print("\nAvailable materials:")
    for mat_name in THRESHOLDS_DB.keys():
        print(f"  - {mat_name}")
    
    print("\n" + "="*70)
    print("SAMPLE THRESHOLDS: Glass")
    print("="*70)
    
    thresh = get_thresholds('Glass')
    print(f"Material: {thresh.material_name}")
    print(f"\nIntensity Thresholds:")
    print(f"  I_modification:   {thresh.I_modification:.2e} W/cm²")
    print(f"  I_ablation:       {thresh.I_ablation:.2e} W/cm²")
    print(f"  I_damage:         {thresh.I_damage:.2e} W/cm²")
    print(f"  I_catastrophic:   {thresh.I_catastrophic:.2e} W/cm²")
    
    print(f"\nTemporal Constraints:")
    print(f"  τ_critical:       {thresh.tau_critical*1e12:.1f} ps")
    print(f"  τ_optimal:        {thresh.tau_optimal[0]*1e15:.0f}-{thresh.tau_optimal[1]*1e15:.0f} fs")
    
    print(f"\nClean Ablation Windows:")
    if thresh.clean_window_fs:
        print(f"  Femtosecond:      {thresh.clean_window_fs[0]:.1e} - {thresh.clean_window_fs[1]:.1e} W/cm²")
    if thresh.clean_window_ps:
        print(f"  Picosecond:       {thresh.clean_window_ps[0]:.1e} - {thresh.clean_window_ps[1]:.1e} W/cm²")
    if thresh.clean_window_ns:
        print(f"  Nanosecond:       {thresh.clean_window_ns[0]:.1e} - {thresh.clean_window_ns[1]:.1e} W/cm²")
    else:
        print(f"  Nanosecond:       NO CLEAN WINDOW")
    
    print(f"\nAdditional Constraints:")
    print(f"  Max rep rate:     {thresh.max_rep_rate/1e3:.0f} kHz")
    print(f"  Recommended λ:    {thresh.recommended_wavelength}")
    
    print(f"\nNotes: {thresh.notes}")
    
    print("\n" + "="*70)
    print("TEMPORAL CONSTRAINT CHECK")
    print("="*70)
    
    test_pulses = [100e-15, 10e-12, 50e-12, 100e-9]
    for tau in test_pulses:
        print(f"\nτ = {tau*1e12:.1f} ps:")
        result = check_temporal_constraint('Glass', tau)
        print(f"  Within optimal: {result['within_optimal']}")
        print(f"  Within critical: {result['within_critical']}")
        for warning in result['warnings']:
            print(f"  ⚠ {warning}")