## File 4: `regime_classification.py`

"""
regime_classification.py

Regime classification logic for laser-material interaction
Determines operating regime based on intensity, pulse duration, and material properties
"""

import numpy as np
from typing import Dict, Tuple
from enum import Enum

class ProcessingRegime(Enum):
    """Processing regime enumeration"""
    SUB_THRESHOLD = "Sub-Threshold"
    MODIFICATION = "Modification"
    PHOTOTHERMAL = "Photothermal"
    STRESS_CONFINED = "Stress-Confined"
    PLASMA_ATHERMAL = "Plasma-Athermal"
    DAMAGE = "Damage"
    CATASTROPHIC = "Catastrophic"


def classify_regime(intensity: float,
                   pulse_duration: float,
                   material_props: dict,
                   thresholds: dict) -> Dict:
    """
    Classify processing regime based on parameters
    
    Args:
        intensity: Peak intensity [W/cm²]
        pulse_duration: Pulse duration [s]
        material_props: Material properties dictionary
        thresholds: Thresholds dictionary
    
    Returns:
        Dictionary with regime classification and quality factor
    """
    result = {
        'regime': None,
        'regime_quality': 0.0,
        'mechanism': '',
        'characteristics': [],
        'warnings': []
    }
    
    # Extract key parameters
    tau_e_phonon = material_props.get('tau_e_phonon', 1e-12)
    sound_velocity = material_props.get('sound_velocity', 5000)
    thermal_diffusivity = material_props.get('thermal_diffusivity', 1e-6)
    
    # Calculate characteristic timescales
    feature_size = 10e-6  # 10 μm typical
    tau_stress = feature_size / sound_velocity
    tau_thermal = feature_size**2 / (4 * thermal_diffusivity)
    
    # ═══════════════════════════════════════════════════════════════════════
    # INTENSITY-BASED CLASSIFICATION
    # ═══════════════════════════════════════════════════════════════════════
    
    if intensity < thresholds['I_modification']:
        # SUB-THRESHOLD
        result['regime'] = ProcessingRegime.SUB_THRESHOLD
        result['regime_quality'] = 0.0
        result['mechanism'] = 'No observable effect'
        result['characteristics'] = [
            'Photon absorption without permanent change',
            'Energy dissipated non-destructively',
            'No processing effect'
        ]
        
    elif intensity < thresholds['I_ablation']:
        # MODIFICATION WITHOUT REMOVAL
        result['regime'] = ProcessingRegime.MODIFICATION
        result['regime_quality'] = 0.2
        result['mechanism'] = 'Structural modification without material removal'
        result['characteristics'] = [
            'Permanent structural/chemical change',
            'No material removed',
            'Defect creation, bond breaking'
        ]
        
    elif intensity < thresholds['I_damage']:
        # ABLATION REGIME (temporal-dependent sub-classification)
        
        if pulse_duration < tau_e_phonon:
            # PLASMA-MEDIATED ATHERMAL
            result['regime'] = ProcessingRegime.PLASMA_ATHERMAL
            result['regime_quality'] = 1.0
            result['mechanism'] = 'Nonlinear ionization → plasma → athermal removal'
            result['characteristics'] = [
                'Electrons hot, lattice cold',
                'Minimal HAZ (<100 nm)',
                'Sharp, clean edges',
                'No phase transitions'
            ]
            
        elif pulse_duration < tau_stress:
            # STRESS-CONFINED
            result['regime'] = ProcessingRegime.STRESS_CONFINED
            result['regime_quality'] = 0.7
            result['mechanism'] = 'Rapid heating → thermoelastic stress → spallation'
            result['characteristics'] = [
                'Mechanical removal (explosive)',
                'Reduced HAZ (1-10 μm)',
                'Solid-phase debris',
                'Mixed thermal-mechanical'
            ]
            
        elif pulse_duration < tau_thermal:
            # THERMAL-CONFINED
            result['regime'] = ProcessingRegime.PHOTOTHERMAL
            result['regime_quality'] = 0.5
            result['mechanism'] = 'Heat confinement → melting/evaporation'
            result['characteristics'] = [
                'Thermal removal',
                'Moderate HAZ (10-50 μm)',
                'Phase transitions occur',
                'Some thermal damage'
            ]
            
        else:
            # THERMAL-DIFFUSIVE
            result['regime'] = ProcessingRegime.PHOTOTHERMAL
            result['regime_quality'] = 0.3
            result['mechanism'] = 'Heat diffusion → thermal ablation'
            result['characteristics'] = [
                'Thermal removal',
                'Large HAZ (>50 μm)',
                'Significant thermal damage',
                'Melt flow, resolidification'
            ]
    
    elif intensity < thresholds['I_catastrophic']:
        # DAMAGE REGIME
        result['regime'] = ProcessingRegime.DAMAGE
        result['regime_quality'] = 0.3
        result['mechanism'] = 'Ablation with damage (over-processing)'
        result['characteristics'] = [
            'Material removed BUT damaged',
            'Plasma shielding or thermal damage',
            'Loss of precision',
            'Mechanical/thermal stress'
        ]
        result['warnings'].append("Intensity too high - reduce to damage threshold")
        
    else:
        # CATASTROPHIC
        result['regime'] = ProcessingRegime.CATASTROPHIC
        result['regime_quality'] = 0.05
        result['mechanism'] = 'Catastrophic breakdown'
        result['characteristics'] = [
            'Complete loss of control',
            'Optical breakdown',
            'Severe damage',
            'Unusable for processing'
        ]
        result['warnings'].append("CRITICAL: Catastrophic intensity - severe damage")
    
    # ═══════════════════════════════════════════════════════════════════════
    # TEMPORAL WARNINGS
    # ═══════════════════════════════════════════════════════════════════════
    
    if pulse_duration > thresholds['tau_critical']:
        result['warnings'].append(
            f"τ = {pulse_duration*1e12:.1f} ps > τ_critical = {thresholds['tau_critical']*1e12:.1f} ps"
        )
        result['warnings'].append("Damage likely - strongly recommend shorter pulses")
        result['regime_quality'] *= 0.5  # Severe penalty
    
    # Material-specific warnings
    material_name = material_props.get('name', 'Unknown')
    if 'Glass' in material_name or 'Sapphire' in material_name:
        if pulse_duration > 10e-12:
            result['warnings'].append(
                "CRITICAL FOR GLASS: τ > 10 ps → cracking inevitable"
            )
            result['regime_quality'] *= 0.3
    
    # ═══════════════════════════════════════════════════════════════════════
    # REGIME-SPECIFIC NOTES
    # ═══════════════════════════════════════════════════════════════════════
    
    result['timescale_info'] = {
        'tau_pulse': pulse_duration,
        'tau_e_phonon': tau_e_phonon,
        'tau_stress': tau_stress,
        'tau_thermal': tau_thermal,
        'pulse_vs_e_phonon': pulse_duration / tau_e_phonon,
        'pulse_vs_stress': pulse_duration / tau_stress,
        'pulse_vs_thermal': pulse_duration / tau_thermal,
    }
    
    return result


def determine_dominant_mechanism(intensity: float,
                                wavelength: float,
                                material_props: dict) -> str:
    """
    Determine dominant physical mechanism
    
    Args:
        intensity: Peak intensity [W/cm²]
        wavelength: Wavelength [m]
        material_props: Material properties
    
    Returns:
        String describing dominant mechanism
    """
    bandgap = material_props.get('bandgap', 0)
    
    if bandgap == 0:  # Metal
        if intensity > 1e12:
            return "Electronic excitation → plasma (ultrafast)"
        elif intensity > 1e8:
            return "Rapid heating → melting/evaporation"
        else:
            return "Conventional heating → thermal diffusion"
    
    else:  # Semiconductor or dielectric
        # Calculate photon energy
        from scipy import constants
        photon_energy_eV = (constants.h * constants.c / wavelength) / constants.e
        
        if photon_energy_eV > bandgap:
            # Linear absorption possible
            if intensity > 1e12:
                return "Linear + nonlinear absorption → plasma"
            else:
                return "Linear absorption → thermal"
        else:
            # Nonlinear absorption required
            if intensity > 1e12:
                return "Multiphoton + avalanche → plasma (unified)"
            elif intensity > 1e8:
                return "Weak nonlinear absorption → thermal"
            else:
                return "Negligible absorption"


def assess_processing_window_position(intensity: float,
                                      thresholds: dict) -> Dict:
    """
    Assess position within processing window
    
    Args:
        intensity: Peak intensity [W/cm²]
        thresholds: Thresholds dictionary
    
    Returns:
        Dictionary with window position analysis
    """
    result = {
        'in_clean_window': False,
        'position': '',
        'distance_to_optimal': 0.0,
        'recommendations': []
    }
    
    I_abl = thresholds['I_ablation']
    I_dam = thresholds['I_damage']
    I_opt = np.sqrt(I_abl * I_dam)  # Geometric mean
    
    if intensity < I_abl:
        result['position'] = 'Below ablation threshold'
        result['distance_to_optimal'] = (I_opt - intensity) / I_opt
        result['recommendations'].append(f"Increase intensity to > {I_abl:.2e} W/cm²")
        
    elif intensity < I_dam:
        result['in_clean_window'] = True
        
        # Calculate normalized position in window
        log_pos = (np.log10(intensity) - np.log10(I_abl)) / (np.log10(I_dam) - np.log10(I_abl))
        
        if 0.4 <= log_pos <= 0.6:
            result['position'] = 'Optimal (center of window)'
            result['distance_to_optimal'] = 0.0
        elif log_pos < 0.4:
            result['position'] = 'Lower part of window'
            result['distance_to_optimal'] = 0.5 - log_pos
            if log_pos < 0.2:
                result['recommendations'].append("Near ablation threshold - consider increasing intensity")
        else:
            result['position'] = 'Upper part of window'
            result['distance_to_optimal'] = log_pos - 0.5
            if log_pos > 0.8:
                result['recommendations'].append("Near damage threshold - consider reducing intensity")
    
    else:
        result['position'] = 'Above damage threshold'
        result['distance_to_optimal'] = (intensity - I_opt) / I_opt
        result['recommendations'].append(f"Reduce intensity to < {I_dam:.2e} W/cm²")
    
    return result


def predict_haz_extent(intensity: float,
                      pulse_duration: float,
                      material_props: dict,
                      thresholds: dict) -> Dict:
    """
    Predict heat-affected zone extent
    
    Args:
        intensity: Peak intensity [W/cm²]
        pulse_duration: Pulse duration [s]
        material_props: Material properties
        thresholds: Thresholds
    
    Returns:
        Dictionary with HAZ predictions
    """
    from physics_calculations import calculate_thermal_diffusion_length
    
    result = {
        'haz_estimate': 0.0,
        'haz_category': '',
        'quality_impact': ''
    }
    
    D = material_props.get('thermal_diffusivity', 1e-6)
    tau_e_phonon = material_props.get('tau_e_phonon', 1e-12)
    
    if pulse_duration < tau_e_phonon:
        # Athermal regime
        l_diff = calculate_thermal_diffusion_length(D, pulse_duration)
        result['haz_estimate'] = l_diff
        result['haz_category'] = 'Minimal (<100 nm)'
        result['quality_impact'] = 'Excellent - athermal processing'
        
    elif pulse_duration < 1e-9:
        # Picosecond regime
        l_diff = calculate_thermal_diffusion_length(D, pulse_duration)
        result['haz_estimate'] = l_diff * 5  # Empirical factor
        result['haz_category'] = 'Small (0.1-1 μm)'
        result['quality_impact'] = 'Good - minimal thermal effects'
        
    else:
        # Nanosecond and longer
        l_diff = calculate_thermal_diffusion_length(D, pulse_duration)
        
        # Additional factor for melting/evaporation
        if intensity > thresholds['I_ablation'] * 5:
            result['haz_estimate'] = l_diff * 20
            result['haz_category'] = 'Large (>10 μm)'
            result['quality_impact'] = 'Poor - significant thermal damage'
        else:
            result['haz_estimate'] = l_diff * 10
            result['haz_category'] = 'Moderate (1-10 μm)'
            result['quality_impact'] = 'Fair - some thermal effects'
    
    return result


def recommend_regime_optimization(current_regime: ProcessingRegime,
                                  target_quality: float,
                                  material_props: dict,
                                  thresholds: dict) -> list:
    """
    Recommend optimizations to reach target quality
    
    Args:
        current_regime: Current processing regime
        target_quality: Target quality (0-1)
        material_props: Material properties
        thresholds: Thresholds
    
    Returns:
        List of recommendations
    """
    recommendations = []
    
    if target_quality > 0.90:
        # High quality required
        if current_regime != ProcessingRegime.PLASMA_ATHERMAL:
            recommendations.append(
                "Target Q > 0.90 requires athermal regime"
            )
            recommendations.append(
                f"→ Use τ < {material_props.get('tau_e_phonon', 1e-12)*1e12:.1f} ps (femtosecond laser)"
            )
            tau_opt = thresholds.get('tau_optimal', (100e-15, 500e-15))
            recommendations.append(
                f"→ Optimal range: {tau_opt[0]*1e15:.0f}-{tau_opt[1]*1e15:.0f} fs"
            )
            
            I_window = thresholds.get('clean_window_fs')
            if I_window:
                I_opt = np.sqrt(I_window[0] * I_window[1])
                recommendations.append(
                    f"→ Intensity: {I_opt:.2e} W/cm² (center of clean window)"
                )
    
    elif target_quality > 0.75:
        # Good quality required
        if current_regime == ProcessingRegime.PHOTOTHERMAL:
            recommendations.append(
                "Target Q > 0.75 requires stress confinement or better"
            )
            recommendations.append(
                "→ Use picosecond pulses (τ = 1-10 ps)"
            )
            recommendations.append(
                "→ Or use femtosecond for higher quality"
            )
    
    elif target_quality > 0.60:
        # Moderate quality acceptable
        if current_regime == ProcessingRegime.SUB_THRESHOLD:
            recommendations.append(
                "No processing occurring - increase intensity"
            )
        elif current_regime == ProcessingRegime.MODIFICATION:
            recommendations.append(
                "Only modification, no removal - increase intensity above ablation threshold"
            )
    
    # Material-specific recommendations
    material_name = material_props.get('name', '')
    
    if 'Glass' in material_name or 'Sapphire' in material_name:
        if current_regime != ProcessingRegime.PLASMA_ATHERMAL:
            recommendations.append(
                "⚠ GLASS REQUIRES femtosecond pulses (τ < 10 ps mandatory)"
            )
            recommendations.append(
                "Longer pulses will cause cracking (thermal stress)"
            )
    
    if 'Polymer' in material_name or 'ABF' in material_name:
        recommendations.append(
            "For polymers: UV wavelength (355 nm) optimal"
        )
        recommendations.append(
            "Use N₂ atmosphere to prevent carbonization"
        )
    
    if material_props.get('category') == 'metal':
        recommendations.append(
            "For metals: Use green (532 nm) or UV (355 nm) wavelength"
        )
        recommendations.append(
            "Avoid NIR (high reflectivity)"
        )
    
    return recommendations


def classify_by_temporal_hierarchy(pulse_duration: float,
                                   material_props: dict,
                                   feature_size: float = 10e-6) -> Dict:
    """
    Classify regime purely by temporal hierarchy
    
    Args:
        pulse_duration: Pulse duration [s]
        material_props: Material properties
        feature_size: Characteristic feature size [m]
    
    Returns:
        Dictionary with temporal classification
    """
    tau_e_phonon = material_props.get('tau_e_phonon', 1e-12)
    sound_velocity = material_props.get('sound_velocity', 5000)
    thermal_diffusivity = material_props.get('thermal_diffusivity', 1e-6)
    
    tau_stress = feature_size / sound_velocity
    tau_thermal = feature_size**2 / (4 * thermal_diffusivity)
    
    result = {
        'temporal_regime': '',
        'confinement_type': '',
        'expected_quality_range': (0, 0),
        'timescales': {
            'tau_pulse': pulse_duration,
            'tau_e_phonon': tau_e_phonon,
            'tau_stress': tau_stress,
            'tau_thermal': tau_thermal,
        }
    }
    
    if pulse_duration < tau_e_phonon:
        result['temporal_regime'] = 'Ultrafast (< τ_e-phonon)'
        result['confinement_type'] = 'Electronic excitation only (athermal)'
        result['expected_quality_range'] = (0.85, 0.95)
        result['description'] = 'Electrons hot, lattice cold - optimal quality'
        
    elif pulse_duration < tau_stress:
        result['temporal_regime'] = 'Fast (τ_e-phonon < τ < τ_stress)'
        result['confinement_type'] = 'Stress confined'
        result['expected_quality_range'] = (0.65, 0.85)
        result['description'] = 'Mechanical removal, reduced HAZ'
        
    elif pulse_duration < tau_thermal:
        result['temporal_regime'] = 'Medium (τ_stress < τ < τ_thermal)'
        result['confinement_type'] = 'Thermally confined'
        result['expected_quality_range'] = (0.45, 0.65)
        result['description'] = 'Local heating, moderate HAZ'
        
    else:
        result['temporal_regime'] = 'Slow (τ > τ_thermal)'
        result['confinement_type'] = 'Thermal diffusion active'
        result['expected_quality_range'] = (0.30, 0.50)
        result['description'] = 'Heat spreads, large HAZ'
    
    return result


if __name__ == "__main__":
    # Test regime classification
    print("="*70)
    print("REGIME CLASSIFICATION TEST")
    print("="*70)
    
    # Mock material properties (Glass)
    material_props = {
        'name': 'Fused Silica',
        'category': 'dielectric',
        'tau_e_phonon': 1e-11,
        'sound_velocity': 5900,
        'thermal_diffusivity': 8e-7,
    }
    
    # Mock thresholds
    thresholds = {
        'I_modification': 1e12,
        'I_ablation': 1e13,
        'I_damage': 1e14,
        'I_catastrophic': 5e14,
        'tau_critical': 10e-12,
        'tau_optimal': (50e-15, 500e-15),
        'clean_window_fs': (3e13, 1e15),
    }
    
    # Test cases
    test_cases = [
        {'I': 1e11, 'tau': 100e-15, 'name': 'Sub-threshold'},
        {'I': 5e12, 'tau': 100e-15, 'name': 'Modification only'},
        {'I': 5e13, 'tau': 100e-15, 'name': 'Athermal plasma (optimal)'},
        {'I': 5e13, 'tau': 10e-12, 'name': 'At critical τ'},
        {'I': 5e13, 'tau': 20e-9, 'name': 'Nanosecond (thermal)'},
        {'I': 5e14, 'tau': 100e-15, 'name': 'Damage regime'},
    ]
    
    for i, case in enumerate(test_cases, 1):
        print(f"\n{'─'*70}")
        print(f"Test Case {i}: {case['name']}")
        print(f"  I = {case['I']:.2e} W/cm², τ = {case['tau']*1e12:.2f} ps")
        print(f"{'─'*70}")
        
        result = classify_regime(
            case['I'], case['tau'], material_props, thresholds
        )
        
        print(f"\nRegime: {result['regime'].value}")
        print(f"Quality: {result['regime_quality']:.2f}")
        print(f"Mechanism: {result['mechanism']}")
        
        print(f"\nCharacteristics:")
        for char in result['characteristics']:
            print(f"  • {char}")
        
        if result['warnings']:
            print(f"\n⚠ Warnings:")
            for warning in result['warnings']:
                print(f"  • {warning}")
        
        # Window position
        window_pos = assess_processing_window_position(case['I'], thresholds)
        print(f"\nWindow Position: {window_pos['position']}")
        if window_pos['recommendations']:
            print(f"Recommendations:")
            for rec in window_pos['recommendations']:
                print(f"  → {rec}")
    
    print("\n" + "="*70)
    print("TEMPORAL HIERARCHY ANALYSIS")
    print("="*70)
    
    test_durations = [100e-15, 1e-12, 10e-12, 1e-9, 10e-9]
    
    for tau in test_durations:
        result = classify_by_temporal_hierarchy(tau, material_props)
        print(f"\nτ = {tau*1e12:.2f} ps:")
        print(f"  Regime: {result['temporal_regime']}")
        print(f"  Confinement: {result['confinement_type']}")
        print(f"  Expected Q: {result['expected_quality_range'][0]:.2f}-{result['expected_quality_range'][1]:.2f}")
        print(f"  {result['description']}")
