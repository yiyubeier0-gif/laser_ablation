"""
quality_assessment.py

Quality assessment and scoring system for laser processing
Calculates comprehensive quality factor Q and provides detailed breakdown
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from regime_classification import classify_regime, assess_processing_window_position
from physics_calculations import (
    calculate_thermal_diffusion_length,
    calculate_keldysh_parameter,
    calculate_multiphoton_order,
    photon_energy_eV
)


def calculate_quality_total(wavelength: float,
                           pulse_duration: float,
                           intensity: float,
                           material_name: str,
                           material_props: dict,
                           thresholds: dict,
                           feature_size: float = 10e-6,
                           num_pulses: int = 1,
                           rep_rate: float = 1000) -> Dict:
    """
    Calculate comprehensive quality factor Q
    
    Args:
        wavelength: Laser wavelength [m]
        pulse_duration: Pulse duration [s]
        intensity: Peak intensity [W/cm²]
        material_name: Material name
        material_props: Material properties dictionary
        thresholds: Thresholds dictionary
        feature_size: Feature size [m]
        num_pulses: Number of pulses
        rep_rate: Repetition rate [Hz]
    
    Returns:
        Dictionary with complete quality assessment
    """
    result = {
        'Q_total': 0.0,
        'components': {},
        'regime': '',
        'warnings': [],
        'recommendations': [],
        'details': {}
    }
    
    # ═══════════════════════════════════════════════════════════════════════
    # COMPONENT 1: REGIME QUALITY (Q_regime)
    # ═══════════════════════════════════════════════════════════════════════
    
    regime_result = classify_regime(
        intensity, pulse_duration, material_props, thresholds
    )
    
    result['regime'] = regime_result['regime'].value
    result['components']['Q_regime'] = regime_result['regime_quality']
    result['warnings'].extend(regime_result['warnings'])
    result['details']['regime_info'] = regime_result
    
    # ═══════════════════════════════════════════════════════════════════════
    # COMPONENT 2: INTENSITY WINDOW QUALITY (Q_intensity)
    # ═══════════════════════════════════════════════════════════════════════
    
    Q_intensity = calculate_intensity_quality(
        intensity, thresholds, pulse_duration
    )
    
    result['components']['Q_intensity'] = Q_intensity['quality']
    result['warnings'].extend(Q_intensity['warnings'])
    result['details']['intensity_info'] = Q_intensity
    
    # ═══════════════════════════════════════════════════════════════════════
    # COMPONENT 3: TEMPORAL OPTIMIZATION (Q_temporal)
    # ═══════════════════════════════════════════════════════════════════════
    
    Q_temporal = calculate_temporal_quality(
        pulse_duration, material_props, thresholds, feature_size
    )
    
    result['components']['Q_temporal'] = Q_temporal['quality']
    result['warnings'].extend(Q_temporal['warnings'])
    result['details']['temporal_info'] = Q_temporal
    
    # ═══════════════════════════════════════════════════════════════════════
    # COMPONENT 4: WAVELENGTH OPTIMIZATION (Q_wavelength)
    # ═══════════════════════════════════════════════════════════════════════
    
    Q_wavelength = calculate_wavelength_quality(
        wavelength, intensity, material_props, material_name, result['regime']
    )
    
    result['components']['Q_wavelength'] = Q_wavelength['quality']
    result['warnings'].extend(Q_wavelength['warnings'])
    result['details']['wavelength_info'] = Q_wavelength
    
    # ═══════════════════════════════════════════════════════════════════════
    # COMPONENT 5: MULTI-PULSE ACCUMULATION (Q_accumulation)
    # ═══════════════════════════════════════════════════════════════════════
    
    if num_pulses > 1:
        Q_accumulation = calculate_accumulation_quality(
            intensity, pulse_duration, num_pulses, rep_rate,
            material_props, feature_size
        )
        result['components']['Q_accumulation'] = Q_accumulation['quality']
        result['warnings'].extend(Q_accumulation['warnings'])
        result['details']['accumulation_info'] = Q_accumulation
    else:
        result['components']['Q_accumulation'] = 1.0
        result['details']['accumulation_info'] = {'note': 'Single pulse - no accumulation'}
    
    # ═══════════════════════════════════════════════════════════════════════
    # CALCULATE TOTAL QUALITY (Multiplicative)
    # ═══════════════════════════════════════════════════════════════════════
    
    Q_total = 1.0
    for component, value in result['components'].items():
        Q_total *= value
    
    result['Q_total'] = Q_total
    
    # ═══════════════════════════════════════════════════════════════════════
    # GENERATE RECOMMENDATIONS
    # ═══════════════════════════════════════════════════════════════════════
    
    result['recommendations'] = generate_recommendations(
        result, material_name, wavelength, pulse_duration, intensity, thresholds
    )
    
    # ═══════════════════════════════════════════════════════════════════════
    # QUALITY INTERPRETATION
    # ═══════════════════════════════════════════════════════════════════════
    
    result['interpretation'] = interpret_quality(Q_total)
    
    return result


def calculate_intensity_quality(intensity: float,
                                thresholds: dict,
                                pulse_duration: float) -> Dict:
    """
    Calculate quality factor based on intensity position in window
    
    Args:
        intensity: Peak intensity [W/cm²]
        thresholds: Thresholds dictionary
        pulse_duration: Pulse duration [s]
    
    Returns:
        Dictionary with intensity quality assessment
    """
    result = {
        'quality': 0.0,
        'position': '',
        'warnings': [],
        'optimal_intensity': 0.0
    }
    
    I_mod = thresholds['I_modification']
    I_abl = thresholds['I_ablation']
    I_dam = thresholds['I_damage']
    I_cat = thresholds['I_catastrophic']
    
    # Calculate optimal intensity (geometric mean)
    I_opt = np.sqrt(I_abl * I_dam)
    result['optimal_intensity'] = I_opt
    
    if intensity < I_mod:
        # Sub-threshold
        result['quality'] = 0.0
        result['position'] = 'Sub-threshold'
        result['warnings'].append("Below modification threshold - no effect")
        
    elif intensity < I_abl:
        # Modification only
        result['quality'] = 0.2
        result['position'] = 'Modification only (no removal)'
        result['warnings'].append("No material removal - increase intensity")
        
    elif intensity < I_dam:
        # CLEAN ABLATION WINDOW
        result['position'] = 'Within clean ablation window'
        
        # Gaussian profile centered at I_opt
        sigma = (np.log10(I_dam) - np.log10(I_abl)) / 4
        Q = np.exp(-((np.log10(intensity) - np.log10(I_opt)) / sigma)**2)
        result['quality'] = Q
        
        # Position-specific warnings
        if intensity < I_abl * 1.5:
            result['warnings'].append("Near ablation threshold - incomplete removal possible")
            result['quality'] *= 0.95
        elif intensity > I_dam * 0.7:
            result['warnings'].append("Approaching damage threshold - reduce intensity")
            result['quality'] *= 0.95
        
    elif intensity < I_cat:
        # DAMAGE REGIME
        result['position'] = 'Damage regime'
        # Exponential decay from damage threshold
        decay_factor = np.exp(-(np.log10(intensity) - np.log10(I_dam)) / 0.5)
        result['quality'] = 0.3 * decay_factor
        result['warnings'].append("DAMAGE REGIME - thermal/mechanical damage present")
        result['warnings'].append(f"Reduce intensity to < {I_dam:.2e} W/cm²")
        
    else:
        # CATASTROPHIC
        result['position'] = 'Catastrophic'
        result['quality'] = 0.05
        result['warnings'].append("CATASTROPHIC - severe damage, process unusable")
        result['warnings'].append("CRITICAL: Reduce intensity immediately")
    
    return result


def calculate_temporal_quality(pulse_duration: float,
                               material_props: dict,
                               thresholds: dict,
                               feature_size: float) -> Dict:
    """
    Calculate quality factor based on temporal confinement
    
    Args:
        pulse_duration: Pulse duration [s]
        material_props: Material properties
        thresholds: Thresholds
        feature_size: Feature size [m]
    
    Returns:
        Dictionary with temporal quality assessment
    """
    result = {
        'quality': 0.0,
        'confinement_type': '',
        'warnings': [],
        'timescales': {}
    }
    
    # Extract timescales
    tau_e_phonon = material_props.get('tau_e_phonon', 1e-12)
    sound_velocity = material_props.get('sound_velocity', 5000)
    thermal_diffusivity = material_props.get('thermal_diffusivity', 1e-6)
    
    tau_stress = feature_size / sound_velocity
    tau_thermal = feature_size**2 / (4 * thermal_diffusivity)
    
    result['timescales'] = {
        'tau_pulse': pulse_duration,
        'tau_e_phonon': tau_e_phonon,
        'tau_stress': tau_stress,
        'tau_thermal': tau_thermal
    }
    
    # Classify and assign quality
    if pulse_duration < tau_e_phonon:
        # ATHERMAL - Optimal
        result['quality'] = 1.0
        result['confinement_type'] = 'Athermal (electronic only)'
        
        if pulse_duration < 10e-15:
            result['warnings'].append("Very short pulse - ensure sufficient energy delivery")
            result['quality'] *= 0.98
            
    elif pulse_duration < tau_stress:
        # STRESS-CONFINED - Good
        result['quality'] = 0.80
        result['confinement_type'] = 'Stress-confined'
        
    elif pulse_duration < tau_thermal:
        # THERMAL-CONFINED - Moderate
        result['quality'] = 0.50
        result['confinement_type'] = 'Thermally-confined'
        result['warnings'].append("Thermal confinement only - some HAZ expected")
        
    else:
        # THERMAL-DIFFUSIVE - Poor
        result['quality'] = 0.30
        result['confinement_type'] = 'Thermal diffusion active'
        result['warnings'].append("Thermal diffusion active - significant HAZ")
    
    # Check critical constraint
    tau_critical = thresholds.get('tau_critical', 1e-9)
    if pulse_duration > tau_critical:
        result['quality'] *= 0.3  # Severe penalty
        result['warnings'].append(
            f"CRITICAL: τ = {pulse_duration*1e12:.1f} ps > τ_critical = {tau_critical*1e12:.1f} ps"
        )
        result['warnings'].append("Damage highly likely - STRONGLY recommend shorter pulses")
    
    # Check optimal range
    if 'tau_optimal' in thresholds:
        tau_min, tau_max = thresholds['tau_optimal']
        
        if tau_min <= pulse_duration <= tau_max:
            result['quality'] = max(result['quality'], 0.95)  # Boost if in optimal range
        elif pulse_duration < tau_min:
            result['warnings'].append(
                f"Below optimal range ({tau_min*1e15:.0f} fs)"
            )
        elif pulse_duration > tau_max:
            result['warnings'].append(
                f"Above optimal range ({tau_max*1e12:.1f} ps)"
            )
    
    # Material-specific checks
    material_name = material_props.get('name', '')
    if 'Glass' in material_name or 'Sapphire' in material_name:
        if pulse_duration > 10e-12:
            result['quality'] *= 0.2
            result['warnings'].append(
                "CRITICAL FOR GLASS: τ > 10 ps → cracking inevitable!"
            )
    
    return result


def calculate_wavelength_quality(wavelength: float,
                                 intensity: float,
                                 material_props: dict,
                                 material_name: str,
                                 regime: str) -> Dict:
    """
    Calculate quality factor based on wavelength optimization
    
    Args:
        wavelength: Wavelength [m]
        intensity: Intensity [W/cm²]
        material_props: Material properties
        material_name: Material name
        regime: Processing regime
    
    Returns:
        Dictionary with wavelength quality assessment
    """
    from materials_database import get_absorption_coefficient
    
    result = {
        'quality': 0.0,
        'absorption_regime': '',
        'warnings': [],
        'recommendations': []
    }
    
    wavelength_nm = wavelength * 1e9
    
    # Get absorption coefficient
    try:
        from materials_database import get_material
        material_obj = get_material(material_name)
        alpha = get_absorption_coefficient(material_obj, wavelength_nm)
    except:
        # Fallback: use dictionary
        alpha = material_props.get('absorption_coeff', 1e3)
    
    if 'Plasma' in regime or 'Athermal' in regime:
        # NONLINEAR REGIME - Check multiphoton efficiency
        
        bandgap = material_props.get('bandgap', 5.0)
        n_photon = calculate_multiphoton_order(bandgap, wavelength)
        
        # Efficiency decreases with multiphoton order
        efficiency = 1.0 / np.sqrt(max(n_photon, 1))
        result['quality'] = min(efficiency, 1.0)
        result['absorption_regime'] = f'Nonlinear ({n_photon}-photon)'
        
        # Check Keldysh parameter (if dielectric)
        if bandgap > 3.0:
            gamma = calculate_keldysh_parameter(wavelength, intensity, bandgap)
            result['keldysh_parameter'] = gamma
            
            if gamma > 1.5:
                result['quality'] *= 0.5
                result['warnings'].append(
                    f"γ = {gamma:.2f} > 1.5: Inefficient ionization"
                )
                result['recommendations'].append(
                    "Consider shorter wavelength for lower multiphoton order"
                )
        
        # Wavelength-specific notes
        if wavelength_nm > 1000:
            result['warnings'].append(
                "Long wavelength (NIR) - high multiphoton order"
            )
        
    else:
        # LINEAR ABSORPTION REGIME
        
        result['absorption_regime'] = 'Linear'
        
        if alpha > 1e4:
            result['quality'] = 1.0
            result['absorption_regime'] = 'Linear (strong absorption)'
        elif alpha > 1e3:
            result['quality'] = 0.8
            result['absorption_regime'] = 'Linear (moderate absorption)'
        elif alpha > 1e2:
            result['quality'] = 0.6
            result['absorption_regime'] = 'Linear (weak absorption)'
            result['warnings'].append("Weak absorption - consider shorter wavelength")
        else:
            result['quality'] = 0.4
            result['absorption_regime'] = 'Linear (very weak absorption)'
            result['warnings'].append("Very weak absorption")
            result['recommendations'].append("Change to wavelength with stronger absorption")
    
    # Material-specific wavelength recommendations
    category = material_props.get('category', 'unknown')
    
    if category == 'metal':
        # Check reflectivity
        reflectivity = material_props.get('reflectivity', {})
        wavelength_str = f'{int(wavelength_nm)}nm'
        R = reflectivity.get(wavelength_str, 0.5)
        
        if R > 0.90:
            result['quality'] *= 0.5
            result['warnings'].append(
                f"High reflectivity (R = {R*100:.0f}%) at {wavelength_nm:.0f} nm"
            )
            result['recommendations'].append(
                "Use green (532 nm) or UV (355 nm) for lower reflectivity"
            )
    
    elif category == 'polymer':
        if wavelength_nm > 400:
            result['recommendations'].append(
                "UV wavelength (355 nm) optimal for polymers - photochemical ablation"
            )
    
    return result


def calculate_accumulation_quality(intensity: float,
                                   pulse_duration: float,
                                   num_pulses: int,
                                   rep_rate: float,
                                   material_props: dict,
                                   feature_size: float) -> Dict:
    """
    Calculate quality factor for multi-pulse heat accumulation
    
    Args:
        intensity: Intensity [W/cm²]
        pulse_duration: Pulse duration [s]
        num_pulses: Number of pulses
        rep_rate: Repetition rate [Hz]
        material_props: Material properties
        feature_size: Feature size [m]
    
    Returns:
        Dictionary with accumulation quality assessment
    """
    from physics_calculations import calculate_heat_accumulation
    
    result = {
        'quality': 1.0,
        'warnings': [],
        'heat_accumulation': False,
        'details': {}
    }
    
    D = material_props.get('thermal_diffusivity', 1e-6)
    time_between_pulses = 1 / rep_rate
    
    # Thermal diffusion length between pulses
    l_diff = calculate_thermal_diffusion_length(D, time_between_pulses)
    
    result['details']['diffusion_length'] = l_diff
    result['details']['feature_size'] = feature_size
    result['details']['rep_rate'] = rep_rate
    
    if l_diff < feature_size * 0.1:
        # Heat accumulates
        result['heat_accumulation'] = True
        
        # Estimate temperature rise
        fluence = intensity * pulse_duration
        pulse_energy_absorbed = fluence * np.pi * (feature_size/2)**2 * 1e4  # J
        
        thermal_conductivity = material_props.get('thermal_conductivity', 1.0)
        delta_T = calculate_heat_accumulation(
            pulse_energy_absorbed, rep_rate, feature_size, thermal_conductivity
        )
        
        result['details']['temperature_rise'] = delta_T
        
        # Quality penalty based on temperature rise
        if delta_T > 500:
            result['quality'] = 0.3
            result['warnings'].append(
                f"Severe heat accumulation: ΔT ≈ {delta_T:.0f} K"
            )
            result['warnings'].append(
                "CRITICAL: Reduce rep rate or increase scan speed"
            )
        elif delta_T > 200:
            result['quality'] = 0.6
            result['warnings'].append(
                f"Moderate heat accumulation: ΔT ≈ {delta_T:.0f} K"
            )
            result['warnings'].append(
                "Consider reducing rep rate"
            )
        elif delta_T > 100:
            result['quality'] = 0.8
            result['warnings'].append(
                f"Minor heat accumulation: ΔT ≈ {delta_T:.0f} K"
            )
        else:
            result['quality'] = 0.9
    else:
        # Heat dissipates between pulses
        result['heat_accumulation'] = False
        result['quality'] = 1.0
    
    # Check rep rate limits
    max_rep_rate = material_props.get('max_rep_rate')
    if max_rep_rate and rep_rate > max_rep_rate:
        result['quality'] *= 0.7
        result['warnings'].append(
            f"Rep rate {rep_rate/1e3:.0f} kHz exceeds recommended {max_rep_rate/1e3:.0f} kHz"
        )
    
    return result


def generate_recommendations(assessment: Dict,
                            material_name: str,
                            wavelength: float,
                            pulse_duration: float,
                            intensity: float,
                            thresholds: dict) -> List[str]:
    """
    Generate actionable recommendations based on quality assessment
    
    Args:
        assessment: Complete quality assessment dictionary
        material_name: Material name
        wavelength: Wavelength [m]
        pulse_duration: Pulse duration [s]
        intensity: Intensity [W/cm²]
        thresholds: Thresholds dictionary
    
    Returns:
        List of recommendation strings
    """
    recommendations = []
    Q = assessment['Q_total']
    components = assessment['components']
    
    # Overall quality assessment
    if Q >= 0.85:
        recommendations.append("✓ GOOD QUALITY ACHIEVED")
        recommendations.append("Parameters well-optimized for production")
        recommendations.append("Monitor process stability and consistency")
        return recommendations
    
    recommendations.append(f"Current Quality: Q = {Q:.3f}")
    
    # Prioritized recommendations (address worst component first)
    component_names = {
        'Q_regime': 'Regime',
        'Q_intensity': 'Intensity',
        'Q_temporal': 'Temporal',
        'Q_wavelength': 'Wavelength',
        'Q_accumulation': 'Accumulation'
    }
    
    # Sort components by quality (lowest first)
    sorted_components = sorted(components.items(), key=lambda x: x[1])
    
    for component, value in sorted_components:
        if value < 0.7:  # Address components below 0.7
            comp_name = component_names.get(component, component)
            recommendations.append(f"\n{comp_name} Quality Low (Q = {value:.2f}):")
            
            if component == 'Q_regime':
                recommendations.extend([
                    "→ Issue: Processing regime not optimal",
                    f"→ Consider reducing pulse duration (currently {pulse_duration*1e12:.1f} ps)",
                    "→ Target athermal regime for highest quality"
                ])
                
            elif component == 'Q_intensity':
                I_opt = assessment['details']['intensity_info'].get('optimal_intensity', 0)
                if I_opt > 0:
                    recommendations.extend([
                        "→ Issue: Intensity not optimal in window",
                        f"→ Adjust intensity toward optimal: {I_opt:.2e} W/cm²",
                        f"→ Clean window: [{thresholds['I_ablation']:.1e}, {thresholds['I_damage']:.1e}] W/cm²"
                    ])
                
            elif component == 'Q_temporal':
                tau_opt = thresholds.get('tau_optimal', (100e-15, 500e-15))
                recommendations.extend([
                    "→ Issue: Pulse duration not optimal",
                    f"→ Optimal range: {tau_opt[0]*1e15:.0f}-{tau_opt[1]*1e15:.0f} fs",
                    "→ May require laser system upgrade for shorter pulses"
                ])
                
            elif component == 'Q_wavelength':
                wl_info = assessment['details']['wavelength_info']
                if wl_info.get('recommendations'):
                    recommendations.append("→ Issue: Wavelength not optimal")
                    recommendations.extend([f"  → {rec}" for rec in wl_info['recommendations']])
                
            elif component == 'Q_accumulation':
                recommendations.extend([
                    "→ Issue: Heat accumulation affecting quality",
                    "→ Reduce repetition rate",
                    "→ Or increase scan speed",
                    "→ Consider burst mode processing"
                ])
    
    # Material-specific recommendations
    if 'Glass' in material_name or 'Sapphire' in material_name:
        if pulse_duration > 10e-12:
            recommendations.append("\n⚠ CRITICAL FOR GLASS:")
            recommendations.append("→ MUST use τ < 10 ps to avoid cracking")
            recommendations.append("→ Femtosecond laser mandatory")
    
    if 'Polymer' in material_name or 'ABF' in material_name:
        recommendations.append("\nPolymer-Specific Tips:")
        recommendations.append("→ Use UV wavelength (355 nm) for photochemical ablation")
        recommendations.append("→ Process in N₂ atmosphere to prevent carbonization")
    
    if Q < 0.5:
        recommendations.append("\n⚠ QUALITY CRITICALLY LOW")
        recommendations.append("→ Current parameters not suitable for production")
        recommendations.append("→ Significant optimization required")
    
    return recommendations


def interpret_quality(Q: float) -> Dict:
    """
    Interpret quality score with descriptive categories
    
    Args:
        Q: Quality score [0-1]
    
    Returns:
        Dictionary with interpretation
    """
    if Q >= 0.90:
        category = "EXCELLENT"
        description = "Optimal processing quality"
        haz = "<100 nm"
        roughness = "<50 nm Ra"
        applications = "Medical devices, optical components, precision micromachining"
        
    elif Q >= 0.70:
        category = "GOOD"
        description = "Acceptable for most applications"
        haz = "0.1-1 μm"
        roughness = "50-200 nm Ra"
        applications = "Microelectronics, general industrial micromachining"
        
    elif Q >= 0.50:
        category = "FAIR"
        description = "Marginal quality, some defects expected"
        haz = "1-10 μm"
        roughness = "200-500 nm Ra"
        applications = "Cutting, separation (edge quality not critical)"
        
    elif Q >= 0.30:
        category = "POOR"
        description = "Significant defects, limited applications"
        haz = "10-50 μm"
        roughness = ">500 nm Ra"
        applications = "Research, proof-of-concept only"
        
    else:
        category = "UNUSABLE"
        description = "Process failure, severe damage"
        haz = ">50 μm"
        roughness = ">1 μm Ra"
        applications = "Not suitable for any practical application"
    
    return {
        'category': category,
        'description': description,
        'expected_haz': haz,
        'expected_roughness': roughness,
        'suitable_applications': applications,
        'production_ready': Q >= 0.70
    }


def compare_parameter_sets(parameter_sets: List[Dict],
                          material_name: str,
                          material_props: dict,
                          thresholds: dict) -> Dict:
    """
    Compare multiple parameter sets and rank by quality
    
    Args:
        parameter_sets: List of parameter dictionaries
        material_name: Material name
        material_props: Material properties
        thresholds: Thresholds
    
    Returns:
        Dictionary with comparison results
    """
    results = []
    
    for i, params in enumerate(parameter_sets):
        assessment = calculate_quality_total(
            wavelength=params['wavelength'],
            pulse_duration=params['pulse_duration'],
            intensity=params['intensity'],
            material_name=material_name,
            material_props=material_props,
            thresholds=thresholds,
            feature_size=params.get('feature_size', 10e-6),
            num_pulses=params.get('num_pulses', 1),
            rep_rate=params.get('rep_rate', 1000)
        )
        
        results.append({
            'set_id': i,
            'parameters': params,
            'Q_total': assessment['Q_total'],
            'assessment': assessment
        })
    
    # Sort by quality (descending)
    results.sort(key=lambda x: x['Q_total'], reverse=True)
    
    return {
        'ranked_results': results,
        'best_set': results[0],
        'worst_set': results[-1],
        'quality_range': (results[-1]['Q_total'], results[0]['Q_total'])
    }


if __name__ == "__main__":
    # Test quality assessment
    print("="*70)
    print("QUALITY ASSESSMENT TEST")
    print("="*70)
    
    # Mock material properties and thresholds
    material_props = {
        'name': 'Fused Silica',
        'category': 'dielectric',
        'bandgap': 9.0,
        'tau_e_phonon': 1e-11,
        'sound_velocity': 5900,
        'thermal_diffusivity': 8e-7,
        'thermal_conductivity': 1.38,
        'refractive_index': 1.45,
        'absorption_coeff': 1e5,
    }
    
    thresholds = {
        'I_modification': 1e12,
        'I_ablation': 1e13,
        'I_damage': 1e14,
        'I_catastrophic': 5e14,
        'tau_critical': 10e-12,
        'tau_optimal': (50e-15, 500e-15),
    }
    
    # Test case: Optimal femtosecond glass processing
    print("\nTest Case: Femtosecond Glass Processing")
    print("-"*70)
    
    assessment = calculate_quality_total(
        wavelength=800e-9,
        pulse_duration=100e-15,
        intensity=5e13,
        material_name='Glass',
        material_props=material_props,
        thresholds=thresholds
    )
    
    print(f"\n{'QUALITY ASSESSMENT RESULTS'}")
    print(f"{'='*70}")
    print(f"Overall Quality: Q = {assessment['Q_total']:.3f}")
    print(f"Category: {assessment['interpretation']['category']}")
    print(f"Description: {assessment['interpretation']['description']}")
    
    print(f"\n{'COMPONENT BREAKDOWN'}")
    print(f"{'-'*70}")
    for component, value in assessment['components'].items():
        bar_length = int(value * 40)
        bar = '█' * bar_length + '░' * (40 - bar_length)
        print(f"{component:20s} {value:.3f} │{bar}│")
    
    print(f"\n{'REGIME INFORMATION'}")
    print(f"{'-'*70}")
    print(f"Regime: {assessment['regime']}")
    
    print(f"\n{'EXPECTED QUALITY CHARACTERISTICS'}")
    print(f"{'-'*70}")
    print(f"HAZ: {assessment['interpretation']['expected_haz']}")
    print(f"Roughness: {assessment['interpretation']['expected_roughness']}")
    print(f"Production Ready: {assessment['interpretation']['production_ready']}")
if assessment['warnings']:
    print(f"\n{'WARNINGS'}")
    print(f"{'-'*70}")
    for warning in assessment['warnings']:
        print(f"  ⚠ {warning}")

if assessment['recommendations']:
    print(f"\n{'RECOMMENDATIONS'}")
    print(f"{'-'*70}")
    for rec in assessment['recommendations']:
        print(f"  {rec}")

print("\n" + "="*70)
print("PARAMETER SET COMPARISON")
print("="*70)

# Compare multiple parameter sets
parameter_sets = [
    {
        'name': 'Femtosecond optimal',
        'wavelength': 800e-9,
        'pulse_duration': 100e-15,
        'intensity': 5e13,
    },
    {
        'name': 'Femtosecond high intensity',
        'wavelength': 800e-9,
        'pulse_duration': 100e-15,
        'intensity': 5e14,
    },
    {
        'name': 'Picosecond',
        'wavelength': 800e-9,
        'pulse_duration': 10e-12,
        'intensity': 5e13,
    },
    {
        'name': 'Nanosecond (unsuitable)',
        'wavelength': 800e-9,
        'pulse_duration': 10e-9,
        'intensity': 5e13,
    },
]

comparison = compare_parameter_sets(
    parameter_sets, 'Glass', material_props, thresholds
)

print(f"\n{'Rank':<6} {'Q':<8} {'Description':<30}")
print("-"*70)
for i, result in enumerate(comparison['ranked_results'], 1):
    params = result['parameters']
    name = params.get('name', f'Set {i}')
    Q = result['Q_total']
    print(f"{i:<6} {Q:.3f}    {name:<30}")

print(f"\nBest: {comparison['best_set']['parameters'].get('name', 'Unknown')}")
print(f"  Q = {comparison['best_set']['Q_total']:.3f}")

print(f"\nWorst: {comparison['worst_set']['parameters'].get('name', 'Unknown')}")
print(f"  Q = {comparison['worst_set']['Q_total']:.3f}")