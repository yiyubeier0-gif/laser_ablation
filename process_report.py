"""
process_report.py

Generate comprehensive processing reports with visualizations
"""

import numpy as np
from typing import Dict, Optional
from datetime import datetime
import json

from materials_database import get_material
from thresholds_database import get_thresholds
from physics_calculations import calculate_all_parameters
from quality_assessment import calculate_quality_total
from regime_classification import classify_regime


def generate_process_report(material_name: str,
                           wavelength: float,
                           pulse_duration: float,
                           intensity: float,
                           feature_size: float = 20e-6,
                           target_depth: float = 10e-6,
                           num_pulses: int = 1,
                           rep_rate: float = 50000,
                           spot_diameter: float = 20e-6,
                           output_format: str = 'text') -> str:
    """
    Generate comprehensive processing report
    
    Args:
        material_name: Material name
        wavelength: Wavelength [m]
        pulse_duration: Pulse duration [s]
        intensity: Peak intensity [W/cm²]
        feature_size: Feature size [m]
        target_depth: Target depth [m]
        num_pulses: Number of pulses
        rep_rate: Repetition rate [Hz]
        spot_diameter: Spot diameter [m]
        output_format: 'text', 'markdown', or 'json'
    
    Returns:
        Formatted report string
    """
    
    # Load material data
    material = get_material(material_name)
    thresholds = get_thresholds(material_name)
    
    material_props = {
        'name': material.name,
        'category': material.category,
        'bandgap': material.bandgap,
        'tau_e_phonon': material.tau_e_phonon,
        'sound_velocity': material.sound_velocity,
        'thermal_diffusivity': material.thermal_diffusivity,
        'thermal_conductivity': material.thermal_conductivity,
        'refractive_index': material.refractive_index,
        'density': material.density,
        'specific_heat': material.specific_heat,
    }
    
    # Calculate all parameters
    from physics_calculations import calculate_fluence, calculate_pulse_energy
    fluence = calculate_fluence(intensity, pulse_duration)
    pulse_energy = calculate_pulse_energy(intensity, spot_diameter, pulse_duration)
    
    # Quality assessment
    quality_result = calculate_quality_total(
        wavelength=wavelength,
        pulse_duration=pulse_duration,
        intensity=intensity,
        material_name=material_name,
        material_props=material_props,
        thresholds=thresholds.__dict__,
        feature_size=feature_size,
        num_pulses=num_pulses,
        rep_rate=rep_rate
    )
    
    # Physics calculations
    physics_params = calculate_all_parameters(
        pulse_energy=pulse_energy,
        spot_diameter=spot_diameter,
        pulse_duration=pulse_duration,
        wavelength=wavelength,
        material_props={
            'bandgap': material.bandgap,
            'refractive_index': material.refractive_index,
            'thermal_diffusivity': material.thermal_diffusivity,
            'sound_velocity': material.sound_velocity,
            'n2': material.n2,
        }
    )
    
    # Ablation depth estimation
    from physics_calculations import calculate_ablation_depth_multipulse
    if thresholds.F_ablation_single_pulse:
        F_th = thresholds.F_ablation_single_pulse
    else:
        F_th = 1.0  # Default
    
    alpha = 1e5  # Default absorption coefficient (effective)
    ablation_depth = calculate_ablation_depth_multipulse(
        fluence, F_th, alpha, num_pulses, material.incubation_parameter
    )
    
    if ablation_depth > 0:
        pulses_needed = int(np.ceil(target_depth * 1e6 / ablation_depth))
    else:
        pulses_needed = 'N/A (below threshold)'
    
    # Generate report based on format
    if output_format == 'json':
        return _generate_json_report(
            material, thresholds, quality_result, physics_params,
            wavelength, pulse_duration, intensity, fluence, pulse_energy,
            spot_diameter, ablation_depth, pulses_needed
        )
    elif output_format == 'markdown':
        return _generate_markdown_report(
            material, thresholds, quality_result, physics_params,
            wavelength, pulse_duration, intensity, fluence, pulse_energy,
            spot_diameter, ablation_depth, pulses_needed, target_depth
        )
    else:  # text
        return _generate_text_report(
            material, thresholds, quality_result, physics_params,
            wavelength, pulse_duration, intensity, fluence, pulse_energy,
            spot_diameter, ablation_depth, pulses_needed, target_depth, rep_rate
        )


def _generate_text_report(material, thresholds, quality_result, physics_params,
                         wavelength, pulse_duration, intensity, fluence, pulse_energy,
                         spot_diameter, ablation_depth, pulses_needed, target_depth, rep_rate):
    """Generate text format report"""
    
    report = []
    
    # Header
    report.append("=" * 80)
    report.append("LASER PROCESSING ANALYSIS REPORT")
    report.append("=" * 80)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")
    
    # Material Information
    report.append("MATERIAL INFORMATION")
    report.append("-" * 80)
    report.append(f"Material:              {material.name}")
    report.append(f"Category:              {material.category}")
    report.append(f"Bandgap:               {material.bandgap:.2f} eV")
    report.append(f"Refractive Index:      {material.refractive_index:.2f}")
    report.append(f"Thermal Diffusivity:   {material.thermal_diffusivity:.2e} m²/s")
    if material.notes:
        report.append(f"Notes:                 {material.notes}")
    report.append("")
    
    # Laser Parameters
    report.append("LASER PARAMETERS")
    report.append("-" * 80)
    report.append(f"Wavelength:            {wavelength*1e9:.0f} nm")
    report.append(f"Pulse Duration:        {pulse_duration*1e12:.2f} ps ({pulse_duration*1e15:.1f} fs)")
    report.append(f"Peak Intensity:        {intensity:.2e} W/cm²")
    report.append(f"Fluence:               {fluence:.2f} J/cm²")
    report.append(f"Pulse Energy:          {pulse_energy*1e6:.1f} μJ")
    report.append(f"Spot Diameter:         {spot_diameter*1e6:.1f} μm")
    report.append(f"Repetition Rate:       {rep_rate/1e3:.0f} kHz")
    report.append("")
    
    # Processing Regime
    report.append("PROCESSING REGIME")
    report.append("-" * 80)
    report.append(f"Regime:                {quality_result['regime']}")
    report.append(f"Regime Quality:        Q_regime = {quality_result['components']['Q_regime']:.3f}")
    report.append("")
    
    # Quality Assessment
    report.append("QUALITY ASSESSMENT")
    report.append("-" * 80)
    report.append(f"Overall Quality:       Q = {quality_result['Q_total']:.3f}")
    report.append(f"Category:              {quality_result['interpretation']['category']}")
    report.append(f"Description:           {quality_result['interpretation']['description']}")
    report.append("")
    
    report.append("Component Breakdown:")
    for component, value in quality_result['components'].items():
        bar_length = int(value * 40)
        bar = '█' * bar_length + '░' * (40 - bar_length)
        report.append(f"  {component:20s} {value:.3f} │{bar}│")
    report.append("")
    
    # Expected Quality Characteristics
    report.append("EXPECTED QUALITY CHARACTERISTICS")
    report.append("-" * 80)
    report.append(f"Heat-Affected Zone:    {quality_result['interpretation']['expected_haz']}")
    report.append(f"Surface Roughness:     {quality_result['interpretation']['expected_roughness']}")
    report.append(f"Production Ready:      {'Yes' if quality_result['interpretation']['production_ready'] else 'No'}")
    report.append(f"Suitable For:          {quality_result['interpretation']['suitable_applications']}")
    report.append("")
    
    # Threshold Comparison
    report.append("THRESHOLD COMPARISON")
    report.append("-" * 80)
    report.append(f"I_modification:        {thresholds.I_modification:.2e} W/cm²")
    report.append(f"I_ablation:            {thresholds.I_ablation:.2e} W/cm²")
    report.append(f"I_damage:              {thresholds.I_damage:.2e} W/cm²")
    report.append(f"I_catastrophic:        {thresholds.I_catastrophic:.2e} W/cm²")
    report.append("")
    
    # Window position
    if thresholds.I_ablation < intensity < thresholds.I_damage:
        I_opt = np.sqrt(thresholds.I_ablation * thresholds.I_damage)
        position_pct = (np.log10(intensity) - np.log10(thresholds.I_ablation)) / \
                       (np.log10(thresholds.I_damage) - np.log10(thresholds.I_ablation)) * 100
        report.append(f"Clean Window:          [{thresholds.I_ablation:.1e}, {thresholds.I_damage:.1e}] W/cm²")
        report.append(f"Optimal Intensity:     {I_opt:.2e} W/cm²")
        report.append(f"Current Position:      ✓ WITHIN WINDOW ({position_pct:.0f}% across)")
    elif intensity < thresholds.I_ablation:
        report.append(f"Position:              ✗ BELOW ablation threshold")
    else:
        report.append(f"Position:              ✗ ABOVE damage threshold")
    report.append("")
    
    # Physics Calculations
    report.append("PHYSICS CALCULATIONS")
    report.append("-" * 80)
    report.append(f"Photon Energy:         {physics_params['photon_energy_eV']:.2f} eV")
    
    if 'multiphoton_order' in physics_params:
        report.append(f"Multiphoton Order:     {physics_params['multiphoton_order']}")
        report.append(f"Keldysh Parameter:     {physics_params['keldysh_parameter']:.2f}")
        report.append(f"Critical Density:      {physics_params['critical_density']:.2e} cm⁻³")
        report.append(f"Est. Electron Density: {physics_params['estimated_electron_density']:.2e} cm⁻³")
    
    report.append("")
    report.append("Temporal Hierarchy:")
    report.append(f"  τ_pulse:             {pulse_duration*1e12:.2f} ps")
    if 'stress_confinement_time' in physics_params:
        report.append(f"  τ_stress:            {physics_params['stress_confinement_time']*1e9:.1f} ns")
    report.append(f"  τ_thermal:           {physics_params['thermal_diffusion_time']*1e6:.1f} μs")
    report.append(f"  l_diff (thermal):    {physics_params['thermal_diffusion_length']*1e9:.1f} nm")
    report.append("")
    
    report.append("Beam Propagation:")
    report.append(f"  Rayleigh Range:      {physics_params['rayleigh_range']*1e6:.1f} μm")
    report.append(f"  Depth of Focus:      {physics_params['depth_of_focus']*1e6:.1f} μm")
    
    if 'self_focusing_warning' in physics_params and physics_params['self_focusing_warning']:
        report.append(f"  ⚠ Self-focusing at:  {physics_params['self_focusing_distance']*1e3:.2f} mm")
    report.append("")
    
    # Ablation Predictions
    report.append("ABLATION PREDICTIONS")
    report.append("-" * 80)
    report.append(f"Ablation per Pulse:    {ablation_depth:.2f} μm")
    report.append(f"Target Depth:          {target_depth*1e6:.0f} μm")
    report.append(f"Pulses Required:       {pulses_needed}")
    
    if isinstance(pulses_needed, int):
        time_per_feature = pulses_needed / rep_rate
        report.append(f"Time per Feature:      {time_per_feature*1e6:.1f} μs")
        features_per_second = 1 / time_per_feature
        report.append(f"Throughput:            {features_per_second:.0f} features/s")
    report.append("")
    
    # Warnings
    if quality_result['warnings']:
        report.append("WARNINGS")
        report.append("-" * 80)
        for warning in quality_result['warnings']:
            report.append(f"  ⚠ {warning}")
        report.append("")
    
    # Recommendations
    if quality_result['recommendations']:
        report.append("RECOMMENDATIONS")
        report.append("-" * 80)
        for rec in quality_result['recommendations']:
            report.append(f"  {rec}")
        report.append("")
    
    # Footer
    report.append("=" * 80)
    report.append("END OF REPORT")
    report.append("=" * 80)
    
    return "\n".join(report)


def _generate_markdown_report(material, thresholds, quality_result, physics_params,
                              wavelength, pulse_duration, intensity, fluence, pulse_energy,
                              spot_diameter, ablation_depth, pulses_needed, target_depth):
    """Generate markdown format report"""
    
    md = []
    
    # Header
    md.append("# Laser Processing Analysis Report")
    md.append("")
    md.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    md.append("")
    
    # Material
    md.append("## Material Information")
    md.append("")
    md.append(f"- **Material:** {material.name}")
    md.append(f"- **Category:** {material.category}")
    md.append(f"- **Bandgap:** {material.bandgap:.2f} eV")
    md.append(f"- **Refractive Index:** {material.refractive_index:.2f}")
    md.append("")
    
    # Laser Parameters
    md.append("## Laser Parameters")
    md.append("")
    md.append("| Parameter | Value |")
    md.append("|-----------|-------|")
    md.append(f"| Wavelength | {wavelength*1e9:.0f} nm |")
    md.append(f"| Pulse Duration | {pulse_duration*1e12:.2f} ps |")
    md.append(f"| Peak Intensity | {intensity:.2e} W/cm² |")
    md.append(f"| Fluence | {fluence:.2f} J/cm² |")
    md.append(f"| Pulse Energy | {pulse_energy*1e6:.1f} μJ |")
    md.append(f"| Spot Diameter | {spot_diameter*1e6:.1f} μm |")
    md.append("")
    
    # Quality Assessment
    md.append("## Quality Assessment")
    md.append("")
    md.append(f"**Overall Quality:** Q = {quality_result['Q_total']:.3f}")
    md.append("")
    md.append(f"**Category:** {quality_result['interpretation']['category']}")
    md.append("")
    md.append(f"*{quality_result['interpretation']['description']}*")
    md.append("")
    
    md.append("### Component Breakdown")
    md.append("")
    for component, value in quality_result['components'].items():
        percentage = int(value * 100)
        md.append(f"- **{component}:** {value:.3f} ({percentage}%)")
    md.append("")
    
    # Threshold Comparison
    md.append("## Threshold Comparison")
    md.append("")
    md.append("| Threshold | Value |")
    md.append("|-----------|-------|")
    md.append(f"| I_modification | {thresholds.I_modification:.2e} W/cm² |")
    md.append(f"| I_ablation | {thresholds.I_ablation:.2e} W/cm² |")
    md.append(f"| I_damage | {thresholds.I_damage:.2e} W/cm² |")
    md.append(f"| I_catastrophic | {thresholds.I_catastrophic:.2e} W/cm² |")
    md.append("")
    
    # Window position
    if thresholds.I_ablation < intensity < thresholds.I_damage:
        md.append("✅ **Within Clean Ablation Window**")
    elif intensity < thresholds.I_ablation:
        md.append("❌ **Below Ablation Threshold**")
    else:
        md.append("⚠️ **Above Damage Threshold**")
    md.append("")
    
    # Warnings
    if quality_result['warnings']:
        md.append("## Warnings")
        md.append("")
        for warning in quality_result['warnings']:
            md.append(f"- ⚠️ {warning}")
        md.append("")
    
    # Recommendations
    if quality_result['recommendations']:
        md.append("## Recommendations")
        md.append("")
        for rec in quality_result['recommendations']:
            md.append(f"- {rec}")
        md.append("")
    
    return "\n".join(md)


def _generate_json_report(material, thresholds, quality_result, physics_params,
                         wavelength, pulse_duration, intensity, fluence, pulse_energy,
                         spot_diameter, ablation_depth, pulses_needed):
    """Generate JSON format report"""
    
    report_dict = {
        'metadata': {
            'generated': datetime.now().isoformat(),
            'version': '1.0'
        },
        'material': {
            'name': material.name,
            'category': material.category,
            'bandgap_eV': material.bandgap,
            'refractive_index': material.refractive_index,
        },
        'laser_parameters': {
            'wavelength_nm': wavelength * 1e9,
            'pulse_duration_ps': pulse_duration * 1e12,
            'pulse_duration_fs': pulse_duration * 1e15,
            'intensity_W_per_cm2': intensity,
            'fluence_J_per_cm2': fluence,
            'pulse_energy_uJ': pulse_energy * 1e6,
            'spot_diameter_um': spot_diameter * 1e6,
        },
        'quality_assessment': {
            'Q_total': quality_result['Q_total'],
            'regime': quality_result['regime'],
            'components': quality_result['components'],
            'interpretation': quality_result['interpretation'],
            'warnings': quality_result['warnings'],
            'recommendations': quality_result['recommendations'],
        },
        'thresholds': {
            'I_modification': thresholds.I_modification,
            'I_ablation': thresholds.I_ablation,
            'I_damage': thresholds.I_damage,
            'I_catastrophic': thresholds.I_catastrophic,
        },
        'physics': physics_params,
        'ablation': {
            'depth_per_pulse_um': ablation_depth,
            'pulses_needed': pulses_needed if isinstance(pulses_needed, int) else None,
        }
    }
    
    return json.dumps(report_dict, indent=2)


def generate_comparison_report(parameter_sets: list,
                              material_name: str,
                              output_format: str = 'text') -> str:
    """
    Generate comparison report for multiple parameter sets
    
    Args:
        parameter_sets: List of parameter dictionaries
        material_name: Material name
        output_format: 'text' or 'markdown'
    
    Returns:
        Formatted comparison report
    """
    
    material = get_material(material_name)
    thresholds = get_thresholds(material_name)
    
    results = []
    
    for params in parameter_sets:
        material_props = {
            'name': material.name,
            'category': material.category,
            'bandgap': material.bandgap,
            'tau_e_phonon': material.tau_e_phonon,
            'sound_velocity': material.sound_velocity,
            'thermal_diffusivity': material.thermal_diffusivity,
            'thermal_conductivity': material.thermal_conductivity,
            'refractive_index': material.refractive_index,
        }
        
        quality_result = calculate_quality_total(
            wavelength=params['wavelength'],
            pulse_duration=params['pulse_duration'],
            intensity=params['intensity'],
            material_name=material_name,
            material_props=material_props,
            thresholds=thresholds.__dict__
        )
        
        results.append({
            'name': params.get('name', 'Unnamed'),
            'params': params,
            'Q': quality_result['Q_total'],
            'regime': quality_result['regime'],
            'quality_result': quality_result
        })
    
    # Sort by quality
    results.sort(key=lambda x: x['Q'], reverse=True)
    
    if output_format == 'markdown':
        return _generate_markdown_comparison(results, material_name)
    else:
        return _generate_text_comparison(results, material_name)


def _generate_text_comparison(results, material_name):
    """Generate text comparison report"""
    
    report = []
    
    report.append("=" * 80)
    report.append(f"PARAMETER COMPARISON REPORT - {material_name}")
    report.append("=" * 80)
    report.append("")
    
    # Summary table
    report.append("SUMMARY")
    report.append("-" * 80)
    report.append(f"{'Rank':<6} {'Name':<25} {'Q':<8} {'Regime':<25}")
    report.append("-" * 80)
    
    for i, result in enumerate(results, 1):
        report.append(f"{i:<6} {result['name']:<25} {result['Q']:.3f}    {result['regime']:<25}")
    
    report.append("")
    
    # Detailed comparison
    report.append("DETAILED COMPARISON")
    report.append("-" * 80)
    
    for i, result in enumerate(results, 1):
        report.append(f"\n{i}. {result['name']} (Q = {result['Q']:.3f})")
        report.append("-" * 40)
        
        params = result['params']
        report.append(f"  Wavelength:      {params['wavelength']*1e9:.0f} nm")
        report.append(f"  Pulse Duration:  {params['pulse_duration']*1e12:.2f} ps")
        report.append(f"  Intensity:       {params['intensity']:.2e} W/cm²")
        report.append(f"  Regime:          {result['regime']}")
        
        components = result['quality_result']['components']
        report.append(f"  Q_regime:        {components['Q_regime']:.3f}")
        report.append(f"  Q_intensity:     {components['Q_intensity']:.3f}")
        report.append(f"  Q_temporal:      {components['Q_temporal']:.3f}")
        report.append(f"  Q_wavelength:    {components['Q_wavelength']:.3f}")
    
    report.append("")
    report.append("=" * 80)
    
    return "\n".join(report)


def _generate_markdown_comparison(results, material_name):
    """Generate markdown comparison report"""
    
    md = []
    
    md.append(f"# Parameter Comparison Report - {material_name}")
    md.append("")
    
    # Summary table
    md.append("## Summary")
    md.append("")
    md.append("| Rank | Name | Quality | Regime |")
    md.append("|------|------|---------|--------|")
    
    for i, result in enumerate(results, 1):
        md.append(f"| {i} | {result['name']} | {result['Q']:.3f} | {result['regime']} |")
    
    md.append("")
    
    # Detailed sections
    md.append("## Detailed Comparison")
    md.append("")
    
    for i, result in enumerate(results, 1):
        md.append(f"### {i}. {result['name']}")
        md.append("")
        md.append(f"**Quality:** Q = {result['Q']:.3f}")
        md.append("")
        
        params = result['params']
        md.append("**Parameters:**")
        md.append(f"- Wavelength: {params['wavelength']*1e9:.0f} nm")
        md.append(f"- Pulse Duration: {params['pulse_duration']*1e12:.2f} ps")
        md.append(f"- Intensity: {params['intensity']:.2e} W/cm²")
        md.append("")
        
        md.append("**Quality Components:**")
        components = result['quality_result']['components']
        for comp, value in components.items():
            md.append(f"- {comp}: {value:.3f}")
        md.append("")
    
    return "\n".join(md)


if __name__ == "__main__":
    # Test report generation
    print("="*80)
    print("PROCESS REPORT GENERATION TEST")
    print("="*80)
    
    # Generate sample report
    report = generate_process_report(
        material_name='Glass',
        wavelength=800e-9,
        pulse_duration=100e-15,
        intensity=5e13,
        feature_size=20e-6,
        target_depth=10e-6,
        num_pulses=1,
        rep_rate=50000,
        spot_diameter=20e-6,
        output_format='text'
    )
    
    print(report)
    
    print("\n" + "="*80)
    print("COMPARISON REPORT TEST")
    print("="*80)
    
    # Generate comparison report
    parameter_sets = [
        {
            'name': 'Femtosecond Optimal',
            'wavelength': 800e-9,
            'pulse_duration': 100e-15,
            'intensity': 5e13,
        },
        {
            'name': 'Picosecond',
            'wavelength': 800e-9,
            'pulse_duration': 10e-12,
            'intensity': 5e13,
        },
        {
            'name': 'Nanosecond',
            'wavelength': 800e-9,
            'pulse_duration': 10e-9,
            'intensity': 5e13,
        },
    ]
    
    comparison = generate_comparison_report(
        parameter_sets,
        'Glass',
        output_format='text'
    )
    
    print(comparison)