"""
examples.py

Example usage of the laser processing framework
Demonstrates various use cases and workflows
"""

import numpy as np
from materials_database import get_material, list_materials
from thresholds_database import get_thresholds
from quality_assessment import calculate_quality_total
from optimization import (
    optimize_parameters_grid_search,
    optimize_for_application,
    suggest_laser_system
)
from process_report import generate_process_report, generate_comparison_report


def example_1_basic_quality_assessment():
    """
    Example 1: Basic quality assessment for femtosecond glass processing
    """
    print("\n" + "="*80)
    print("EXAMPLE 1: Basic Quality Assessment")
    print("="*80)
    
    # Define parameters
    material_name = 'Glass'
    wavelength = 800e-9  # 800 nm
    pulse_duration = 100e-15  # 100 fs
    intensity = 5e13  # 50 TW/cm²
    
    # Load material
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
    }
    
    # Calculate quality
    result = calculate_quality_total(
        wavelength=wavelength,
        pulse_duration=pulse_duration,
        intensity=intensity,
        material_name=material_name,
        material_props=material_props,
        thresholds=thresholds.__dict__
    )
    
    # Display results
    print(f"\nMaterial: {material_name}")
    print(f"Parameters: λ={wavelength*1e9:.0f}nm, τ={pulse_duration*1e15:.0f}fs, I={intensity:.1e} W/cm²")
    print(f"\nOverall Quality: Q = {result['Q_total']:.3f}")
    print(f"Category: {result['interpretation']['category']}")
    print(f"Regime: {result['regime']}")
    
    print(f"\nComponent Breakdown:")
    for component, value in result['components'].items():
        print(f"  {component}: {value:.3f}")
    
    if result['warnings']:
        print(f"\nWarnings:")
        for warning in result['warnings']:
            print(f"  ⚠ {warning}")


def example_2_parameter_optimization():
    """
    Example 2: Automated parameter optimization for polymer processing
    """
    print("\n" + "="*80)
    print("EXAMPLE 2: Parameter Optimization")
    print("="*80)
    
    # Optimize for ABF polymer processing
    result = optimize_parameters_grid_search(
        material_name='ABF_Polymer',
        quality_target=0.75,
        wavelength_options=[355e-9, 532e-9],
        pulse_duration_options=[10e-9, 20e-9, 50e-9])
    
    print(f"\n{result['message']}")
    print(f"\nBest Parameters:")
    params = result['best_parameters']
    print(f"  Wavelength:      {params['wavelength']*1e9:.0f} nm")
    print(f"  Pulse Duration:  {params['pulse_duration']*1e9:.0f} ns")
    print(f"  Intensity:       {params['intensity']:.2e} W/cm²")
    print(f"  Achieved Q:      {result['best_Q']:.3f}")
    
    print(f"\nTop 3 Solutions:")
    for i, sol in enumerate(result['top_5_results'][:3], 1):
        print(f"  {i}. Q={sol['Q']:.3f}: "
              f"λ={sol['wavelength']*1e9:.0f}nm, "
              f"τ={sol['pulse_duration']*1e9:.0f}ns, "
              f"I={sol['intensity']:.1e}")


def example_3_application_specific():
    """
    Example 3: Application-specific optimization (microvia drilling)
    """
    print("\n" + "="*80)
    print("EXAMPLE 3: Application-Specific Optimization")
    print("="*80)
    
    # Optimize for microvia drilling
    result = optimize_for_application(
        application='microvia',
        material_name='ABF_Polymer',
        feature_size=20e-6,
        throughput_target=1000  # 1000 vias/second
    )
    
    print(f"\nApplication: Microvia Drilling in ABF")
    print(f"Target Quality: {result['application_specs']['quality_target']}")
    print(f"Achieved Quality: {result['best_Q']:.3f}")
    print(f"Budget Level: {result['budget_level']}")
    
    if result['success']:
        params = result['best_parameters']
        print(f"\nRecommended Parameters:")
        print(f"  Wavelength:      {params['wavelength']*1e9:.0f} nm")
        print(f"  Pulse Duration:  {params['pulse_duration']*1e9:.0f} ns")
        print(f"  Intensity:       {params['intensity']:.2e} W/cm²")
        
        if 'throughput_estimate' in result:
            throughput = result['throughput_estimate']
            print(f"\nThroughput Estimate:")
            print(f"  Ablation/pulse:  {throughput['depth_per_pulse_um']:.2f} μm")
            print(f"  Pulses needed:   {throughput['pulses_needed']}")
            print(f"  Features/second: {throughput['features_per_second']:.0f}")
            print(f"  Target met:      {'Yes' if throughput['meets_target'] else 'No'}")


def example_4_laser_system_selection():
    """
    Example 4: Laser system recommendation based on requirements
    """
    print("\n" + "="*80)
    print("EXAMPLE 4: Laser System Selection")
    print("="*80)
    
    # Case 1: High-quality glass processing
    print("\nCase 1: High-Quality Glass Processing")
    print("-" * 40)
    
    result1 = suggest_laser_system(
        material_name='Glass',
        quality_target=0.90,
        budget=350000
    )
    
    if result1['success']:
        system = result1['recommended_system']
        print(f"Recommended: {system['name']}")
        print(f"  Cost:            ${system['cost']/1000:.0f}k")
        print(f"  Pulse Duration:  {system['pulse_duration']*1e15:.0f} fs")
        print(f"  Wavelength:      {system['wavelength']*1e9:.0f} nm")
        print(f"  Max Quality:     {system['max_quality']:.2f}")
        print(f"  Notes:           {system['notes']}")
    else:
        print(f"No suitable system found")
        for rec in result1.get('recommendations', []):
            print(f"  → {rec}")
    
    # Case 2: Cost-effective polymer processing
    print("\nCase 2: Cost-Effective Polymer Processing")
    print("-" * 40)
    
    result2 = suggest_laser_system(
        material_name='ABF_Polymer',
        quality_target=0.70,
        budget=100000
    )
    
    if result2['success']:
        system = result2['recommended_system']
        print(f"Recommended: {system['name']}")
        print(f"  Cost:            ${system['cost']/1000:.0f}k")
        print(f"  Pulse Duration:  {system['pulse_duration']*1e9:.0f} ns")
        print(f"  Max Quality:     {system['max_quality']:.2f}")


def example_5_process_report():
    """
    Example 5: Generate comprehensive process report
    """
    print("\n" + "="*80)
    print("EXAMPLE 5: Process Report Generation")
    print("="*80)
    
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


def example_6_parameter_comparison():
    """
    Example 6: Compare multiple parameter sets
    """
    print("\n" + "="*80)
    print("EXAMPLE 6: Parameter Comparison")
    print("="*80)
    
    # Define parameter sets to compare
    parameter_sets = [
        {
            'name': 'Femtosecond - Optimal',
            'wavelength': 800e-9,
            'pulse_duration': 100e-15,
            'intensity': 5e13,
        },
        {
            'name': 'Femtosecond - High Intensity',
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
    
    comparison = generate_comparison_report(
        parameter_sets,
        'Glass',
        output_format='text'
    )
    
    print(comparison)


def example_7_regime_boundaries():
    """
    Example 7: Explore regime boundaries
    """
    print("\n" + "="*80)
    print("EXAMPLE 7: Regime Boundary Exploration")
    print("="*80)
    
    material_name = 'Glass'
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
    }
    
    # Explore intensity at fixed pulse duration
    print("\nIntensity Scan (τ = 100 fs):")
    print("-" * 60)
    print(f"{'Intensity':<20} {'Regime':<25} {'Q':<10}")
    print("-" * 60)
    
    intensities = [1e11, 1e12, 1e13, 5e13, 1e14, 5e14]
    pulse_duration = 100e-15
    
    for I in intensities:
        result = calculate_quality_total(
            wavelength=800e-9,
            pulse_duration=pulse_duration,
            intensity=I,
            material_name=material_name,
            material_props=material_props,
            thresholds=thresholds.__dict__
        )
        
        print(f"{I:.1e} W/cm²      {result['regime']:<25} {result['Q_total']:.3f}")
    
    # Explore pulse duration at fixed intensity
    print("\n\nPulse Duration Scan (I = 5×10¹³ W/cm²):")
    print("-" * 60)
    print(f"{'Pulse Duration':<20} {'Regime':<25} {'Q':<10}")
    print("-" * 60)
    
    pulse_durations = [100e-15, 500e-15, 1e-12, 10e-12, 100e-12, 1e-9, 10e-9]
    intensity = 5e13
    
    for tau in pulse_durations:
        result = calculate_quality_total(
            wavelength=800e-9,
            pulse_duration=tau,
            intensity=intensity,
            material_name=material_name,
            material_props=material_props,
            thresholds=thresholds.__dict__
        )
        
        tau_str = f"{tau*1e15:.0f} fs" if tau < 1e-12 else f"{tau*1e12:.1f} ps" if tau < 1e-9 else f"{tau*1e9:.0f} ns"
        print(f"{tau_str:<20} {result['regime']:<25} {result['Q_total']:.3f}")


def example_8_material_survey():
    """
    Example 8: Survey quality across different materials
    """
    print("\n" + "="*80)
    print("EXAMPLE 8: Material Survey")
    print("="*80)
    
    # Fixed laser parameters
    wavelength = 800e-9
    pulse_duration = 100e-15
    
    print("\nFixed Parameters:")
    print(f"  Wavelength:      {wavelength*1e9:.0f} nm")
    print(f"  Pulse Duration:  {pulse_duration*1e15:.0f} fs")
    
    print("\n" + "-" * 80)
    print(f"{'Material':<20} {'Optimal I (W/cm²)':<25} {'Max Q':<10} {'Notes':<30}")
    print("-" * 80)
    
    materials_to_test = ['Glass', 'ABF_Polymer', 'Copper', 'Silicon']
    
    for mat_name in materials_to_test:
        try:
            material = get_material(mat_name)
            thresholds = get_thresholds(mat_name)
            
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
            
            # Find optimal intensity
            I_opt = np.sqrt(thresholds.I_ablation * thresholds.I_damage)
            
            # Calculate quality at optimal
            result = calculate_quality_total(
                wavelength=wavelength,
                pulse_duration=pulse_duration,
                intensity=I_opt,
                material_name=mat_name,
                material_props=material_props,
                thresholds=thresholds.__dict__
            )
            
            # Get main limitation
            components = result['components']
            limiting_factor = min(components.items(), key=lambda x: x[1])
            
            print(f"{mat_name:<20} {I_opt:.2e}            {result['Q_total']:.3f}      "
                  f"Limited by {limiting_factor[0]}")
            
        except Exception as e:
            print(f"{mat_name:<20} ERROR: {str(e)}")


def example_9_cost_quality_tradeoff():
    """
    Example 9: Analyze cost-quality trade-offs
    """
    print("\n" + "="*80)
    print("EXAMPLE 9: Cost-Quality Trade-off Analysis")
    print("="*80)
    
    material_name = 'ABF_Polymer'
    quality_target = 0.75
    
    print(f"\nMaterial: {material_name}")
    print(f"Quality Target: Q ≥ {quality_target}")
    print("\n" + "-" * 80)
    print(f"{'Budget Level':<15} {'Cost':<15} {'Achieved Q':<12} {'Success':<10}")
    print("-" * 80)
    
    from optimization import optimize_parameters_cost_constrained
    
    for budget_level in ['low', 'moderate', 'high']:
        result = optimize_parameters_cost_constrained(
            material_name=material_name,
            quality_target=quality_target,
            budget_level=budget_level
        )
        
        success = "Yes" if result['success'] else "No"
        print(f"{budget_level:<15} ${result['estimated_laser_cost']/1000:.0f}k         "
              f"{result['best_Q']:.3f}        {success:<10}")
        
        if result['success']:
            params = result['best_parameters']
            print(f"  → λ={params['wavelength']*1e9:.0f}nm, "
                  f"τ={params['pulse_duration']*1e12:.2f}ps, "
                  f"I={params['intensity']:.1e} W/cm²")


def example_10_troubleshooting():
    """
    Example 10: Troubleshooting problematic parameters
    """
    print("\n" + "="*80)
    print("EXAMPLE 10: Troubleshooting Guide")
    print("="*80)
    
    # Problematic case: Glass with nanosecond pulse
    print("\nProblematic Case: Glass with Nanosecond Pulse")
    print("-" * 80)
    
    material = get_material('Glass')
    thresholds = get_thresholds('Glass')
    
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
    
    result = calculate_quality_total(
        wavelength=800e-9,
        pulse_duration=10e-9,  # 10 ns - TOO LONG
        intensity=5e13,
        material_name='Glass',
        material_props=material_props,
        thresholds=thresholds.__dict__
    )
    
    print(f"Quality: Q = {result['Q_total']:.3f} (POOR)")
    
    print(f"\nIdentified Problems:")
    for warning in result['warnings']:
        print(f"  ⚠ {warning}")
    
    print(f"\nRecommendations:")
    for rec in result['recommendations']:
        print(f"  {rec}")
    
    # Show corrected parameters
    print("\n\nCorrected Parameters:")
    print("-" * 80)
    
    result_corrected = calculate_quality_total(
        wavelength=800e-9,
        pulse_duration=100e-15,  # 100 fs - CORRECTED
        intensity=5e13,
        material_name='Glass',
        material_props=material_props,
        thresholds=thresholds.__dict__
    )
    
    print(f"Quality: Q = {result_corrected['Q_total']:.3f} (EXCELLENT)")
    print(f"Regime: {result_corrected['regime']}")
    print(f"\nImprovement: ΔQ = +{result_corrected['Q_total'] - result['Q_total']:.3f}")


def run_all_examples():
    """Run all examples"""
    
    examples = [
        example_1_basic_quality_assessment,
        example_2_parameter_optimization,
        example_3_application_specific,
        example_4_laser_system_selection,
        example_5_process_report,
        example_6_parameter_comparison,
        example_7_regime_boundaries,
        example_8_material_survey,
        example_9_cost_quality_tradeoff,
        example_10_troubleshooting,
    ]
    
    for example in examples:
        try:
            example()
        except Exception as e:
            print(f"\nERROR in {example.__name__}: {str(e)}")
        
        input("\nPress Enter to continue to next example...")


if __name__ == "__main__":
    print("="*80)
    print("LASER PROCESSING FRAMEWORK - EXAMPLE DEMONSTRATIONS")
    print("="*80)
    print("\nThis script demonstrates various use cases of the framework.")
    print("Examples will run sequentially. Press Enter after each to continue.")
    
    input("\nPress Enter to start...")
    
    run_all_examples()
    
    print("\n" + "="*80)
    print("ALL EXAMPLES COMPLETED")
    print("="*80)
