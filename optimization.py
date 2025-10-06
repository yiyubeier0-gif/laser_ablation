## File 6: `optimization.py`

"""
optimization.py

Parameter optimization routines for laser processing
Automated parameter selection based on requirements
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from scipy.optimize import minimize, differential_evolution
from quality_assessment import calculate_quality_total
from materials_database import get_material, MATERIALS_DB
from thresholds_database import get_thresholds, THRESHOLDS_DB


def optimize_parameters_grid_search(material_name: str,
                                    quality_target: float = 0.85,
                                    wavelength_options: Optional[List[float]] = None,
                                    pulse_duration_options: Optional[List[float]] = None,
                                    feature_size: float = 10e-6) -> Dict:
    """
    Grid search optimization over parameter space
    
    Args:
        material_name: Material to process
        quality_target: Target quality (0-1)
        wavelength_options: List of available wavelengths [m]
        pulse_duration_options: List of available pulse durations [s]
        feature_size: Feature size [m]
    
    Returns:
        Dictionary with optimized parameters
    """
    
    # Load material properties and thresholds
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
    
    # Default wavelength options
    if wavelength_options is None:
        wavelength_options = [355e-9, 532e-9, 800e-9, 1064e-9]
    
    # Default pulse duration options based on material
    if pulse_duration_options is None:
        if material.category == 'dielectric':
            pulse_duration_options = [100e-15, 500e-15, 10e-12]
        else:
            pulse_duration_options = [100e-15, 1e-12, 10e-12, 10e-9]
    
    # Intensity range within clean window
    I_abl = thresholds.I_ablation
    I_dam = thresholds.I_damage
    I_options = np.logspace(np.log10(I_abl * 1.2), np.log10(I_dam * 0.8), 5)
    
    print(f"Grid Search Optimization for {material_name}")
    print(f"Target Quality: Q ≥ {quality_target}")
    print(f"Exploring {len(wavelength_options)} wavelengths × {len(pulse_duration_options)} durations × {len(I_options)} intensities")
    print(f"Total combinations: {len(wavelength_options) * len(pulse_duration_options) * len(I_options)}")
    
    best_result = None
    best_Q = 0
    all_results = []
    
    # Grid search
    for wavelength in wavelength_options:
        for tau in pulse_duration_options:
            for intensity in I_options:
                
                assessment = calculate_quality_total(
                    wavelength=wavelength,
                    pulse_duration=tau,
                    intensity=intensity,
                    material_name=material_name,
                    material_props=material_props,
                    thresholds=thresholds.__dict__,
                    feature_size=feature_size
                )
                
                Q = assessment['Q_total']
                
                result = {
                    'wavelength': wavelength,
                    'pulse_duration': tau,
                    'intensity': intensity,
                    'Q': Q,
                    'assessment': assessment
                }
                
                all_results.append(result)
                
                if Q > best_Q:
                    best_Q = Q
                    best_result = result
    
    # Sort results by quality
    all_results.sort(key=lambda x: x['Q'], reverse=True)
    
    optimization_result = {
        'success': best_Q >= quality_target,
        'best_parameters': {
            'wavelength': best_result['wavelength'],
            'pulse_duration': best_result['pulse_duration'],
            'intensity': best_result['intensity'],
        },
        'best_Q': best_Q,
        'target_Q': quality_target,
        'full_assessment': best_result['assessment'],
        'top_5_results': all_results[:5],
        'all_results': all_results,
    }
    
    if best_Q >= quality_target:
        optimization_result['message'] = f"✓ Target achieved: Q = {best_Q:.3f} ≥ {quality_target}"
    else:
        optimization_result['message'] = f"✗ Target not achieved: Q = {best_Q:.3f} < {quality_target}"
        optimization_result['recommendation'] = "Consider upgrading laser system or adjusting quality target"
    
    return optimization_result


def optimize_parameters_cost_constrained(material_name: str,
                                         quality_target: float,
                                         budget_level: str = 'moderate') -> Dict:
    """
    Optimize parameters within cost constraints
    
    Args:
        material_name: Material to process
        quality_target: Target quality
        budget_level: 'low', 'moderate', 'high'
    
    Returns:
        Dictionary with cost-optimized parameters
    """
    
    # Define cost-constrained options
    if budget_level == 'low':
        wavelength_options = [532e-9, 1064e-9]  # Frequency-doubled, fundamental
        pulse_duration_options = [10e-9, 20e-9]  # Nanosecond only
        laser_cost = 50000  # USD
        
    elif budget_level == 'moderate':
        wavelength_options = [355e-9, 532e-9, 1064e-9]
        pulse_duration_options = [1e-12, 10e-12, 100e-12]  # Picosecond
        laser_cost = 150000
        
    else:  # high
        wavelength_options = [343e-9, 515e-9, 800e-9, 1030e-9]
        pulse_duration_options = [100e-15, 300e-15, 500e-15, 1e-12]  # Femtosecond
        laser_cost = 300000
    
    result = optimize_parameters_grid_search(
        material_name=material_name,
        quality_target=quality_target,
        wavelength_options=wavelength_options,
        pulse_duration_options=pulse_duration_options
    )
    
    result['budget_level'] = budget_level
    result['estimated_laser_cost'] = laser_cost
    result['cost_per_quality_point'] = laser_cost / result['best_Q'] if result['best_Q'] > 0 else float('inf')
    
    return result


def optimize_for_application(application: str,
                             material_name: str,
                             feature_size: float,
                             throughput_target: Optional[float] = None) -> Dict:
    """
    Optimize parameters for specific application
    
    Args:
        application: Application type (e.g., 'microvia', 'cutting', 'marking', 'optical')
        material_name: Material name
        feature_size: Feature size [m]
        throughput_target: Target throughput [features/s] (optional)
    
    Returns:
        Dictionary with application-optimized parameters
    """
    
    # Application-specific requirements
    application_specs = {
        'microvia': {
            'quality_target': 0.70,
            'diameter': 20e-6,
            'depth': 10e-6,
            'throughput_importance': 'high',
        },
        'cutting': {
            'quality_target': 0.65,
            'feature_size': 50e-6,
            'throughput_importance': 'high',
        },
        'marking': {
            'quality_target': 0.60,
            'feature_size': 100e-6,
            'throughput_importance': 'very_high',
        },
        'optical': {
            'quality_target': 0.90,
            'feature_size': 10e-6,
            'throughput_importance': 'low',
        },
        'medical': {
            'quality_target': 0.95,
            'feature_size': 50e-6,
            'throughput_importance': 'low',
        },
        'microfluidics': {
            'quality_target': 0.90,
            'feature_size': 20e-6,
            'throughput_importance': 'moderate',
        },
    }
    
    if application not in application_specs:
        raise ValueError(f"Unknown application: {application}. Available: {list(application_specs.keys())}")
    
    specs = application_specs[application]
    quality_target = specs['quality_target']
    
    # Determine budget level based on quality requirement
    if quality_target >= 0.90:
        budget_level = 'high'
    elif quality_target >= 0.75:
        budget_level = 'moderate'
    else:
        budget_level = 'low'
    
    result = optimize_parameters_cost_constrained(
        material_name=material_name,
        quality_target=quality_target,
        budget_level=budget_level
    )
    
    result['application'] = application
    result['application_specs'] = specs
    
    # Throughput estimation
    if result['success']:
        params = result['best_parameters']
        
        # Estimate ablation depth per pulse (simplified)
        material = get_material(material_name)
        thresholds = get_thresholds(material_name)
        
        from physics_calculations import calculate_fluence
        fluence = calculate_fluence(params['intensity'], params['pulse_duration'])
        
        # Rough depth per pulse (logarithmic model)
        F_th = thresholds.F_ablation_single_pulse if thresholds.F_ablation_single_pulse else 1.0
        alpha = 1e5  # Rough estimate
        
        if fluence > F_th:
            depth_per_pulse = (1 / alpha) * np.log(fluence / F_th) * 1e4  # μm
        else:
            depth_per_pulse = 0
        
        # Pulses needed
        target_depth = specs.get('depth', feature_size)
        if depth_per_pulse > 0:
            pulses_needed = int(np.ceil(target_depth / (depth_per_pulse * 1e-6)))
        else:
            pulses_needed = 100  # Default estimate
        
        # Throughput calculation
        rep_rate = 50000  # 50 kHz typical
        features_per_second = rep_rate / pulses_needed
        
        result['throughput_estimate'] = {
            'depth_per_pulse_um': depth_per_pulse,
            'pulses_needed': pulses_needed,
            'rep_rate_Hz': rep_rate,
            'features_per_second': features_per_second,
            'meets_target': throughput_target is None or features_per_second >= throughput_target
        }
    
    return result


def suggest_laser_system(material_name: str,
                        quality_target: float,
                        budget: float,
                        throughput_requirement: Optional[float] = None) -> Dict:
    """
    Suggest appropriate laser system based on requirements
    
    Args:
        material_name: Material to process
        quality_target: Target quality [0-1]
        budget: Budget in USD
        throughput_requirement: Throughput in features/s (optional)
    
    Returns:
        Dictionary with laser system recommendations
    """
    
    # Laser system database
    laser_systems = {
        'fs_tisapphire': {
            'name': 'Femtosecond Ti:Sapphire',
            'cost': 800000,
            'pulse_duration': 100e-15,
            'rep_rate': 1000,
            'wavelength': 800e-9,
            'max_quality': 0.95,
            'pulse_energy': 1e-3,
            'notes': 'Ultimate quality, research-grade'
        },
        'fs_fiber': {
            'name': 'Femtosecond Fiber (Yb)',
            'cost': 300000,
            'pulse_duration': 300e-15,
            'rep_rate': 1e6,
            'wavelength': 1030e-9,
            'max_quality': 0.92,
            'pulse_energy': 50e-6,
            'notes': 'High reliability, industrial'
        },
        'ps_dpss': {
            'name': 'Picosecond DPSS',
            'cost': 150000,
            'pulse_duration': 10e-12,
            'rep_rate': 100e3,
            'wavelength': 1064e-9,
            'max_quality': 0.85,
            'pulse_energy': 100e-6,
            'notes': 'Balanced cost/quality'
        },
        'ps_fiber': {
            'name': 'Picosecond Fiber',
            'cost': 180000,
            'pulse_duration': 10e-12,
            'rep_rate': 500e3,
            'wavelength': 1030e-9,
            'max_quality': 0.82,
            'pulse_energy': 50e-6,
            'notes': 'High throughput'
        },
        'ns_dpss_uv': {
            'name': 'Nanosecond DPSS (UV)',
            'cost': 50000,
            'pulse_duration': 10e-9,
            'rep_rate': 50e3,
            'wavelength': 355e-9,
            'max_quality': 0.75,
            'pulse_energy': 500e-6,
            'notes': 'Cost-effective, polymers'
        },
        'ns_fiber': {
            'name': 'Nanosecond Fiber',
            'cost': 70000,
            'pulse_duration': 20e-9,
            'rep_rate': 100e3,
            'wavelength': 1064e-9,
            'max_quality': 0.70,
            'pulse_energy': 1e-3,
            'notes': 'High power, low cost'
        },
    }
    
    # Material-specific constraints
    material = get_material(material_name)
    
    # Filter systems by constraints
    suitable_systems = []
    
    for system_id, system in laser_systems.items():
        # Check budget
        if system['cost'] > budget:
            continue
        
        # Check quality capability
        if system['max_quality'] < quality_target:
            continue
        
        # Material-specific checks
        if material.category == 'dielectric':
            # Glass requires femtosecond
            if system['pulse_duration'] > 10e-12:
                continue
        
        # Throughput check (if specified)
        if throughput_requirement:
            estimated_throughput = system['rep_rate'] / 50  # Rough estimate
            if estimated_throughput < throughput_requirement:
                continue
        
        suitable_systems.append({
            'system_id': system_id,
            'system': system,
            'cost_efficiency': system['max_quality'] / system['cost'] * 1e6,
        })
    
    # Sort by cost efficiency
    suitable_systems.sort(key=lambda x: x['cost_efficiency'], reverse=True)
    
    result = {
        'material': material_name,
        'quality_target': quality_target,
        'budget': budget,
        'suitable_systems': suitable_systems,
    }
    
    if suitable_systems:
        result['recommended_system'] = suitable_systems[0]['system']
        result['success'] = True
        result['message'] = f"Found {len(suitable_systems)} suitable system(s)"
    else:
        result['success'] = False
        result['message'] = "No suitable system within constraints"
        result['recommendations'] = []
        
        # Provide recommendations
        if quality_target > 0.85:
            result['recommendations'].append(
                f"Quality target Q={quality_target} requires femtosecond laser (budget ~$300k+)"
            )
        
        if budget < 150000:
            result['recommendations'].append(
                f"Budget ${budget/1000:.0f}k limits to nanosecond systems (Q < 0.75)"
            )
        
        if material.category == 'dielectric':
            result['recommendations'].append(
                f"{material_name} requires femtosecond pulses (not achievable with current budget)"
            )
    
    return result


def multi_objective_optimization(material_name: str,
                                 quality_weight: float = 0.5,
                                 cost_weight: float = 0.3,
                                 throughput_weight: float = 0.2) -> Dict:
    """
    Multi-objective optimization balancing quality, cost, and throughput
    
    Args:
        material_name: Material to process
        quality_weight: Weight for quality (0-1)
        cost_weight: Weight for cost (0-1)
        throughput_weight: Weight for throughput (0-1)
    
    Returns:
        Pareto-optimal solutions
    """
    
    # Normalize weights
    total_weight = quality_weight + cost_weight + throughput_weight
    quality_weight /= total_weight
    cost_weight /= total_weight
    throughput_weight /= total_weight
    
    # Generate candidate solutions
    budget_levels = ['low', 'moderate', 'high']
    quality_targets = [0.65, 0.75, 0.85, 0.95]
    
    solutions = []
    
    for budget_level in budget_levels:
        for q_target in quality_targets:
            try:
                result = optimize_parameters_cost_constrained(
                    material_name, q_target, budget_level
                )
                
                if result['success']:
                    # Calculate multi-objective score
                    Q_norm = result['best_Q']
                    cost_norm = 1 - (result['estimated_laser_cost'] / 800000)  # Normalize by max cost
                    throughput_norm = 0.5  # Placeholder (would need actual calculation)
                    
                    score = (quality_weight * Q_norm + 
                            cost_weight * cost_norm + 
                            throughput_weight * throughput_norm)
                    
                    solutions.append({
                        'parameters': result['best_parameters'],
                        'quality': result['best_Q'],
                        'cost': result['estimated_laser_cost'],
                        'budget_level': budget_level,
                        'score': score,
                        'full_result': result
                    })
            except:
                pass
    
    # Sort by score
    solutions.sort(key=lambda x: x['score'], reverse=True)
    
    return {
        'solutions': solutions,
        'best_solution': solutions[0] if solutions else None,
        'pareto_front': solutions[:5],  # Top 5
        'weights': {
            'quality': quality_weight,
            'cost': cost_weight,
            'throughput': throughput_weight
        }
    }


if __name__ == "__main__":
    # Test optimization
    print("="*70)
    print("PARAMETER OPTIMIZATION TEST")
    print("="*70)
    
    # Test 1: Grid search for glass
    print("\nTest 1: Grid Search Optimization (Glass)")
    print("-"*70)
    
    result1 = optimize_parameters_grid_search(
        material_name='Glass',
        quality_target=0.90,
        wavelength_options=[800e-9],
        pulse_duration_options=[100e-15, 500e-15, 10e-12]
    )
    
    print(f"\n{result1['message']}")
    print(f"\nBest Parameters:")
    print(f"  Wavelength: {result1['best_parameters']['wavelength']*1e9:.0f} nm")
    print(f"  Pulse Duration: {result1['best_parameters']['pulse_duration']*1e15:.0f} fs")
    print(f"  Intensity: {result1['best_parameters']['intensity']:.2e} W/cm²")
    print(f"  Achieved Q: {result1['best_Q']:.3f}")
    
    print(f"\nTop 3 Solutions:")
    for i, sol in enumerate(result1['top_5_results'][:3], 1):
        print(f"  {i}. Q={sol['Q']:.3f}: λ={sol['wavelength']*1e9:.0f}nm, "
              f"τ={sol['pulse_duration']*1e12:.1f}ps, I={sol['intensity']:.1e}")
    
    # Test 2: Cost-constrained optimization
    print("\n" + "="*70)
    print("Test 2: Cost-Constrained Optimization (Polymer)")
    print("-"*70)
    
    for budget_level in ['low', 'moderate', 'high']:
        result2 = optimize_parameters_cost_constrained(
            material_name='ABF_Polymer',
            quality_target=0.75,
            budget_level=budget_level
        )
        
        print(f"\nBudget Level: {budget_level.upper()}")
        print(f"  Estimated Cost: ${result2['estimated_laser_cost']/1000:.0f}k")
        print(f"  Achieved Q: {result2['best_Q']:.3f}")
        print(f"  Success: {result2['success']}")
    
    # Test 3: Application-specific optimization
    print("\n" + "="*70)
    print("Test 3: Application-Specific Optimization")
    print("-"*70)
    
    applications = ['microvia', 'optical', 'medical']
    
    for app in applications:
        result3 = optimize_for_application(
            application=app,
            material_name='ABF_Polymer' if app == 'microvia' else 'Glass',
            feature_size=20e-6
        )
        
        print(f"\nApplication: {app.upper()}")
        print(f"  Target Q: {result3['application_specs']['quality_target']}")
        print(f"  Achieved Q: {result3['best_Q']:.3f}")
        print(f"  Budget Level: {result3['budget_level']}")
        if 'throughput_estimate' in result3:
            print(f"  Throughput: {result3['throughput_estimate']['features_per_second']:.0f} features/s")
    
    # Test 4: Laser system recommendation
    print("\n" + "="*70)
    print("Test 4: Laser System Recommendation")
    print("-"*70)
    
    result4 = suggest_laser_system(
        material_name='Glass',
        quality_target=0.90,
        budget=350000
    )
    
    print(f"\nQuery: Glass processing, Q≥0.90, Budget=$350k")
    print(f"Result: {result4['message']}")
    
    if result4['success']:
        system = result4['recommended_system']
        print(f"\nRecommended System: {system['name']}")
        print(f"  Cost: ${system['cost']/1000:.0f}k")
        print(f"  Pulse Duration: {system['pulse_duration']*1e15:.0f} fs")
        print(f"  Max Quality: {system['max_quality']:.2f}")
        print(f"  Notes: {system['notes']}")