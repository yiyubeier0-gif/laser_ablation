
## File 9: `main.py`

"""
main.py

Main entry point for the laser processing framework
Provides interactive command-line interface
"""

import sys
import argparse
from materials_database import list_materials, get_material
from thresholds_database import get_thresholds
from quality_assessment import calculate_quality_total
from optimization import (
    optimize_parameters_grid_search,
    optimize_for_application,
    suggest_laser_system
)
from process_report import generate_process_report


def main():
    """Main entry point"""
    
    parser = argparse.ArgumentParser(
        description='Laser-Material Interaction Analysis Framework',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze specific parameters
  python main.py analyze --material Glass --wavelength 800 --pulse-duration 100 --intensity 5e13
  
  # Optimize parameters
  python main.py optimize --material Glass --quality-target 0.90
  
  # Application-specific optimization
  python main.py application --app microvia --material ABF_Polymer
  
  # Laser system recommendation
  python main.py recommend --material Glass --quality 0.90 --budget 350000
  
  # List available materials
  python main.py list-materials
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze specific parameters')
    analyze_parser.add_argument('--material', required=True, help='Material name')
    analyze_parser.add_argument('--wavelength', type=float, required=True, help='Wavelength [nm]')
    analyze_parser.add_argument('--pulse-duration', type=float, required=True, help='Pulse duration [fs or ps]')
    analyze_parser.add_argument('--intensity', type=float, required=True, help='Intensity [W/cm²]')
    analyze_parser.add_argument('--spot-diameter', type=float, default=20, help='Spot diameter [μm]')
    analyze_parser.add_argument('--feature-size', type=float, default=20, help='Feature size [μm]')
    analyze_parser.add_argument('--target-depth', type=float, default=10, help='Target depth [μm]')
    analyze_parser.add_argument('--output', choices=['text', 'markdown', 'json'], default='text', help='Output format')
    
    # Optimize command
    optimize_parser = subparsers.add_parser('optimize', help='Optimize parameters')
    optimize_parser.add_argument('--material', required=True, help='Material name')
    optimize_parser.add_argument('--quality-target', type=float, default=0.85, help='Target quality [0-1]')
    optimize_parser.add_argument('--budget', choices=['low', 'moderate', 'high'], default='moderate', help='Budget level')
    
    # Application command
    app_parser = subparsers.add_parser('application', help='Application-specific optimization')
    app_parser.add_argument('--app', required=True, 
                           choices=['microvia', 'cutting', 'marking', 'optical', 'medical', 'microfluidics'],
                           help='Application type')
    app_parser.add_argument('--material', required=True, help='Material name')
    app_parser.add_argument('--feature-size', type=float, default=20, help='Feature size [μm]')
    app_parser.add_argument('--throughput', type=float, help='Throughput target [features/s]')
    
    # Recommend command
    recommend_parser = subparsers.add_parser('recommend', help='Recommend laser system')
    recommend_parser.add_argument('--material', required=True, help='Material name')
    recommend_parser.add_argument('--quality', type=float, required=True, help='Quality target [0-1]')
    recommend_parser.add_argument('--budget', type=float, required=True, help='Budget [USD]')
    recommend_parser.add_argument('--throughput', type=float, help='Throughput requirement [features/s]')
    
    # List materials command
    list_parser = subparsers.add_parser('list-materials', help='List available materials')
    list_parser.add_argument('--category', choices=['metal', 'semiconductor', 'dielectric', 'polymer'],
                            help='Filter by category')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
    
    # Execute command
    if args.command == 'analyze':
        cmd_analyze(args)
    elif args.command == 'optimize':
        cmd_optimize(args)
    elif args.command == 'application':
        cmd_application(args)
    elif args.command == 'recommend':
        cmd_recommend(args)
    elif args.command == 'list-materials':
        cmd_list_materials(args)


def cmd_analyze(args):
    """Execute analyze command"""
    
    # Convert units
    wavelength = args.wavelength * 1e-9  # nm to m
    
    # Determine pulse duration unit (fs or ps)
    if args.pulse_duration < 1000:
        pulse_duration = args.pulse_duration * 1e-15  # fs to s
    else:
        pulse_duration = args.pulse_duration * 1e-12  # ps to s
    
    spot_diameter = args.spot_diameter * 1e-6  # μm to m
    feature_size = args.feature_size * 1e-6
    target_depth = args.target_depth * 1e-6
    
    # Generate report
    report = generate_process_report(
        material_name=args.material,
        wavelength=wavelength,
        pulse_duration=pulse_duration,
        intensity=args.intensity,
        feature_size=feature_size,
        target_depth=target_depth,
        spot_diameter=spot_diameter,
        output_format=args.output
    )
    
    print(report)


def cmd_optimize(args):
    """Execute optimize command"""
    
    from optimization import optimize_parameters_cost_constrained
    
    result = optimize_parameters_cost_constrained(
        material_name=args.material,
        quality_target=args.quality_target,
        budget_level=args.budget
    )
    
    print("\n" + "="*80)
    print("PARAMETER OPTIMIZATION RESULTS")
    print("="*80)
    print(f"\n{result['message']}")
    
    if result['success']:
        params = result['best_parameters']
        print(f"\nOptimized Parameters:")
        print(f"  Wavelength:      {params['wavelength']*1e9:.0f} nm")
        print(f"  Pulse Duration:  {params['pulse_duration']*1e12:.2f} ps")
        print(f"  Intensity:       {params['intensity']:.2e} W/cm²")
        print(f"\nAchieved Quality: Q = {result['best_Q']:.3f}")
        print(f"Budget Level:     {result['budget_level']}")
        print(f"Estimated Cost:   ${result['estimated_laser_cost']/1000:.0f}k")
    else:
        if 'recommendation' in result:
            print(f"\nRecommendation: {result['recommendation']}")


def cmd_application(args):
    """Execute application command"""
    
    result = optimize_for_application(
        application=args.app,
        material_name=args.material,
        feature_size=args.feature_size * 1e-6,
        throughput_target=args.throughput
    )
    
    print("\n" + "="*80)
    print(f"APPLICATION-SPECIFIC OPTIMIZATION: {args.app.upper()}")
    print("="*80)
    
    print(f"\nApplication Requirements:")
    print(f"  Quality Target:  {result['application_specs']['quality_target']}")
    
    if result['success']:
        params = result['best_parameters']
        print(f"\nRecommended Parameters:")
        print(f"  Wavelength:      {params['wavelength']*1e9:.0f} nm")
        print(f"  Pulse Duration:  {params['pulse_duration']*1e12:.2f} ps")
        print(f"  Intensity:       {params['intensity']:.2e} W/cm²")
        print(f"\nAchieved Quality: Q = {result['best_Q']:.3f}")
        print(f"Budget Level:     {result['budget_level']}")
        
        if 'throughput_estimate' in result:
            throughput = result['throughput_estimate']
            print(f"\nThroughput Estimate:")
            print(f"  Features/second: {throughput['features_per_second']:.0f}")
            print(f"  Pulses needed:   {throughput['pulses_needed']}")
            print(f"  Target met:      {'Yes' if throughput['meets_target'] else 'No'}")


def cmd_recommend(args):
    """Execute recommend command"""
    
    result = suggest_laser_system(
        material_name=args.material,
        quality_target=args.quality,
        budget=args.budget,
        throughput_requirement=args.throughput
    )
    
    print("\n" + "="*80)
    print("LASER SYSTEM RECOMMENDATION")
    print("="*80)
    
    print(f"\nRequirements:")
    print(f"  Material:        {args.material}")
    print(f"  Quality Target:  Q ≥ {args.quality}")
    print(f"  Budget:          ${args.budget/1000:.0f}k")
    if args.throughput:
        print(f"  Throughput:      {args.throughput} features/s")
    
    print(f"\n{result['message']}")
    
    if result['success']:
        system = result['recommended_system']
        print(f"\nRecommended System: {system['name']}")
        print(f"  Cost:            ${system['cost']/1000:.0f}k")
        print(f"  Pulse Duration:  {system['pulse_duration']*1e15:.0f} fs")
        print(f"  Wavelength:      {system['wavelength']*1e9:.0f} nm")
        print(f"  Max Quality:     {system['max_quality']:.2f}")
        print(f"  Repetition Rate: {system['rep_rate']/1e3:.0f} kHz")
        print(f"  Notes:           {system['notes']}")
        
        if len(result['suitable_systems']) > 1:
            print(f"\nAlternative Systems ({len(result['suitable_systems'])-1}):")
            for alt in result['suitable_systems'][1:3]:
                s = alt['system']
                print(f"  • {s['name']} (${s['cost']/1000:.0f}k, Q_max={s['max_quality']:.2f})")
    else:
        print(f"\nRecommendations:")
        for rec in result.get('recommendations', []):
            print(f"  • {rec}")


def cmd_list_materials(args):
    """Execute list-materials command"""
    
    materials = list_materials(args.category)
    
    if args.category:
        print(f"\n{args.category.upper()} MATERIALS:")
    else:
        print("\nALL AVAILABLE MATERIALS:")
    
    print("-" * 80)
    
    for mat_name in materials:
        try:
            material = get_material(mat_name)
            thresholds = get_thresholds(mat_name)
            print(f"\n{mat_name}:")
            print(f"  Category:        {material.category}")
            print(f"  Bandgap:         {material.bandgap:.2f} eV")
            print(f"  I_ablation:      {thresholds.I_ablation:.1e} W/cm²")
            print(f"  Recommended λ:   {thresholds.recommended_wavelength}")
            if material.notes:
                print(f"  Notes:           {material.notes[:60]}...")
        except:
            print(f"\n{mat_name}: (No thresholds available)")


if __name__ == "__main__":
    main()

