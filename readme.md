
```markdown
# Laser-Material Interaction Framework

A comprehensive, physics-based framework for laser processing analysis, optimization, and prediction.

## Overview

This framework provides:
- **Rigorous regime classification** (photothermal, stress-confined, plasma-mediated)
- **Quality assessment** with quantitative Q metric
- **Parameter optimization** for specific applications
- **Material-specific databases** with validated thresholds
- **Process reports** with detailed analysis

## Installation

```bash
# Clone repository
git clone <[repository-url](https://github.com/yiyubeier0-gif/laser_ablation.git)>
cd laser-processing-framework

# Install dependencies
pip install numpy scipy
```

## Quick Start

### 1. Basic Quality Assessment

```python
from materials_database import get_material
from thresholds_database import get_thresholds
from quality_assessment import calculate_quality_total

# Load material
material = get_material('Glass')
thresholds = get_thresholds('Glass')

# Define parameters
wavelength = 800e-9  # 800 nm
pulse_duration = 100e-15  # 100 fs
intensity = 5e13  # 50 TW/cm²

# Calculate quality
result = calculate_quality_total(
    wavelength=wavelength,
    pulse_duration=pulse_duration,
    intensity=intensity,
    material_name='Glass',
    material_props={...},  # See documentation
    thresholds=thresholds.__dict__
)

print(f"Quality: Q = {result['Q_total']:.3f}")
print(f"Regime: {result['regime']}")
```

### 2. Parameter Optimization

```python
from optimization import optimize_parameters_grid_search

# Optimize for target quality
result = optimize_parameters_grid_search(
    material_name='Glass',
    quality_target=0.90
)

print(f"Best Parameters:")
print(f"  Wavelength: {result['best_parameters']['wavelength']*1e9:.0f} nm")
print(f"  Pulse Duration: {result['best_parameters']['pulse_duration']*1e15:.0f} fs")
print(f"  Intensity: {result['best_parameters']['intensity']:.2e} W/cm²")
```

### 3. Application-Specific Optimization

```python
from optimization import optimize_for_application

# Optimize for microvia drilling
result = optimize_for_application(
    application='microvia',
    material_name='ABF_Polymer',
    feature_size=20e-6
)

print(f"Recommended system: {result['budget_level']}")
print(f"Quality: Q = {result['best_Q']:.3f}")
```

### 4. Generate Process Report

```python
from process_report import generate_process_report

report = generate_process_report(
    material_name='Glass',
    wavelength=800e-9,
    pulse_duration=100e-15,
    intensity=5e13,
    output_format='text'
)

print(report)
```

## Command-Line Interface

```bash
# Analyze specific parameters
python main.py analyze --material Glass --wavelength 800 --pulse-duration 100 --intensity 5e13

# Optimize parameters
python main.py optimize --material Glass --quality-target 0.90 --budget moderate

# Application-specific optimization
python main.py application --app microvia --material ABF_Polymer

# Laser system recommendation
python main.py recommend --material Glass --quality 0.90 --budget 350000

# List available materials
python main.py list-materials
```

## Framework Architecture

### Core Modules

1. **materials_database.py** - Material properties database
   - Electronic properties (bandgap, refractive index)
   - Thermal properties (conductivity, diffusivity)
   - Mechanical properties (sound velocity, Young's modulus)

2. **thresholds_database.py** - Processing thresholds
   - I_modification, I_ablation, I_damage, I_catastrophic
   - Temporal constraints (τ_critical, τ_optimal)
   - Clean ablation windows

3. **physics_calculations.py** - Core physics
   - Intensity/fluence calculations
   - Plasma parameters (critical density, Keldysh parameter)
   - Thermal diffusion, ablation depth models
   - Nonlinear optics (self-focusing)

4. **regime_classification.py** - Regime analysis
   - Three-regime framework classification
   - Temporal hierarchy analysis
   - Mechanism identification

5. **quality_assessment.py** - Quality metrics
   - Comprehensive Q calculation
   - Component breakdown (Q_regime, Q_intensity, Q_temporal, Q_wavelength)
   - Recommendations generation

6. **optimization.py** - Parameter optimization
   - Grid search optimization
   - Cost-constrained optimization
   - Application-specific optimization
   - Laser system recommendation

7. **process_report.py** - Report generation
   - Text, Markdown, JSON formats
   - Comprehensive analysis
   - Comparison reports

## Key Concepts

### Three-Regime Framework

```
1. PHOTOTHERMAL (τ > τ_thermal)
   - Linear absorption, heat diffusion
   - Large HAZ (10-100 μm)
   - Q = 0.3-0.6

2. STRESS-CONFINED (τ_stress < τ < τ_thermal)
   - Thermomechanical spallation
   - Reduced HAZ (1-10 μm)
   - Q = 0.6-0.8

3. PLASMA-MEDIATED (τ < τ_e-phonon)
   - Athermal ablation
   - Minimal HAZ (<100 nm)
   - Q = 0.85-0.95
```

### Four-Level Intensity Thresholds

```
I_modification < I_ablation < I_damage < I_catastrophic

Clean Ablation Window: [I_ablation, I_damage]
  → Material removed WITHOUT damage
  → Target processing regime
```

### Quality Factor Q

```
Q = Q_regime × Q_intensity × Q_temporal × Q_wavelength × Q_accumulation

Q = 0.9-1.0:  Excellent (optimal)
Q = 0.7-0.9:  Good (acceptable)
Q = 0.5-0.7:  Fair (marginal)
Q < 0.5:      Poor (problematic)
```

## Material Database

Currently supported materials:

**Dielectrics:**
- Fused Silica (Glass)
- BK7 Glass
- Sapphire

**Polymers:**
- ABF (Ajinomoto Build-up Film)
- Polyimide (Kapton)
- PMMA

**Metals:**
- Copper
- Aluminum
- Stainless Steel
- Titanium
- Gold

**Semiconductors:**
- Silicon
- GaAs
- GaN
- SiC

## Examples

See `examples.py` for comprehensive demonstrations:

1. Basic quality assessment
2. Parameter optimization
3. Application-specific optimization
4. Laser system selection
5. Process report generation
6. Parameter comparison
7. Regime boundary exploration
8. Material survey
9. Cost-quality trade-offs
10. Troubleshooting guide

Run all examples:
```bash
python examples.py
```

## Critical Material-Specific Guidelines

### Glass/Dielectrics
- **τ < 10 ps MANDATORY** (thermal stress cracking above)
- Femtosecond strongly recommended
- I = 10¹³ - 10¹⁴ W/cm² (clean window)
- λ = 800 nm standard

### Polymers
- UV wavelength (355 nm) optimal
- Nanosecond acceptable for cost
- N₂ atmosphere prevents carbonization
- I = 10⁷ - 10⁸ W/cm²

### Metals
- Green (532 nm) or UV (355 nm) preferred
- Avoid NIR (high reflectivity)
- τ < 1 ps for melt-free
- I = 10¹² - 10¹³ W/cm²

### Semiconductors
- λ near/above bandgap for linear absorption
- τ < 10 ps to avoid amorphization
- I = 10¹¹ - 10¹² W/cm²

## Validation and Testing

Run test suite:
```bash
python -m pytest tests/
```

Individual module tests:
```bash
python materials_database.py
python thresholds_database.py
python physics_calculations.py
python quality_assessment.py
```

## API Documentation

### calculate_quality_total()

```python
def calculate_quality_total(
    wavelength: float,          # [m]
    pulse_duration: float,      # [s]
    intensity: float,           # [W/cm²]
    material_name: str,
    material_props: dict,
    thresholds: dict,
    feature_size: float = 10e-6,  # [m]
    num_pulses: int = 1,
    rep_rate: float = 1000       # [Hz]
) -> Dict
```

Returns dictionary with:
- `Q_total`: Overall quality [0-1]
- `components`: Component breakdown
- `regime`: Processing regime
- `warnings`: List of warnings
- `recommendations`: Actionable recommendations
- `interpretation`: Quality interpretation

### optimize_parameters_grid_search()

```python
def optimize_parameters_grid_search(
    material_name: str,
    quality_target: float = 0.85,
    wavelength_options: Optional[List[float]] = None,
    pulse_duration_options: Optional[List[float]] = None,
    feature_size: float = 10e-6
) -> Dict
```

Returns dictionary with:
- `success`: Whether target achieved
- `best_parameters`: Optimized parameters
- `best_Q`: Achieved quality
- `top_5_results`: Top 5 parameter sets

### generate_process_report()

```python
def generate_process_report(
    material_name: str,
    wavelength: float,
    pulse_duration: float,
    intensity: float,
    output_format: str = 'text'  # 'text', 'markdown', 'json'
) -> str
```

Returns formatted report string.

## Contributing

To add a new material:

1. Add entry to `materials_database.py`:
```python
MATERIALS_DB['NewMaterial'] = MaterialProperties(
    name='New Material',
    category='metal',  # or 'semiconductor', 'dielectric', 'polymer'
    bandgap=0.0,
    # ... other properties
)
```

2. Add thresholds to `thresholds_database.py`:
```python
THRESHOLDS_DB['NewMaterial'] = ProcessingThresholds(
    material_name='NewMaterial',
    I_modification=1e6,
    I_ablation=1e7,
    # ... other thresholds
)
```

3. Validate with experimental data

## References

Key papers incorporated into framework:

1. Chichkov et al. (1996) - "Femtosecond, picosecond and nanosecond laser ablation of solids"
2. Stuart et al. (1996) - "Optical breakdown in dielectrics with nanosecond to sub-picosecond pulses"
3. Gattass & Mazur (2008) - "Femtosecond laser micromachining in transparent materials"
4. Bulgakova et al. (2010) - "A general continuum approach to describe fast electronic transport in pulsed laser irradiated materials"

## License

[Specify license]

## Contact

[Contact information]

## Changelog

### Version 1.0 (2025-10-06)
- Initial release
- 14 materials supported
- Complete quality assessment framework
- Parameter optimization
- Process reporting
- Command-line interface

## Known Limitations

1. **Simplified thermal models**: Use 1D heat diffusion (adequate for most cases)
2. **Plasma dynamics**: Simplified rate equations (not full PIC simulation)
3. **Multi-pulse accumulation**: Empirical models (material-dependent validation needed)
4. **Wavelength dependence**: Limited data for some materials at all wavelengths
5. **Nonlinear effects**: Self-focusing approximations (Marburger formula)

## Future Development

Planned features:
- [ ] Time-resolved dynamics visualization
- [ ] 3D thermal simulation integration
- [ ] Machine learning parameter optimization
- [ ] Expanded material database (ceramics, composites)
- [ ] Experimental validation database
- [ ] GUI interface
- [ ] Real-time process monitoring integration

## Troubleshooting

### Common Issues

**Issue: Low quality (Q < 0.5)**
- Check temporal constraint (τ vs τ_critical)
- Verify intensity within clean window
- Consider wavelength absorption

**Issue: "Material not found"**
- Check spelling (case-sensitive)
- Use `list-materials` to see available materials
- Add material if not in database

**Issue: "No clean window for nanosecond"**
- Expected for glass/brittle materials
- Must use femtosecond for these materials
- See material-specific guidelines

**Issue: Self-focusing warning**
- Reduce pulse energy
- Increase spot size
- Use shorter wavelength

## Performance

Typical calculation times (on standard laptop):
- Single quality assessment: <0.1 s
- Grid search (50 points): ~5 s
- Full optimization: ~10-30 s
- Report generation: <1 s

## Citation

If you use this framework in research, please cite:

```
[Citation format to be added]
```

---

**Last Updated:** 2025-10-06  
**Framework Version:** 1.0  
**Python Version:** 3.7+  
**Dependencies:** numpy, scipy
```

---
