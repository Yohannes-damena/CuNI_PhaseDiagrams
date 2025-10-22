# Cu-Ni Binary Phase Diagram Analysis Project

A comprehensive computational materials science project for analyzing the Copper-Nickel (Cu-Ni) binary phase diagram using CALPHAD (CALculation of PHAse Diagrams) methodology. This project demonstrates complete solid solubility behavior, phase transformation kinetics, and practical alloy design applications.

## 📋 Table of Contents

- [Project Overview](#project-overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation & Setup](#installation--setup)
- [Project Structure](#project-structure)
- [Usage](#usage)
- [Diagrams Generated](#diagrams-generated)
- [Packages Used](#packages-used)
- [Scientific Background](#scientific-background)
- [Results & Analysis](#results--analysis)
- [Applications](#applications)
- [Contributing](#contributing)

## 🔬 Project Overview

This project provides a complete computational analysis of the Cu-Ni binary system, demonstrating:

- **Complete solid solubility** between copper and nickel
- **Phase diagram visualization** with liquidus and solidus curves
- **Equilibrium phase fraction calculations** for specific alloy compositions
- **Thermodynamic analysis** using CALPHAD databases
- **Practical applications** in alloy design and processing

The analysis focuses on the Cu-40Ni composition, which represents a typical commercial alloy with excellent corrosion resistance and mechanical properties.

## ✨ Features

- 🧪 **Binary Phase Diagram Generation**: High-quality phase diagrams with detailed annotations
- 📊 **Phase Fraction Analysis**: Temperature-dependent phase evolution calculations
- 🔬 **Thermodynamic Modeling**: CALPHAD-based equilibrium calculations
- 📈 **Data Visualization**: Professional scientific plots with comprehensive labeling
- 📚 **Educational Content**: Detailed explanations of materials science concepts
- 🛠️ **Practical Applications**: Real-world alloy design and processing insights

## 🔧 Prerequisites

- Python 3.11 or higher
- Conda package manager
- Basic understanding of materials science and phase diagrams

## 🚀 Installation & Setup

### Method 1: Using Conda Environment (Recommended)

1. **Clone the repository**:
   ```bash
   git clone <https://github.com/Yohannes-damena/CuNI_PhaseDiagrams>
   cd "Engineering Material Assignment"
   ```

2. **Create and activate the conda environment**:
   ```bash
   cd CuNi_PhaseDiagram
   conda env create -f environment.yml
   conda activate cu_ni_env
   ```

3. **Verify installation**:
   ```bash
   python -c "import pycalphad, matplotlib, numpy; print('All packages installed successfully!')"
   ```

### Method 2: Manual Installation

1. **Create a new conda environment**:
   ```bash
   conda create -n cu_ni_env python=3.11
   conda activate cu_ni_env
   ```

2. **Install required packages**:
   ```bash
   conda install -c conda-forge pycalphad matplotlib numpy scipy jupyter
   ```

3. **Navigate to the project directory**:
   ```bash
   cd CuNi_PhaseDiagram
   ```

## 📁 Project Structure

```
Engineering Material Assignment/
├── README.md                           # This file
├── plot_CuNi.py                       # Alternative plotting script
└── CuNi_PhaseDiagram/
    ├── CuNi_Binary_PhaseDiagram.py    # Main analysis script
    ├── CuNi_Project_Guide.md          # Project guide (empty)
    ├── environment.yml                # Conda environment specification
    └── data/
        ├── CuNi_db.tdb               # CALPHAD thermodynamic database
        ├── 1_solidification_temperatures.txt
        ├── 2_phase_transformation_path.txt
        ├── 3_complete_solid_solubility.txt
        └── 4_alloy_design_applications.txt
```

## 🎯 Usage

### Running the Main Analysis

1. **Activate the environment**:
   ```bash
   conda activate cu_ni_env
   ```

2. **Run the main analysis script**:
   ```bash
   python CuNi_Binary_PhaseDiagram.py
   ```

This will generate:
- Binary phase diagram plot
- Phase fraction analysis for Cu-40Ni alloy
- Detailed temperature calculations
- Professional scientific visualizations

### Running the Alternative Script

```bash
python plot_CuNi.py
```

This script provides additional analysis including:
- Energy surface calculations
- Gibbs energy analysis
- Equilibrium composition calculations

## 📊 Diagrams Generated

### 1. Binary Phase Diagram
- **X-axis**: Mole fraction of Ni (0 to 1)
- **Y-axis**: Temperature (400K to 1800K)
- **Features**:
  - Liquidus and solidus curves
  - Phase region labels
  - Melting point annotations
  - Professional scientific formatting

### 2. Phase Fraction Analysis
- **Composition**: Cu-40Ni (X(Ni) = 0.40)
- **Temperature Range**: 400K to 1800K
- **Features**:
  - Liquid and solid phase fractions vs. temperature
  - Solidus and liquidus temperature identification
  - Melting range calculation
  - Two-phase region visualization

### 3. Energy Surface Analysis (plot_CuNi.py)
- **Gibbs Energy Surfaces**: All phases at 1550K
- **Equilibrium Compositions**: Calculated using Newton's method
- **Phase Stability**: Visual comparison of phase energies

## 📦 Packages Used

### Core Scientific Computing
- **Python 3.11**: Programming language
- **NumPy**: Numerical computations and array operations
- **SciPy**: Scientific computing and optimization algorithms

### CALPHAD and Thermodynamics
- **pycalphad**: CALPHAD methodology implementation
  - Database management
  - Equilibrium calculations
  - Phase diagram plotting
  - Thermodynamic property calculations

### Visualization
- **Matplotlib**: Scientific plotting and visualization
  - High-quality figure generation
  - Professional formatting
  - Publication-ready plots

### Development Environment
- **Jupyter**: Interactive development and analysis
- **Conda**: Package and environment management

## 🔬 Scientific Background

### Why Cu-Ni Exhibits Complete Solid Solubility

1. **Crystal Structure Compatibility**:
   - Both Cu and Ni have FCC crystal structure
   - Identical atomic packing arrangement
   - Similar atomic radii (Cu: 1.28 Å, Ni: 1.24 Å)

2. **Electronic Structure Similarity**:
   - Both are transition metals
   - Similar electronegativities (Cu: 1.90, Ni: 1.91)
   - Compatible electronic configurations

3. **Thermodynamic Factors**:
   - Negative enthalpy of mixing (exothermic)
   - Favorable atomic interactions
   - Low strain energy due to similar atomic sizes

### CALPHAD Methodology

CALPHAD (CALculation of PHAse Diagrams) is a computational approach that:
- Uses thermodynamic databases to predict phase equilibria
- Combines experimental data with theoretical models
- Enables prediction of phase diagrams for complex systems
- Provides accurate property calculations for alloy design

## 📈 Results & Analysis

### Key Findings for Cu-40Ni Alloy

- **Solidus Temperature**: 1231.9°C (1504.9 K)
- **Liquidus Temperature**: 1526.8°C (1799.8 K)
- **Melting Range**: 295.0°C
- **Phase Behavior**:
  - Below 1231.9°C: 100% FCC-A1 solid solution
  - 1231.9°C to 1526.8°C: Two-phase region (Liquid + FCC-A1)
  - Above 1526.8°C: 100% Liquid phase

### Scientific Significance

- Demonstrates complete solid solubility principles
- Shows importance of atomic size and structure compatibility
- Illustrates thermodynamic principles in alloy design
- Provides model system for understanding solid solutions

## 🏭 Applications

### Industrial Applications
- **Marine Engineering**: Seawater corrosion resistance
- **Chemical Processing**: Acid resistance and high-temperature stability
- **Electrical Engineering**: Good conductivity with corrosion resistance
- **Aerospace**: Lightweight design with excellent properties

### Processing Applications
- **Casting**: Wide melting range enables good castability
- **Welding**: Gradual solidification reduces cracking
- **Heat Treatment**: Homogenization and property optimization
- **Powder Metallurgy**: Sintering in two-phase region

### Alloy Design
- **Property Tailoring**: Continuous property variation with composition
- **Cost Optimization**: Balance between performance and economics
- **Quality Control**: Predictable behavior and consistent properties

## 🤝 Contributing

Contributions are welcome! Please feel free to:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

### Areas for Contribution
- Additional alloy compositions
- Enhanced visualization features
- Extended thermodynamic analysis
- Educational content improvements
- Documentation enhancements


## 📚 Additional Resources

- [CALPHAD Methodology](https://www.calphad.org/)
- [pycalphad Documentation](https://pycalphad.org/)
- [Materials Science Education](https://www.materialseducation.org/)
- [Phase Diagram Resources](https://www.phase-trans.msm.cam.ac.uk/)

## 📞 Support

For questions, issues, or contributions, please:
- Open an issue on GitHub
- Contact the project maintainers
- Refer to the documentation in the `data/` directory

---

**Note**: This project is designed for educational and research purposes. The thermodynamic database and calculations are based on established CALPHAD models and experimental data.
