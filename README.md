# Jumping Bubbles

A high-performance computational framework for studying bubble coalescence and jumping phenomena on substrates using [Basilisk C](http://basilisk.fr/). This repository provides a complete suite of simulation tools, post-processing utilities, and validated test cases for investigating two-phase free-surface flows with adaptive mesh refinement.

## Key Features

### **Advanced Two-Phase Flow Modeling**
- Volume of Fluid (VOF) method for sharp interface tracking
- Conservative momentum advection coupled with VOF field
- Accurate surface tension implementation using the Brackbill method
- Reduced gravity approach for specific physical scenarios. This keeps the system of equations well-balanced.

### **High-Performance Computing**
- Adaptive octree mesh refinement based on multiple criteria:
  - Interface fraction gradients
  - Local curvature
  - Velocity field variations
- Parallel computation support:
  - OpenMP for shared-memory systems
  - MPI for distributed computing (HPC clusters)

### **Comprehensive Analysis Tools**
- Real-time monitoring of physical quantities
- Advanced post-processing capabilities:
  - 2D/3D visualization scripts
  - Interface reconstruction
  - Energy and momentum analysis
  - Trajectory tracking

## Project Structure

```plaintext
├── basilisk/src/               # Core Basilisk C framework
│   ├── navier-stokes/         # Flow solvers
│   │   ├── centered.h         # Main centered NS solver
│   │   └── conserving.h       # Conservative form solver
│   ├── two-phase*.h           # Two-phase flow implementations
│   ├── vof.h                  # Volume of Fluid method
│   ├── curvature.h           # Interface property calculations
│   ├── tension.h             # Surface tension (Brackbill)
│   ├── integral.h            # Surface tension integral form
│   ├── viscosity.h           # Implicit viscous stress solver
│   └── axi.h                 # Axisymmetric computations
├── src-local/                 # Project-specific modifications
├── postProcess/               # Analysis and visualization
│   ├── Video2DSlice.py       # 2D slice animations
│   ├── Video3D.py            # 3D visualization
│   ├── Visualization3D.ipynb  # Interactive 3D notebook
│   └── getFacets3D.c         # Interface extraction
└── testCases/                # Validation cases
    ├── JumpingBubbles.c      # Standard test case
    ├── JumpingBubbles-hydrophilic.c  # Hydrophilic substrate
    └── runCases*.sh          # Execution scripts
```

## Installation

### Prerequisites
- GCC compiler (version 7.0+)
- MPI implementation (OpenMPI or MPICH)
- Python 3.7+ (for post-processing)
- Required Python packages: numpy, matplotlib, scipy
- Xcode Command Line Tools (for MacOS users)

### Setting up Basilisk
1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/jumping-bubbles.git
   cd jumping-bubbles
   ```

2. Run the installation script:
   ```bash
   ./reset_install_requirements.sh
   ```
   This script will:
   - Install Basilisk if not present (use `--hard` flag for fresh installation)
   - Configure the environment automatically
   - Create `.project_config` with required paths
   - Verify the installation

3. Verify installation:
   ```bash
   source .project_config
   qcc --version
   ```

## Usage

### Running Simulations

#### Local Development (OpenMP)
```bash
# Compile
qcc -O2 -Wall -disable-dimensions -fopenmp -I$(PWD)/src-local testCases/JumpingBubbles.c -o JumpingBubbles -lm

# Run with 4 threads
export OMP_NUM_THREADS=4
./JumpingBubbles
```

#### HPC Deployment (MPI)
1. Generate initial condition locally first
2. Compile with MPI support:
   ```bash
   CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions \
   -I$(PWD)/src-local testCases/JumpingBubbles.c -o JumpingBubbles -lm
   ```
3. Use provided Slurm script: `testCases/runSnellius.sbatch`

### Post-Processing

#### Visualization
1. 2D Slice Videos:
   ```bash
   python postProcess/Video2DSlice.py <simulation_directory>
   ```

2. 3D Visualization:
   ```bash
   python postProcess/Video3D.py <simulation_directory>
   ```

3. Interactive Analysis:
   - Open `postProcess/Visualization3D.ipynb` in Jupyter

#### Data Analysis
- Use `getFacets3D.c` to extract interface geometry
- Various slice extraction tools available in `postProcess/`

## Contributing
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with:
   - Clear description of changes
   - Updated test cases if applicable
   - Modified documentation as needed

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Citation
If you use this code in your research, please cite:
```bibtex
```

## Acknowledgments
- Based on the [Basilisk C](http://basilisk.fr/) framework.