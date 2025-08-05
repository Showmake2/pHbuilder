#!/bin/bash

# Constant pH Molecular Dynamics Simulation Script (Selective GPU Optimization + Multi-task CPU Core Management)
# Author: [Your Name]
# Date: $(date)

# ====================================
# Configuration Parameters (can be modified as needed)
# ====================================

# Basic parameters - now set as defaults, will be overridden by command line args
PDB_FILE=""                   # Will be set by command line argument
PH_VALUE=""                   # Will be set by command line argument  
BOX_DISTANCE=""               # Will be set by command line argument

# GPU settings (only for MD simulation)
GPU_ID="0"                    # GPU device ID, can be set to "0,1" for multiple GPUs
USE_GPU_FOR_MD=true          # Whether to use GPU acceleration in MD stage
GPU_TASKS="nb,pme,bonded"    # GPU acceleration task types

# Parallel settings and CPU core management
NPROCS=$(nproc)               # Auto-detect CPU core count (current: 54 cores)
NPME=0                        # PME thread count, 0 for auto
CPU_OFFSET=0                  # CPU core offset, used to avoid conflicts in multi-tasking

# Optimization configuration for 54-core CPU + RTX4090
OPTIMAL_NTOMP=12              # Empirically optimal OpenMP thread count (for MD stage)
OPTIMAL_NTMPI=1               # Use 1 MPI rank for single GPU

# CPU optimization configuration (for EM/NVT/NPT stages)
CPU_NTOMP=$NPROCS            # CPU mode uses all cores

# Logging settings
LOG_DIR="logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# ====================================
# Function Definitions
# ====================================

# Create log directory
setup_logging() {
    mkdir -p $LOG_DIR
    exec > >(tee -a "$LOG_DIR/simulation_${TIMESTAMP}.log")
    exec 2>&1
}

# Error handling function
handle_error() {
    echo "Error: $1"
    echo "Script failed at $(date)"
    exit 1
}

# Check if file exists
check_file() {
    if [ ! -f "$1" ]; then
        handle_error "File $1 does not exist"
    fi
}

# Generate CPU affinity list
generate_cpu_affinity() {
    local offset=$1
    local nthreads=$2
    local max_cores=$3
    
    local cpu_list=""
    for ((i=0; i<nthreads; i++)); do
        local core_id=$(((offset + i) % max_cores))
        if [ $i -eq 0 ]; then
            cpu_list="$core_id"
        else
            cpu_list="$cpu_list,$core_id"
        fi
    done
    echo $cpu_list
}

# Check CPU core availability
check_cpu_availability() {
    local offset=$1
    local nthreads=$2
    
    echo "Checking CPU core availability..."
    echo "Requested cores: $CPU_OFFSET to $((CPU_OFFSET + OPTIMAL_NTOMP - 1))"
    
    # Show current CPU usage
    echo "Current CPU usage (top 5 most CPU-intensive processes):"
    ps -eo pid,ppid,cmd,pcpu --sort=-pcpu | head -6
    
    # Check if other GROMACS processes are running
    local gromacs_procs=$(pgrep -f "gmx mdrun" | wc -l)
    if [ $gromacs_procs -gt 0 ]; then
        echo "Warning: Detected $gromacs_procs GROMACS processes running"
        echo "Recommend using different CPU offset to avoid conflicts"
        
        # Show running GROMACS processes
        echo "Running GROMACS processes:"
        ps -f -C mdrun 2>/dev/null || echo "  (Unable to get detailed information)"
    fi
}

# CPU mode mdrun command (for EM/NVT/NPT) - strictly following original commands
run_mdrun_cpu() {
    local deffnm=$1
    local input_structure=$2
    local additional_opts="$3"
    
    echo "Running $deffnm using CPU (original command mode)..."
    
    # If CPU offset is set, use taskset to limit CPU usage
    local taskset_cmd=""
    if [ $CPU_OFFSET -gt 0 ]; then
        local cpu_list=$(generate_cpu_affinity $CPU_OFFSET $CPU_NTOMP $NPROCS)
        taskset_cmd="taskset -c $cpu_list"
        echo "  Using CPU cores: $cpu_list"
    fi
    
    if [ "$deffnm" = "EM" ]; then
        # EM: use -dd 13 1 1 parameter
        $taskset_cmd gmx mdrun -v -deffnm $deffnm -c ${deffnm}.pdb -npme 0 -dd 13 1 1 || handle_error "$deffnm simulation failed"
    else
        # NVT/NPT: only use -npme 0 parameter
        $taskset_cmd gmx mdrun -v -deffnm $deffnm -c ${deffnm}.pdb -npme 0 || handle_error "$deffnm simulation failed"
    fi
              
    echo "  ✓ $deffnm completed"
}

# GPU-optimized mdrun command (only for MD)
run_mdrun_gpu() {
    local deffnm=$1
    local input_structure=$2
    local additional_opts="$3"
    
    if [ "$USE_GPU_FOR_MD" = true ] && [ "$deffnm" = "MD" ]; then
        echo "Running $deffnm with GPU acceleration..."
        
        # Use optimized thread configuration
        local ntmpi=$OPTIMAL_NTMPI
        local ntomp=$OPTIMAL_NTOMP
        
        echo "  GPU configuration: RTX 4090 (24GB)"
        echo "  CPU configuration: $NPROCS total cores"
        echo "  Thread count: $ntomp"
        echo "  CPU offset: $CPU_OFFSET"
        echo "  MPI ranks: $ntmpi"
        echo "  GPU tasks: Non-bonded interactions + PME + Bonded interactions"
        
        # Check CPU availability
        check_cpu_availability $CPU_OFFSET $ntomp
        
        # Build CPU affinity parameters
        local pin_options=""
        local taskset_cmd=""
        
        if [ $CPU_OFFSET -gt 0 ]; then
            # Use taskset to limit CPU range, disable GROMACS built-in CPU binding
            local cpu_list=$(generate_cpu_affinity $CPU_OFFSET $ntomp $NPROCS)
            taskset_cmd="taskset -c $cpu_list"
            pin_options="-pin off"  # Disable GROMACS built-in binding, use taskset control
            echo "  Actually using CPU cores: $cpu_list"
        else
            # Use GROMACS built-in CPU binding
            pin_options="-pin on -pinstride 1"
            echo "  Using GROMACS default CPU binding (first $ntomp cores)"
        fi
        
        # MD simulation: use GPU acceleration + intelligent CPU allocation
        $taskset_cmd gmx mdrun -v -deffnm $deffnm \
                  -c ${deffnm}.pdb \
                  -gpu_id $GPU_ID \
                  -nb gpu -pme gpu -bonded gpu \
                  -ntmpi $ntmpi \
                  -ntomp $ntomp \
                  -npme 0 \
                  $pin_options \
                  -resethway \
                  $additional_opts || handle_error "$deffnm simulation failed"
                  
        echo "  ✓ $deffnm completed (GPU acceleration + CPU cores $CPU_OFFSET-$((CPU_OFFSET + ntomp - 1)))"
    else
        # For non-MD steps or when GPU is disabled, use CPU
        run_mdrun_cpu $deffnm $input_structure "$additional_opts"
    fi
}

# Intelligent GPU settings configuration
configure_gpu_settings() {
    if [ "$USE_GPU_FOR_MD" = false ]; then
        echo "GPU acceleration disabled, all steps will use CPU"
        return
    fi
    
    local ngpus=$(nvidia-smi -L | wc -l)
    local gpu_memory=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
    local gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
    
    echo "Detected GPU: $gpu_name ($gpu_memory MB)"
    echo "GPU will only be used for MD simulation stage"
    
    # Special optimization for RTX 4090 + 54-core CPU
    if [[ "$gpu_name" == *"4090"* ]] && [ $NPROCS -ge 50 ]; then
        echo "Detected RTX 4090 + high core count CPU, applying special optimization..."
        GPU_ID="0"
        OPTIMAL_NTMPI=1
        OPTIMAL_NTOMP=12  # For 4090, 12 threads is usually optimal
        
        echo "MD stage GPU optimization configuration:"
        echo "  - Single GPU maximum utilization"
        echo "  - 12 OpenMP threads"
        echo "  - Support CPU core offset to avoid multi-task conflicts"
        echo "  - Enable all GPU acceleration tasks"
        
    elif [ $ngpus -eq 1 ]; then
        GPU_ID="0"
        OPTIMAL_NTMPI=1
        # For other single GPU configurations, use more conservative settings
        if [ $NPROCS -gt 32 ]; then
            OPTIMAL_NTOMP=16
        elif [ $NPROCS -gt 16 ]; then
            OPTIMAL_NTOMP=8
        else
            OPTIMAL_NTOMP=$NPROCS
        fi
        echo "Single GPU configuration optimized (only for MD)"
        
    elif [ $ngpus -eq 2 ]; then
        GPU_ID="0,1"
        OPTIMAL_NTMPI=2
        OPTIMAL_NTOMP=$((NPROCS / 2))
        echo "Dual GPU configuration optimized (only for MD)"
        
    elif [ $ngpus -ge 4 ]; then
        GPU_ID="0,1,2,3"
        OPTIMAL_NTMPI=4
        OPTIMAL_NTOMP=$((NPROCS / 4))
        echo "Multi-GPU configuration optimized (only for MD)"
    fi
    
    # Validate CPU offset reasonableness
    if [ $((CPU_OFFSET + OPTIMAL_NTOMP)) -gt $NPROCS ]; then
        echo "Warning: CPU offset ($CPU_OFFSET) + thread count ($OPTIMAL_NTOMP) exceeds total cores ($NPROCS)"
        echo "Will use wrap-around method for CPU core allocation"
    fi
    
    # 4090's 24GB memory is sufficient for most simulations
    if [ $gpu_memory -gt 20000 ]; then
        echo "GPU memory sufficient, can handle large systems"
    fi
}

# Check GROMACS environment
check_environment() {
    echo "Checking GROMACS environment..."
    if ! command -v gmx &> /dev/null; then
        handle_error "GROMACS not found, please check environment setup"
    fi
    
    # Check taskset command
    if ! command -v taskset &> /dev/null; then
        echo "Warning: taskset command not found, CPU affinity control may not be available"
    fi
    
    if [ "$USE_GPU_FOR_MD" = true ]; then
        # Check GPU availability
        if command -v nvidia-smi &> /dev/null; then
            echo "GPU information:"
            nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv
            
            # Intelligently configure GPU parameters
            configure_gpu_settings
        else
            echo "Warning: No NVIDIA GPU detected, MD stage will use CPU mode"
            USE_GPU_FOR_MD=false
        fi
    fi
    
    echo "Current configuration:"
    echo "  Total CPU cores: $NPROCS"
    echo "  CPU offset: $CPU_OFFSET"
    echo "  MD stage thread count: $OPTIMAL_NTOMP"
    if [ $CPU_OFFSET -gt 0 ]; then
        echo "  MD stage will use cores: $CPU_OFFSET to $((CPU_OFFSET + OPTIMAL_NTOMP - 1))"
    else
        echo "  MD stage will use cores: 0 to $((OPTIMAL_NTOMP - 1))"
    fi
    echo "  EM/NVT/NPT: strictly follow original commands"
    echo "  - EM: gmx mdrun -v -deffnm EM -c EM.pdb -npme 0 -dd 13 1 1"
    echo "  - NVT/NPT: gmx mdrun -v -deffnm [step] -c [step].pdb -npme 0"
    echo "  MD stage GPU acceleration: $USE_GPU_FOR_MD"
    if [ "$USE_GPU_FOR_MD" = true ]; then
        echo "  MD stage GPU device: $GPU_ID"
    fi
}

# Recommend CPU offset
suggest_cpu_offset() {
    echo "Multi-task running suggestions:"
    echo "For 54-core system, recommended CPU offsets:"
    echo "  Task1: --cpu-offset 0   (use cores 0-11)"
    echo "  Task2: --cpu-offset 12  (use cores 12-23)"
    echo "  Task3: --cpu-offset 24  (use cores 24-35)"
    echo "  Task4: --cpu-offset 36  (use cores 36-47)"
    echo "  Task5: --cpu-offset 48  (use cores 48-53, then 0-5)"
    echo ""
    echo "Example commands:"
    echo "  # Terminal 1"
    echo "  $0 -p protein1.pdb --ph 6.0 --box 1.5 --cpu-offset 0"
    echo "  # Terminal 2"  
    echo "  $0 -p protein2.pdb --ph 6.0 --box 1.5 --cpu-offset 12"
    echo "  # Terminal 3"
    echo "  $0 -p protein3.pdb --ph 6.0 --box 1.5 --cpu-offset 24"
}

# ====================================
# Main Simulation Pipeline
# ====================================

main_simulation() {
    echo "======================================"
    echo "Starting Constant pH Molecular Dynamics Simulation"
    echo "Time: $(date)"
    echo "pH value: $PH_VALUE"
    echo "Input file: $PDB_FILE"
    echo "Box distance: $BOX_DISTANCE"
    echo "CPU offset: $CPU_OFFSET"
    echo "Strategy: EM/NVT/NPT strictly follow original commands, MD uses GPU+intelligent CPU allocation"
    echo "======================================"
    
    # Check input file
    check_file $PDB_FILE
    
    # Step 1: Setup GROMACS environment
    echo "1. Setting up GROMACS environment..."
    source /opt/software/gromacs-constantph-gpu/env.sh || handle_error "Environment setup failed"
    
    # Step 2: Generate topology file
    echo "2. Generating topology file..."
    echo "1" | phbuilder gentopol -f $PDB_FILE -ph $PH_VALUE || handle_error "Topology generation failed"
    check_file "phprocessed.pdb"
    
    # Step 3: Create simulation box
    echo "3. Creating simulation box..."
    gmx editconf -f phprocessed.pdb -o box.pdb -bt cubic -d $BOX_DISTANCE || handle_error "Box creation failed"
    
    # Step 4: Add solvent
    echo "4. Adding solvent..."
    gmx solvate -cp box.pdb -p topol.top -o solvated.pdb || handle_error "Solvation failed"
    
    # Step 5: Neutralize system
    echo "5. Neutralizing system..."
    phbuilder neutralize -f solvated.pdb || handle_error "Neutralization failed"
    
    # Step 6: Generate parameter files
    echo "6. Generating constant pH parameters..."
    phbuilder genparams -f phneutral.pdb -ph $PH_VALUE || handle_error "Parameter generation failed"
    
    # Step 7: Energy minimization (strictly follow original command)
    echo "7. Energy minimization (original command)..."
    gmx grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -maxwarn 1 || handle_error "EM preprocessing failed"
    run_mdrun_cpu "EM" "phneutral.pdb"
    
    # Step 8: NVT equilibration (strictly follow original command)
    echo "8. NVT equilibration (original command)..."
    gmx grompp -f NVT.mdp -c EM.pdb -p topol.top -n index.ndx -o NVT.tpr || handle_error "NVT preprocessing failed"
    run_mdrun_cpu "NVT" "EM.pdb"
    
    # Step 9: NPT equilibration (strictly follow original command)
    echo "9. NPT equilibration (original command)..."
    gmx grompp -f NPT.mdp -c NVT.pdb -p topol.top -n index.ndx -o NPT.tpr || handle_error "NPT preprocessing failed"
    run_mdrun_cpu "NPT" "NVT.pdb"
    
    # Step 10: Molecular dynamics simulation (use GPU or CPU)
    echo "10. Molecular dynamics simulation ($([ "$USE_GPU_FOR_MD" = true ] && echo "GPU" || echo "CPU") mode)..."
    gmx grompp -f MD.mdp -c NPT.pdb -p topol.top -n index.ndx -o MD.tpr || handle_error "MD preprocessing failed"
    run_mdrun_gpu "MD" "NPT.pdb" "-x MD.xtc"
    
    # Step 11: Result analysis
    echo "11. Analyzing results..."
    echo "0 0" | gmx rms -s MD.tpr -f MD.xtc -o rmsd.xvg -tu ns || handle_error "RMSD analysis failed"
    gmx rmsf -s MD.tpr -f MD.xtc -o rmsf.xvg -res || handle_error "RMSF analysis failed"
    
    echo "======================================"
    echo "Simulation completed! Time: $(date)"
    echo "Result files:"
    echo "  - Final structure: MD.pdb"
    echo "  - Trajectory file: MD.xtc"
    echo "  - RMSD analysis: rmsd.xvg"
    echo "  - RMSF analysis: rmsf.xvg"
    echo "Performance summary:"
    echo "  - EM: original command"
    echo "  - NVT/NPT: original command"
    echo "  - MD: $([ "$USE_GPU_FOR_MD" = true ] && echo "GPU acceleration + CPU cores$CPU_OFFSET-$((CPU_OFFSET + OPTIMAL_NTOMP - 1))" || echo "CPU mode")"
    echo "======================================"
}

# ====================================
# Batch Processing Functionality (supports automatic CPU offset)
# ====================================

batch_simulation() {
    local ph_values=("3.0" "4.0" "5.0" "6.0" "7.0")  # pH value list
    local pdb_files=("protein1.pdb" "protein2.pdb")   # Protein file list
    
    echo "Starting batch simulation..."
    
    local task_count=0
    for pdb in "${pdb_files[@]}"; do
        for ph in "${ph_values[@]}"; do
            if [ -f "$pdb" ]; then
                echo "Processing $pdb at pH=$ph conditions..."
                
                # Automatically calculate CPU offset
                local auto_offset=$((task_count * OPTIMAL_NTOMP))
                if [ $auto_offset -ge $NPROCS ]; then
                    auto_offset=$((auto_offset % NPROCS))
                fi
                
                # Create work directory
                work_dir="${pdb%.pdb}_pH${ph}_cpu${auto_offset}"
                mkdir -p $work_dir
                cd $work_dir
                
                # Copy necessary files
                cp ../$pdb .
                cp ../*.mdp .
                
                # Set parameters
                PDB_FILE=$pdb
                PH_VALUE=$ph
                CPU_OFFSET=$auto_offset
                
                echo "Auto-assigned CPU offset: $auto_offset"
                
                # Run simulation
                main_simulation
                
                cd ..
                task_count=$((task_count + 1))
            else
                echo "Warning: File $pdb does not exist, skipping..."
            fi
        done
    done
}

# ====================================
# Command Line Argument Processing
# ====================================

show_usage() {
    echo "Usage: $0 -p <PDB_FILE> --ph <PH_VALUE> --box <BOX_DISTANCE> [options]"
    echo "Required arguments:"
    echo "  -p, --pdb FILE         Input PDB file (required)"
    echo "  --ph VALUE             pH value (required)"
    echo "  --box DISTANCE         Box distance (required)"
    echo ""
    echo "Optional arguments:"
    echo "  -g, --gpu ID           GPU device ID (default: $GPU_ID)"
    echo "  --cpu-offset N         CPU core offset (default: $CPU_OFFSET)"
    echo "  --no-gpu-md            Disable GPU acceleration for MD stage"
    echo "  -b, --batch            Batch processing mode"
    echo "  --suggest-cpu          Show multi-task CPU allocation suggestions"
    echo "  -h, --help             Show help information"
    echo ""
    echo "Note: EM/NVT/NPT stages always use CPU, only MD stage can optionally use GPU acceleration"
    echo ""
    echo "Multi-task examples:"
    echo "  # Task 1 (use cores 0-11)"
    echo "  $0 -p protein1.pdb --ph 6.0 --box 1.5 --cpu-offset 0"
    echo "  # Task 2 (use cores 12-23)"  
    echo "  $0 -p protein2.pdb --ph 6.0 --box 1.5 --cpu-offset 12"
    echo "  # Task 3 (use cores 24-35)"
    echo "  $0 -p protein3.pdb --ph 6.0 --box 1.5 --cpu-offset 24"
    echo ""
    echo "Other examples:"
    echo "  $0 -p protein.pdb --ph 5.0 --box 1.2     # Specify file, pH and box distance"
    echo "  $0 -p protein.pdb --ph 6.0 --box 1.5 --no-gpu-md    # Disable GPU for MD stage"
    echo "  $0 -b                                     # Batch processing (auto CPU allocation)"
}

# Validate required arguments
validate_required_args() {
    local missing_args=()
    
    if [ -z "$PDB_FILE" ]; then
        missing_args+=("PDB file (-p/--pdb)")
    fi
    
    if [ -z "$PH_VALUE" ]; then
        missing_args+=("pH value (--ph)")
    fi
    
    if [ -z "$BOX_DISTANCE" ]; then
        missing_args+=("box distance (--box)")
    fi
    
    if [ ${#missing_args[@]} -gt 0 ]; then
        echo "Error: Missing required arguments:"
        printf '  %s\n' "${missing_args[@]}"
        echo ""
        show_usage
        exit 1
    fi
    
    # Validate numeric values
    if ! [[ "$PH_VALUE" =~ ^[0-9]+\.?[0-9]*$ ]]; then
        echo "Error: pH value must be a number"
        exit 1
    fi
    
    if ! [[ "$BOX_DISTANCE" =~ ^[0-9]+\.?[0-9]*$ ]]; then
        echo "Error: Box distance must be a number"
        exit 1
    fi
    
    if ! [[ "$CPU_OFFSET" =~ ^[0-9]+$ ]]; then
        echo "Error: CPU offset must be a non-negative integer"
        exit 1
    fi
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--pdb)
            PDB_FILE="$2"
            shift 2
            ;;
        --ph)
            PH_VALUE="$2"
            shift 2
            ;;
        --box)
            BOX_DISTANCE="$2"
            shift 2
            ;;
        -g|--gpu)
            GPU_ID="$2"
            shift 2
            ;;
        --cpu-offset)
            CPU_OFFSET="$2"
            shift 2
            ;;
        --no-gpu-md)
            USE_GPU_FOR_MD=false
            shift
            ;;
        -b|--batch)
            BATCH_MODE=true
            shift
            ;;
        --suggest-cpu)
            suggest_cpu_offset
            exit 0
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# ====================================
# Main Program Entry Point
# ====================================

# For batch mode, skip argument validation (uses predefined values)
if [ "$BATCH_MODE" != true ]; then
    validate_required_args
fi

# Setup logging
setup_logging

# Check environment
check_environment

# Run simulation
if [ "$BATCH_MODE" = true ]; then
    batch_simulation
else
    main_simulation
fi

echo "Script execution completed: $(date)"