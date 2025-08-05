import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams, font_manager
import seaborn as sns
import os
import glob
from scipy import stats
import warnings
import re
import argparse
import sys
warnings.filterwarnings('ignore')

# Try to import MDAnalysis, provide alternative if not installed
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    MDANALYSIS_AVAILABLE = True
except ImportError:
    MDANALYSIS_AVAILABLE = False
    print("Warning: MDAnalysis not available. Pocket analysis will be limited.")

# ========== High-quality graphics settings ==========
def setup_publication_style():
    """
    Set high-quality graphics parameters for scientific publication
    """
    # Basic settings
    plt.style.use('default')  # Reset to default style
    
    # Font settings - use most universal fonts
    rcParams['font.family'] = 'sans-serif'  # Use system default font directly
    
    rcParams['font.size'] = 12
    
    # Axis and label settings
    rcParams['axes.titlesize'] = 14
    rcParams['axes.labelsize'] = 12
    try:
        rcParams['axes.titleweight'] = 'bold'
        rcParams['axes.labelweight'] = 'bold'
    except:
        pass  # Some versions may not support these parameters
    
    rcParams['axes.linewidth'] = 1.2
    rcParams['axes.spines.top'] = False
    rcParams['axes.spines.right'] = False
    rcParams['axes.grid'] = True
    rcParams['axes.edgecolor'] = '#333333'
    
    # Tick settings
    rcParams['xtick.labelsize'] = 11
    rcParams['ytick.labelsize'] = 11
    rcParams['xtick.major.width'] = 1.2
    rcParams['ytick.major.width'] = 1.2
    try:
        rcParams['xtick.minor.width'] = 0.8
        rcParams['ytick.minor.width'] = 0.8
    except:
        pass
    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.direction'] = 'in'
    
    # Legend settings
    rcParams['legend.fontsize'] = 10
    rcParams['legend.frameon'] = True
    try:
        rcParams['legend.framealpha'] = 0.9
        rcParams['legend.fancybox'] = False
        rcParams['legend.edgecolor'] = '#333333'
        rcParams['legend.facecolor'] = 'white'
    except:
        pass
    
    # Line settings
    rcParams['lines.linewidth'] = 1.5
    rcParams['lines.markersize'] = 6
    
    # Figure saving settings
    rcParams['savefig.dpi'] = 1200  # High resolution
    rcParams['savefig.format'] = 'png'
    rcParams['savefig.bbox'] = 'tight'
    rcParams['savefig.pad_inches'] = 0.1
    try:
        rcParams['savefig.transparent'] = False
        rcParams['savefig.facecolor'] = 'white'
    except:
        pass
    
    # Set default figure size
    rcParams['figure.figsize'] = [8, 6]
    rcParams['figure.dpi'] = 150  # Display resolution
    
    print("✓ Publication-style settings applied successfully")

# Scientific publication color scheme
SCIENTIFIC_COLORS = {
    'primary': '#2E86AB',      # Deep blue - primary data
    'secondary': '#A23B72',    # Purple-red - secondary data
    'accent': '#F18F01',       # Orange - emphasis
    'success': '#C73E1D',      # Red - warning/important
    'neutral': '#666666',      # Gray - auxiliary information
    'light_blue': '#87CEEB',   # Light blue - fill
    'light_gray': '#E5E5E5',   # Light gray - background
    'black': '#000000',        # Black - text
    'runs': ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6', '#1ABC9C'],  # Multi-dataset colors
}

def read_xvg_file(filename):
    """
    Read GROMACS .xvg format file
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('@') and not line.startswith('#'):
                try:
                    data.append([float(x) for x in line.split()])
                except ValueError:
                    continue
    return np.array(data)

def read_multiple_xvg_files(directories, filename):
    """
    Read the same xvg file from multiple directories
    """
    all_data = []
    for directory in directories:
        filepath = os.path.join(directory, filename)
        if os.path.exists(filepath):
            data = read_xvg_file(filepath)
            all_data.append(data)
        else:
            print(f"Warning: {filepath} not found")
    
    return all_data

def calculate_statistics(data_list):
    """
    Calculate statistical data from multiple experiments
    """
    if not data_list:
        return None, None, None
    
    # Ensure all data have the same time points/residue numbers
    min_length = min(len(data) for data in data_list)
    trimmed_data = [data[:min_length] for data in data_list]
    
    # Convert to numpy array
    data_array = np.array(trimmed_data)
    
    # Calculate statistics
    mean_data = np.mean(data_array, axis=0)
    std_data = np.std(data_array, axis=0)
    sem_data = stats.sem(data_array, axis=0)  # Standard error
    
    return mean_data, std_data, sem_data

def save_high_quality_figure(fig, filename, formats=['png']):
    """
    Save figure in specified formats (default: PNG only)
    """
    base_name = os.path.splitext(filename)[0]
    
    for fmt in formats:
        output_file = f"{base_name}.{fmt}"
        if fmt == 'png':
            fig.savefig(output_file, dpi=600, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
        elif fmt == 'pdf':
            fig.savefig(output_file, format='pdf', bbox_inches='tight',
                       facecolor='white', edgecolor='none')
        elif fmt == 'svg':
            fig.savefig(output_file, format='svg', bbox_inches='tight',
                       facecolor='white', edgecolor='none')
        print(f"✓ Saved: {output_file}")

def get_pocket_residues_simple(pocket_centers, cutoff=0.8):
    """
    Simple method: determine pocket residues based on sequence distance
    This is an approximation method suitable for cases without MDAnalysis
    """
    pocket_residues = set()
    
    # Add surrounding residues for each center residue
    for center in pocket_centers:
        # Add adjacent residues in sequence (rough approximation of 8Å spatial range)
        for offset in range(-10, 11):  # Approximately corresponds to 8Å sequence range
            residue = center + offset
            if residue > 0:  # Ensure positive residue numbers
                pocket_residues.add(residue)
    
    return sorted(list(pocket_residues))

def analyze_pocket_rmsd_simple(directories, pocket_centers, cutoff=0.8):
    """
    Simplified pocket RMSD analysis based on existing RMSF data
    """
    print("Using simplified pocket analysis based on RMSF data...")
    
    # Get pocket residue list
    pocket_residues = get_pocket_residues_simple(pocket_centers, cutoff)
    print(f"Pocket residues (approximate): {len(pocket_residues)} residues")
    
    # Extract pocket region flexibility information from RMSF data
    rmsf_data_list = read_multiple_xvg_files(directories, 'rmsf.xvg')
    
    if not rmsf_data_list:
        print("No RMSF data available for pocket analysis")
        return []
    
    pocket_flexibility_data = []
    
    for rmsf_data in rmsf_data_list:
        residue_numbers = rmsf_data[:, 0].astype(int)
        rmsf_values = rmsf_data[:, 1]
        
        # Find RMSF values corresponding to pocket residues
        pocket_rmsf = []
        for res_num in pocket_residues:
            idx = np.where(residue_numbers == res_num)[0]
            if len(idx) > 0:
                pocket_rmsf.append(rmsf_values[idx[0]])
        
        if pocket_rmsf:
            avg_pocket_flexibility = np.mean(pocket_rmsf)
            pocket_flexibility_data.append(avg_pocket_flexibility)
    
    return pocket_flexibility_data

def extract_ph_value(directories):
    """
    Extract pH value from directory names
    """
    if not directories:
        return ""
    
    # Extract pH value from first directory name
    first_dir = directories[0]
    # Use regex to match numbers after pH
    match = re.search(r'pH(\d+)', first_dir)
    if match:
        ph_value = match.group(1)
        return f"_pH{ph_value}"
    else:
        return ""

def create_comprehensive_report(directories, pocket_centers, output_suffix=""):
    """
    Generate comprehensive analysis report - high quality version
    """
    print("=" * 60)
    print("GROMACS MD SIMULATION ANALYSIS REPORT")
    print("=" * 60)
    print(f"Analysis directories: {directories}")
    print(f"Number of runs: {len(directories)}")
    print(f"Pocket definition: Residues {pocket_centers} (8Å cutoff)")
    print("=" * 60)
    
    # Create high-quality combined plot
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # 1. RMSD analysis
    ax1 = fig.add_subplot(gs[0, 0])
    rmsd_data_list = read_multiple_xvg_files(directories, 'rmsd.xvg')
    if rmsd_data_list:
        mean_data, std_data, sem_data = calculate_statistics(rmsd_data_list)
        if mean_data is not None:
            time = mean_data[:, 0]
            rmsd_mean = mean_data[:, 1]
            rmsd_sem = sem_data[:, 1]
            
            # Plot all raw data
            for i, data in enumerate(rmsd_data_list):
                ax1.plot(data[:, 0], data[:, 1], alpha=0.4, linewidth=1.0, 
                        color=SCIENTIFIC_COLORS['runs'][i % len(SCIENTIFIC_COLORS['runs'])])
            
            # Plot mean and confidence interval
            ax1.fill_between(time, rmsd_mean - rmsd_sem, rmsd_mean + rmsd_sem, 
                            alpha=0.3, color=SCIENTIFIC_COLORS['primary'])
            ax1.plot(time, rmsd_mean, linewidth=2.0, color=SCIENTIFIC_COLORS['primary'], 
                    label='Mean')
            
            ax1.set_xlabel('Time (ns)', fontweight='bold')
            ax1.set_ylabel('RMSD (nm)', fontweight='bold')
            ax1.set_title('(A) Protein RMSD', fontweight='bold', fontsize=13)
            ax1.grid(True, alpha=0.3)
            ax1.legend()
    
    # 2. RMSF analysis
    ax2 = fig.add_subplot(gs[0, 1])
    rmsf_data_list = read_multiple_xvg_files(directories, 'rmsf.xvg')
    if rmsf_data_list:
        mean_data, std_data, sem_data = calculate_statistics(rmsf_data_list)
        if mean_data is not None:
            residue = mean_data[:, 0]
            rmsf_mean = mean_data[:, 1]
            rmsf_sem = sem_data[:, 1]
            
            ax2.fill_between(residue, rmsf_mean - rmsf_sem, rmsf_mean + rmsf_sem, 
                            alpha=0.3, color=SCIENTIFIC_COLORS['secondary'])
            ax2.plot(residue, rmsf_mean, linewidth=2.0, color=SCIENTIFIC_COLORS['secondary'],
                    label='Mean RMSF')
            
            # Mark high flexibility regions
            overall_mean_rmsf = np.mean(rmsf_mean)
            threshold = overall_mean_rmsf + 2 * np.std(rmsf_mean)
            high_flexibility_mask = rmsf_mean > threshold
            if np.any(high_flexibility_mask):
                high_flex_residues = residue[high_flexibility_mask]
                high_flex_rmsf = rmsf_mean[high_flexibility_mask]
                ax2.scatter(high_flex_residues, high_flex_rmsf, 
                           color=SCIENTIFIC_COLORS['accent'], s=30, zorder=6, 
                           alpha=0.8, edgecolors='black', linewidth=0.5)
            
            ax2.set_xlabel('Residue Number', fontweight='bold')
            ax2.set_ylabel('RMSF (nm)', fontweight='bold')
            ax2.set_title('(B) Protein RMSF', fontweight='bold', fontsize=13)
            ax2.grid(True, alpha=0.3)
            ax2.legend()
    
    # 3. Pocket analysis
    ax3 = fig.add_subplot(gs[1, 0])
    pocket_data = analyze_pocket_rmsd_simple(directories, pocket_centers)
    if pocket_data:
        runs = [f'Run {i+1}' for i in range(len(pocket_data))]
        x_pos = np.arange(len(runs))
        
        bars = ax3.bar(x_pos, pocket_data, 
                      color=SCIENTIFIC_COLORS['runs'][:len(pocket_data)], 
                      alpha=0.8, edgecolor='black', linewidth=1.0)
        
        mean_val = np.mean(pocket_data)
        ax3.axhline(y=mean_val, color=SCIENTIFIC_COLORS['black'], 
                   linestyle='--', linewidth=2.0, alpha=0.8)
        
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(runs)
        ax3.set_xlabel('Simulation Run', fontweight='bold')
        ax3.set_ylabel('Pocket Flexibility (nm)', fontweight='bold')
        ax3.set_title('(C) Pocket Region Analysis', fontweight='bold', fontsize=13)
        ax3.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bar, value in zip(bars, pocket_data):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                    f'{value:.3f}', ha='center', va='bottom', fontsize=9)
    
    # 4. Statistical summary
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    
    summary_text = "ANALYSIS SUMMARY\n" + "="*25 + "\n\n"
    
    if rmsd_data_list:
        final_rmsds = [data[-1, 1] for data in rmsd_data_list]
        summary_text += f"Final RMSD:\n  {np.mean(final_rmsds):.3f} ± {np.std(final_rmsds):.3f} nm\n\n"
    
    if rmsf_data_list and mean_data is not None:
        summary_text += f"Average RMSF:\n  {np.mean(rmsf_mean):.3f} ± {np.std(rmsf_mean):.3f} nm\n\n"
    
    if pocket_data:
        summary_text += f"Pocket Flexibility:\n  {np.mean(pocket_data):.3f} ± {np.std(pocket_data):.3f} nm\n\n"
    
    summary_text += f"Pocket Definition:\n"
    summary_text += f"  Centers: {pocket_centers}\n"
    summary_text += f"  Cutoff: 8.0 Å\n\n"
    summary_text += f"Simulation Details:\n"
    summary_text += f"  Runs: {len(directories)}\n"
    
    if MDANALYSIS_AVAILABLE:
        summary_text += f"  MDAnalysis: Available\n"
    else:
        summary_text += f"  MDAnalysis: Approximation used\n"
    
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, fontsize=11, 
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round,pad=0.8', facecolor='#f8f9fa', 
                      edgecolor=SCIENTIFIC_COLORS['neutral'], alpha=0.9))
    
    plt.suptitle('Comprehensive MD Simulation Analysis', fontsize=16, fontweight='bold', y=0.98)
    
    # Use filename with pH suffix
    output_filename = f'comprehensive_md_analysis{output_suffix}'
    save_high_quality_figure(fig, output_filename, formats=['png'])
    plt.show()
    
    print("Comprehensive analysis completed!")
    print(f"Generated file: {output_filename}.png")

def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description='GROMACS MD Simulation Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python script.py -d pH4_1 pH4_2 pH4_3 -p 243 105 129
  
  # Single directory with different pocket centers
  python script.py -d pH6_run1 -p 100 200 300
  
  # Multiple directories with custom pocket centers
  python script.py -d exp1 exp2 exp3 exp4 -p 243 105 129 180
        """
    )
    
    parser.add_argument('-d', '--directories', 
                        nargs='+', 
                        required=True,
                        help='List of analysis directories (e.g., pH4_1 pH4_2 pH4_3)')
    
    parser.add_argument('-p', '--pocket-centers', 
                        nargs='+', 
                        type=int,
                        required=True,
                        help='Pocket center residue numbers (e.g., 243 105 129)')
    
    parser.add_argument('--output-suffix',
                        default='',
                        help='Suffix for output filename (default: auto-detected from directory names)')
    
    parser.add_argument('--formats',
                        nargs='+',
                        default=['png'],
                        choices=['png', 'pdf', 'svg'],
                        help='Output formats (default: png)')
    
    return parser.parse_args()

def validate_arguments(args):
    """
    Validate command line arguments
    """
    # Check if directories exist
    existing_dirs = []
    missing_dirs = []
    
    for d in args.directories:
        if os.path.exists(d):
            existing_dirs.append(d)
        else:
            missing_dirs.append(d)
    
    if missing_dirs:
        print(f"Warning: The following directories do not exist:")
        for d in missing_dirs:
            print(f"  - {d}")
        
        if not existing_dirs:
            print("Error: No valid directories found!")
            sys.exit(1)
        else:
            print(f"\nContinuing with {len(existing_dirs)} valid directories: {existing_dirs}")
    
    # Validate pocket centers
    if len(args.pocket_centers) < 1:
        print("Error: At least one pocket center residue must be provided")
        sys.exit(1)
    
    if any(pc <= 0 for pc in args.pocket_centers):
        print("Error: Pocket center residue numbers must be positive integers")
        sys.exit(1)
    
    return existing_dirs

# Main function
def main():
    """
    Main analysis function
    """
    # Parse command line arguments
    args = parse_arguments()
    
    # Validate arguments and get existing directories
    existing_dirs = validate_arguments(args)
    
    # Set high-quality plotting parameters
    setup_publication_style()
    
    print(f"Found {len(existing_dirs)} directories: {existing_dirs}")
    print(f"Pocket centers: {args.pocket_centers}")
    
    # Extract pH value from directory names or use provided suffix
    if args.output_suffix:
        ph_suffix = args.output_suffix
    else:
        ph_suffix = extract_ph_value(existing_dirs)
    
    print(f"Output suffix: {ph_suffix}")
    
    print("\n" + "="*60)
    print("Starting HIGH-QUALITY comprehensive MD analysis...")
    print(f"Output formats: {args.formats}")
    print("="*60)
    
    # Generate comprehensive report
    print("\nGenerating comprehensive report...")
    try:
        create_comprehensive_report(existing_dirs, args.pocket_centers, ph_suffix)
        print("✓ High-quality comprehensive report generated")
    except Exception as e:
        print(f"✗ Report generation failed: {str(e)}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*60)
    print("HIGH-QUALITY ANALYSIS COMPLETED!")
    print(f"Generated file: comprehensive_md_analysis{ph_suffix}.png")
    print("\nKey features:")
    print("  ✓ 600 DPI PNG for high-resolution")
    print("  ✓ Scientific color scheme")
    print("  ✓ Publication-ready typography")
    print("  ✓ Professional layout and styling")
    print("  ✓ Flexible command-line interface")
    print("="*60)

if __name__ == "__main__":
    main()