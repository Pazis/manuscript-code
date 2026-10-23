"""antiSMASH vs GECCO cluster comparison per genome: counts, overlap, statistics and plots."""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import glob

def parse_gff_attributes(attribute_string):
    attrs = {}
    for item in attribute_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs

def parse_gff_file(filepath):
    """(antiSMASH clusters, GECCO clusters) from biosynthetic-gene-cluster features."""
    as_clusters = []
    gecco_clusters = []

    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) != 9:
                    continue

                tool_source = parts[1]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                attributes = parse_gff_attributes(parts[8])

                if tool_source == 'antiSMASH' and feature_type == 'biosynthetic-gene-cluster':
                    as_clusters.append({
                        'start': start, 'end': end, 'length': end - start,
                        'product': attributes.get('product', 'unknown'),
                        'id': attributes.get('ID', 'unknown')
                    })
                elif tool_source == 'GECCO' and feature_type == 'biosynthetic-gene-cluster':
                    gecco_clusters.append({
                        'start': start, 'end': end, 'length': end - start,
                        'product': attributes.get('product', 'unknown'),
                        'id': attributes.get('ID', 'unknown')
                    })
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")

    return as_clusters, gecco_clusters

def calculate_overlap_metrics(as_clusters, gecco_clusters):
    """(antiSMASH clusters overlapping GECCO, overlap bp, antiSMASH bp covered, GECCO bp covered)."""
    overlapping_as_count = 0
    total_overlap_bp = 0
    total_as_bp_covered = 0
    total_gecco_bp_covered = 0

    for as_c in as_clusters:
        is_overlapping = False
        bp_covered_in_this_cluster = 0
        for gc_c in gecco_clusters:
            overlap_start = max(as_c['start'], gc_c['start'])
            overlap_end = min(as_c['end'], gc_c['end'])
            if overlap_start < overlap_end:
                is_overlapping = True
                overlap_len = overlap_end - overlap_start
                total_overlap_bp += overlap_len
                bp_covered_in_this_cluster += overlap_len

        if is_overlapping:
            overlapping_as_count += 1
            total_as_bp_covered += min(bp_covered_in_this_cluster, as_c['length'])

    for gc_c in gecco_clusters:
        bp_covered_in_this_cluster = 0
        for as_c in as_clusters:
            overlap_start = max(as_c['start'], gc_c['start'])
            overlap_end = min(as_c['end'], gc_c['end'])
            if overlap_start < overlap_end:
                overlap_len = overlap_end - overlap_start
                bp_covered_in_this_cluster += overlap_len
        total_gecco_bp_covered += min(bp_covered_in_this_cluster, gc_c['length'])

    return overlapping_as_count, total_overlap_bp, total_as_bp_covered, total_gecco_bp_covered

def process_files(file_list, OUTPUT_DIR):
    results = []
    print(f"Found {len(file_list)} GFF files. Processing...")

    for filepath in file_list:
        filename = os.path.basename(filepath)
        locus_tag = os.path.splitext(filename)[0] # Assuming filename is the Locus Tag

        as_clusters, gecco_clusters = parse_gff_file(filepath)

        n_as = len(as_clusters)
        n_gecco = len(gecco_clusters)

        n_overlap, overlap_bp, as_bp_cov, gecco_bp_cov = calculate_overlap_metrics(as_clusters, gecco_clusters)

        total_as_length = sum(c['length'] for c in as_clusters) if n_as > 0 else 0
        total_gecco_length = sum(c['length'] for c in gecco_clusters) if n_gecco > 0 else 0

        pct_as_covered = (as_bp_cov / total_as_length * 100) if total_as_length > 0 else 0
        pct_gecco_covered = (gecco_bp_cov / total_gecco_length * 100) if total_gecco_length > 0 else 0

        # Inter-tool consensus: overlapping antiSMASH clusters / all antiSMASH clusters
        consensus = (n_overlap / n_as) if n_as > 0 else 0

        # Overlap efficiency (Dice-like): 2 * overlap / (antiSMASH length + GECCO length)
        efficiency = (2 * overlap_bp / (total_as_length + total_gecco_length)) if (total_as_length + total_gecco_length) > 0 else 0

        # Consolidation ratio: >1 means antiSMASH splits loci that GECCO merges
        consolidation = (n_as / n_gecco) if n_gecco > 0 else 0

        gecco_products = ", ".join(list(set([c['product'] for c in gecco_clusters])))

        results.append({
            'Locus Tag': locus_tag,
            'AntiSMASH Count': n_as,
            'GECCO Count': n_gecco,
            'Overlapping Clusters': n_overlap,
            'GECCO Products': gecco_products,
            'Total Overlap (bp)': overlap_bp,
            '% AntiSMASH Covered': round(pct_as_covered, 2),
            '% GECCO Covered': round(pct_gecco_covered, 2),
            'Inter-tool Consensus': round(consensus, 2),
            'Avg Overlap Efficiency': round(efficiency, 2),
            'Consolidation Ratio': round(consolidation, 2)
        })

    df = pd.DataFrame(results)

    means = df.mean(numeric_only=True)
    means['Locus Tag'] = 'AVERAGE'
    df = pd.concat([df, pd.DataFrame([means])], ignore_index=True)

    df.to_excel(f'{OUTPUT_DIR}/bgc_analysis_results.xlsx', index=False)
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/bgc_analysis_results.xlsx")
    return f'{OUTPUT_DIR}/bgc_analysis_results.xlsx'


def analyze_bgc_data(file_path, OUTPUT_DIR):
    """Paired antiSMASH vs GECCO count test and Figures 1-4, on genomes where both tools found overlapping clusters."""
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        return

    df = pd.read_excel(file_path)
    df_clean = df[df['Locus Tag'] != 'AVERAGE'].copy()
    print(f"\nLoaded {len(df_clean)} genomes for statistical analysis.")

    df_clean = df_clean[
        (df_clean['AntiSMASH Count'] > 0) &
        (df_clean['GECCO Count'] > 0) &
        (df_clean['Overlapping Clusters'] > 0)
    ].copy()

    print("--- STATISTICAL ANALYSIS ---")

    differences = df_clean['AntiSMASH Count'] - df_clean['GECCO Count']
    shapiro_stat, shapiro_p = stats.shapiro(differences)
    print(f"\n1. Normality Check (Shapiro-Wilk): P={shapiro_p:.5f}")

    if shapiro_p > 0.05:
        test_name = "Paired T-Test"
        stat, p_value = stats.ttest_rel(df_clean['AntiSMASH Count'], df_clean['GECCO Count'])
    else:
        test_name = "Wilcoxon Signed-Rank"
        stat, p_value = stats.wilcoxon(df_clean['AntiSMASH Count'], df_clean['GECCO Count'])

    print(f"\n2. Detection Power ({test_name}): P={p_value:.5e}")

    print("\n--- GENERATING PLOTS ---")
    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 12})

    # PLOT 1: Detection Sensitivity (Box Plot)
    plt.figure(figsize=(8, 6))
    plot_data = df_clean.melt(value_vars=['AntiSMASH Count', 'GECCO Count'],
                              var_name='Tool', value_name='Count')

    sns.boxplot(x='Tool', y='Count', data=plot_data, palette=['#3498db', '#e74c3c'], width=0.5)
    sns.stripplot(x='Tool', y='Count', data=plot_data, color='black', alpha=0.3, jitter=True)

    plt.title('BGC Detection Sensitivity:\nAntiSMASH vs GECCO', fontsize=14, fontweight='bold')
    plt.ylabel('Number of BGCs per Genome')
    plt.xlabel('')
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure_1_Detection_Sensitivity.png', dpi=300)
    plt.close()
    print("-> Saved 'Figure_1_Detection_Sensitivity.png'")

    # PLOT 2: Fragmentation Analysis (Histogram)
    plt.figure(figsize=(8, 6))
    sns.histplot(df_clean['Consolidation Ratio'], bins=10, kde=True, color='#9b59b6', edgecolor='black')

    mean_ratio = df_clean['Consolidation Ratio'].mean()
    plt.axvline(mean_ratio, color='red', linestyle='--', label=f'Mean: {mean_ratio:.2f}')

    plt.title('Locus Fragmentation Analysis\n(Consolidation Ratio Distribution)', fontsize=14, fontweight='bold')
    plt.xlabel('Ratio (AntiSMASH Clusters / GECCO Clusters)')
    plt.ylabel('Frequency (Genomes)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure_2_Consolidation_Ratio.png', dpi=300)
    plt.close()
    print("-> Saved 'Figure_2_Consolidation_Ratio.png'")

    # PLOT 3: Boundary Precision (Scatter Plot)
    plt.figure(figsize=(9, 7))

    sns.scatterplot(
        data=df_clean,
        x='% AntiSMASH Covered',
        y='% GECCO Covered',
        color="red",
        edgecolor='black',
        alpha=0.6,
        s=80
    )

    plt.plot([0, 100], [0, 100], color='black', linestyle='--', linewidth=2, label='Perfect Match (1:1)')

    plt.title('Boundary Definition Precision (Density)', fontsize=14, fontweight='bold')
    plt.xlabel('% AntiSMASH Locus Covered by GECCO')
    plt.ylabel('% GECCO Locus Covered by AntiSMASH')

    plt.xlim(0, 105)
    plt.ylim(0, 105)

    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure_3_Boundary_Precision.png', dpi=300)
    plt.close()
    print("-> Saved 'Figure_3_Boundary_Precision.png'")

    # PLOT 4: Performance Metrics (Bar Chart)
    plt.figure(figsize=(10, 6))

    avg_consensus = df_clean['Inter-tool Consensus'].mean()
    avg_efficiency = df_clean['Avg Overlap Efficiency'].mean()
    avg_ratio = df_clean['Consolidation Ratio'].mean()

    metrics = ['Inter-tool Consensus', 'Overlap Efficiency', 'Consolidation Ratio']
    values = [avg_consensus, avg_efficiency, avg_ratio]

    bars = plt.bar(metrics, values, color="blue", edgecolor='black', alpha=0.8)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height:.2f}', ha='center', va='bottom', fontweight='bold')

    plt.title('Average Pipeline Performance Metrics', fontsize=14, fontweight='bold')
    plt.ylabel('Metric Value')
    plt.ylim(0, max(values) * 1.2)

    plt.figtext(0.5, 0.01,
                "Note: Consensus & Efficiency range [0-1]. Ratio > 1 indicates fragmentation.",
                ha="center", fontsize=10, style='italic')

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Figure_4_Performance_Metrics.png', dpi=300)
    plt.close()
    print("-> Saved 'Figure_4_Performance_Metrics.png'")

    stats_summary = {
        'Metric': ['Shapiro_P_Value', 'Test_Used', 'Statistical_P_Value'],
        'Value': [shapiro_p, test_name, p_value]
    }
    stats_df = pd.DataFrame(stats_summary)

    with pd.ExcelWriter(f'{OUTPUT_DIR}/bgc_statistical_summary.xlsx') as writer:
        df_clean.to_excel(writer, sheet_name='Raw_Data', index=False)
        stats_df.to_excel(writer, sheet_name='Stats_Summary', index=False)

    print(f"-> Stats saved to '{OUTPUT_DIR}/bgc_statistical_summary.xlsx'")

    print("\nAnalysis Complete.")

def cmd_plot(args):
    input_dir = args.dir
    if not os.path.isdir(input_dir):
        print(f"Error: The path '{input_dir}' is not a valid directory.")
        return

    files_to_process = glob.glob(os.path.join(input_dir, "*.gff")) + glob.glob(os.path.join(input_dir, "*.gff3"))
    if not files_to_process:
        print(f"No .gff or .gff3 files found in '{input_dir}'")
        return

    bgc_excel = process_files(files_to_process, args.output)
    if not args.skip_plots:
        analyze_bgc_data(bgc_excel, args.output)
