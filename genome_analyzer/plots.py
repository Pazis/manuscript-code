"""Per-genome statistics and Figures 1-7 from the summary CSV: Bakta (Initial) vs pipeline (Enhanced), with Prokka where present."""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import glob


def count_unique_contigs(gff_path):
    """Counts unique sequence regions (contigs) in a GFF file."""
    unique_contigs = set()
    try:
        with open(gff_path, 'r') as file:
            for line in file:
                if line.startswith('#') or not line.strip():
                    continue
                columns = line.split('\t')
                if len(columns) > 0:
                    unique_contigs.add(columns[0])
        return len(unique_contigs)
    except Exception:
        return 0


def calculate_genome_size(gff_path):
    """Total genome size (bp): ##sequence-region lengths, else the largest feature end per contig."""
    contig_sizes = {}
    try:
        with open(gff_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line: continue

                if line.startswith('##sequence-region'):
                    parts = line.split()
                    if len(parts) >= 4:
                        contig_sizes[parts[1]] = int(parts[3])
                    continue

                if line.startswith('#'): continue

                cols = line.split('\t')
                if len(cols) >= 5:
                    seqid = cols[0]
                    try:
                        end = int(cols[4])
                        if end > contig_sizes.get(seqid, 0):
                            contig_sizes[seqid] = end
                    except ValueError:
                        continue
        return sum(contig_sizes.values())
    except Exception:
        return 0


def load_data(file_path):
    """Loads a CSV or Excel file and converts the *_Undesc_PCT columns ("12.3%") to floats."""
    ext = os.path.splitext(file_path)[-1].lower()

    if ext in ['.xlsx', '.xls']:
        df = pd.read_excel(file_path)
    else:
        try:
            df = pd.read_csv(file_path)
        except UnicodeDecodeError:
            df = pd.read_csv(file_path, encoding='latin1')

    def clean_pct(val):
        if pd.isna(val): return 0.0
        if isinstance(val, str):
            return float(val.replace('%', '').strip())
        return float(val)

    for col in ['Initial_Undesc_PCT', 'Enhanced_Undesc_PCT', 'Prokka_Undesc_PCT']:
        if col in df.columns:
            df[col] = df[col].apply(clean_pct)

    return df


def enrich_data(df, gff_dir):
    """Adds 'Contigs' and 'Genome_Size' columns from GFF files matched on the 'File' column (0 if unmatched)."""
    if not gff_dir or not os.path.isdir(gff_dir):
        print(f"Error: GFF directory '{gff_dir}' not found.")
        return df

    gff_files = [f for f in glob.glob(os.path.join(gff_dir, "*")) if f.endswith('.gff') or f.endswith('.gff3')]
    print(f"Found {len(gff_files)} GFF files in '{gff_dir}'.")

    # Match on full name, name without extension, or name without the "_enhanced" suffix
    gff_map = {}
    for p in gff_files:
        base = os.path.basename(p)
        gff_map[base] = p
        gff_map[os.path.splitext(base)[0]] = p
        clean = base.replace(".gff3", "").replace(".gff", "")
        if clean.endswith("_enhanced"):
            clean = clean.replace("_enhanced", "")
        gff_map[clean] = p

    contigs = []
    sizes = []
    matches = 0

    for index, row in df.iterrows():
        csv_name = str(row['File']).strip()
        target = gff_map.get(csv_name)

        if target:
            contigs.append(count_unique_contigs(target))
            sizes.append(calculate_genome_size(target))
            matches += 1
        else:
            if matches == 0 and index < 3:
                print(f"   No GFF file matches '{csv_name}'.")
            contigs.append(0)
            sizes.append(0)

    df['Contigs'] = contigs
    df['Genome_Size'] = sizes
    df['Genome_Size_MB'] = df['Genome_Size'].apply(lambda x: x / 1_000_000 if x > 0 else 0)

    print(f"Enriched {matches}/{len(df)} genomes.")
    return df


def correlate_size(df):
    """Genome size vs Enhanced hypothetical %: Pearson if both pass Shapiro-Wilk, else Spearman."""
    print("\n--- Genome Size Correlation ---")
    if 'Genome_Size' not in df.columns or df['Genome_Size'].sum() == 0:
        print("Skipping: No size data.")
        return

    stat_size, p_size = stats.shapiro(df['Genome_Size'])
    stat_hypo, p_hypo = stats.shapiro(df['Enhanced_Undesc_PCT'])

    print(f"   Normality (Genome Size): p={p_size:.5f}")
    print(f"   Normality (Hypothetical %): p={p_hypo:.5f}")

    if p_size > 0.05 and p_hypo > 0.05:
        print("-> Both Normal. Using Pearson Correlation.")
        test_type = "Pearson (r)"
        corr, p_val = stats.pearsonr(df['Genome_Size'], df['Enhanced_Undesc_PCT'])
    else:
        print("-> Non-Normal Data. Using Spearman Correlation.")
        test_type = "Spearman (rho)"
        corr, p_val = stats.spearmanr(df['Genome_Size'], df['Enhanced_Undesc_PCT'])

    print(f"   Correlation: {test_type}={corr:.3f}, p={p_val:.5e}")

    if p_val < 0.05:
        strength = "Strong" if abs(corr) > 0.66 else "Moderate" if abs(corr) > 0.33 else "Weak"
        direction = "Positive" if corr > 0 else "Negative"
        print(f"-> Significant {strength} {direction} correlation detected.")
    else:
        print("-> No significant correlation.")


def hypothetical_reduction(df):
    """Initial vs Enhanced hypothetical %: paired t-test if the differences pass Shapiro-Wilk, else Wilcoxon."""
    print("\n--- Hypothetical Reduction (Initial vs Enhanced) ---")

    diff = df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']

    if len(diff) < 3:
        print("Warning: Not enough data points for normality test (need >= 3).")
        shapiro_p = 1.0
    else:
        shapiro_stat, shapiro_p = stats.shapiro(diff)

    print(f"Normality Test (p-value): {shapiro_p:.5f}")

    if shapiro_p > 0.05:
        print("-> Data distribution is Normal. Using Paired t-test.")
        stat, p_val = stats.ttest_rel(df['Initial_Undesc_PCT'], df['Enhanced_Undesc_PCT'])
        test_name = "Paired t-test"
    else:
        print("-> Data distribution is Non-Normal. Using Wilcoxon Signed-Rank Test.")
        stat, p_val = stats.wilcoxon(df['Initial_Undesc_PCT'], df['Enhanced_Undesc_PCT'])
        test_name = "Wilcoxon Test"

    print(f"Test Used: {test_name}")
    print(f"Statistic: {stat:.2f}")
    print(f"P-Value:   {p_val:.5e}")
    print(f"Average Hypothetical Reduction: {diff.mean():.2f}%")

    return test_name, p_val, diff.mean()


def per_database_increase(df):
    """Percentage increase in total GO, COG, KEGG and EC entries (Initial -> Enhanced)."""
    print("\n" + "="*60)
    print("FUNCTIONAL GAINS: PERCENTAGE INCREASE PER DB (Initial -> Enhanced)")
    print("="*60)

    metrics = ['GO_entries', 'COG_entries', 'KEGG_entries', 'EC_entries']
    pct_increases = {}

    for metric in metrics:
        col_init = f'Initial_{metric}'
        col_enh = f'Enhanced_{metric}'

        if col_init in df.columns and col_enh in df.columns:
            total_init = pd.to_numeric(df[col_init], errors='coerce').sum()
            total_enh = pd.to_numeric(df[col_enh], errors='coerce').sum()

            if total_init > 0:
                pct_increase = ((total_enh - total_init) / total_init) * 100
            else:
                pct_increase = 0

            pct_increases[metric.replace('_entries', '')] = pct_increase

            print(f"{metric.replace('_entries', ''):<10} | Initial: {total_init:<6} -> Enhanced: {total_enh:<6} | Increase: +{pct_increase:.2f}%")

    print("="*60)
    return pct_increases


def total_improvement(df):
    """Mean relative reduction in hypothetical % (Initial -> Enhanced); adds 'Rel_Reduction_PCT'."""
    raw_improvement = (
        (df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']) /
        df['Initial_Undesc_PCT'].replace(0, 1)
    ) * 100

    df['Rel_Reduction_PCT'] = raw_improvement.abs()

    avg_init = df['Initial_Undesc_PCT'].mean()
    avg_enh = df['Enhanced_Undesc_PCT'].mean()
    avg_total_improvement = df['Rel_Reduction_PCT'].mean()

    print(f"\nAverage Initial Hypotheticals:  {avg_init:.2f}%")
    print(f"Average Enhanced Hypotheticals: {avg_enh:.2f}%")
    print(f"TOTAL HYPOTHETICAL REDUCTION EFFICIENCY: {avg_total_improvement:.2f}%")

    return avg_total_improvement


def generate_plots(df, test_name, p_val, pct_increases, avg_total_improvement, OUTPUT_DIR):
    """Figures 1-7. Statistics compare Initial vs Enhanced only; Prokka is plotted where its columns exist."""
    metrics = ['GO_entries', 'COG_entries', 'KEGG_entries', 'EC_entries']

    TOOL_COLORS = {
        'Prokka': '#0072B2',
        'Initial': '#999999',
        'Enhanced': '#D55E00'
    }

    # Figure 1: hypothetical % per tool (boxplot)
    print("\n--- Generating Figure 1 (Comparison Boxplot) ---")
    melt_vars = ['Initial_Undesc_PCT', 'Enhanced_Undesc_PCT']
    if 'Prokka_Undesc_PCT' in df.columns:
        melt_vars.insert(0, 'Prokka_Undesc_PCT')

    df_long = pd.melt(df,
                      id_vars=['File'],
                      value_vars=melt_vars,
                      var_name='Condition', value_name='Hypothetical_Percentage')

    df_long['Condition'] = df_long['Condition'].replace({
        'Prokka_Undesc_PCT': 'Prokka',
        'Initial_Undesc_PCT': 'Initial (Bakta)',
        'Enhanced_Undesc_PCT': 'Enhanced (Pipeline)'
    })

    plt.figure(figsize=(9, 6))
    sns.boxplot(x='Condition', y='Hypothetical_Percentage', data=df_long, palette=[TOOL_COLORS['Prokka'], TOOL_COLORS['Initial'], TOOL_COLORS['Enhanced']])
    plt.title(f"Comparison of Undescribed Coding Space\n(Init vs Enhanced: {test_name}, p < {p_val:.1e})")
    plt.ylabel("Undescribed Coding Space (%)")
    plt.savefig(f"{OUTPUT_DIR}/Figure1_Hypothetical_Reduction.png", dpi=300)

    # Figure 2: total Dbxref entries per database and tool
    print("\n--- Generating Figure 2 (Functional Gains Counts) ---")
    plot_data = []
    prefixes = [('Prokka', 'Prokka'), ('Initial', 'Initial'), ('Enhanced', 'Enhanced')]

    for metric in metrics:
        clean_metric_name = metric.replace('_entries', '')

        for prefix, label in prefixes:
            col_name = f'{prefix}_{metric}'
            if col_name in df.columns:
                total_sum = pd.to_numeric(df[col_name], errors='coerce').sum()
                plot_data.append({'Metric': clean_metric_name, 'Condition': label, 'Count': total_sum})

    if plot_data:
        df_counts = pd.DataFrame(plot_data)
        plt.figure(figsize=(12, 6))
        sns.barplot(x='Metric', y='Count', hue='Condition', data=df_counts, palette=[TOOL_COLORS['Prokka'],TOOL_COLORS['Initial'],TOOL_COLORS['Enhanced']])
        plt.title("Total Functional Annotations Recovered")
        plt.ylabel("Total Count")
        plt.legend(title="Tool")
        plt.savefig(f"{OUTPUT_DIR}/Figure2_Functional_Gains.jpeg", dpi=300)

    # Figure 3: contigs vs hypothetical %
    print("\n--- Generating Figure 3 (Contigs vs Hypotheticals) ---")
    if 'Contigs' in df.columns and df['Contigs'].sum() > 0:
        df_sorted = df.sort_values(by='Contigs')

        plt.figure(figsize=(10, 6))

        if 'Prokka_Undesc_PCT' in df.columns:
            plt.plot(df_sorted['Contigs'], df_sorted['Prokka_Undesc_PCT'],
                     marker='^', label='Prokka', linestyle='--', color=TOOL_COLORS['Prokka'], alpha=0.6)

        plt.plot(df_sorted['Contigs'], df_sorted['Initial_Undesc_PCT'],
                 marker='o', label='Initial (Bakta)', linestyle='-', color=TOOL_COLORS['Initial'], alpha=0.7)

        plt.plot(df_sorted['Contigs'], df_sorted['Enhanced_Undesc_PCT'],
                 marker='o', label='Enhanced', linestyle='-', color=TOOL_COLORS['Enhanced'], alpha=0.8)

        plt.title("Impact of Assembly Quality on Annotations")
        plt.xlabel("Number of Contigs")
        plt.ylabel("Undescribed Coding Space (%)")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(f"{OUTPUT_DIR}/Figure3_Contigs_vs_Hypotheticals.png", dpi=300)
    else:
        print("Skipping Figure 3: No GFF directory provided or no contigs found.")

    # Figure 4: pseudogene candidates per genome and tool
    print("\n--- Generating Figure 4 (Pseudogenes) ---")
    pseudo_cols = [
        ('Prokka_Pseudogene_candidates', 'Prokka'),
        ('Initial_Pseudogene_candidates', 'Initial'),
        ('Enhanced_Pseudogene_candidates', 'Enhanced')
    ]

    pseudo_data = []
    for _, row in df.iterrows():
        for col, label in pseudo_cols:
            if col in df.columns:
                pseudo_data.append({'File': row['File'], 'Condition': label, 'Count': row[col]})

    if pseudo_data:
        df_pseudo = pd.DataFrame(pseudo_data)
        plt.figure(figsize=(12, 6))
        sns.barplot(x='File', y='Count', hue='Condition', data=df_pseudo, palette="viridis")
        plt.title("Comparison of Detected Pseudogenes")
        plt.xticks([])
        plt.ylabel("Number of Pseudogenes")
        plt.tight_layout()
        plt.savefig(f"{OUTPUT_DIR}/Figure4_Pseudogene_Addition.png", dpi=300)

    # Figure 5: relative hypothetical reduction per genome (Initial -> Enhanced)
    print("\n--- Generating Figure 5 (Total Improvement - Init vs Enh) ---")
    if 'Rel_Reduction_PCT' not in df.columns:
        raw_improvement = ((df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']) / df['Initial_Undesc_PCT'].replace(0, 1)) * 100
        df['Rel_Reduction_PCT'] = raw_improvement.abs()

    plt.figure(figsize=(10, 6))
    sns.barplot(x='File', y='Rel_Reduction_PCT', data=df, palette="viridis")
    plt.axhline(avg_total_improvement, color='red', linestyle='--', label=f'Mean Improvement ({avg_total_improvement:.1f}%)')
    plt.title("Pipeline Efficiency: % of Hypothetical Proteins Annotated\n(Initial -> Enhanced)")
    plt.ylabel("Relative Reduction (%)")
    plt.xlabel("Genomes")
    plt.xticks([])
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/Figure5_Total_Improvement.png", dpi=300)

    # Figure 6: % increase per database (Initial -> Enhanced)
    print("\n--- Generating Figure 6 (DB % Increase - Init vs Enh) ---")
    if pct_increases:
        db_df = pd.DataFrame(list(pct_increases.items()), columns=['Database', 'Percent_Increase'])
        plt.figure(figsize=(8, 6))
        ax = sns.barplot(x='Database', y='Percent_Increase', data=db_df, palette="colorblind")
        for i in ax.containers:
            ax.bar_label(i, fmt='%.1f%%', padding=3)
        plt.title("Percentage Increase in Annotations (Initial -> Enhanced)")
        plt.ylabel("Increase (%)")
        plt.ylim(0, db_df['Percent_Increase'].max() * 1.15)
        plt.tight_layout()
        plt.savefig(f"{OUTPUT_DIR}/Figure6_DB_Percentage_Increase.png", dpi=300)

    # Figure 7: genome size vs Enhanced hypothetical %
    if 'Genome_Size_MB' in df.columns and df['Genome_Size_MB'].sum() > 0:
        print("\n--- Generating Figure 7 (Genome Size Correlation) ---")
        plt.figure(figsize=(10, 6))
        sns.regplot(x='Genome_Size_MB', y='Enhanced_Undesc_PCT', data=df,
                    scatter_kws={'color': '#D55E00', 'label': 'Enhanced'},
                    line_kws={'color': 'darkred'})

        plt.title("Genome Size vs. Remaining Hypothetical Proteins")
        plt.xlabel("Genome Size (Mb)")
        plt.ylabel("Enhanced Undescribed Coding Space (%)")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.savefig(f"{OUTPUT_DIR}/Figure7_Size_vs_Hypotheticals.png", dpi=300)
        plt.close()

    print(f"\nDone! Check the '{OUTPUT_DIR}' folder.")


def cmd_plot(args):
    """Entry point for `stats`: statistics, Figures 1-7 and pipeline_stats_summary.xlsx."""
    input_file = args.input
    output_dir = args.output
    gff_dir_arg = getattr(args, 'gff_dir', None)

    if not os.path.exists(input_file):
        print(f"Error: File {input_file} not found.")
        return

    df = load_data(input_file)

    if gff_dir_arg:
        df = enrich_data(df, gff_dir_arg)
    else:
        print("No GFF directory given (-g): skipping contig counts and genome size.")

    test_name, p_val, _ = hypothetical_reduction(df)
    pct_increases = per_database_increase(df)
    avg_total_improvement = total_improvement(df)
    correlate_size(df)

    os.makedirs(output_dir, exist_ok=True)
    generate_plots(df, test_name, p_val, pct_increases, avg_total_improvement, output_dir)

    gains_data = [{'Metric': f'Increase_{db}', 'Value': inc} for db, inc in pct_increases.items()]
    summary_stats = [
        {'Metric': 'Hypothetical_Reduction_P_Value', 'Value': p_val},
        {'Metric': 'Test_Name', 'Value': test_name},
        {'Metric': 'Avg_Total_Improvement_PCT', 'Value': avg_total_improvement},
    ]
    df_stats = pd.DataFrame(summary_stats + gains_data)

    output_path = os.path.join(output_dir, "pipeline_stats_summary.xlsx")
    with pd.ExcelWriter(output_path) as writer:
        df.to_excel(writer, sheet_name='Genomes_Data', index=False)
        df_stats.to_excel(writer, sheet_name='Summary_Metrics', index=False)

    print(f"Full statistics exported to {output_path}")
