#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ZA Primary Enzyme Index & MycoDeg Phylogeny Index Calculator
=============================================================
Based on GAI index method (Lai et al, 2023)

Formula (2-1) ZA Primary Enzyme Index:
  1. Convert to relative abundance: rel = data / rowSums(data)
  2. p_num = number of present (>0) MZ genera in sample
     n_num = number of present (>0) MM genera in sample
  3. numerator   = sum(rel[j] * (p_num / |MZ|)) for j in MZ + beta
     denominator = sum(rel[j] * (n_num / |MM|)) for j in MM + beta
  4. ZA_score = log10(numerator / denominator)

Formula (2-2) MycoDeg Phylogeny Index:
  Same logic with MA (positive) and MS (negative) groups
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import logging
from datetime import datetime

# ============================================================
# Configuration
# ============================================================
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

# *_t families (not single enzymes, need sub-family aggregation)
TOTAL_FAMILIES = {'GH55_t', 'PL1_t', 'AA3_t', 'GH51_t', 'GH43_t', 'GH13_t', 'GH16_t'}

BETA = 1e-3  # smoothing factor to avoid log(0)

# ============================================================
# Data Loading
# ============================================================
def load_weights(filepath):
    """Load gene family weight file, return {family_name: estimate_value}"""
    for enc in ['utf-8', 'gbk', 'gb18030', 'latin1']:
        try:
            df = pd.read_csv(filepath, sep='\t', encoding=enc, header=None)
            break
        except UnicodeDecodeError:
            continue
    else:
        df = pd.read_csv(filepath, sep='\t', encoding='utf-8', errors='replace', header=None)

    # Find the row containing 'genus' and 'Estimate' as header
    header_row = None
    for i, row in df.iterrows():
        row_str = ' '.join(str(v) for v in row.values)
        if 'genus' in row_str.lower() and 'estimate' in row_str.lower():
            header_row = i
            break

    if header_row is not None:
        df.columns = df.iloc[header_row].astype(str).str.strip().values
        df = df.iloc[header_row + 1:].reset_index(drop=True)
    else:
        df = pd.read_csv(filepath, sep='\t', encoding=enc, header=0)

    col_map = {}
    for col in df.columns:
        cl = str(col).lower().strip()
        if 'genus' in cl:
            col_map[col] = 'genus'
        elif 'estimate' in cl:
            col_map[col] = 'Estimate'
    df = df.rename(columns=col_map)

    df = df.dropna(subset=['genus', 'Estimate'])
    df['Estimate'] = pd.to_numeric(df['Estimate'], errors='coerce')
    df = df.dropna(subset=['Estimate'])
    df['genus'] = df['genus'].astype(str).str.strip()
    return dict(zip(df['genus'], df['Estimate']))


def load_strain_data(filepath, preserve_class=True):
    """
    Load strain data, return (DataFrame, class_col_name_or_None, class_data_or_None)
    DataFrame: index = strain name, columns = gene families (numeric)
    """
    for enc in ['utf-8', 'gbk', 'gb18030', 'latin1']:
        try:
            df = pd.read_csv(filepath, sep='\t', encoding=enc)
            break
        except UnicodeDecodeError:
            continue
    else:
        df = pd.read_csv(filepath, sep='\t', encoding='utf-8', errors='replace')

    name_col = df.columns[0]
    df = df.set_index(name_col)
    df.index = df.index.astype(str).str.strip()

    # Detect class column (second column with values A/B/C/D)
    class_col = None
    class_data = None
    if preserve_class and len(df.columns) >= 1:
        first_col = df.columns[0]
        vals = df[first_col].astype(str).str.strip().str.upper().unique()
        if len(vals) <= 6 and all(v in ('A', 'B', 'C', 'D', 'E', 'F', '') for v in vals) and len(vals) >= 2:
            class_col = first_col
            class_data = df[class_col].copy()
            df = df.drop(columns=[class_col])

    # Convert all columns to numeric
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

    return df, class_col, class_data


# ============================================================
# *_t Family Handling
# ============================================================
def resolve_total_family_value(family_name, strain_row, all_columns):
    """
    Handle *_t families:
    1. If *_t column exists and value > 0, use it directly
    2. Otherwise, sum all sub-family columns matching the prefix
    Returns (value, method)
    """
    if family_name in strain_row.index:
        val = strain_row[family_name]
        if pd.notna(val) and int(val) > 0:
            return int(val), 'direct'

    prefix = family_name.replace('_t', '')
    matching_cols = [c for c in all_columns if c.startswith(prefix) and c != family_name]
    if matching_cols:
        total = sum(int(strain_row.get(c, 0)) for c in matching_cols if c in strain_row.index)
        if total > 0:
            return int(total), 'computed'
    return 0, 'not_found'


# ============================================================
# Index Calculation (GAI-based, corrected)
# ============================================================
def compute_index(strain_row_raw, weights, all_columns, index_name='Index', verbose=True):
    """
    Compute GAI-based index following the R script:

    1. Convert to relative abundance
    2. For positive group (Estimate > 0, i.e. MZ or MA):
       - p_num = count of positive-group families present (>0) in this sample
       - numerator = sum( rel[j] * (p_num / |positive_total|) ) for j in positive + beta
    3. For negative group (Estimate < 0, i.e. MM or MS):
       - n_num = count of negative-group families present (>0) in this sample
       - denominator = sum( rel[j] * (n_num / |negative_total|) ) for j in negative + beta
    4. Score = log10(numerator / denominator)

    Parameters:
        strain_row_raw: raw abundance Series (not yet relative)
        weights: {family_name: estimate_value}
        all_columns: list of all column names
        index_name: label for logging
        verbose: whether to log details

    Returns:
        (score, matched_dict, unmatched_list)
    """
    # Step 1: Convert to relative abundance
    raw_values = {}
    for col in all_columns:
        if col in strain_row_raw.index:
            raw_values[col] = int(strain_row_raw[col]) if pd.notna(strain_row_raw[col]) else 0
        else:
            raw_values[col] = 0

    # Handle *_t families - compute their values first
    resolved_values = dict(raw_values)
    for family in weights:
        if family in TOTAL_FAMILIES:
            val, method = resolve_total_family_value(family, strain_row_raw, all_columns)
            if val > 0:
                resolved_values[family] = val
                if method == 'computed' and verbose:
                    logger.info(f"  {family}: computed from sub-families = {val}")

    # Compute row sum for relative abundance
    row_sum = sum(resolved_values.values())
    if row_sum == 0:
        logger.warning(f"  {index_name}: all values are zero, cannot compute relative abundance")
        return 0.0, {}, list(weights.keys())

    # Relative abundance
    rel_values = {k: v / row_sum for k, v in resolved_values.items()}

    # Separate positive/negative groups
    positive_families = [f for f, e in weights.items() if e > 0]
    negative_families = [f for f, e in weights.items() if e < 0]

    matched = {}
    unmatched = []

    # Step 2: Positive group
    p_num = 0  # count of present positive families
    positive_sum_weighted = 0.0
    for family in positive_families:
        if family in rel_values and rel_values[family] > 0:
            p_num += 1
            positive_sum_weighted += rel_values[family]
            matched[family] = resolved_values[family]
        elif family not in rel_values or resolved_values.get(family, 0) == 0:
            # Check if it's a *_t family that might have been resolved
            if family in resolved_values and resolved_values[family] > 0:
                p_num += 1
                positive_sum_weighted += rel_values[family]
                matched[family] = resolved_values[family]
            else:
                unmatched.append(family)
                if verbose:
                    logger.warning(f"  {index_name}: '{family}' not found or zero")

    # Weight by presence ratio
    numerator = positive_sum_weighted * (p_num / len(positive_families)) + BETA if len(positive_families) > 0 else BETA

    # Step 3: Negative group
    n_num = 0  # count of present negative families
    negative_sum_weighted = 0.0
    for family in negative_families:
        if family in rel_values and rel_values[family] > 0:
            n_num += 1
            negative_sum_weighted += rel_values[family]
            matched[family] = resolved_values[family]
        else:
            if family in resolved_values and resolved_values[family] > 0:
                n_num += 1
                negative_sum_weighted += rel_values[family]
                matched[family] = resolved_values[family]
            else:
                unmatched.append(family)
                if verbose:
                    logger.warning(f"  {index_name}: '{family}' not found or zero")

    # Weight by presence ratio
    denominator = negative_sum_weighted * (n_num / len(negative_families)) + BETA if len(negative_families) > 0 else BETA

    # Step 4: log10 ratio
    score = np.log10(numerator / denominator)

    if verbose:
        logger.info(f"  {index_name}: p_num={p_num}/{len(positive_families)}, n_num={n_num}/{len(negative_families)}")
        logger.info(f"  {index_name}: numerator={numerator:.6f}, denominator={denominator:.6f}")
        logger.info(f"  {index_name}: score = log10({numerator:.6f}/{denominator:.6f}) = {score:.4f}")

    return score, matched, unmatched


def compute_index_for_dataframe(df, weights, verbose=False):
    """Compute index for all strains in a DataFrame"""
    all_columns = df.columns.tolist()
    results = {}
    for strain_name in df.index:
        strain_row = df.loc[strain_name]
        score, matched, unmatched = compute_index(
            strain_row, weights, all_columns,
            index_name=strain_name, verbose=verbose
        )
        results[strain_name] = {
            'score': score,
            'matched': matched,
            'unmatched': unmatched,
            'matched_count': len(matched),
            'total_families': len(weights),
            'matched_ratio': len(matched) / len(weights) if len(weights) > 0 else 0
        }
    return results


# ============================================================
# Visualization
# ============================================================
def plot_formula_2_1(input_name, input_score, dikaryon_results, non_dikaryon_results,
                     save_path):
    """Formula (2-1) visualization: input strain vs all dikaryon + non-dikaryon"""
    all_scores = {}
    for name, res in non_dikaryon_results.items():
        all_scores[name] = res['score']

    dikaryon_scores = {}
    for name, res in dikaryon_results.items():
        all_scores[name] = res['score']
        dikaryon_scores[name] = res['score']

    all_scores[input_name] = input_score

    sorted_strains = sorted(all_scores.items(), key=lambda x: x[1], reverse=True)
    rank = next(i + 1 for i, (name, _) in enumerate(sorted_strains) if name == input_name)

    dikaryon_min_score = min(dikaryon_scores.values())
    is_at_bottom = (input_score <= dikaryon_min_score)

    fig, ax = plt.subplots(figsize=(14, max(8, len(sorted_strains) * 0.16)))

    names = [s[0] for s in sorted_strains]
    scores = [s[1] for s in sorted_strains]

    colors = []
    for name in names:
        if name == input_name:
            colors.append('#E74C3C')
        elif name in dikaryon_scores:
            colors.append('#3498DB')
        else:
            colors.append('#95A5A6')

    y_pos = range(len(names))
    ax.barh(y_pos, scores, color=colors, height=0.7)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel('ZA Primary Enzyme Index Score', fontsize=12)
    ax.set_title(
        f'Formula (2-1) ZA Primary Enzyme Index Ranking\n'
        f'{input_name} Score: {input_score:.4f} | Rank: {rank}/{len(sorted_strains)}',
        fontsize=14, fontweight='bold')

    legend_patches = [
        mpatches.Patch(color='#E74C3C', label=f'{input_name} (Input)'),
        mpatches.Patch(color='#3498DB', label='Dikaryon'),
        mpatches.Patch(color='#95A5A6', label='Non-dikaryon'),
    ]
    ax.legend(handles=legend_patches, loc='lower right', fontsize=10)

    input_idx = names.index(input_name)
    ax.annotate(f'>> {input_score:.4f}', xy=(scores[input_idx], input_idx),
                xytext=(10, 0), textcoords='offset points', fontsize=9,
                fontweight='bold', color='#E74C3C')

    if is_at_bottom:
        ax.text(0.5, 0.02, '[!] Possibly lacking primary enzymes', transform=ax.transAxes,
                fontsize=16, color='red', fontweight='bold',
                ha='center', va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#FEE', edgecolor='red', alpha=0.8))

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

    return rank, is_at_bottom


def plot_formula_2_2(input_name, input_score, dikaryon_results, dikaryon_df,
                     class_col, class_data, save_path):
    """Formula (2-2) visualization: input strain vs dikaryon (with class labels)"""
    strain_scores = {}
    strain_classes = {}
    for name, res in dikaryon_results.items():
        strain_scores[name] = res['score']
        if class_col and class_data is not None and name in class_data.index:
            strain_classes[name] = str(class_data[name]).strip().upper()

    strain_scores[input_name] = input_score

    sorted_strains = sorted(strain_scores.items(), key=lambda x: x[1], reverse=True)
    rank = next(i + 1 for i, (name, _) in enumerate(sorted_strains) if name == input_name)

    class_means = {}
    for cls in ['A', 'B', 'C', 'D']:
        cls_scores = [res['score'] for n, res in dikaryon_results.items()
                      if strain_classes.get(n) == cls]
        if cls_scores:
            class_means[cls] = np.mean(cls_scores)

    closest_class = min(class_means.keys(), key=lambda c: abs(class_means[c] - input_score)) if class_means else 'N/A'

    fig, ax = plt.subplots(figsize=(14, max(8, len(sorted_strains) * 0.16)))

    names = [s[0] for s in sorted_strains]
    scores = [s[1] for s in sorted_strains]

    class_colors = {'A': '#E74C3C', 'B': '#2ECC71', 'C': '#F39C12', 'D': '#9B59B6'}

    colors = []
    for name in names:
        if name == input_name:
            colors.append('#E74C3C')
        else:
            cls = strain_classes.get(name, '')
            colors.append(class_colors.get(cls, '#95A5A6'))

    y_pos = range(len(names))
    bars = ax.barh(y_pos, scores, color=colors, height=0.7, alpha=0.8)

    input_idx = names.index(input_name)
    bars[input_idx].set_alpha(1.0)
    bars[input_idx].set_edgecolor('black')
    bars[input_idx].set_linewidth(2)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel('MycoDeg Phylogeny Index Score', fontsize=12)
    ax.set_title(
        f'Formula (2-2) MycoDeg Phylogeny Index Ranking\n'
        f'{input_name} Score: {input_score:.4f} | Rank: {rank}/{len(sorted_strains)} | Closest: Class {closest_class}',
        fontsize=14, fontweight='bold')

    legend_patches = [mpatches.Patch(color='#E74C3C', label=f'{input_name} (Input)')]
    for cls in sorted(class_colors.keys()):
        count = sum(1 for v in strain_classes.values() if v == cls)
        if count > 0:
            legend_patches.append(mpatches.Patch(color=class_colors[cls], label=f'Class {cls} ({count} strains)'))
    ax.legend(handles=legend_patches, loc='lower right', fontsize=10)

    ax.annotate(f'>> {input_score:.4f}', xy=(scores[input_idx], input_idx),
                xytext=(10, 0), textcoords='offset points', fontsize=9,
                fontweight='bold', color='#E74C3C')

    for cls, mean_score in class_means.items():
        ax.axvline(x=mean_score, color=class_colors.get(cls, 'gray'),
                   linestyle='--', alpha=0.5, linewidth=1)
        ax.text(mean_score, -1.5, f'Class {cls} mean\n{mean_score:.3f}', fontsize=7,
                color=class_colors.get(cls, 'gray'), ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

    return rank, closest_class


# ============================================================
# Main
# ============================================================
def main():
    import argparse
    parser = argparse.ArgumentParser(description='ZA Primary Enzyme Index & MycoDeg Phylogeny Index Calculator')
    parser.add_argument('input_files', nargs='+', help='Input file path(s)')
    parser.add_argument('--data-dir', default=None, help='Reference data directory (default: data/ next to script)')
    parser.add_argument('--output-dir', default=None, help='Output directory (default: result/ next to script)')
    args = parser.parse_args()

    # Default directories relative to script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = args.data_dir or os.path.join(script_dir, 'data')
    output_dir = args.output_dir or os.path.join(script_dir, 'result')
    os.makedirs(output_dir, exist_ok=True)

    # Load weight files
    logger.info("=" * 60)
    logger.info("Loading reference data...")
    weights_za = load_weights(os.path.join(data_dir, '\u4e3b\u9176\u8d4b\u503c.txt'))
    weights_phy = load_weights(os.path.join(data_dir, '\u7cfb\u7edf\u53d1\u80b2\u8d4b\u503c.txt'))
    logger.info(f"  Formula (2-1) ZA Primary Enzyme Index: {len(weights_za)} gene families")
    logger.info(f"    MZ group (Estimate>0): {sum(1 for v in weights_za.values() if v > 0)} families")
    logger.info(f"    MM group (Estimate<0): {sum(1 for v in weights_za.values() if v < 0)} families")
    logger.info(f"  Formula (2-2) MycoDeg Phylogeny Index: {len(weights_phy)} gene families")
    logger.info(f"    MA group (Estimate>0): {sum(1 for v in weights_phy.values() if v > 0)} families")
    logger.info(f"    MS group (Estimate<0): {sum(1 for v in weights_phy.values() if v < 0)} families")
    logger.info(f"  Beta (smoothing): {BETA}")

    # Load reference strain data
    non_dikaryon_df, _, _ = load_strain_data(os.path.join(data_dir, '\u975e\u53cc\u6838.txt'), preserve_class=False)
    dikaryon_df, class_col, class_data = load_strain_data(os.path.join(data_dir, '\u53cc\u6838.txt'), preserve_class=True)

    logger.info(f"  Non-dikaryon strains: {len(non_dikaryon_df)}")
    logger.info(f"  Dikaryon strains: {len(dikaryon_df)}" + (f" (class column: {class_col})" if class_col else " (no class column!)"))

    if class_data is not None:
        logger.info(f"  Dikaryon class distribution: {dict(class_data.value_counts())}")

    # Pre-compute reference strain indices
    logger.info("Pre-computing reference strain indices...")
    non_dikaryon_za = compute_index_for_dataframe(non_dikaryon_df, weights_za, verbose=False)
    non_dikaryon_phy = compute_index_for_dataframe(non_dikaryon_df, weights_phy, verbose=False)
    dikaryon_za = compute_index_for_dataframe(dikaryon_df, weights_za, verbose=False)
    dikaryon_phy = compute_index_for_dataframe(dikaryon_df, weights_phy, verbose=False)
    logger.info("  Pre-computation done")

    # Process each input file
    overall_summary = []

    for input_file in args.input_files:
        input_basename = os.path.splitext(os.path.basename(input_file))[0]
        logger.info("=" * 60)
        logger.info(f"Processing input file: {input_file}")

        # Create output subdirectory
        strain_output_dir = os.path.join(output_dir, input_basename)
        figures_dir = os.path.join(strain_output_dir, 'figures')
        os.makedirs(figures_dir, exist_ok=True)

        # Load input data
        input_df, _, _ = load_strain_data(input_file, preserve_class=False)
        if len(input_df) == 0:
            logger.error(f"Input file {input_file} has no data!")
            continue

        strain_name = input_df.index[0]
        strain_row = input_df.iloc[0]
        all_columns = input_df.columns.tolist()

        summary_lines = []
        summary_lines.append(f"{'=' * 60}")
        summary_lines.append(f"Strain: {strain_name}")
        summary_lines.append(f"Input file: {os.path.basename(input_file)}")
        summary_lines.append(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        summary_lines.append(f"Beta: {BETA}")
        summary_lines.append(f"")

        # ========================================
        # Formula (2-1) ZA Primary Enzyme Index
        # ========================================
        logger.info("-" * 40)
        logger.info("Computing Formula (2-1) ZA Primary Enzyme Index...")

        za_score, za_matched, za_unmatched = compute_index(
            strain_row, weights_za, all_columns, index_name='Formula(2-1)', verbose=True
        )

        logger.info(f"  ZA Primary Enzyme Index score: {za_score:.4f}")
        logger.info(f"  Matched gene families: {len(za_matched)}/{len(weights_za)} ({len(za_matched)/len(weights_za)*100:.1f}%)")
        if za_unmatched:
            logger.info(f"  Unmatched gene families ({len(za_unmatched)}): {', '.join(za_unmatched)}")

        # Visualization
        za_rank, is_at_bottom = plot_formula_2_1(
            strain_name, za_score, dikaryon_za, non_dikaryon_za,
            os.path.join(figures_dir, 'formula_2_1_ranking.png')
        )
        total_strains = len(non_dikaryon_za) + len(dikaryon_za) + 1
        logger.info(f"  Rank: {za_rank}/{total_strains}")
        if is_at_bottom:
            logger.info(f"  [!] Score at bottom of dikaryon strains - possibly lacking primary enzymes!")

        summary_lines.append(f"[Formula (2-1) ZA Primary Enzyme Index]")
        summary_lines.append(f"  Score: {za_score:.4f}")
        summary_lines.append(f"  Matched families: {len(za_matched)}/{len(weights_za)} ({len(za_matched)/len(weights_za)*100:.1f}%)")
        if za_unmatched:
            summary_lines.append(f"  Unmatched families ({len(za_unmatched)}): {', '.join(za_unmatched)}")
        else:
            summary_lines.append(f"  Unmatched families: none")
        summary_lines.append(f"  Overall rank: {za_rank}/{total_strains}")
        if is_at_bottom:
            summary_lines.append(f"  ** Possibly lacking primary enzymes **")
        else:
            summary_lines.append(f"  Not at bottom of dikaryon strains")
        summary_lines.append(f"")

        # ========================================
        # Formula (2-2) MycoDeg Phylogeny Index
        # ========================================
        logger.info("-" * 40)
        logger.info("Computing Formula (2-2) MycoDeg Phylogeny Index...")

        phy_score, phy_matched, phy_unmatched = compute_index(
            strain_row, weights_phy, all_columns, index_name='Formula(2-2)', verbose=True
        )

        logger.info(f"  MycoDeg Phylogeny Index score: {phy_score:.4f}")
        logger.info(f"  Matched gene families: {len(phy_matched)}/{len(weights_phy)} ({len(phy_matched)/len(weights_phy)*100:.1f}%)")
        if phy_unmatched:
            logger.info(f"  Unmatched gene families ({len(phy_unmatched)}): {', '.join(phy_unmatched)}")

        # Visualization
        phy_rank, closest_class = plot_formula_2_2(
            strain_name, phy_score, dikaryon_phy, dikaryon_df,
            class_col, class_data, os.path.join(figures_dir, 'formula_2_2_ranking.png')
        )
        dikaryon_total = len(dikaryon_phy) + 1
        logger.info(f"  Dikaryon rank: {phy_rank}/{dikaryon_total}")
        logger.info(f"  Closest class: {closest_class}")

        summary_lines.append(f"[Formula (2-2) MycoDeg Phylogeny Index]")
        summary_lines.append(f"  Score: {phy_score:.4f}")
        summary_lines.append(f"  Matched families: {len(phy_matched)}/{len(weights_phy)} ({len(phy_matched)/len(weights_phy)*100:.1f}%)")
        if phy_unmatched:
            summary_lines.append(f"  Unmatched families ({len(phy_unmatched)}): {', '.join(phy_unmatched)}")
        else:
            summary_lines.append(f"  Unmatched families: none")
        summary_lines.append(f"  Dikaryon rank: {phy_rank}/{dikaryon_total}")
        summary_lines.append(f"  Closest class: {closest_class}")
        summary_lines.append(f"")

        # Save per-strain summary
        summary_path = os.path.join(strain_output_dir, f'{input_basename}_summary.txt')
        with open(summary_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(summary_lines))
        logger.info(f"Summary saved: {summary_path}")

        # Collect overall summary
        overall_summary.append({
            'Strain': strain_name,
            'F2-1_Matched': len(za_matched),
            'F2-1_Total': len(weights_za),
            'F2-1_Ratio': f"{len(za_matched)/len(weights_za)*100:.1f}%",
            'F2-1_Unmatched': ', '.join(za_unmatched) if za_unmatched else 'none',
            'F2-1_Score': f"{za_score:.4f}",
            'F2-1_Rank': f"{za_rank}/{total_strains}",
            'F2-1_LackingPrimary': 'Yes' if is_at_bottom else 'No',
            'F2-2_Matched': len(phy_matched),
            'F2-2_Total': len(weights_phy),
            'F2-2_Ratio': f"{len(phy_matched)/len(weights_phy)*100:.1f}%",
            'F2-2_Unmatched': ', '.join(phy_unmatched) if phy_unmatched else 'none',
            'F2-2_Score': f"{phy_score:.4f}",
            'F2-2_DikaryonRank': f"{phy_rank}/{dikaryon_total}",
            'F2-2_ClosestClass': f"{closest_class}"
        })

    # Save overall summary
    if overall_summary:
        overall_df = pd.DataFrame(overall_summary)
        overall_path = os.path.join(output_dir, 'overall_summary.txt')
        with open(overall_path, 'w', encoding='utf-8') as f:
            f.write(f"ZA Primary Enzyme Index & MycoDeg Phylogeny Index - Calculation Summary\n")
            f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Beta: {BETA}\n")
            f.write(f"{'=' * 80}\n\n")
            for _, row in overall_df.iterrows():
                f.write(f"Strain: {row['Strain']}\n")
                f.write(f"  Formula (2-1) ZA Primary Enzyme Index:\n")
                f.write(f"    Matched: {row['F2-1_Matched']}/{row['F2-1_Total']} ({row['F2-1_Ratio']})\n")
                f.write(f"    Unmatched: {row['F2-1_Unmatched']}\n")
                f.write(f"    Score: {row['F2-1_Score']}\n")
                f.write(f"    Rank: {row['F2-1_Rank']}\n")
                f.write(f"    Lacking primary enzymes: {row['F2-1_LackingPrimary']}\n")
                f.write(f"  Formula (2-2) MycoDeg Phylogeny Index:\n")
                f.write(f"    Matched: {row['F2-2_Matched']}/{row['F2-2_Total']} ({row['F2-2_Ratio']})\n")
                f.write(f"    Unmatched: {row['F2-2_Unmatched']}\n")
                f.write(f"    Score: {row['F2-2_Score']}\n")
                f.write(f"    Dikaryon rank: {row['F2-2_DikaryonRank']}\n")
                f.write(f"    Closest class: {row['F2-2_ClosestClass']}\n")
                f.write(f"\n")

        csv_path = os.path.join(output_dir, 'overall_summary.csv')
        overall_df.to_csv(csv_path, index=False, encoding='utf-8-sig')

        logger.info("=" * 60)
        logger.info(f"Overall summary saved: {overall_path}")
        logger.info(f"CSV summary saved: {csv_path}")

    logger.info("=" * 60)
    logger.info("Calculation complete!")


if __name__ == '__main__':
    main()
