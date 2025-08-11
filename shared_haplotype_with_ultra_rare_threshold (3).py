from __future__ import annotations

import argparse
from typing import List, Set
import pandas as pd
from cyvcf2 import VCF
from tqdm import tqdm

# --------------------------------------------------------------------------- #
# Helper functions
# --------------------------------------------------------------------------- #

def variant_carriers(v) -> Set[int]:
    return {i for i, gts in enumerate(v.genotypes) if any(allele > 0 for allele in gts[:2])}

def all_shared_homozygous(v, shared_inds: Set[int]) -> bool:
    gt = v.gt_types
    return all(gt[i] == 3 for i in shared_inds)

def overlaps_gene_region(chrom: str, start: int, end: int, regions: list[tuple[str, int, int]]) -> bool:
    chrom = chrom.lstrip("chr")
    for region_chrom, region_start, region_end in regions:
        if chrom == region_chrom and not (end < region_start or start > region_end):
            return True
    return False

def is_ultra_rare_variant(variant, gnomad_idx: int, af_threshold:float) -> bool:
    csq_entries = variant.INFO.get('CSQ', '')
    if not csq_entries:
        return False

    for entry in csq_entries.split(','):
        fields = entry.split('|')
        if len(fields) <= gnomad_idx:
            continue
        val = fields[gnomad_idx]
        if val in ('', '.', None):
            return True
        try:
            af = float(val)
            if af <= af_threshold:
                return True
        except ValueError:
            continue
    return False

# --------------------------------------------------------------------------- #
# Core routine to flush (sub-)clusters
# --------------------------------------------------------------------------- #

def flush_cluster(
    variants: List,
    prev_nonshared_pos: int | None,
    next_nonshared_pos: int | None,
    *,
    min_individuals: int,
    min_variants: int,
    max_individuals: int,
    min_segment_length: int,
    max_segment_length: int,
    sample_ids: List[str],
    gene_regions: list[tuple[str, int, int]],
    gnomad_idx: int,
    min_ultra_rare: int,
    af_threshold: float
):
    if not variants:
        return

    individual_sets = [variant_carriers(v) for v in variants]
    total_variants = len(individual_sets)

    from collections import Counter
    flat_list = [ind for s in individual_sets for ind in s]
    count = Counter(flat_list)
    threshold = int(0.8 * total_variants)
    shared_inds = {i for i, c in count.items() if c >= threshold}
    n_shared = len(shared_inds)

    if n_shared < min_individuals:
        return
    if max_individuals is not None and n_shared > max_individuals:
        return

    chrom = variants[0].CHROM
    inds_str = ",".join(sample_ids[i] for i in sorted(shared_inds))

    # --- Segment by genotype mode ----------------------------------------- #
    blocks: List[tuple[int, int, str]] = []
    cur_mode: str | None = None
    block_start = 0

    for idx, v in enumerate(variants):
        mode = "double" if all_shared_homozygous(v, shared_inds) else "single"
        if cur_mode is None:
            cur_mode = mode
        elif mode != cur_mode:
            blocks.append((block_start, idx - 1, cur_mode))
            block_start = idx
            cur_mode = mode
    blocks.append((block_start, len(variants) - 1, cur_mode))

    for b_start, b_end, mode in blocks:
        block_vars = variants[b_start : b_end + 1]
        segment_length = block_vars[-1].POS - block_vars[0].POS

        if (len(block_vars) < min_variants or
            (min_segment_length is not None and segment_length < min_segment_length) or
            (max_segment_length is not None and segment_length > max_segment_length)):
            continue

        start = block_vars[0].POS
        end = block_vars[-1].POS

        if not overlaps_gene_region(chrom, start, end, gene_regions):
            continue

        n_ultra_rare = sum(is_ultra_rare_variant(v, gnomad_idx, af_threshold) for v in block_vars)
        if n_ultra_rare < min_ultra_rare:
            continue

        left_flank = variants[b_start - 1].POS if b_start > 0 else prev_nonshared_pos
        right_flank = variants[b_end + 1].POS if b_end + 1 < len(variants) else next_nonshared_pos
        err_left = None if left_flank is None else start - left_flank
        err_right = None if right_flank is None else right_flank - end

        print(
            f"{chrom}\t{start}\t{end}\t{len(block_vars)}\t{len(shared_inds)}\t{inds_str}"
            f"\t{mode}\t{err_left}\t{err_right}\t{n_ultra_rare}\tULTRA_RARE"
        )

# --------------------------------------------------------------------------- #
# Main driver
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        description="Identify shared ultra-rare variant segments from a pre-filtered pVCF."
    )
    parser.add_argument("vcf", help="Input pVCF (bgzipped + indexed)")
    parser.add_argument("--cluster-distance", type=int, default=2000)
    parser.add_argument("--min-individuals", type=int, default=2)
    parser.add_argument("--max-individuals", type=int, default=None)
    parser.add_argument("--min-variants", type=int, default=8)
    parser.add_argument("--min-segment-length", type=int, default=0)
    parser.add_argument("--max-segment-length", type=int, default=None)
    parser.add_argument("--min-ultra-rare", type=int, default=1,
                        help="Minimum number of ultra-rare variants required in segment")
    parser.add_argument("--ultra-rare-threshold", type=float, default=0.0001,
                    help="Max GnomAD AF to consider a variant ultra-rare (default: 0.0001)")
    parser.add_argument("--regions", type=str, default="gene_data.csv")

    args = parser.parse_args()

    gene_df = pd.read_csv(args.regions)
    BUFFER = 10000  # Add 10kb on either side
    gene_regions = [
    (
        str(row["chromosome_name"]).lstrip("chr"),
        max(0, int(row["start_position"]) - BUFFER),
        int(row["end_position"]) + BUFFER
    )
    for _, row in gene_df.iterrows()
]


    vcf = VCF(args.vcf)
    samples = vcf.samples

    csq_header = vcf.get_header_type('CSQ')['Description']
    csq_fields = csq_header.split('Format: ')[1].strip('"').split('|')

    try:
        gnomad_idx = csq_fields.index('GnomAD_v4_1_AF_popmax')
    except ValueError as e:
        raise RuntimeError(f"Required CSQ field 'GnomAD_v4_1_AF_popmax' not found in VCF header.") from e

    current_cluster: List = []
    prev_nonshared_pos: int | None = None
    last_pos: int | None = None

    pbar = tqdm(vcf, desc="Scanning variants", unit="variants")
    for v in pbar:
        if last_pos is None:
            current_cluster.append(v)
            last_pos = v.POS
            continue

        if v.CHROM != current_cluster[0].CHROM:
            flush_cluster(
                current_cluster,
                prev_nonshared_pos,
                None,
                min_individuals=args.min_individuals,
                max_individuals=args.max_individuals,
                min_variants=args.min_variants,
                min_segment_length=args.min_segment_length,
                max_segment_length=args.max_segment_length,
                sample_ids=samples,
                gene_regions=gene_regions,
                gnomad_idx=gnomad_idx,
                min_ultra_rare=args.min_ultra_rare, 
                af_threshold=args.ultra_rare_threshold
            )
            current_cluster = [v]
            last_pos = v.POS
            prev_nonshared_pos = None
            continue

        distance = v.POS - last_pos
        if distance <= args.cluster_distance:
            current_cluster.append(v)
            last_pos = v.POS
        else:
            flush_cluster(
                current_cluster,
                prev_nonshared_pos,
                v.POS,
                min_individuals=args.min_individuals,
                max_individuals=args.max_individuals,
                min_variants=args.min_variants,
                min_segment_length=args.min_segment_length,
                max_segment_length=args.max_segment_length,
                sample_ids=samples,
                gene_regions=gene_regions,
                gnomad_idx=gnomad_idx,
                min_ultra_rare=args.min_ultra_rare,
                af_threshold=args.ultra_rare_threshold
            )
            prev_nonshared_pos = current_cluster[-1].POS
            current_cluster = [v]
            last_pos = v.POS

    flush_cluster(
        current_cluster,
        prev_nonshared_pos,
        None,
        min_individuals=args.min_individuals,
        max_individuals=args.max_individuals,
        min_variants=args.min_variants,
        min_segment_length=args.min_segment_length,
        max_segment_length=args.max_segment_length,
        sample_ids=samples,
        gene_regions=gene_regions,
        gnomad_idx=gnomad_idx,
        min_ultra_rare=args.min_ultra_rare,
        af_threshold=args.ultra_rare_threshold
    )

if __name__ == "__main__":
    main()
