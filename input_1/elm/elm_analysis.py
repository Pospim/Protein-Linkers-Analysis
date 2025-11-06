#!/usr/bin/env python3
"""
ELM Analysis Script with Checkpointing
Run this in a tmux session for long-running ELM analysis that can survive laptop closure.

Usage:
    python run_elm_analysis.py [OPTIONS]

Options:
    --resume                    Resume from checkpoint
    --workdir PATH              Working directory (default: ~/Desktop/work/protein_linkers/input_2)
    --max-workers N             Number of parallel workers (default: 4)
    --checkpoint-freq N         Save checkpoint every N regions (default: 10)
    --rate-limit SECONDS        Seconds between ELM requests (default: 60)
    --filter-confidence         Only return high-confidence motif classes (LIG, DOC, MOD, TRG)
    --help                      Show this help message
"""

import requests
import time
import json
import os
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict, Counter
import threading
from datetime import datetime

# ================================
# Configuration (will be set from args)
# ================================
WORKDIR = None
JSON_FILE = None
FASTA_FILE = None
CHECKPOINT_FILE = None
FINAL_OUTPUT_JSON = None
FINAL_OUTPUT_CSV = None
LOG_FILE = None

# ELM Configuration (will be set from args)
ELM_API_URL = "http://elm.eu.org/start_search/"
MAX_SEQUENCE_LENGTH = 2000
ELM_RATE_LIMIT = 60  # seconds between requests (configurable)
MAX_WORKERS = 4  # Conservative for rate-limited API (configurable)

# High confidence classes
ALLOWED_CLASSES = {"LIG", "DOC", "MOD", "TRG"}

# Thread-safe rate limiting
_elm_lock = threading.Lock()
_last_elm_request = 0

# Checkpoint saving frequency (configurable)
CHECKPOINT_FREQUENCY = 10


# ================================
# Logging
# ================================
def log(message):
    """Log message to both console and file with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_line = f"[{timestamp}] {message}"
    print(log_line)
    with open(LOG_FILE, 'a') as f:
        f.write(log_line + "\n")


# ================================
# ELM Analysis Functions
# ================================
def _rate_limited_elm_request(sequence):
    """Submit ELM request with thread-safe rate limiting."""
    global _last_elm_request

    with _elm_lock:
        current_time = time.time()
        time_since_last = current_time - _last_elm_request

        if time_since_last < ELM_RATE_LIMIT:
            wait_time = ELM_RATE_LIMIT - time_since_last
            # Don't log every wait, too verbose
            time.sleep(wait_time)

        _last_elm_request = time.time()

    # Make the actual request outside the lock
    url = f"{ELM_API_URL}{sequence}"
    headers = {"Accept": "text/tab-separated-values"}

    try:
        response = requests.get(url, headers=headers, timeout=30)
        if response.status_code == 200:
            return response.text
        elif response.status_code == 429:
            log(f"Rate limit hit (429), waiting 60s and retrying...")
            time.sleep(60)
            response = requests.get(url, headers=headers, timeout=30)
            if response.status_code == 200:
                return response.text
        response.raise_for_status()
    except Exception as e:
        log(f"ELM request failed: {e}")
        return None


def _parse_elm_tsv(tsv_text, filter_for_confidence=False):
    """Parse ELM TSV response and return motif information."""
    if not tsv_text:
        return {}

    motifs = {}

    for line in tsv_text.strip().split("\n")[1:]:  # Skip header
        if not line.strip():
            continue

        cols = line.split('\t')
        if len(cols) < 9:
            continue

        elm_id = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        is_filtered = cols[5].lower() == "true"

        # Skip filtered motifs
        if is_filtered:
            continue

        # Filter by allowed classes if requested
        if filter_for_confidence and (elm_id.split("_")[0] not in ALLOWED_CLASSES):
            continue

        motifs[elm_id] = {
            'start': start,
            'end': end
        }

    return motifs


def analyze_region_elm(region_id, sequence, region_type, filter_for_confidence=False):
    """Analyze a single region (domain or linker) with ELM."""
    if len(sequence) > MAX_SEQUENCE_LENGTH:
        log(f"Skipping {region_id}: sequence too long ({len(sequence)} > {MAX_SEQUENCE_LENGTH})")
        return None

    tsv_result = _rate_limited_elm_request(sequence)
    if not tsv_result:
        return None

    motifs = _parse_elm_tsv(tsv_result, filter_for_confidence=filter_for_confidence)

    return {
        'region_id': region_id,
        'region_type': region_type,
        'sequence_length': len(sequence),
        'motifs': motifs,
        'motif_count': len(motifs)
    }


# ================================
# Checkpointing
# ================================
def save_checkpoint(results, completed_regions):
    """Save current progress to checkpoint file."""
    checkpoint_data = {
        'timestamp': datetime.now().isoformat(),
        'completed_regions': list(completed_regions),
        'results': results
    }

    # Write to temporary file first, then rename (atomic operation)
    temp_file = CHECKPOINT_FILE + ".tmp"
    with open(temp_file, 'w') as f:
        json.dump(checkpoint_data, f, indent=2)
    os.rename(temp_file, CHECKPOINT_FILE)

    log(f"Checkpoint saved: {len(completed_regions)} regions completed")


def load_checkpoint():
    """Load checkpoint if it exists."""
    if not os.path.exists(CHECKPOINT_FILE):
        return None, set()

    try:
        with open(CHECKPOINT_FILE, 'r') as f:
            checkpoint_data = json.load(f)

        results = checkpoint_data.get('results', {})
        completed_regions = set(checkpoint_data.get('completed_regions', []))

        log(f"Checkpoint loaded: {len(completed_regions)} regions already completed")
        log(f"Checkpoint timestamp: {checkpoint_data.get('timestamp', 'unknown')}")

        return results, completed_regions
    except Exception as e:
        log(f"Error loading checkpoint: {e}")
        return None, set()


# ================================
# Main Analysis
# ================================
def analyze_domains_and_linkers_elm(proteins_dict, sequences_dict, resume=False, filter_for_confidence=False):
    """
    Perform ELM analysis on all domains and linkers with checkpointing.

    Parameters:
    -----------
    proteins_dict : dict
        Dictionary with protein accessions and domain/linker info
    sequences_dict : dict
        Dictionary with protein sequences
    resume : bool
        If True, resume from checkpoint
    filter_for_confidence : bool
        If True, only return high-confidence motif classes
    """
    # Load checkpoint if resuming
    if resume:
        results, completed_regions = load_checkpoint()
        if results is None:
            log("No checkpoint found, starting from scratch")
            results = {
                'domains': {},
                'n-terminus': {},
                'c-terminus': {},
                'inner': {}
            }
            completed_regions = set()
    else:
        results = {
            'domains': {},
            'n-terminus': {},
            'c-terminus': {},
            'inner': {}
        }
        completed_regions = set()

    # Prepare all tasks
    tasks = []

    for accession, data in proteins_dict.items():
        seq = sequences_dict.get(accession)
        if not seq:
            log(f"Warning: No sequence found for {accession}")
            continue

        # Add domain tasks
        for idx, (domain_name, start, end) in enumerate(data['domains'], 1):
            domain_seq = seq[start-1:end]  # Convert to 0-indexed
            region_id = f"{accession}_domain_{idx}_{domain_name}"

            # Skip if already completed
            if region_id not in completed_regions:
                tasks.append((region_id, domain_seq, 'domain'))

        # Add linker tasks
        for idx, (linker_type, start, end) in enumerate(data['linkers'], 1):
            linker_seq = seq[start-1:end]  # Convert to 0-indexed
            region_id = f"{accession}_linker_{idx}_{linker_type}"

            # Skip if already completed
            if region_id not in completed_regions:
                tasks.append((region_id, linker_seq, linker_type))

    log(f"\n{'='*60}")
    log(f"ELM ANALYSIS: {len(tasks)} regions to analyze (out of {len(tasks) + len(completed_regions)} total)")
    if completed_regions:
        log(f"Resuming: {len(completed_regions)} regions already completed")
    if filter_for_confidence:
        log("Filtering for high-confidence motifs only (LIG, DOC, MOD, TRG)")
    log(f"{'='*60}\n")

    if len(tasks) == 0:
        log("All regions already completed!")
        return results

    # Execute tasks with thread pool
    completed = len(completed_regions)
    total = len(tasks) + len(completed_regions)
    new_completions = 0

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(analyze_region_elm, region_id, seq, region_type, filter_for_confidence): (region_id, region_type)
            for region_id, seq, region_type in tasks
        }

        # Process completed tasks
        for future in as_completed(future_to_task):
            region_id, region_type = future_to_task[future]
            completed += 1
            new_completions += 1

            try:
                result = future.result()
                if result:
                    # Categorize result
                    if region_type == 'domain':
                        results['domains'][region_id] = result
                    else:
                        results[region_type][region_id] = result

                    # Mark as completed
                    completed_regions.add(region_id)

                    log(f"[{completed}/{total}] ✓ {region_id}: {result['motif_count']} motifs found")
                else:
                    log(f"[{completed}/{total}] ✗ {region_id}: failed")

            except Exception as e:
                log(f"[{completed}/{total}] ✗ {region_id}: error - {e}")

            # Save checkpoint periodically
            if new_completions > 0 and new_completions % CHECKPOINT_FREQUENCY == 0:
                save_checkpoint(results, completed_regions)

    # Final checkpoint save
    save_checkpoint(results, completed_regions)

    # Print summary
    log(f"\n{'='*60}")
    log("ELM ANALYSIS SUMMARY")
    log(f"{'='*60}")
    log(f"Domains analyzed:      {len(results['domains'])}")
    log(f"N-terminus linkers:    {len(results['n-terminus'])}")
    log(f"C-terminus linkers:    {len(results['c-terminus'])}")
    log(f"Inner linkers:         {len(results['inner'])}")

    # Count total motifs
    total_motifs = sum(
        r['motif_count']
        for category in results.values()
        for r in category.values()
    )
    log(f"Total motifs found:    {total_motifs}")
    log(f"{'='*60}\n")

    return results


def save_elm_results(elm_results, output_json, output_csv=None):
    """Save ELM results to files."""
    # Save JSON
    with open(output_json, 'w') as f:
        json.dump(elm_results, f, indent=2)
    log(f"✓ Saved ELM results to {output_json}")

    # Save CSV if requested
    if output_csv:
        with open(output_csv, 'w') as f:
            f.write("region_id,region_type,sequence_length,motif_id,motif_start,motif_end\n")

            for region_type, regions in elm_results.items():
                for region_id, region_data in regions.items():
                    for motif_id, motif_info in region_data['motifs'].items():
                        f.write(f"{region_id},{region_type},{region_data['sequence_length']},"
                               f"{motif_id},{motif_info['start']},{motif_info['end']}\n")

        log(f"✓ Saved ELM results to {output_csv}")


# ================================
# Main
# ================================
def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='ELM Analysis Script with Checkpointing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic run
  python run_elm_analysis.py

  # Resume from checkpoint
  python run_elm_analysis.py --resume

  # Custom settings
  python run_elm_analysis.py --max-workers 8 --checkpoint-freq 20

  # Different working directory
  python run_elm_analysis.py --workdir ~/my_proteins

  # High confidence motifs only with faster checkpoints
  python run_elm_analysis.py --filter-confidence --checkpoint-freq 5
        '''
    )

    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume from checkpoint if available'
    )

    parser.add_argument(
        '--workdir',
        type=str,
        default='~/Desktop/work/protein_linkers/input_2',
        help='Working directory containing input files (default: ~/Desktop/work/protein_linkers/input_2)'
    )

    parser.add_argument(
        '--max-workers',
        type=int,
        default=4,
        help='Maximum number of parallel workers (default: 4)'
    )

    parser.add_argument(
        '--checkpoint-freq',
        type=int,
        default=10,
        help='Save checkpoint every N completed regions (default: 10)'
    )

    parser.add_argument(
        '--rate-limit',
        type=int,
        default=60,
        help='Seconds to wait between ELM requests (default: 60)'
    )

    parser.add_argument(
        '--filter-confidence',
        action='store_true',
        help='Only return high-confidence motif classes (LIG, DOC, MOD, TRG)'
    )

    return parser.parse_args()


def setup_config(args):
    """Setup global configuration from arguments."""
    global WORKDIR, JSON_FILE, FASTA_FILE, CHECKPOINT_FILE
    global FINAL_OUTPUT_JSON, FINAL_OUTPUT_CSV, LOG_FILE
    global ELM_RATE_LIMIT, MAX_WORKERS, CHECKPOINT_FREQUENCY

    WORKDIR = os.path.expanduser(args.workdir)
    JSON_FILE = os.path.join(WORKDIR, "formatted_proteins.json")
    FASTA_FILE = os.path.join(WORKDIR, "proteins.fa")
    CHECKPOINT_FILE = os.path.join(WORKDIR, "elm_checkpoint.json")
    FINAL_OUTPUT_JSON = os.path.join(WORKDIR, "elm_results.json")
    FINAL_OUTPUT_CSV = os.path.join(WORKDIR, "elm_results.csv")
    LOG_FILE = os.path.join(WORKDIR, "elm_analysis.log")

    ELM_RATE_LIMIT = args.rate_limit
    MAX_WORKERS = args.max_workers
    CHECKPOINT_FREQUENCY = args.checkpoint_freq

    # Create working directory if it doesn't exist
    os.makedirs(WORKDIR, exist_ok=True)


def main():
    args = parse_args()
    setup_config(args)

    log(f"\n{'='*60}")
    log("ELM ANALYSIS SCRIPT - STARTING")
    log(f"Working directory: {WORKDIR}")
    log(f"Resume mode: {args.resume}")
    log(f"Max workers: {MAX_WORKERS}")
    log(f"Checkpoint frequency: every {CHECKPOINT_FREQUENCY} regions")
    log(f"Rate limit: {ELM_RATE_LIMIT} seconds between requests")
    log(f"Filter confidence: {args.filter_confidence}")
    log(f"{'='*60}\n")    # Load formatted proteins
    log(f"Loading proteins from {JSON_FILE}")
    with open(JSON_FILE, 'r') as f:
        proteins_dict = json.load(f)
    log(f"✓ Loaded {len(proteins_dict)} proteins")

    # Load FASTA sequences
    log(f"Loading sequences from {FASTA_FILE}")
    sequences_dict = {}
    with open(FASTA_FILE, 'r') as f:
        current_acc = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_acc:
                    sequences_dict[current_acc] = ''.join(current_seq)
                current_acc = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_acc:
            sequences_dict[current_acc] = ''.join(current_seq)
    log(f"✓ Loaded {len(sequences_dict)} sequences")

    # Run ELM analysis
    elm_results = analyze_domains_and_linkers_elm(
        proteins_dict,
        sequences_dict,
        resume=args.resume,
        filter_for_confidence=args.filter_confidence
    )

    # Save final results
    save_elm_results(elm_results, FINAL_OUTPUT_JSON, FINAL_OUTPUT_CSV)

    log(f"\n{'='*60}")
    log("ELM ANALYSIS COMPLETE!")
    log(f"Results saved to:")
    log(f"  - {FINAL_OUTPUT_JSON}")
    log(f"  - {FINAL_OUTPUT_CSV}")
    log(f"  - {LOG_FILE}")
    log(f"{'='*60}\n")


if __name__ == "__main__":
    main()
