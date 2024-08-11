import argparse
import subprocess
import os
import logging

def run_command(command):
    """Run a shell command and handle errors."""
    logging.info(f"Running command: {command}")
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.returncode != 0:
        logging.error(f"Error running command: {command}\n{result.stderr}")
    return result

def main(args):
    sample = args.sample
    human_db = args.human_db
    germs_db = args.germs_db
    output_dir = args.output_dir
    fastq1 = args.fastq1
    fastq2 = args.fastq2
    read_len = args.read_len
    min_qual = args.min_qual
    min_hit = args.min_hit
    thread_num = args.thread_num

    os.makedirs(output_dir, exist_ok=True)

    # BWA MEM alignment
    logging.info("Starting BWA MEM alignment...")
    bwa_command = (
        f"bwa mem {human_db} {fastq1} {fastq2} "
        f"-M -t {thread_num} > {output_dir}/{sample}.sam 2> {output_dir}/{sample}.human_align.log"
    )
    run_command(bwa_command)

    # Extract unmapped reads from SAM to BAM
    logging.info("Extracting unmapped reads from SAM to BAM...")
    samtools_view_command = (
        f"samtools view -@ {thread_num} -f 12 -F 256 "
        f"-o {output_dir}/{sample}.unmapped.bam "
        f"-bhS {output_dir}/{sample}.sam 2> {output_dir}/{sample}.human_filter.log"
    )
    run_command(samtools_view_command)

    # Remove raw SAM file
    logging.info("Removing raw SAM file...")
    os.remove(f"{output_dir}/{sample}.sam")

    # Query sort BAM file
    logging.info("Query sorting BAM file...")
    samtools_sort_command = (
        f"samtools sort --threads {thread_num} -n "
        f"-o {output_dir}/{sample}.unmapped.sorted.bam --output-fmt BAM "
        f"{output_dir}/{sample}.unmapped.bam 2>> {output_dir}/{sample}.human_filter.log"
    )
    run_command(samtools_sort_command)

    # Convert BAM to FASTQ
    logging.info("Converting BAM to FASTQ...")
    bedtools_command = (
        f"bedtools bamtofastq -i {output_dir}/{sample}.unmapped.sorted.bam "
        f"-fq {output_dir}/{sample}.unmapped_R1.fastq -fq2 {output_dir}/{sample}.unmapped_R2.fastq"
    )
    run_command(bedtools_command)

    # Log read counts
    logging.info("Logging read counts...")
    with open(f"{output_dir}/{sample}.human_filter.log", "a") as log_file:
        log_file.write(f"ORIGINAL RAW READ: {run_command(f'wc -l {fastq1}').stdout}")
        log_file.write(f"ORIGINAL RAW READ: {run_command(f'wc -l {output_dir}/{sample}.unmapped_R1.fastq').stdout}")

    # Kraken bacterial genome mapping
    logging.info("Starting Kraken bacterial genome mapping...")
    kraken_command = (
        f"kraken2 --paired --db {germs_db} --threads {thread_num} "
        f"--minimum-base-quality {min_qual} --report {output_dir}/{sample}.k2report "
        f"--report-minimizer-data --minimum-hit-groups {min_hit} "
        f"{output_dir}/{sample}.unmapped_R1.fastq {output_dir}/{sample}.unmapped_R2.fastq "
        f"> {output_dir}/{sample}.kraken2 2> {output_dir}/{sample}.germs_align.log"
    )
    run_command(kraken_command)

    # Abundance profiling using Bracken
    logging.info("Starting abundance profiling using Bracken...")
    bracken_command = (
        f"bracken -d {germs_db} -i {output_dir}/{sample}.k2report -r {read_len} -l S "
        f"-t {thread_num} -o {output_dir}/{sample}.braken -w {output_dir}/{sample}.breport "
        f"2> {output_dir}/{sample}.germs_annot.log"
    )
    run_command(bracken_command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove human read from orignal NGS dataset and alignment to bacterial genome pipeline.",
        epilog="Example: python dambs-align.py -s SAMPLE_PREFIX -hd HUMAN_GENOME_DIR/human_genome.fa.gz -gd GERMS_GENOME_DIR -o OUTPUT_DIR -f1 test_R1.fastq -f2 test_R2.fastq -rl 100 -mq 10 -mh 3 -t 8"
    )
    parser.add_argument('--sample', '-s', required=True, help="Sample prefiex")
    parser.add_argument('--human_db', '-hd', required=True, help="Directory for human reference database")
    parser.add_argument('--germs_db', '-gd', required=True, help="Directory for germs reference database")
    parser.add_argument('--output_dir', '-o', required=True, help="Output directory")
    parser.add_argument('--fastq1', '-f1', required=True, help="Path to raw FASTQ file 1")
    parser.add_argument('--fastq2', '-f2', required=True, help="Path to raw FASTQ file 2")
    parser.add_argument('--read_len', '-rl', type=int, required=True, help="Read length")
    parser.add_argument('--min_qual', '-mq', type=int, required=True, help="Minimum quality for Kraken")
    parser.add_argument('--min_hit', '-mh', type=int, required=True, help="Minimum hit groups for Kraken")
    parser.add_argument('--thread_num', '-t', type=int, required=True, help="Number of threads")

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    args = parser.parse_args()
    main(args)

