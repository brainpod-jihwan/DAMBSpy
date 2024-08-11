import argparse
import subprocess
import os
import logging
import gzip
from Bio import SeqIO
import pandas as pd

def run_command(command):
    """Run a shell command and handle errors."""
    logging.info(f"Running command: {command}")
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.returncode != 0:
        logging.error(f"Error running command: {command}\n{result.stderr}")
    return result

def parse_input_file(input_file):
    tax_data = []
    with open(input_file, 'r') as file:
        # Skip the header line
        next(file)
        for line in file:
            columns = line.strip().split('\t')
            name = columns[0]
            tax_id = columns[1]
            if tax_id != '9606':  # Exclude taxonomy ID 9606
                tax_data.append((name, tax_id))
    return tax_data

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def generate_and_run_commands(tax_data, output_dir, human_db, sample):
    for name, tax_id in tax_data:
        output_path = os.path.join(output_dir, 'cds', 'temp_cds', tax_id)
        target_dir = os.path.join(output_dir, 'cds', 'target_dir')
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        command = (
            f"ncbi-genome-download --taxids {tax_id} --formats cds-fasta --type-materials type "
            f"--output-folder {output_path} --verbose bacteria"
        )

        logging.info(f"Downloading for {name} with Taxonomy ID {tax_id}")

        try:
            run_command(command)
            logging.info(f"Successfully downloaded for {name} with Taxonomy ID {tax_id}")
            process_downloaded_files(name, tax_id, output_path, target_dir)
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to download for {name} with Taxonomy ID {tax_id}: {e}")

    align_command = (
        f"puppy-align -pr {output_dir}/cds/target_dir -nt {human_db} -o {output_dir}/MMseq2 "
        f"> {output_dir}/{sample}.germs_uniq.log"
    )
    run_command(align_command)

    primer_command = (
        f"puppy-primers -pr {output_dir}/cds/target_dir -i {output_dir}/MMseq2/ResultDB.tsv "
        f"--primers_type unique -o {output_dir}/{sample} 2> {output_dir}/{sample}.germs_primer.log"
    )
    run_command(primer_command)

def process_downloaded_files(name, tax_id, output_path, target_dir):
    seq_records = []
    for root, dirs, files in os.walk(output_path):
        for file in files:
            if file.endswith('.fna.gz'):
                gz_file_path = os.path.join(root, file)
                with gzip.open(gz_file_path, 'rt') as f_in:
                    for record in SeqIO.parse(f_in, 'fasta'):
                        seq_records.append([record.description, str(record.seq), gz_file_path])
    
    df = pd.DataFrame(seq_records, columns=['header', 'sequence', 'original_file_path'])
    original_count = len(df)
    df_unique = df.drop_duplicates(subset='sequence')
    unique_count = len(df_unique)
    duplicate_count = original_count - unique_count

    for _, row in df_unique.iterrows():
        cleaned_name = name.replace(' ', '_').replace('.', '_')
        cleaned_file_name = os.path.basename(row['original_file_path'])[:-3].replace('.', '_')
        output_file_name = f"{cleaned_name}_{tax_id}_{cleaned_file_name}.fna"
        output_file_path = os.path.join(target_dir, output_file_name)

        with open(output_file_path, 'a') as f_out:
            f_out.write(f">{row['header']}\n{row['sequence']}\n")

    logging.info(f"Processed {original_count} sequences for Taxonomy ID {tax_id}")
    logging.info(f"Found {duplicate_count} duplicate sequences for Taxonomy ID {tax_id}")
    logging.info(f"Removed {duplicate_count} duplicates, resulting in {unique_count} unique sequences for Taxonomy ID {tax_id}")

def main():
    parser = argparse.ArgumentParser(
        description='Generate and execute genome download and processing commands.'
    )
    parser.add_argument('--input_braken', '-i', required=True, help='Input file containing taxonomy IDs.')
    parser.add_argument('--output_dir', '-o', required=True, help='Output directory for the downloads.')
    parser.add_argument('--human_db', '-hd', required=True, help='Path to the human database CDS.')
    parser.add_argument('--sample', '-s', required=True, help='Sample name.')

    args = parser.parse_args()

    setup_logging()

    tax_data = parse_input_file(args.input_braken)
    generate_and_run_commands(tax_data, args.output_dir, args.human_db, args.sample)

if __name__ == '__main__':
    main()
