###P. Raimondeau _ Generates a bed file of up/dowstream region for each 'feature' entry in a gff file
import pandas as pd
import pybedtools
import re
import sys

def parse_gff(gff_file):
    gff_data = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            if feature_type == 'mRNA':
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                # Extract gene name or ID from the attributes field
                attributes = fields[8]
                match = re.search(r'ID=([^;]+)', attributes)
                gene_id = match.group(1) if match else "unknown"
                gff_data.append([chrom, start, end, strand, gene_id])
    return pd.DataFrame(gff_data, columns=['chrom', 'start', 'end', 'strand', 'gene_id'])

def create_bed(df, upstream=2000, downstream=2000):
    bed_entries = []
    for _, row in df.iterrows():
        if row['strand'] == '+':
            upstream_region_start = max(0, row['start'] - upstream)
            upstream_region_end = row['start']
            downstream_region_start = row['end']
            downstream_region_end = row['end'] + downstream
        else:
            upstream_region_start = row['end']
            upstream_region_end = row['end'] + upstream
            downstream_region_start = max(0, row['start'] - downstream)
            downstream_region_end = row['start']

        bed_entries.append([row['chrom'], upstream_region_start, upstream_region_end, f"{row['gene_id']}_upstream", '.', row['strand']])
        bed_entries.append([row['chrom'], downstream_region_start, downstream_region_end, f"{row['gene_id']}_downstream", '.', row['strand']])
    
    return pd.DataFrame(bed_entries, columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])

def main(gff_file, output_bed):
    gff_df = parse_gff(gff_file)
    bed_df = create_bed(gff_df)
    # Convert DataFrame to BedTool object and save to a BED file
    bed = pybedtools.BedTool.from_dataframe(bed_df)
    bed.saveas(output_bed)
    print(f"BED file saved to {output_bed}")

if __name__ == "__main__":
    gff_file = sys.argv[1]  # Replace with your input GFF file path
    output_bed = sys.argv[1]+'.bed'  # Replace with your desired output BED file path
    main(gff_file, output_bed)
