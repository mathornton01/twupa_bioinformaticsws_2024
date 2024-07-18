import argparse
import requests

def query_ensembl(snp_id):
    base_url = f"https://rest.ensembl.org/variation/human/{snp_id}?"
    headers = {"Content-Type": "application/json"}
    response = requests.get(base_url, headers=headers)
    response.raise_for_status()
    return response.json()

def parse_snp_info(snp_data):
    location = snp_data['mappings'][0]['location']
    chromosome, position_range = location.split(":")
    position = position_range.split("-")[0]
    reference_base = snp_data['mappings'][0]['ancestral_allele']
    alleles = snp_data['mappings'][0]['allele_string'].split("/")
    possible_alleles = alleles
    
    return {
        "chromosome": chromosome,
        "position": position,
        "reference_base": reference_base,
        "alleles": possible_alleles
    }

def main():
    parser = argparse.ArgumentParser(description="Retrieve SNP information from Ensembl")
    parser.add_argument("snp_id", help="The SNP ID (e.g., rs12203592)")
    args = parser.parse_args()

    try:
        snp_data = query_ensembl(args.snp_id)
        snp_info = parse_snp_info(snp_data)
        print(f"SNP ID: {args.snp_id}")
        print(f"Chromosome: {snp_info['chromosome']}")
        print(f"Position (GRCh38): {snp_info['position']}")
        print(f"Reference Base: {snp_info['reference_base']}")
        print(f"Possible Alleles: {', '.join(snp_info['alleles'])}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()

