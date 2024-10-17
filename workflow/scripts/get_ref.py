import requests
from bs4 import BeautifulSoup
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--organism', required=True)
    parser.add_argument('--outfile', required=True)
    return parser.parse_args()


def get_assembly_summary(organism_name):
    """Download and parse the assembly summary for RefSeq genomes."""
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    response = requests.get(url)

    if response.status_code == 200:
        lines = response.text.splitlines()
        headers = lines[1].split("\t")
        organism_col = headers.index("organism_name")
        ftp_col = headers.index("ftp_path")

        # Search for organism_name in the assembly summary
        for line in lines[2:]:
            columns = line.split("\t")
            if organism_name.lower() in columns[organism_col].lower():
                return columns[ftp_col]
    else:
        raise Exception(f"Failed to download assembly summary, status code: {response.status_code}")


def download_genome_file(ftp_url, outfile):
    """Download the full genomic .fna.gz file from the NCBI FTP directory."""
    response = requests.get(ftp_url)

    if response.status_code == 200:
        # Parse the FTP directory HTML to find .fna.gz files that contain the full genome
        soup = BeautifulSoup(response.text, 'html.parser')
        file_link = None

        # Look for the link to the full genome .fna.gz file
        for link in soup.find_all('a'):
            href = link.get('href')
            # Ensure that it ends with "_genomic.fna.gz" and excludes "_cds_from_genomic"
            if href.endswith("_genomic.fna.gz") and "cds_from_genomic" not in href:
                file_link = href
                break
        
        if file_link:
            full_file_url = f"{ftp_url}/{file_link}"

            # Download the full genomic .fna.gz file
            print(f"Downloading {outfile} from {full_file_url}...")
            with open(outfile, "wb") as f:
                file_response = requests.get(full_file_url, stream=True)
                if file_response.status_code == 200:
                    for chunk in file_response.iter_content(chunk_size=1024):
                        if chunk:  # Filter out keep-alive chunks
                            f.write(chunk)
                else:
                    raise Exception(f"Failed to download file, status code: {file_response.status_code}")
        else:
            raise Exception(f"No valid _genomic.fna.gz file found in the FTP directory.")
    else:
        raise Exception(f"Failed to access the FTP directory, status code: {response.status_code}")


def main():
    args = parse_args()
    print(f"Searching for the reference genome of {args.organism}...")
    ftp_url = get_assembly_summary(args.organism)
    print(f"Found FTP path: {ftp_url}")
    download_genome_file(ftp_url, args.outfile)


if __name__ == "__main__":
    main()