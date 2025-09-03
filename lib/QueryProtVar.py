import requests
import csv
import sys

letterDict = {
    "Ala" : "A",
    "Arg" : "R",
    "Asn" : "N",
    "Asp" : "D",
    "Cys" : "C",
    "Glu" : "E",
    "Gln" : "Q",
    "Gly" : "G",
    "His" : "H",
    "Ile" : "I",
    "Leu" : "L",
    "Lys" : "K",
    "Met" : "M",
    "Phe" : "F",
    "Pro" : "P",
    "Ser" : "S",
    "Thr" : "T",
    "Trp" : "W",
    "Tyr" : "Y",
    "Val" : "V"
}

def encodeAAByLETTER(AA:str):
    vec = letterDict[AA]
    return vec

def query_protvar_api(vcf_file_path, output_file_path, assembly="AUTO"):
    """
    Query the ProtVar API with a VCF file and save the response as a 4-column TSV,
    limited to entries where "canonical": true.

    Args:
        vcf_file_path (str): Path to the VCF file to upload.
        output_file_path (str): Path to save the API response as TSV.
        assembly (str): Genome assembly (default "AUTO").
    """
    url = f"https://www.ebi.ac.uk/ProtVar/api/mapping/input?assembly={assembly}&idOnly=false"
    headers = {
        "accept": "application/json"
    }
    files = {
        "file": (vcf_file_path, open(vcf_file_path, "rb"), "text/vcard")
    }
    response = requests.post(url, headers=headers, files=files)
    response.raise_for_status()
    data = response.json()

    out_f = open(output_file_path, "w", newline="") if output_file_path else sys.stdout
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["Protein", "AA_Ref", "AA_Pos", "AA_Var"])
    for input_entry in data.get("content", {}).get("inputs", []):
        for mapping in input_entry.get("mappings", []):
            for gene in mapping.get("genes", []):
                for isoform in gene.get("isoforms", []):
                    if not isoform.get("canonical", False):
                        continue
                    protein = isoform.get("accession")
                    aa_ref = isoform.get("refAA")
                    aa_ref = encodeAAByLETTER(aa_ref)
                    aa_pos = isoform.get("isoformPosition")
                    aa_var = isoform.get("variantAA")
                    aa_var = encodeAAByLETTER(aa_var)
                    if all([protein, aa_ref, aa_pos, aa_var]):
                        writer.writerow([protein, aa_ref, aa_pos, aa_var])
    if output_file_path:
        out_f.close()