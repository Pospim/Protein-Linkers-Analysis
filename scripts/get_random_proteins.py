import requests
import random
import os

workdir = "~/Desktop/work/protein_linkers"

# Expand the ~ in the path
workdir = os.path.expanduser(workdir)

# Ensure the input_2 directory exists
output_dir = os.path.join(workdir, "input_2")
os.makedirs(output_dir, exist_ok=True)

# 1) Fetch a big list of accessions from UniProt
url = "https://rest.uniprot.org/uniprotkb/stream"
params = {
    "compressed": "false",
    "format": "list",                # one accession per line
    "query": "(* AND reviewed:true)",# all reviewed proteins
    "size": 5000
}

resp = requests.get(url, params=params)
resp.raise_for_status()

accessions = resp.text.strip().splitlines()
accessions = list(dict.fromkeys(accessions))  # remove duplicates, keep order

# 2) Shuffle and take about 1000 random proteins
random.shuffle(accessions)
subset = accessions[:1000]

# 3) Fetch organism for each accession in a single batch query
info_url = "https://rest.uniprot.org/uniprotkb/search"
query = " OR ".join(f"accession:{a}" for a in subset)

info_params = {
    "query": query,
    "format": "tsv",
    "fields": "accession,organism_name",
    "size": 1000
}

info_resp = requests.get(info_url, params=info_params)
info_resp.raise_for_status()

lines = info_resp.text.strip().splitlines()
header = lines[0].split("\t")
organism_index = header.index("Organism")

# Map accession → organism
acc_to_org = {}
for line in lines[1:]:
    parts = line.split("\t")
    if len(parts) >= 2:
        acc_to_org[parts[0]] = parts[1]

# 4) Write output
output_file = os.path.join(output_dir, "random_uniprot_accessions.tsv")

with open(output_file, "w") as f:
    f.write("accession\torganism\n")
    for acc in subset:
        org = acc_to_org.get(acc, "UNKNOWN")
        f.write(f"{acc}\t{org}\n")

print(f"Wrote {len(subset)} accessions to {output_file}")
