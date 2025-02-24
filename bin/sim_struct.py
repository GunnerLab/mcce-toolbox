#!/usr/bin/env python
import requests
import csv
import os
import time
import json

def get_script_directory():
    print(os.path.dirname(os.path.abspath(__file__)))
    return os.path.dirname(os.path.abspath(__file__))

def fetch_protein_sequence(pdb_code):
    url = "https://data.rcsb.org/graphql"
    
    query = {
        "query": f"""
        {{
            entry(entry_id: "{pdb_code}") {{
                polymer_entities {{
                    entity_poly {{
                        pdbx_seq_one_letter_code
                    }}
                }}
            }}
        }}
        """
    }
    
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=query, headers=headers)
    
    if response.status_code == 200:
        data = response.json()
        try:
            sequence = data["data"]["entry"]["polymer_entities"][0]["entity_poly"]["pdbx_seq_one_letter_code"]
            return sequence
        except (KeyError, TypeError, IndexError):
            return f"No sequence found for {pdb_code}."
    else:
        return f"Error {response.status_code}: {response.text}"

def fetch_similar_structures(protein_code):
    url = f'https://search.rcsb.org/rcsbsearch/v2/query'
    seq = fetch_protein_sequence(protein_code)
    print(seq)
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 1.0,
                "sequence_type": "protein",
                "value": seq
            }
        },
        "request_options": {
            "scoring_strategy": "sequence",
            "paginate": {
                "start": 0,
                "rows": 10000  # Requesting up to 10,000 results
            }
        },
        "return_type": "polymer_entity"
        }
    
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=query, headers=headers)
    if response.status_code == 200:
        data = response.json()
        return [entry["identifier"][:4] for entry in data.get("result_set", [])]
    else:
        print(f"Error fetching data for {protein_code}: {response.status_code}")
        return []

def read_existing_results(filename):
    if not os.path.exists(filename):
        return {}
    existing_data = {}
    with open(filename, mode='r', newline='') as file:
        reader = csv.reader(file)
        next(reader, None)  # Skip header
        for row in reader:
            existing_data[row[0]] = row[1:]
    return existing_data

def append_results(filename, json_filename, protein_code, similar_codes):
    existing_data = read_existing_results(filename)
    if protein_code in existing_data:
        print(f"{protein_code} already processed. Skipping.")
        return
    
    with open(filename, mode='a', newline='') as file:
        writer = csv.writer(file)
        if os.stat(filename).st_size == 0:
            writer.writerow(["Protein Code", "Similar Structures"])
        writer.writerow([protein_code] + similar_codes)

    json_data = {}
    if os.path.exists(json_filename):
        with open(json_filename, 'r') as json_file:
            try:
                json_data = json.load(json_file)
            except json.JSONDecodeError:
                json_data = {}
    
    json_data[protein_code] = similar_codes
    with open(json_filename, 'w') as json_file:
        json.dump(json_data, json_file, indent=4)    
    
def main():
    script_dir = get_script_directory()
    csv_filename = os.path.join(script_dir, "sim_structs.csv")
    json_filename = os.path.join(script_dir, "sim_structs.json")
    protein_code = input("Enter a 4-letter protein code: ").strip().upper()
    print(f"Fetching similar structures for {protein_code}...")
    
    similar_codes = fetch_similar_structures(protein_code)
    if similar_codes:
        append_results(csv_filename, json_filename, protein_code, similar_codes)
        print(f"Results saved for {protein_code}.")
    else:
        print(f"No matches found for {protein_code}.")
    
    time.sleep(1)  # Avoid overwhelming the server

if __name__ == "__main__":
    main()

