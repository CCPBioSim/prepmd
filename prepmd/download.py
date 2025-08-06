#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from Bio.PDB import PDBList
import urllib.request
from os.path import sep
import requests


def get_structure(pdb_id, directory, file_format="mmCif"):
    if file_format == "mmCif":
        format_str = "cif"
    if file_format == "pdb":
        format_str = "pdb"
    urllib.request.urlretrieve("https://files.rcsb.org/download/" +
                               pdb_id+"."+format_str, directory+sep+pdb_id+"."+format_str)
    return directory+sep+pdb_id+"."+format_str


def get_uniprot_sequence(pdb_id, merge_sequence=True, write_to_file=None,
                         verbose=False):
    query = '''{
      entries(entry_ids:["'''+pdb_id+'''"]){
        polymer_entities {
          rcsb_id
          rcsb_polymer_entity_container_identifiers {
            reference_sequence_identifiers {
              database_accession
              database_name
            }
          }
        }
      }
    }'''
    data = requests.get('https://data.rcsb.org/graphql?query='+query)
    accessions = {}
    for item in data.json()['data']['entries']:
        entities = item['polymer_entities']
        for entity in entities:
            if entity['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers']:
                for identifier in entity['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers']:
                    accessions[int(entity["rcsb_id"].split("_")[1])
                               ] = identifier['database_accession']
            else:
                accessions[int(entity["rcsb_id"].split("_")[1])] = None  # rna?
    if verbose:
        print("Accessions: "+str(accessions))
    sequences = {}
    for key, value in accessions.items():
        if value:
            url = "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=(accession:" + \
                value+")"
            seq = requests.get(url)
            sequences[key] = seq.text
        else:
            sequences[key] = None

    if not merge_sequence:
        return sequences
    names = []
    combined_sequences = []
    for key, sequence in sequences.items():
        if sequence is None:
            continue
        names.append(sequence.split("\n")[0].replace(
            "sp|", "").replace(">", ""))
        seq_noheader = sequence.split("\n")[1:]
        combined_sequences.append("".join(seq_noheader).replace("\n", ""))

    fasta = ">P1;"+pdb_id+"\n" + \
        ";".join(names)+"\n"+"/".join(combined_sequences)
    if write_to_file:
        with open(write_to_file, "w") as file:
            file.write(fasta)
    return fasta
