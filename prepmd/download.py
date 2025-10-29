#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download data from the PDB and UNIPROT
"""

#from Bio.PDB import PDBList
import urllib.request
from os.path import sep
import requests


def get_structure(pdb_id, directory, file_format="mmCif"):
    """
    Download a structure from the PDB.
    Args:
        pdb_id: id of the pdb to download, a string
        directory: directory to download the file into, a string
        file_format: mmCif or pdb, a string
    returns:
        path to the downloaded file.
    """
    
    if file_format == "mmCif" or file_format == "cif":
        format_str = "cif"
    if file_format == "pdb":
        format_str = "pdb"
    try:
        url = "https://files.rcsb.org/download/"+pdb_id+"."+format_str
        destination = directory+sep+pdb_id+"."+format_str
        urllib.request.urlretrieve(url, destination)
    except urllib.error.HTTPError as e:
        r = requests.get(url.replace(".pdb", ".cif"))
        if r.status_code == 200:
            msg = "No PDB for "+pdb_id+" exists (but an mmcif structure does). "
            "Run with --fmt cif to use it."
            raise IOError(msg)
        else:
            raise e
    return directory+sep+pdb_id+"."+format_str


def get_uniprot_sequence(pdb_id, merge_sequence=True, write_to_file=None,
                         verbose=False):
    """
    For a given pdb id, find the fasta sequence for all chains from UNIPROT.
    Args:
        pdb_id: the id of the pdb, a string
        merge_sequence: whether to merge the sequences together into a single
        fasta file, a bool
        write_to_file: path of file to write to, a string
        verbose: a bool, whether to write debug info out
    Returns:
        the fasta sequences as a dictionry keyed by pdb id, or a string
        with the fasta sequence.
    """
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
