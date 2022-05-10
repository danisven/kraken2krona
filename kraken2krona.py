#!/usr/bin/env python3

__version__ = "0.1.0"

import argparse
import logging
import pandas as pd
import gzip
from stringmeup.taxonomy import TaxonomyTree
from dataclasses import dataclass
from io import TextIOWrapper
from os import path

# The different taxonomical levels that we are interested in. Values are corresponding attributes in the TaxonInfo dataclass.
taxonomy_ranks = {
    'superkingdom': 't_superkingdom', 
    'kingdom': 't_kingdom', 
    'phylum': 't_phylum', 
    'class': 't_class', 
    'order': 't_order', 
    'family': 't_family', 
    'genus': 't_genus', 
    'species': 't_species'}

@dataclass
class TaxonInfo:
    # Number of classified reads
    hits_at_clade: int
    hits_at_taxon: int
    
    # Taxonomy
    t_id: int
    t_superkingdom: str
    t_kingdom: str
    t_phylum: str
    t_class: str
    t_order: str
    t_family: str
    t_genus: str
    t_species: str

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d [%H:%M:%S]')
log = logging.getLogger(path.basename(__file__))

def get_taxonomy(tax_id: int, tt: TaxonomyTree) -> dict[str: str]:
    """
    For a given tax ID, find its lineage. Saves the lineage information in 'taxonomy_dict' which
    is also returned.

    taxonomy_dict can be used when creating an object of class TaxonInfo.
    """
    # Set up the taxonomy_dict
    taxonomy_dict = {taxonomic_rank: None for taxonomic_rank in taxonomy_ranks.values()}
    
    # The unclassified node (special case since TaxonomyTree can't process tax ID 0)
    if tax_id == 0:
        return taxonomy_dict

    lineage = tt.get_lineage([tax_id])[tax_id]

    # Get the rank of each tax_id in lineage, and if that rank is of interest we save the scientific name of the tax_id
    # in taxonomy_dict
    for lineage_tax_id in lineage:
        lineage_tax_id_rank = tt.get_rank([lineage_tax_id])[lineage_tax_id]  # Taxonomic rank of current tax ID

        if lineage_tax_id_rank in taxonomy_ranks:
            taxon_info_rank_name = tt.get_name([lineage_tax_id])[lineage_tax_id]  # Scientific name of current tax ID
            taxon_info_rank = taxonomy_ranks[lineage_tax_id_rank]  # Taxonomic rank (in terms of TaxonInfo attribute names)
            taxonomy_dict[taxon_info_rank] = taxon_info_rank_name

    return taxonomy_dict
    
def parse_report(report_file_object: TextIOWrapper, taxon_info_dict: dict[int: TaxonInfo], tt: TaxonomyTree) -> dict[int: TaxonInfo]:
    """
    Parses a Kraken 2 style report line by line. Adds the number of reads classified to 
    each tax_id to a dictionary where keys are tax IDs and values are number of classified reads.

    Accepts a pre-filled dictionary. If there's a record of a tax ID in the report file currently 
    being parsed that is already in the dictionary, their values will be summed. This 
    is for when summing many report files instead of just a single one.
    """
    
    log.info('Reading report file  "{file}".'.format(file=report_file_object.name))

    for line in report_file_object:
        line = line.split()
        hits_at_clade = int(line[1])
        hits_at_taxon = int(line[2])
        tax_id = int(line[4])

        if tax_id in taxon_info_dict:
            taxon_info_dict[tax_id].hits_at_clade += hits_at_clade
            taxon_info_dict[tax_id].hits_at_taxon += hits_at_taxon
        else:
            # Fetch the taxonomy for the taxon
            taxon_dict = get_taxonomy(tax_id, tt)

            # Add tax_id and info on number of hits to the dict
            taxon_dict['t_id'] = tax_id
            taxon_dict['hits_at_clade'] = hits_at_clade
            taxon_dict['hits_at_taxon'] = hits_at_taxon

            # Create a TaxonInfo object with data in taxon_dict, and att it to the datastructure
            taxon_info_dict[tax_id] = TaxonInfo(**taxon_dict)
        
    return taxon_info_dict

def write_output(krona_output_filename: str, result_dict: dict) -> None:
    """
    Writes the contents of result_dict to file.
    """

    log.info('Writing output to "{file}".'.format(file=krona_output_filename))


def read_file(filename):
    """
    Wrapper to read either gzipped or ordinary text file input.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def get_arguments() -> argparse.Namespace:
    """
    Wrapper function to get the command line arguments. Inserting this piece of code
    into its own function for conda compatibility.
    """

    parser = argparse.ArgumentParser(
        prog='Kraken 2 Krona',
        usage='kraken2krona --names <FILE> --nodes <FILE> --output_file <FILE> [--include_unclassified] kraken_report <FILE(s)>',
        description='''
            Takes one or more Kraken 2 style reports and creates an output file compatible with Krona. 
            If multiple report files are specified, the read classifications across the reports will be summed.''')
    parser.add_argument(
        'kraken_report',
        metavar='KRAKEN_REPORT(S)',
        type=str,
        nargs='+',
        help='The Kraken 2 style report. Can be many files separated with spaces. Can be gzipped.')
    parser.add_argument(
        '--output_file',
        metavar='FILE',
        type=str,
        required=True,
        help='The output file.')
    parser.add_argument(
        '--names',
        metavar='FILE',
        required=True,
        help='Taxonomy names dump file (e.g. names.dmp)')
    parser.add_argument(
        '--nodes',
        metavar='FILE',
        required=True,
        help='Taxonomy nodes dump file (e.g. nodes.dmp)')
    parser.add_argument(
        '--include_unclassified',
        action='store_true',
        help='Set this flag to include the unclassified reads in the output.')

    args = parser.parse_args()

    return args

def kraken2krona():

    # Get the CL arguments
    args = get_arguments()
    
    # Create taxonomy tree (tt)
    tt = TaxonomyTree(names_filename=args.names, nodes_filename=args.nodes)

    # Main loop. Parse the report files and create TaxonInfo objects that keep track of read counts and taxonomy.
    # Storing TaxonInfo objects in taxon_info_dict with tax ID as keys.
    taxon_info_dict = {}
    for report_filename in args.kraken_report:
        with read_file(report_filename) as report_file_object:
            taxon_info_dict = parse_report(report_file_object, taxon_info_dict, tt)
    
    # Create a pandas DataFrame from the TaxonInfo objects
    taxon_info_list = [taxon_info for taxon_info in taxon_info_dict.values()]
    taxon_info_df = pd.DataFrame(taxon_info_list).set_index('t_id')
    
    # The columns of the dataframe that we want to output to file
    taxonomy_cols = [taxonomy_rank for taxonomy_rank in taxonomy_ranks.values()]
    cols = ['hits_at_taxon'] + taxonomy_cols

    # Drop unclassified row if user specified so
    if not args.include_unclassified:
        if 0 in taxon_info_df.index.values:
            log.info('Removing row for unclassified reads from output file.')
            taxon_info_df.drop(0, inplace=True)
    
    # Set names for the rows without any values in the taxonomic levels
    taxonomically_empty_rows = taxon_info_df[taxon_info_df[taxonomy_cols].isnull().all(axis=1)]
    for tax_id in taxonomically_empty_rows.index.values:
        if tax_id==0:
            tax_id_name = 'Unclassified'
        else:
            tax_id_name = tt.get_name([tax_id])[tax_id]
        taxon_info_df.at[tax_id, 't_superkingdom'] = tax_id_name

    # Write to file
    log.info('Saving output to {file}.'.format(file=args.output_file))
    taxon_info_df[cols].to_csv(args.output_file, index=False, header=False, sep='\t')

if __name__ == '__main__':
    kraken2krona()
