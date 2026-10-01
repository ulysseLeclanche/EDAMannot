"""
Configuration centrale for EDAMannot.

This module contains paths, SPARQL configuration, RDF prefixes, and URI
constants currently defined in EDAMannot.py.

The historical variable names are intentionally preserved for compatibility
with the existing code during the refactoring:
    endpointURL
    rdfFormat
    prefixes
    biotoolsURI
    biotoolsOntologyURI
    edamURI
    bioschemas_file
    edam_file
    current_dir
    neighbor_dir
    notebooks_dir
    edam_neighbors_file

Only configuration belongs here. No pandas DataFrame is loaded at import time.
"""

import os
from typing import Optional


# ---------------------------------------------------------------------------
# PROJECT DIRECTORIES
# ---------------------------------------------------------------------------

current_dir = os.getcwd()

neighbor_dir = os.path.join(current_dir, "edam")
notebooks_dir = os.path.join(current_dir, "Notebooks")


# ---------------------------------------------------------------------------
# INPUT / RESOURCE FILES
# ---------------------------------------------------------------------------

bioschemas_file = os.path.join(
    neighbor_dir,
    "bioschemas-dump_05_01_2025.ttl",
)

edam_file = os.path.join(
    neighbor_dir,
    "EDAM_1.25.owl",
)

edam_neighbors_file = os.path.join(
    notebooks_dir,
    "EDAM_neighbors_result.ttl",
)


# ---------------------------------------------------------------------------
# SPARQL CONFIGURATION
# ---------------------------------------------------------------------------

endpointURL = "http://localhost:3030/sharefair/query"
rdfFormat = "turtle"


# ---------------------------------------------------------------------------
# SPARQL PREFIXES
# ---------------------------------------------------------------------------

prefixes = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX foaf: <http://xmlns.com/foaf/0.1/>
PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

PREFIX bt: <https://bio.tools/>
PREFIX biotools: <https://bio.tools/ontology/>
PREFIX bsc: <http://bioschemas.org/>
PREFIX bsct: <http://bioschemas.org/types/>
PREFIX edam: <http://edamontology.org/>
PREFIX sc: <http://schema.org/>
PREFIX schema: <https://schema.org/>
"""


# ---------------------------------------------------------------------------
# COMMON URIs
# ---------------------------------------------------------------------------

biotoolsURI = "https://bio.tools/"
biotoolsOntologyURI = "https://bio.tools/ontology/"
edamURI = "http://edamontology.org/"


# ---------------------------------------------------------------------------
# COMPATIBILITY / CONFIGURATION HELPERS
# ---------------------------------------------------------------------------

def set_file_paths(
    new_bioschemas_file: Optional[str] = None,
    new_edam_file: Optional[str] = None,
) -> None:
    """
    Update the configured Bioschemas and EDAM ontology file paths.

    This preserves the behavior of the original EDAMannot.set_file_paths()
    function: paths are changed only when a new value is provided, and the
    resulting paths are converted to absolute paths.

    Parameters
    ----------
    new_bioschemas_file:
        Optional path to the Bioschemas Turtle dump.
    new_edam_file:
        Optional path to the EDAM ontology file.
    """
    global bioschemas_file, edam_file

    if new_bioschemas_file:
        bioschemas_file = os.path.abspath(new_bioschemas_file)

    if new_edam_file:
        edam_file = os.path.abspath(new_edam_file)


__all__ = [
    "current_dir",
    "neighbor_dir",
    "notebooks_dir",
    "bioschemas_file",
    "edam_file",
    "edam_neighbors_file",
    "endpointURL",
    "rdfFormat",
    "prefixes",
    "biotoolsURI",
    "biotoolsOntologyURI",
    "edamURI",
    "set_file_paths",
]
