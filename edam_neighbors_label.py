from __future__ import annotations

import os

from SPARQLWrapper import JSON, SPARQLWrapper


# Fuseki endpoint used by the pipeline.
# Can be overridden with:
#   export FUSEKI_ENDPOINT=http://localhost:3030/biotoolsEdam/query
endpointURL = os.getenv(
    "FUSEKI_ENDPOINT",
    "http://localhost:3030/biotoolsEdam/query",
)

EDAM_URI = "http://edamontology.org/"


PREFIXES = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX edam: <http://edamontology.org/>
"""


def get_edam_category(uri: str) -> str:
    """Return topic/operation/data/format from an EDAM URI."""
    local_name = str(uri).rstrip("/").rsplit("/", 1)[-1]
    return local_name.split("_", 1)[0]


def get_edam_neighbors(uri_list, infer_labels: bool = True):
    """
    Retrieve EDAM neighbors for a list of EDAM class URIs.

    Parameters
    ----------
    infer_labels:
        If True, retrieve each neighbor rdfs:label from the ontology. If False,
        do not query labels and return None for neighbor labels.

    Returns
    -------
    dict
        Mapping:
            parent EDAM URI -> {neighbor EDAM URI -> preferred label}

        The parent URI is mapped to all neighbors found for that class.
    """
    # Do not contact Fuseki when neighbor-label retrieval is disabled.
    # The caller can use this flag to disable the complete neighbor step.
    if not infer_labels:
        return {}

    annotation_neighbors = {}

    for uri in uri_list:
        uri = str(uri).strip()
        if not uri:
            continue

        if not uri.startswith(EDAM_URI):
            raise ValueError(f"Invalid EDAM URI: {uri}")

        current_annotation = uri.replace(EDAM_URI, "edam:", 1)

        neighbor_label_select = " ?neighborLabel" if infer_labels else ""
        neighbor_label_clause = (
            "  OPTIONAL {{\n    ?neighbor rdfs:label ?neighborLabel .\n  }}"
            if infer_labels
            else ""
        )

        query = f"""
SELECT DISTINCT ?conceptLabel ?neighbor{neighbor_label_select}
WHERE {{
  VALUES ?conceptRoot {{
    edam:topic_0003
    edam:operation_0004
    edam:data_0006
    edam:format_1915
  }}

  VALUES ?neighborRoot {{
    edam:topic_0003
    edam:operation_0004
    edam:data_0006
    edam:format_1915
  }}

  VALUES ?concept {{ {current_annotation} }}

  ?concept rdfs:subClassOf* ?conceptRoot .
  ?concept rdf:type owl:Class .

  FILTER NOT EXISTS {{
    ?concept rdfs:subClassOf owl:DeprecatedClass
  }}

  OPTIONAL {{
    ?concept rdfs:label ?conceptLabel .
  }}

  {{
    ?concept (rdfs:subClassOf|owl:someValuesFrom)* [
      rdf:type owl:Restriction ;
      owl:onProperty ?relation ;
      owl:someValuesFrom ?neighborDescendant
    ] .
  }}
  UNION
  {{
    ?concept rdfs:subClassOf ?neighborDescendant .
    ?concept ?relation ?neighborDescendant .
  }}

  ?neighborDescendant rdfs:subClassOf* ?neighbor .

  ?neighbor rdfs:subClassOf* ?neighborRoot .
  ?neighbor rdf:type owl:Class .

  FILTER NOT EXISTS {{
    ?neighbor rdfs:subClassOf owl:DeprecatedClass
  }}

{neighbor_label_clause}
}}
"""

        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(PREFIXES + query)
        sparql.setReturnFormat(JSON)

        query_results = sparql.queryAndConvert()

        neighbors_for_uri = {}

        for result in query_results["results"]["bindings"]:
            neighbor = result["neighbor"]["value"]
            neighbor_label = (
                result["neighborLabel"]["value"]
                if "neighborLabel" in result
                else None
            )
            neighbors_for_uri[neighbor] = neighbor_label

        annotation_neighbors[uri] = neighbors_for_uri

    return annotation_neighbors


if __name__ == "__main__":
    # Simple standalone smoke test.
    test_uri = os.getenv(
        "EDAM_TEST_URI",
        "http://edamontology.org/operation_2403",
    )

    print(f"Fuseki endpoint: {endpointURL}")
    print(f"Querying neighbors of: {test_uri}")

    neighbors = get_edam_neighbors([test_uri])

    for parent_uri, parent_neighbors in sorted(neighbors.items()):
        print(f"\n{parent_uri}")
        for neighbor_uri, label in sorted(parent_neighbors.items()):
            print(f"  {neighbor_uri}\t{label or '-'}")
