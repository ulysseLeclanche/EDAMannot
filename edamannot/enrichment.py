"""
EDAM annotation enrichment for EDAMannot.

Architecture:
    - Existing Fuseki dataset: http://localhost:3030/sharefair
    - Existing EDAM / ShareFAIR data: default graph
    - Input Bioschemas TTL: urn:edamannot:input
    - Enriched Bioschemas TTL: urn:edamannot:enriched
    - Temporary EDAM neighbor index: urn:edamannot:neighbors

Python only orchestrates Fuseki HTTP operations.
The complete input TTL is NOT parsed with rdflib.

Enrichment is done in two stages:

    1. Build a temporary neighbor index in Fuseki:
           EDAM annotation -> EDAM neighbor

       This is calculated only for EDAM annotations present in the input
       graph, rather than for every EDAM class. In --tool mode it is limited
       to that tool.

    2. Join the input annotations against that index and insert inferred
       annotations into urn:edamannot:enriched.

Supported annotation categories:

    Topic:
        sc:applicationSubCategory

    Operation:
        sc:featureList

    Data:
        sc:additionalType on bsct:FormalParameter under bsc:input/bsc:output

    Format:
        sc:encodingFormat on bsct:FormalParameter under bsc:input/bsc:output

The neighbor logic follows the existing EDAMannot
edam_neighbors_label(7).py implementation:
    - rdfs:subClassOf*
    - owl:Restriction
    - owl:onProperty
    - owl:someValuesFrom
    - owl:DeprecatedClass exclusion

Neighbors are constrained to the same EDAM branch:
    Topic     -> Topic
    Operation -> Operation
    Data      -> Data
    Format    -> Format
"""

from __future__ import annotations

import argparse
from pathlib import Path

try:
    # Normal package import when called through CLI.py at repository root.
    from .fuseki_client import (
        DEFAULT_FUSEKI_BASE,
        ENRICHED_GRAPH,
        INPUT_GRAPH,
        FusekiClient,
        cleanup_enrichment_graphs,
        prepare_enrichment_graphs,
    )
except ImportError:
    # Direct execution from edamannot/, e.g.:
    # python enrichment.py --input ... --output ...
    from fuseki_client import (
        DEFAULT_FUSEKI_BASE,
        ENRICHED_GRAPH,
        INPUT_GRAPH,
        FusekiClient,
        cleanup_enrichment_graphs,
        prepare_enrichment_graphs,
    )


EDAM_PREFIXES = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX bsc: <https://bioschemas.org/>
PREFIX bsct: <https://bioschemas.org/types/>
PREFIX edam: <http://edamontology.org/>
PREFIX sc: <http://schema.org/>
PREFIX edamannot: <urn:edamannot:>
"""

NEIGHBORS_GRAPH = "urn:edamannot:neighbors"
NEIGHBOR_PREDICATE = "edamannot:neighbor"


EDAM_ROOTS = {
    "topic": "http://edamontology.org/topic_0003",
    "operation": "http://edamontology.org/operation_0004",
    "data": "http://edamontology.org/data_0006",
    "format": "http://edamontology.org/format_1915",
}


def _tool_filter(tool_uri: str | None) -> str:
    """Restrict enrichment to one tool when --tool is supplied."""
    if not tool_uri:
        return ""

    if "<" in tool_uri or ">" in tool_uri or any(
        char.isspace() for char in tool_uri
    ):
        raise ValueError(f"Invalid tool URI: {tool_uri}")

    return f"""
        VALUES ?tool {{ <{tool_uri}> }}
"""


def _neighbor_pattern(root: str) -> str:
    """
    Return the EDAM neighbor pattern from the existing EDAMannot logic,
    restricted to one EDAM branch.
    """
    return f"""
        ?concept rdfs:subClassOf* <{root}> .
        ?concept rdf:type owl:Class .

        FILTER NOT EXISTS {{
            ?concept rdfs:subClassOf owl:DeprecatedClass
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

        ?neighbor rdfs:subClassOf* <{root}> .
        ?neighbor rdf:type owl:Class .

        FILTER NOT EXISTS {{
            ?neighbor rdfs:subClassOf owl:DeprecatedClass
        }}

        FILTER (?neighbor != ?concept)
    """


def build_neighbor_index(
    client: FusekiClient,
    root: str,
    tool_uri: str | None = None,
    category: str | None = None,
) -> None:
    """
    Materialize only the EDAM annotation -> neighbor pairs that are needed.

    Without --tool:
        concepts are restricted to EDAM annotations actually present in
        the input graph, which avoids computing neighbors for the whole
        EDAM ontology.

    With --tool:
        concepts are restricted to annotations attached to that tool.
        This makes single-tool validation much faster.

    The EDAM ontology is read from the existing /sharefair default graph.
    Results are stored in a temporary named graph.
    """
    if category not in EDAM_ROOTS:
        raise ValueError(f"Unknown EDAM category: {category}")

    root_filter = f"""
        VALUES ?tool {{
            <{tool_uri}>
        }}
""" if tool_uri else ""

    if category == "topic":
        annotation_pattern = f"""
        GRAPH <{INPUT_GRAPH}> {{
            ?tool a sc:SoftwareApplication ;
                  sc:applicationSubCategory ?concept .
            {root_filter}
        }}
"""
    elif category == "operation":
        annotation_pattern = f"""
        GRAPH <{INPUT_GRAPH}> {{
            ?tool a sc:SoftwareApplication ;
                  sc:featureList ?concept .
            {root_filter}
        }}
"""
    elif category == "data":
        annotation_pattern = f"""
        GRAPH <{INPUT_GRAPH}> {{
            ?tool a sc:SoftwareApplication ;
                  ?parameterProperty ?parameter .

            VALUES ?parameterProperty {{
                bsc:input
                bsc:output
            }}

            ?parameter a bsct:FormalParameter ;
                       sc:additionalType ?concept .
            {root_filter}
        }}
"""
    else:  # format
        annotation_pattern = f"""
        GRAPH <{INPUT_GRAPH}> {{
            ?tool a sc:SoftwareApplication ;
                  ?parameterProperty ?parameter .

            VALUES ?parameterProperty {{
                bsc:input
                bsc:output
            }}

            ?parameter a bsct:FormalParameter ;
                       sc:encodingFormat ?concept .
            {root_filter}
        }}
"""

    update = EDAM_PREFIXES + f"""
    INSERT {{
        GRAPH <{NEIGHBORS_GRAPH}> {{
            ?concept {NEIGHBOR_PREDICATE} ?neighbor
        }}
    }}
    WHERE {{
        {annotation_pattern}

        {_neighbor_pattern(root)}
    }}
    """

    client.run_update(update)


def clear_neighbor_index(client: FusekiClient) -> None:
    """Remove the temporary annotation -> neighbor index."""
    client.clear_graph(NEIGHBORS_GRAPH)


def enrich_topics(
    client: FusekiClient,
    tool_uri: str | None = None,
) -> None:
    """Add inferred Topic neighbors to sc:applicationSubCategory."""
    update = EDAM_PREFIXES + f"""
    INSERT {{
        GRAPH <{ENRICHED_GRAPH}> {{
            ?tool sc:applicationSubCategory ?neighbor
        }}
    }}
    WHERE {{
        GRAPH <{INPUT_GRAPH}> {{
            ?tool a sc:SoftwareApplication ;
                  sc:applicationSubCategory ?concept .
        }}

        GRAPH <{NEIGHBORS_GRAPH}> {{
            ?concept {NEIGHBOR_PREDICATE} ?neighbor .
        }}

        {_tool_filter(tool_uri)}

        FILTER NOT EXISTS {{
            GRAPH <{ENRICHED_GRAPH}> {{
                ?tool sc:applicationSubCategory ?neighbor
            }}
        }}
    }}
    """

    client.run_update(update)


def enrich_operations(
    client: FusekiClient,
    tool_uri: str | None = None,
) -> None:
    """Add inferred Operation neighbors to sc:featureList."""
    update = EDAM_PREFIXES + f"""
    INSERT {{
        GRAPH <{ENRICHED_GRAPH}> {{
            ?tool sc:featureList ?neighbor
        }}
    }}
    WHERE {{
        GRAPH <{INPUT_GRAPH}> {{
            ?tool a sc:SoftwareApplication ;
                  sc:featureList ?concept .
        }}

        GRAPH <{NEIGHBORS_GRAPH}> {{
            ?concept {NEIGHBOR_PREDICATE} ?neighbor .
        }}

        {_tool_filter(tool_uri)}

        FILTER NOT EXISTS {{
            GRAPH <{ENRICHED_GRAPH}> {{
                ?tool sc:featureList ?neighbor
            }}
        }}
    }}
    """

    client.run_update(update)


def enrich_data(
    client: FusekiClient,
    tool_uri: str | None = None,
) -> None:
    """
    Add inferred Data neighbors to the same FormalParameter as the original
    data annotation.
    """
    update = EDAM_PREFIXES + f"""
    INSERT {{
        GRAPH <{ENRICHED_GRAPH}> {{
            ?parameter sc:additionalType ?neighbor
        }}
    }}
    WHERE {{
        GRAPH <{INPUT_GRAPH}> {{
            ?tool a sc:SoftwareApplication ;
                  ?parameterProperty ?parameter .

            VALUES ?parameterProperty {{
                bsc:input
                bsc:output
            }}

            ?parameter a bsct:FormalParameter ;
                       sc:additionalType ?concept .
        }}

        GRAPH <{NEIGHBORS_GRAPH}> {{
            ?concept {NEIGHBOR_PREDICATE} ?neighbor .
        }}

        {_tool_filter(tool_uri)}

        FILTER NOT EXISTS {{
            GRAPH <{ENRICHED_GRAPH}> {{
                ?parameter sc:additionalType ?neighbor
            }}
        }}
    }}
    """

    client.run_update(update)


def enrich_formats(
    client: FusekiClient,
    tool_uri: str | None = None,
) -> None:
    """
    Add inferred Format neighbors to the same FormalParameter as the original
    format annotation.
    """
    update = EDAM_PREFIXES + f"""
    INSERT {{
        GRAPH <{ENRICHED_GRAPH}> {{
            ?parameter sc:encodingFormat ?neighbor
        }}
    }}
    WHERE {{
        GRAPH <{INPUT_GRAPH}> {{
            ?tool a sc:SoftwareApplication ;
                  ?parameterProperty ?parameter .

            VALUES ?parameterProperty {{
                bsc:input
                bsc:output
            }}

            ?parameter a bsct:FormalParameter ;
                       sc:encodingFormat ?concept .
        }}

        GRAPH <{NEIGHBORS_GRAPH}> {{
            ?concept {NEIGHBOR_PREDICATE} ?neighbor .
        }}

        {_tool_filter(tool_uri)}

        FILTER NOT EXISTS {{
            GRAPH <{ENRICHED_GRAPH}> {{
                ?parameter sc:encodingFormat ?neighbor
            }}
        }}
    }}
    """

    client.run_update(update)


def enrich_all(
    client: FusekiClient,
    tool_uri: str | None = None,
) -> None:
    """
    Build the EDAM neighbor index once, then enrich all four categories.
    """
    print("Building EDAM Topic neighbor index...")
    build_neighbor_index(
        client,
        EDAM_ROOTS["topic"],
        tool_uri=tool_uri,
        category="topic",
    )
    print("Topic neighbor index ready.")

    print("Building EDAM Operation neighbor index...")
    build_neighbor_index(
        client,
        EDAM_ROOTS["operation"],
        tool_uri=tool_uri,
        category="operation",
    )
    print("Operation neighbor index ready.")

    print("Building EDAM Data neighbor index...")
    build_neighbor_index(
        client,
        EDAM_ROOTS["data"],
        tool_uri=tool_uri,
        category="data",
    )
    print("Data neighbor index ready.")

    print("Building EDAM Format neighbor index...")
    build_neighbor_index(
        client,
        EDAM_ROOTS["format"],
        tool_uri=tool_uri,
        category="format",
    )
    print("Format neighbor index ready.")

    print("Enriching Topics...")
    enrich_topics(client, tool_uri=tool_uri)
    print("Topics enriched.")

    print("Enriching Operations...")
    enrich_operations(client, tool_uri=tool_uri)
    print("Operations enriched.")

    print("Enriching Data...")
    enrich_data(client, tool_uri=tool_uri)
    print("Data enriched.")

    print("Enriching Formats...")
    enrich_formats(client, tool_uri=tool_uri)
    print("Formats enriched.")


def enriched_annotation(
    input_path: str | Path,
    output_path: str | Path,
    fuseki_url: str = DEFAULT_FUSEKI_BASE,
    keep_graphs: bool = False,
    tool_uri: str | None = None,
) -> None:
    """
    Run the complete enriched annotation workflow.
    """
    client = FusekiClient(fuseki_url)

    try:
        print(f"Fuseki dataset: {client.base_url}")
        print(f"Input graph: {INPUT_GRAPH}")
        print(f"Enriched graph: {ENRICHED_GRAPH}")
        print(f"Neighbor graph: {NEIGHBORS_GRAPH}")

        prepare_enrichment_graphs(
            client,
            input_path,
        )

        clear_neighbor_index(client)

        print("Input TTL: loaded.")
        print("Enriched graph: initialized.")

        enrich_all(
            client,
            tool_uri=tool_uri,
        )

        client.export_ttl(
            output_path,
            ENRICHED_GRAPH,
        )

        print(f"Enriched TTL exported to: {output_path}")

    finally:
        if not keep_graphs:
            try:
                clear_neighbor_index(client)
            finally:
                cleanup_enrichment_graphs(client)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Enrich a Bioschemas Turtle file with EDAM neighbor annotations "
            "using the existing Fuseki dataset."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input Bioschemas Turtle file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output enriched Turtle file.",
    )

    parser.add_argument(
        "--fuseki",
        default=DEFAULT_FUSEKI_BASE,
        help=f"Fuseki dataset URL (default: {DEFAULT_FUSEKI_BASE}).",
    )

    parser.add_argument(
        "--tool",
        dest="tool_uri",
        help=(
            "Optional Bio.tools URI to restrict enrichment to one tool for "
            "testing, e.g. https://bio.tools/2dkd."
        ),
    )

    parser.add_argument(
        "--keep-graphs",
        action="store_true",
        help="Keep temporary named graphs in Fuseki after the run.",
    )

    args = parser.parse_args()

    enriched_annotation(
        input_path=args.input,
        output_path=args.output,
        fuseki_url=args.fuseki,
        keep_graphs=args.keep_graphs,
        tool_uri=args.tool_uri,
    )


if __name__ == "__main__":
    main()
