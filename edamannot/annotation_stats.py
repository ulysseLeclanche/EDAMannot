"""
Compare EDAM annotation counts before/after enrichment using Fuseki.

No rdflib is used.

The script temporarily uploads:
    --before -> urn:edamannot:stats_before
    --after  -> urn:edamannot:stats_after

Then it counts, for each EDAM category:
    - annotation assignments
    - unique EDAM terms
    - annotated tools (Topic / Operation)
    - annotated parameters (Data / Format)

Data and Format are counted on the FormalParameter level:
    Data   -> sc:additionalType
    Format -> sc:encodingFormat

Only EDAM URI objects are counted. Literal values are ignored.

Example:
    python annotation_stats.py \
        --before bioschemas-dump.ttl \
        --after bioschemas_enriched.ttl

Optional:
    --keep-graphs
"""

from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request

from fuseki_client import DEFAULT_FUSEKI_BASE, FusekiClient


BEFORE_GRAPH = "urn:edamannot:stats_before"
AFTER_GRAPH = "urn:edamannot:stats_after"

PREFIXES = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX bsc: <https://bioschemas.org/>
PREFIX bsct: <https://bioschemas.org/types/>
PREFIX edam: <http://edamontology.org/>
PREFIX sc: <http://schema.org/>
"""


CATEGORY_QUERIES = {
    "Topic": """
        ?tool a sc:SoftwareApplication ;
              sc:applicationSubCategory ?annotation .
        FILTER(STRSTARTS(STR(?annotation), "http://edamontology.org/topic_"))
    """,
    "Operation": """
        ?tool a sc:SoftwareApplication ;
              sc:featureList ?annotation .
        FILTER(STRSTARTS(STR(?annotation), "http://edamontology.org/operation_"))
    """,
    "Data": """
        ?tool a sc:SoftwareApplication ;
              ?parameter_property ?parameter .

        VALUES ?parameter_property {
            bsc:input
            bsc:output
        }

        ?parameter a bsct:FormalParameter ;
                   sc:additionalType ?annotation .

        FILTER(STRSTARTS(STR(?annotation), "http://edamontology.org/data_"))
    """,
    "Format": """
        ?tool a sc:SoftwareApplication ;
              ?parameter_property ?parameter .

        VALUES ?parameter_property {
            bsc:input
            bsc:output
        }

        ?parameter a bsct:FormalParameter ;
                   sc:encodingFormat ?annotation .

        FILTER(STRSTARTS(STR(?annotation), "http://edamontology.org/format_"))
    """,
}


def query_json(client: FusekiClient, query: str) -> dict:
    """Execute a SPARQL SELECT query and return its JSON result."""
    params = urllib.parse.urlencode({"query": query})
    url = f"{client.query_url}?{params}"

    request = urllib.request.Request(
        url,
        headers={"Accept": "application/sparql-results+json"},
        method="GET",
    )

    with urllib.request.urlopen(request, timeout=client.timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def count_category(
    client: FusekiClient,
    graph: str,
    category: str,
) -> dict[str, int]:
    """Count assignments, unique terms, and annotated entities."""
    pattern = CATEGORY_QUERIES[category]

    entity_var = "?tool" if category in {"Topic", "Operation"} else "?parameter"

    query = PREFIXES + f"""
    SELECT
        (COUNT(*) AS ?assignments)
        (COUNT(DISTINCT ?annotation) AS ?unique_terms)
        (COUNT(DISTINCT {entity_var}) AS ?entities)
    WHERE {{
        GRAPH <{graph}> {{
            {pattern}
        }}
    }}
    """

    result = query_json(client, query)
    binding = result["results"]["bindings"][0]

    return {
        "assignments": int(binding["assignments"]["value"]),
        "unique_terms": int(binding["unique_terms"]["value"]),
        "entities": int(binding["entities"]["value"]),
    }


def print_comparison(
    before: dict[str, dict[str, int]],
    after: dict[str, dict[str, int]],
) -> None:
    """Print a compact before/after comparison."""
    headers = [
        "Category",
        "Assignments before",
        "Assignments after",
        "Added",
        "Unique terms before",
        "Unique terms after",
        "Entities before",
        "Entities after",
    ]

    rows = []
    for category in ["Topic", "Operation", "Data", "Format"]:
        b = before[category]
        a = after[category]
        rows.append([
            category,
            b["assignments"],
            a["assignments"],
            a["assignments"] - b["assignments"],
            b["unique_terms"],
            a["unique_terms"],
            b["entities"],
            a["entities"],
        ])

    widths = [
        max(len(str(row[i])) for row in [headers] + rows)
        for i in range(len(headers))
    ]

    print("\nAnnotation statistics")
    print("=" * sum(widths + [3] * (len(widths) - 1)))

    print(" | ".join(
        str(value).ljust(widths[i])
        for i, value in enumerate(headers)
    ))

    print("-+-".join("-" * width for width in widths))

    for row in rows:
        print(" | ".join(
            str(value).ljust(widths[i])
            for i, value in enumerate(row)
        ))

    total_before = sum(before[c]["assignments"] for c in before)
    total_after = sum(after[c]["assignments"] for c in after)

    print("\nTotal assignments:")
    print(f"  Before: {total_before}")
    print(f"  After : {total_after}")
    print(f"  Added : {total_after - total_before}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare EDAM annotations before and after enrichment."
    )
    parser.add_argument(
        "--before",
        required=True,
        help="Original Bioschemas Turtle file.",
    )
    parser.add_argument(
        "--after",
        required=True,
        help="Enriched Bioschemas Turtle file.",
    )
    parser.add_argument(
        "--fuseki",
        default=DEFAULT_FUSEKI_BASE,
        help=f"Fuseki dataset URL (default: {DEFAULT_FUSEKI_BASE}).",
    )
    parser.add_argument(
        "--keep-graphs",
        action="store_true",
        help="Keep temporary statistics graphs after the run.",
    )

    args = parser.parse_args()

    client = FusekiClient(args.fuseki)

    try:
        print(f"Fuseki dataset: {client.base_url}")
        print(f"Loading original file: {args.before}")
        client.clear_graph(BEFORE_GRAPH)
        client.upload_ttl(args.before, BEFORE_GRAPH)

        print(f"Loading enriched file: {args.after}")
        client.clear_graph(AFTER_GRAPH)
        client.upload_ttl(args.after, AFTER_GRAPH)

        before = {}
        after = {}

        for category in ["Topic", "Operation", "Data", "Format"]:
            print(f"Counting {category}...")
            before[category] = count_category(
                client,
                BEFORE_GRAPH,
                category,
            )
            after[category] = count_category(
                client,
                AFTER_GRAPH,
                category,
            )

        print_comparison(before, after)

    finally:
        if not args.keep_graphs:
            print("Removing temporary statistics graphs...")
            try:
                client.clear_graph(BEFORE_GRAPH)
            finally:
                client.clear_graph(AFTER_GRAPH)
            print("Temporary statistics graphs removed.")


if __name__ == "__main__":
    main()
