from pathlib import Path

from rdflib import Graph, Namespace, URIRef
from SPARQLWrapper import SPARQLWrapper, JSON


DATA_DIR = Path(__file__).resolve().parent.parent / "data"

INPUT_FILE = DATA_DIR / "bioschemas-dump.ttl"
OUTPUT_FILE = DATA_DIR / "bioschemas-dump_enriched.ttl"

ENDPOINT = "http://localhost:3031/edam/query"

EDAM = Namespace("http://edamontology.org/")
SC = Namespace("http://schema.org/")

EDAM_PREDICATES = {
    "topic": SC.applicationSubCategory,
    "operation": SC.featureList,
    "data": SC.additionalType,
    "format": SC.encodingFormat,
}


def edam_type(uri):
    """
    Extract the EDAM type from the URI.
    """
    local_name = str(uri).rsplit("/", 1)[-1]
    return local_name.split("_", 1)[0]


def enriched_annotation(input_file, output_file):
    graph = Graph()
    graph.parse(input_file, format="turtle")  # Parse bioschemas-dump.ttl

    annotation_tools = {}

    for predicate in EDAM_PREDICATES.values():
        for tool, _, annotation in graph.triples((None, predicate, None)):
            if not isinstance(annotation, URIRef):
                continue
            if not str(annotation).startswith(str(EDAM)):
                continue
            annotation_tools.setdefault(str(annotation), set()).add(tool)

    # annotation_tools dict contains the EDAM-tools annotation association

    if annotation_tools:
        values = " ".join(f"<{uri}>" for uri in annotation_tools)

        # Query edam_neighbors.ttl to determine the relationships between annotations
        query = f"""
        SELECT DISTINCT ?concept ?neighbor
        WHERE {{
            VALUES ?concept {{
                {values}
            }}

            ?concept ?relation ?neighbor .
            FILTER(?relation != rdf:type)
            FILTER(isIRI(?neighbor))
        }}
        """

        sparql = SPARQLWrapper(ENDPOINT)
        sparql.setQuery(
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n" + query
        )
        sparql.setReturnFormat(JSON)

        results = sparql.query().convert()

        for result in results["results"]["bindings"]:
            concept = result["concept"]["value"]
            neighbor = result["neighbor"]["value"]

            predicate = EDAM_PREDICATES.get(edam_type(neighbor))
            if predicate is None:
                continue

            # Adds the "neighbor" annotation to tools
            for tool in annotation_tools.get(concept, []):
                graph.add((tool, predicate, URIRef(neighbor)))

    graph.serialize(destination=output_file, format="turtle")

    print(f"Generated file: {output_file}")


if __name__ == "__main__":
    enriched_annotation(INPUT_FILE, OUTPUT_FILE)
