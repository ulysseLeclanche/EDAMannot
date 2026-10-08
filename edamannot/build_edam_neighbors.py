from pathlib import Path

from rdflib import Graph
from SPARQLWrapper import SPARQLWrapper, TURTLE


DATA_DIR = Path(__file__).resolve().parent.parent / "data"

OUTPUT_FILE = DATA_DIR / "edam_neighbors.ttl"
ENDPOINT = "http://localhost:3031/edam/query"


PREFIXES = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX edam: <http://edamontology.org/>
"""

QUERY = """
CONSTRUCT {
    ?concept ?relation ?neighbor .
    ?concept rdf:type owl:Class .
    ?neighbor rdf:type owl:Class .
}
WHERE {
    VALUES ?conceptRoot {
        edam:topic_0003
        edam:operation_0004
        edam:data_0006
        edam:format_1915
    }

    VALUES ?neighborRoot {
        edam:topic_0003
        edam:operation_0004
        edam:data_0006
        edam:format_1915
    }

    ?concept rdfs:subClassOf* ?conceptRoot .
    ?concept rdf:type owl:Class .

    FILTER NOT EXISTS {
        ?concept rdfs:subClassOf owl:DeprecatedClass
    }

    {
        ?concept
            (rdfs:subClassOf|owl:someValuesFrom)*
            [
                rdf:type owl:Restriction ;
                owl:onProperty ?relation ;
                owl:someValuesFrom ?neighborDescendant
            ] .

        ?neighborDescendant rdfs:subClassOf* ?neighbor .
    }

    UNION

    {
        ?concept rdfs:subClassOf+ ?neighbor .
        BIND(rdfs:subClassOf AS ?relation)
    }

    ?neighbor rdfs:subClassOf* ?neighborRoot .
    ?neighbor rdf:type owl:Class .

    FILTER isIRI(?neighbor)

    FILTER NOT EXISTS {
        ?neighbor rdfs:subClassOf owl:DeprecatedClass
    }
}
"""


def build_edam_neighbors():
    sparql = SPARQLWrapper(ENDPOINT)
    sparql.setQuery(PREFIXES + QUERY)
    sparql.setReturnFormat(TURTLE)

    result = sparql.query().convert()

    graph = Graph()
    graph.parse(data=result, format="turtle")
    graph.serialize(destination=OUTPUT_FILE, format="turtle")

    print(f"Generated: {OUTPUT_FILE}")
    print(f"Number of triples: {len(graph)}")


if __name__ == "__main__":
    build_edam_neighbors()
