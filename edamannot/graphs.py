"""
Graph and visualization utilities for EDAMannot.

SPARQL access is delegated to edamannot.sparql, while tool annotation helpers
are delegated to edamannot.annotations. Metric DataFrames are loaded lazily
when graph colouring is requested, avoiding import-time file I/O.
"""

import pandas as pd
import pygraphviz as pgv

from .config import biotoolsURI, edamURI
from .sparql import execute_query, sparql_results_to_dataframe
from .annotations import (
    getToolOperations,
    getToolTopics,
    getToolsCommonOperations,
    getToolsCommonTopics,
)


def _create_graph():
    """Create the Graphviz graph used throughout EDAMannot."""
    return pgv.AGraph(directed=True, rankdir="BT")


def _tool_id(tool_uri: str) -> str:
    return tool_uri.replace(biotoolsURI, "")


def _edam_id(uri: str) -> str:
    return uri.replace(edamURI, "").replace("edam:", "")


def getHierarchyGraph(
    entityURI,
    graph=None,
    direction="ancestors",
    displayIdentifier=False,
    highlightEntity=False,
):
    """Return a graph representing the hierarchy of (in)direct superclasses or subclasses for an entity.

    Keyword arguments:
    entityURI -- the URI for the entity
    graph -- the graph in which the hierarchy is added. A new graph is created if the value is None. (default: None)
    direction -- should the hierarchy concern the ancestors and/or the descendants of the entity. Possible values: "ancestors", "descendants", "both" (default: "ancestors")
    displayIdentifier -- should the nodes also display their URI (default: False)
    highlightEntity -- should the entity be highlighted (default:False)
    """
    if graph is None:
        graph = _create_graph()

    entityIdent = _edam_id(entityURI)

    entityType = "Class"
    if entityIdent.startswith("topic_"):
        entityType = "Topic"
    elif entityIdent.startswith("operation_"):
        entityType = "Operation"

    conceptStyle = {}
    conceptStyle["Class"] = "filled"
    conceptStyle["Tool"] = "filled"
    conceptStyle["Topic"] = "filled"
    conceptStyle["Operation"] = "rounded,filled"

    if entityURI.startswith("http"):
        entityURI = "<" + entityURI + ">"

    if (direction == "ancestors") or (direction == "both"):
        query = (
            """
SELECT DISTINCT ?subConcept ?subConceptLabel ?superConcept ?superConceptLabel
WHERE {
  VALUES ?concept { """
            + entityURI
            + """ }
 
  ?concept rdfs:subClassOf* ?subConcept .
  ?subConcept rdf:type owl:Class .
  FILTER NOT EXISTS { ?subConcept rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?subConcept rdfs:label ?subLabel }
  ?subConcept rdfs:subClassOf ?superConcept .
  ?superConcept rdf:type owl:Class .
  FILTER NOT EXISTS { ?superConcept rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?superConcept rdfs:label ?supLabel }
  BIND(COALESCE(?subLabel, "") AS ?subConceptLabel)
  BIND(COALESCE(?supLabel, "") AS ?superConceptLabel)
}
"""
        )
        results = execute_query(query)
        for result in results["results"]["bindings"]:
            startClassIdent = result["subConcept"]["value"].replace(edamURI, "")
            startClassLabel = result["subConceptLabel"]["value"] + (
                "\n(" + startClassIdent + ")" if displayIdentifier else ""
            )
            endClassIdent = result["superConcept"]["value"].replace(edamURI, "")
            endClassLabel = result["superConceptLabel"]["value"] + (
                "\n(" + endClassIdent + ")" if displayIdentifier else ""
            )

            graph.add_node(
                startClassIdent,
                label=startClassLabel,
                shape="box",
                color="black",
                nodeType=entityType,
                style=conceptStyle[entityType],
                fillcolor="#ffffff",
            )
            graph.add_node(
                endClassIdent,
                label=endClassLabel,
                shape="box",
                color="black",
                nodeType=entityType,
                style=conceptStyle[entityType],
                fillcolor="#ffffff",
            )
            graph.add_edge(startClassIdent, endClassIdent, arrowhead="onormal")

    if (direction == "descendants") or (direction == "both"):
        query = (
            """
SELECT DISTINCT ?subConcept ?subConceptLabel ?superConcept ?superConceptLabel
WHERE {
  VALUES ?concept { """
            + entityURI
            + """ }
  
  ?superConcept rdfs:subClassOf* ?concept .
  ?superConcept rdf:type owl:Class .
  FILTER NOT EXISTS { ?superConcept rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?superConcept rdfs:label ?supLabel }
  ?subConcept rdfs:subClassOf ?superConcept .
  ?subConcept rdf:type owl:Class .
  FILTER NOT EXISTS { ?subConcept rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?subConcept rdfs:label ?subLabel }

  BIND(COALESCE(?subLabel, "") AS ?subConceptLabel)
  BIND(COALESCE(?supLabel, "") AS ?superConceptLabel)
}
"""
        )
        results = execute_query(query)
        for result in results["results"]["bindings"]:
            startClassIdent = result["subConcept"]["value"].replace(edamURI, "")
            startClassLabel = result["subConceptLabel"]["value"] + (
                "\n(" + startClassIdent + ")" if displayIdentifier else ""
            )
            endClassIdent = result["superConcept"]["value"].replace(edamURI, "")
            endClassLabel = result["superConceptLabel"]["value"] + (
                "\n(" + endClassIdent + ")" if displayIdentifier else ""
            )

            graph.add_node(
                startClassIdent,
                label=startClassLabel,
                shape="box",
                color="black",
                nodeType=entityType,
                style=conceptStyle[entityType],
                fillcolor="#ffffff",
            )
            graph.add_node(
                endClassIdent,
                label=endClassLabel,
                shape="box",
                color="black",
                nodeType=entityType,
                style=conceptStyle[entityType],
                fillcolor="#ffffff",
            )
            graph.add_edge(startClassIdent, endClassIdent, arrowhead="onormal")

    if highlightEntity:
        if graph.has_node(entityIdent):
            graph.get_node(entityIdent).attr["color"] = "red"
        else:
            graph.add_node(entityIdent, color="red")
    return graph


def get_edam_neighbors_dataframe(endpointURL, prefixes):
    """Return the EDAM neighbor relationships dataframe."""
    query = """
    # Find EDAM concept with neighbors that inherit other neighbors from one of its ancestors
    SELECT DISTINCT ?concept ?conceptLabel ?neighborRelation ?conceptNeighbor ?conceptNeighborLabel ?conceptAncestor ?conceptAncestorLabel ?ancestorNeighborRelation ?conceptAncestorNeighbor ?conceptAncestorNeighborLabel 
    WHERE {
      ?concept rdf:type owl:Class .
      OPTIONAL { ?concept rdfs:label ?conceptLabel . }
      FILTER NOT EXISTS { ?concept rdfs:subClassOf? owl:DeprecatedClass }
      
      ?concept rdfs:subClassOf [
        rdf:type owl:Restriction ;
        owl:onProperty ?neighborRelation ;
        owl:someValuesFrom ?conceptNeighbor
      ] .
      ?conceptNeighbor rdf:type owl:Class .
      FILTER NOT EXISTS { ?conceptNeighbor rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?conceptNeighbor rdfs:label ?conceptNeighborLabel . }
      
      ?concept rdfs:subClassOf+ ?conceptAncestor .
      OPTIONAL { ?conceptAncestor rdfs:label ?conceptAncestorLabel . }
      FILTER NOT EXISTS { ?conceptAncestor rdfs:subClassOf? owl:DeprecatedClass }
      
      ?conceptAncestor rdfs:subClassOf [
        rdf:type owl:Restriction ;
        owl:onProperty ?ancestorNeighborRelation ;
        owl:someValuesFrom ?conceptAncestorNeighbor
      ] .
      ?conceptAncestorNeighbor rdf:type owl:Class .
      FILTER NOT EXISTS { ?conceptAncestorNeighbor rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?conceptAncestorNeighbor rdfs:label ?conceptAncestorNeighborLabel . }
    }
    """
    results = execute_query(query, endpoint_url=endpointURL, query_prefixes=prefixes)
    return sparql_results_to_dataframe(results)


def get_edam_chained_neighbors_dataframe(endpointURL, prefixes):
    """Return the EDAM concepts dataframe containing chained neighbors."""
    query = """
    # Chained neighbors: Find EDAM concept with neighbors that have neighbors of their own
    SELECT DISTINCT ?concept ?conceptLabel ?neighborRelation ?conceptNeighbor ?conceptNeighborLabel ?neighborNeighborRelation ?neighborNeighbor ?neighborNeighborLabel 
    WHERE {
      ?concept rdf:type owl:Class .
      OPTIONAL { ?concept rdfs:label ?conceptLabel . }
      FILTER NOT EXISTS { ?concept rdfs:subClassOf? owl:DeprecatedClass }
      
      ?concept rdfs:subClassOf [
        rdf:type owl:Restriction ;
        owl:onProperty ?neighborRelation ;
        owl:someValuesFrom ?conceptNeighbor
      ] .
      ?conceptNeighbor rdf:type owl:Class .
      FILTER NOT EXISTS { ?conceptNeighbor rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?conceptNeighbor rdfs:label ?conceptNeighborLabel . }
      
      ?conceptNeighbor rdfs:subClassOf [
        rdf:type owl:Restriction ;
        owl:onProperty ?neighborNeighborRelation ;
        owl:someValuesFrom ?neighborNeighbor
      ] .
      ?neighborNeighbor rdf:type owl:Class .
      FILTER NOT EXISTS { ?neighborNeighbor rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?neighborNeighbor rdfs:label ?neighborNeighborLabel . }
    }
    """
    results = execute_query(query, endpoint_url=endpointURL, query_prefixes=prefixes)
    return sparql_results_to_dataframe(results)


def getEntityDescriptionGraph(
    entityURI,
    graph=None,
    direction="ancestors",
    displayIdentifier=False,
    highlightEntity=False,
):
    """Return a graph representing the neighbors of an entity."""
    entityIdent = _edam_id(entityURI)
    graph = getHierarchyGraph(
        entityURI,
        graph=graph,
        direction=direction,
        displayIdentifier=displayIdentifier,
        highlightEntity=highlightEntity,
    )
    if entityURI.startswith("http"):
        entityURI = "<" + entityURI + ">"
    query = (
        """
SELECT DISTINCT ?concept ?relation ?neighbor ?neighborLabel
WHERE {
  VALUES ?concept { """
        + entityURI
        + """ }
 
  ?concept rdfs:subClassOf* ?conceptAncestor .
  ?conceptAncestor rdf:type owl:Class .
  FILTER NOT EXISTS { ?conceptAncestor rdfs:subClassOf? owl:DeprecatedClass }
  
  ?conceptAncestor rdfs:subClassOf ?restriction .
  ?restriction rdf:type owl:Restriction .
  ?restriction owl:onProperty ?relation .
  ?restriction owl:someValuesFrom ?neighbor .
  OPTIONAL { ?neighbor rdfs:label ?neighborConceptLabel }
  BIND(COALESCE(?neighborConceptLabel, "") AS ?neighborLabel)
}
"""
    )

    conceptStyle = {
        "Class": "filled",
        "Tool": "filled",
        "Topic": "filled",
        "Operation": "rounded,filled",
        "Data": "filled",
        "Format": "rounded,filled",
    }
    conceptShape = {
        "Class": "box",
        "Tool": "oval",
        "Topic": "box",
        "Operation": "box",
        "Data": "hexagon",
        "Format": "parallelogram",
    }
    conceptFillColor = {
        "Class": "white",
        "Tool": "white",
        "Topic": "#ccebc5",
        "Operation": "#b3cde3",
        "Data": "#fbb4ae",
        "Format": "#fed9a6",
    }

    results = execute_query(query)
    for result in results["results"]["bindings"]:
        relationIdent = result["relation"]["value"].replace(edamURI, "")
        neighborIdent = result["neighbor"]["value"].replace(edamURI, "")
        neighborLabel = result["neighborLabel"]["value"] + (
            "\n(" + neighborIdent + ")" if displayIdentifier else ""
        )
        neighborType = "Class"
        if neighborIdent.startswith("topic_"):
            neighborType = "Topic"
        elif neighborIdent.startswith("operation_"):
            neighborType = "Operation"
        elif neighborIdent.startswith("data_"):
            neighborType = "Data"
        elif neighborIdent.startswith("format_"):
            neighborType = "Format"
        graph.add_node(
            neighborIdent,
            label=neighborLabel,
            shape=conceptShape[neighborType],
            color="black",
            nodeType=neighborType,
            style=conceptStyle[neighborType],
            fillcolor=conceptFillColor[neighborType],
        )
        graph.add_edge(
            entityIdent,
            neighborIdent,
            arrowhead="open",
            label=relationIdent,
            color="purple",
            fontcolor="purple",
        )
    return graph


def addToolAndAnnotationsToGraph(
    toolURI,
    graph=None,
    showTopics=True,
    showOperations=True,
    showDeprecatedAnnotations=False,
    highlightDirectAnnotations=False,
):
    """Return a graph representing a tool and its EDAM annotations."""
    if graph is None:
        graph = _create_graph()

    toolIdent = _tool_id(toolURI)

    if toolURI.startswith("http"):
        toolURI = "<" + toolURI + ">"

    directAnnotationColor = "red" if highlightDirectAnnotations else "black"

    conceptStyle = {
        "Tool": "filled",
        "Topic": "filled",
        "TopicDeprecated": "filled,dashed",
        "TopicAlternative": "filled,dotted",
        "Operation": "rounded,filled",
        "OperationDeprecated": "rounded,filled,dashed",
        "OperationAlternative": "rounded,filled,dotted",
    }

    # Tool node
    query = (
        """
SELECT DISTINCT ?toolLabel 
WHERE {
  VALUES ?tool { """
        + toolURI
        + """ }

  ?tool rdf:type sc:SoftwareApplication .
  OPTIONAL { ?tool sc:name ?tLabel }
  BIND(COALESCE(?tLabel, "") AS ?toolLabel)
}
"""
    )
    results = execute_query(query)
    for result in results["results"]["bindings"]:
        if not graph.has_node(toolIdent):
            clusterTools = graph.get_subgraph(name="cluster_tools")
            if clusterTools is None:
                clusterTools = graph.add_subgraph(
                    name="cluster_tools", rankdir="same", style="invis"
                )
            clusterTools.add_node(
                toolIdent,
                label="{}".format(result["toolLabel"]["value"]),
                shape="ellipse",
                color="blue",
                nodeType="Tool",
                style=conceptStyle["Tool"],
                fillcolor="#ffffff",
            )

    # Direct topics and their hierarchy
    if showTopics:
        directConcepts = []
        conceptType = "Topic"
        query = (
            """
SELECT DISTINCT ?conceptURI ?conceptLabel
WHERE {
  VALUES ?tool { """
            + toolURI
            + """ }

  ?tool sc:applicationSubCategory ?conceptURI .
  ?conceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?conceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?conceptURI rdfs:label ?cLabel }
  BIND(COALESCE(?cLabel, "") AS ?conceptLabel)
}
"""
        )
        results = execute_query(query)
        for result in results["results"]["bindings"]:
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            if not graph.has_node(conceptIdent):
                graph.add_node(
                    conceptIdent,
                    label="{}\n({})".format(
                        result["conceptLabel"]["value"], conceptIdent
                    ),
                    nodeType=conceptType,
                    shape="box",
                    color=directAnnotationColor,
                    style=conceptStyle[conceptType],
                    fillcolor="#ffffff",
                )
            graph.add_edge(
                toolIdent,
                conceptIdent,
                arrowhead="vee",
                color="blue",
                fontcolor="blue",
                style="dashed",
            )

        for conceptURI in directConcepts:
            hierarchy_query = (
                """
SELECT DISTINCT ?subConceptURI ?subConceptLabel ?superConceptURI ?superConceptLabel 
WHERE {
  VALUES ?conceptURI { <"""
                + conceptURI
                + """> }
  ?conceptURI rdfs:subClassOf* ?subConceptURI .
  ?subConceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?subConceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?subConceptURI rdfs:label ?subLabel }
  ?subConceptURI rdfs:subClassOf ?superConceptURI .
  ?superConceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?superConceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?superConceptURI rdfs:label ?supLabel }
  BIND(COALESCE(?subLabel, "") AS ?subConceptLabel)
  BIND(COALESCE(?supLabel, "") AS ?superConceptLabel)
}
"""
            )
            results = execute_query(hierarchy_query)
            for result in results["results"]["bindings"]:
                subConceptIdent = result["subConceptURI"]["value"].replace(edamURI, "")
                superConceptIdent = result["superConceptURI"]["value"].replace(edamURI, "")
                if not graph.has_node(subConceptIdent):
                    graph.add_node(
                        subConceptIdent,
                        label="{}\n{}".format(
                            result["subConceptLabel"]["value"], subConceptIdent
                        ),
                        shape="box",
                        color="black",
                        nodeType=conceptType,
                        style=conceptStyle[conceptType],
                        fillcolor="#ffffff",
                    )
                if not graph.has_node(superConceptIdent):
                    graph.add_node(
                        superConceptIdent,
                        label="{}\n{}".format(
                            result["superConceptLabel"]["value"], superConceptIdent
                        ),
                        shape="box",
                        color="black",
                        nodeType=conceptType,
                        style=conceptStyle[conceptType],
                        fillcolor="#ffffff",
                    )
                graph.add_edge(subConceptIdent, superConceptIdent, arrowhead="onormal")

    # Direct operations and their hierarchy
    if showOperations:
        directConcepts = []
        conceptType = "Operation"
        query = (
            """
SELECT DISTINCT ?conceptURI ?conceptLabel
WHERE {
  VALUES ?tool { """
            + toolURI
            + """ }

  ?tool sc:featureList ?conceptURI .
  ?conceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?conceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?conceptURI rdfs:label ?cLabel }
  BIND(COALESCE(?cLabel, "") AS ?conceptLabel)
}
"""
        )
        results = execute_query(query)
        for result in results["results"]["bindings"]:
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            if not graph.has_node(conceptIdent):
                graph.add_node(
                    conceptIdent,
                    label="{}\n({})".format(
                        result["conceptLabel"]["value"], conceptIdent
                    ),
                    nodeType=conceptType,
                    shape="box",
                    color=directAnnotationColor,
                    style=conceptStyle[conceptType],
                    fillcolor="#ffffff",
                )
            graph.add_edge(
                toolIdent,
                conceptIdent,
                arrowhead="vee",
                color="blue",
                fontcolor="blue",
                style="dashed",
            )

        for conceptURI in directConcepts:
            hierarchy_query = (
                """
SELECT DISTINCT ?subConceptURI ?subConceptLabel ?superConceptURI ?superConceptLabel 
WHERE {
  VALUES ?conceptURI { <"""
                + conceptURI
                + """> }
  ?conceptURI rdfs:subClassOf* ?subConceptURI .
  ?subConceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?subConceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?subConceptURI rdfs:label ?subLabel }
  ?subConceptURI rdfs:subClassOf ?superConceptURI .
  ?superConceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?superConceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?superConceptURI rdfs:label ?supLabel }
  BIND(COALESCE(?subLabel, "") AS ?subConceptLabel)
  BIND(COALESCE(?supLabel, "") AS ?superConceptLabel)
}
"""
            )
            results = execute_query(hierarchy_query)
            for result in results["results"]["bindings"]:
                subConceptIdent = result["subConceptURI"]["value"].replace(edamURI, "")
                superConceptIdent = result["superConceptURI"]["value"].replace(edamURI, "")
                if not graph.has_node(subConceptIdent):
                    graph.add_node(
                        subConceptIdent,
                        label="{}\n{}".format(
                            result["subConceptLabel"]["value"], subConceptIdent
                        ),
                        shape="box",
                        color="black",
                        nodeType=conceptType,
                        style=conceptStyle[conceptType],
                        fillcolor="#ffffff",
                    )
                if not graph.has_node(superConceptIdent):
                    graph.add_node(
                        superConceptIdent,
                        label="{}\n{}".format(
                            result["superConceptLabel"]["value"], superConceptIdent
                        ),
                        shape="box",
                        color="black",
                        nodeType=conceptType,
                        style=conceptStyle[conceptType],
                        fillcolor="#ffffff",
                    )
                graph.add_edge(subConceptIdent, superConceptIdent, arrowhead="onormal")

    if showDeprecatedAnnotations:
        # Deprecated topics
        query = (
            """
SELECT DISTINCT ?conceptURI ?conceptLabel ?conceptAlternative ?conceptAlternativeLabel
WHERE {
  VALUES ?tool { """
            + toolURI
            + """ }

  ?tool sc:applicationSubCategory ?conceptURI .
  { ?conceptURI rdfs:subClassOf? owl:DeprecatedClass }
  UNION
  { ?conceptURI owl:deprecated true }
  UNION
  { ?conceptURI owl:deprecated "true" }
  UNION
  { ?conceptURI owl:deprecated "True" }
  OPTIONAL { ?conceptURI rdfs:label ?cLabel }
  BIND(COALESCE(?cLabel, "") AS ?conceptLabel)
  
  OPTIONAL {
    ?conceptURI oboInOwl:consider ?conceptAlternative .
    OPTIONAL { ?conceptAlternative rdfs:label ?caLabel }
    BIND(COALESCE(?caLabel, "") AS ?conceptAlternativeLabel)
  }
}
"""
        )
        results = execute_query(query)
        for result in results["results"]["bindings"]:
            conceptType = "TopicDeprecated"
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            if not graph.has_node(conceptIdent):
                graph.add_node(
                    conceptIdent,
                    label="{}\n({})".format(
                        result["conceptLabel"]["value"], conceptIdent
                    ),
                    nodeType=conceptType,
                    shape="box",
                    color="grey",
                    style=conceptStyle[conceptType],
                    fillcolor="#ffffff",
                )
            graph.add_edge(
                toolIdent,
                conceptIdent,
                arrowhead="vee",
                color="grey",
                fontcolor="grey",
                style="dashed",
            )
            if "conceptAlternative" in result:
                conceptType = "TopicAlternative"
                alternativeURI = result["conceptAlternative"]["value"]
                alternativeIdent = alternativeURI.replace(edamURI, "")
                if not graph.has_node(alternativeIdent):
                    graph.add_node(
                        alternativeIdent,
                        label="{}\n({})".format(
                            result["conceptAlternativeLabel"]["value"], alternativeIdent
                        ),
                        nodeType=conceptType,
                        shape="box",
                        color="grey",
                        style=conceptStyle[conceptType],
                        fillcolor="#ffffff",
                    )
                    queryHierarchy = (
                        """
SELECT DISTINCT ?subConceptURI ?subConceptLabel ?superConceptURI ?superConceptLabel 
WHERE {
  VALUES ?conceptURI { <"""
                        + alternativeURI
                        + """> }
  ?conceptURI rdfs:subClassOf* ?subConceptURI .
  ?subConceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?subConceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?subConceptURI rdfs:label ?subLabel }
  ?subConceptURI rdfs:subClassOf ?superConceptURI .
  ?superConceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?superConceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?superConceptURI rdfs:label ?supLabel }
  BIND(COALESCE(?subLabel, "") AS ?subConceptLabel)
  BIND(COALESCE(?supLabel, "") AS ?superConceptLabel)
}
"""
                    )
                    resultsHierarchy = execute_query(queryHierarchy)
                    for resultHierarchy in resultsHierarchy["results"]["bindings"]:
                        subConceptIdent = resultHierarchy["subConceptURI"]["value"].replace(edamURI, "")
                        superConceptIdent = resultHierarchy["superConceptURI"]["value"].replace(edamURI, "")
                        if not graph.has_node(subConceptIdent):
                            graph.add_node(
                                subConceptIdent,
                                label="{}\n{}".format(
                                    resultHierarchy["subConceptLabel"]["value"],
                                    subConceptIdent,
                                ),
                                shape="box",
                                color="grey",
                                nodeType=conceptType,
                                style=conceptStyle[conceptType],
                                fillcolor="#ffffff",
                            )
                        if not graph.has_node(superConceptIdent):
                            graph.add_node(
                                superConceptIdent,
                                label="{}\n{}".format(
                                    resultHierarchy["superConceptLabel"]["value"],
                                    superConceptIdent,
                                ),
                                shape="box",
                                color="grey",
                                nodeType=conceptType,
                                style=conceptStyle[conceptType],
                                fillcolor="#ffffff",
                            )
                        if not graph.has_edge(subConceptIdent, superConceptIdent):
                            graph.add_edge(
                                subConceptIdent,
                                superConceptIdent,
                                arrowhead="onormal",
                                color="grey",
                                style="dotted",
                            )

                graph.add_edge(
                    conceptIdent,
                    alternativeIdent,
                    arrowhead="vee",
                    color="grey",
                    fontcolor="grey",
                    style="dotted",
                )

        # Deprecated operations
        query = (
            """
SELECT DISTINCT ?conceptURI ?conceptLabel ?conceptAlternative ?conceptAlternativeLabel
WHERE {
  VALUES ?tool { """
            + toolURI
            + """ }

  ?tool sc:featureList ?conceptURI .
  { ?conceptURI rdfs:subClassOf? owl:DeprecatedClass }
  UNION
  { ?conceptURI owl:deprecated true }
  UNION
  { ?conceptURI owl:deprecated "true" }
  UNION
  { ?conceptURI owl:deprecated "True" }
  OPTIONAL { ?conceptURI rdfs:label ?cLabel }
  BIND(COALESCE(?cLabel, "") AS ?conceptLabel)
  
  OPTIONAL {
    ?conceptURI oboInOwl:consider ?conceptAlternative .
    OPTIONAL { ?conceptAlternative rdfs:label ?caLabel }
    BIND(COALESCE(?caLabel, "") AS ?conceptAlternativeLabel)
  }
}
"""
        )
        results = execute_query(query)
        for result in results["results"]["bindings"]:
            conceptType = "OperationDeprecated"
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            if not graph.has_node(conceptIdent):
                graph.add_node(
                    conceptIdent,
                    label="{}\n({})".format(
                        result["conceptLabel"]["value"], conceptIdent
                    ),
                    nodeType=conceptType,
                    shape="box",
                    color="grey",
                    style=conceptStyle[conceptType],
                    fillcolor="#ffffff",
                )
            graph.add_edge(
                toolIdent,
                conceptIdent,
                arrowhead="vee",
                color="grey",
                fontcolor="grey",
                style="dashed",
            )
            if "conceptAlternative" in result:
                conceptType = "OperationAlternative"
                alternativeURI = result["conceptAlternative"]["value"]
                alternativeIdent = alternativeURI.replace(edamURI, "")
                if not graph.has_node(alternativeIdent):
                    graph.add_node(
                        alternativeIdent,
                        label="{}\n({})".format(
                            result["conceptAlternativeLabel"]["value"], alternativeIdent
                        ),
                        nodeType=conceptType,
                        shape="box",
                        color="grey",
                        style=conceptStyle[conceptType],
                        fillcolor="#ffffff",
                    )
                    queryHierarchy = (
                        """
SELECT DISTINCT ?subConceptURI ?subConceptLabel ?superConceptURI ?superConceptLabel 
WHERE {
  VALUES ?conceptURI { <"""
                        + alternativeURI
                        + """> }
  ?conceptURI rdfs:subClassOf* ?subConceptURI .
  ?subConceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?subConceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?subConceptURI rdfs:label ?subLabel }
  ?subConceptURI rdfs:subClassOf ?superConceptURI .
  ?superConceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?superConceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?superConceptURI rdfs:label ?supLabel }
  BIND(COALESCE(?subLabel, "") AS ?subConceptLabel)
  BIND(COALESCE(?supLabel, "") AS ?superConceptLabel)
}
"""
                    )
                    resultsHierarchy = execute_query(queryHierarchy)
                    for resultHierarchy in resultsHierarchy["results"]["bindings"]:
                        subConceptIdent = resultHierarchy["subConceptURI"]["value"].replace(edamURI, "")
                        superConceptIdent = resultHierarchy["superConceptURI"]["value"].replace(edamURI, "")
                        if not graph.has_node(subConceptIdent):
                            graph.add_node(
                                subConceptIdent,
                                label="{}\n{}".format(
                                    resultHierarchy["subConceptLabel"]["value"],
                                    subConceptIdent,
                                ),
                                shape="box",
                                color="grey",
                                nodeType=conceptType,
                                style=conceptStyle[conceptType],
                                fillcolor="#ffffff",
                            )
                        if not graph.has_node(superConceptIdent):
                            graph.add_node(
                                superConceptIdent,
                                label="{}\n{}".format(
                                    resultHierarchy["superConceptLabel"]["value"],
                                    superConceptIdent,
                                ),
                                shape="box",
                                color="grey",
                                nodeType=conceptType,
                                style=conceptStyle[conceptType],
                                fillcolor="#ffffff",
                            )
                        if not graph.has_edge(subConceptIdent, superConceptIdent):
                            graph.add_edge(
                                subConceptIdent,
                                superConceptIdent,
                                arrowhead="onormal",
                                color="grey",
                                style="dotted",
                            )
                graph.add_edge(
                    conceptIdent,
                    alternativeIdent,
                    arrowhead="vee",
                    color="grey",
                    fontcolor="grey",
                    style="dotted",
                )

    return graph


def addToolsAndAnnotationsToGraph(
    listToolURI,
    graph=None,
    showTopics=True,
    showOperations=True,
    showDeprecatedAnnotations=False,
    highlightDirectAnnotations=False,
    highlightIntersection=False,
):
    """Return a graph representing tools and their EDAM annotations."""
    if graph is None:
        graph = _create_graph()

    for toolURI in listToolURI:
        addToolAndAnnotationsToGraph(
            toolURI,
            graph=graph,
            showTopics=showTopics,
            showOperations=showOperations,
            showDeprecatedAnnotations=showDeprecatedAnnotations,
            highlightDirectAnnotations=highlightDirectAnnotations,
        )

    if highlightIntersection:
        commonConcepts = []
        if showTopics:
            commonConcepts += getToolsCommonTopics(listToolURI, transitive=True)
        if showOperations:
            commonConcepts += getToolsCommonOperations(listToolURI, transitive=True)

        for currentConcept, currentLabel in commonConcepts:
            currentConceptIdent = currentConcept.replace(edamURI, "")
            node = graph.get_node(currentConceptIdent)
            node.attr["color"] = "red"
            node.attr["penwidth"] = "3"
            node.attr["style"] = "filled,bold"
            node.attr["fillcolor"] = "red"
            node.attr["highlightIntersection"] = "True"

    return graph


def getScoreColorRGB(scoreValue, scoreMaxValue, color="red"):
    """Return the RGB color (in hex) associated to a score."""
    scoreMaxValue = max(1, scoreMaxValue)
    fraction = scoreValue / scoreMaxValue

    if color == "red":
        target = (255, 0, 0)
    elif color == "green":
        target = (0, 255, 0)
    elif color == "blue":
        target = (0, 0, 255)
    elif color == "orange":
        target = (255, 165, 0)
    elif color == "yellow":
        target = (255, 255, 0)
    elif color == "pink":
        target = (255, 192, 203)
    elif color == "grey":
        target = (128, 128, 128)
    else:
        return "#ffffff"

    r = int(255 - fraction * (255 - target[0]))
    g = int(255 - fraction * (255 - target[1]))
    b = int(255 - fraction * (255 - target[2]))

    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def dictTopic(dfTopicmetrics, metric_col, edamURI=edamURI):
    """Extract a dictionary for a given topic metric column."""
    raw_dict = (
        dfTopicmetrics[["topic", metric_col]].set_index("topic").to_dict()[metric_col]
    )
    return {topic.replace(edamURI, ""): val for topic, val in raw_dict.items()}


def dictOperation(dfOperationmetrics, metric_col, edamURI=edamURI):
    """Extract a dictionary for a given operation metric column."""
    raw_dict = (
        dfOperationmetrics[["operation", metric_col]]
        .set_index("operation")
        .to_dict()[metric_col]
    )
    return {op.replace(edamURI, ""): val for op, val in raw_dict.items()}


def _load_graph_metric_dataframes():
    """Load the metric DataFrames required for graph colouring."""
    dfTopicmetrics = pd.read_csv(
        "Dataframe/dfTopicmetrics.tsv.bz2", sep="\t", compression="bz2"
    )
    dfOperationmetrics = pd.read_csv(
        "Dataframe/dfOperationmetrics.tsv.bz2", sep="\t", compression="bz2"
    )
    return dfTopicmetrics, dfOperationmetrics


def buildTopicOperationDicts(metric_col):
    """Build the topic and operation score dictionaries for a metric column."""
    dfTopicmetrics, dfOperationmetrics = _load_graph_metric_dataframes()
    dictTopicScore = dictTopic(dfTopicmetrics, metric_col)
    dictOperationScore = dictOperation(dfOperationmetrics, metric_col)
    return dictTopicScore, dictOperationScore


def colorGraphNodesAccordingToScore(
    graph, dictTopicScore, dictOperationScore, color="red"
):
    """Modify node fill colours according to associated topic/operation scores."""
    topicMaxValue = max(dictTopicScore.values(), default=0)
    operationMaxValue = max(dictOperationScore.values(), default=0)
    for currentNode in graph.nodes_iter():
        if "nodeType" in currentNode.attr.keys():
            currentNodeType = currentNode.attr["nodeType"]
            if currentNodeType == "Topic":
                if currentNode in dictTopicScore.keys():
                    graph.get_node(currentNode).attr["fillcolor"] = getScoreColorRGB(
                        dictTopicScore[currentNode], topicMaxValue, color=color
                    )
            elif currentNodeType == "Operation":
                if currentNode in dictOperationScore.keys():
                    graph.get_node(currentNode).attr["fillcolor"] = getScoreColorRGB(
                        dictOperationScore[currentNode], operationMaxValue, color=color
                    )


def getToolScore(toolURI, transitive=False, dictTopicScore={}, dictOperationScore={}):
    """Return the sum of the scores of the tool's topic and operation annotations."""
    toolScore = 0
    for currentTopic in getToolTopics(toolURI, transitive=transitive):
        currentTopicIdent = currentTopic[0].replace(edamURI, "")
        if currentTopicIdent in dictTopicScore.keys():
            toolScore += dictTopicScore[currentTopicIdent]
    for currentOperation in getToolOperations(toolURI, transitive=transitive):
        currentOperationIdent = currentOperation[0].replace(edamURI, "")
        if currentOperationIdent in dictOperationScore.keys():
            toolScore += dictOperationScore[currentOperationIdent]
    return toolScore


__all__ = [
    "getHierarchyGraph",
    "get_edam_neighbors_dataframe",
    "get_edam_chained_neighbors_dataframe",
    "getEntityDescriptionGraph",
    "addToolAndAnnotationsToGraph",
    "addToolsAndAnnotationsToGraph",
    "getScoreColorRGB",
    "dictTopic",
    "dictOperation",
    "buildTopicOperationDicts",
    "colorGraphNodesAccordingToScore",
    "getToolScore",
]
