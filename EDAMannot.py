import os
import pyarrow
import pandas as pd
from rdflib import Dataset
import collections
from graphviz import Digraph
import IPython
import json
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pygraphviz as pgv
import rdflib
import rdflib.namespace
import scipy.stats as stats
import sparqldataframe
from SPARQLWrapper import SPARQLWrapper, JSON
import time
from typing import Union, List, Dict, Optional, Literal, Tuple


# === Global variables ===

current_dir = os.getcwd()
neighbor_dir = os.path.join(current_dir, "edam")

bioschemas_file = os.path.join(neighbor_dir, "bioschemas-dump_05_01_2025.ttl")
edam_file = os.path.join(neighbor_dir, "EDAM_1.25.owl")


dfTool = pd.read_csv("Dataframe/dfTool.tsv.bz2", sep="\t")  # all Tool and toolLabel
dfToolTopic = pd.read_csv(
    "Dataframe/dfToolTopic.tsv.bz2", sep="\t"
)  # Tool, topic, topicLabel no transitive
dfToolOperation = pd.read_csv(
    "Dataframe/dfToolOperation.tsv.bz2", sep="\t"
)  # tool, operation, operationLabel no transitive


dfToolallmetrics = pd.read_csv("Dataframe/dfToolallmetrics.tsv.bz2", sep="\t")#All tool metrics on transitive Topic and Operation
dfToolTopic = pd.read_csv("Dataframe/dfToolTopic.tsv.bz2", sep="\t")#Tool, topic, topicLabel no transitive
dfToolTopicTransitive = pd.read_csv("Dataframe/dfToolTopicTransitive.tsv.bz2", sep="\t")#Tool, topic, topicLabel transitive
dfToolOperation = pd.read_csv("Dataframe/dfToolOperation.tsv.bz2", sep="\t")#tool, operation, operationLabel no transitive
dfToolOperationTransitive = pd.read_csv("Dataframe/dfToolOperationTransitive.tsv.bz2", sep="\t")# tool, operation, operationLabel transitive
df_redundancy_topic = pd.read_csv("Dataframe/dfToolTopic_redundancy.tsv.bz2", sep="\t")#Identification of redundancy topic
df_redundancy_operation = pd.read_csv("Dataframe/dfToolOperation_redundancy.tsv.bz2", sep="\t")#Identification of redundancy operation
df_topic_no_redundancy = pd.read_csv("Dataframe/df_topic_no_redundancy.tsv.bz2", sep="\t")#tool, topic and topicLabel with no redundancy and no transitive
df_operation_no_redundancy = pd.read_csv("Dataframe/df_operation_no_redundancy.tsv.bz2", sep="\t")#tool, operation, operationLabel with no redundancy and no transitive
dfToolallmetrics_NT = pd.read_csv("Dataframe/dfToolallmetrics_NT.tsv.bz2", sep="\t")



dfTopicmetrics = pd.read_csv("Dataframe/dfTopicmetrics.tsv.bz2", sep="\t")#frequence, IC and entroypy of topics unique metric inherited  
dfOperationmetrics = pd.read_csv("Dataframe/dfOperationmetrics.tsv.bz2", sep="\t")#frequence, IC and entroypy of operations unique metric inherited 
dfOperationmetrics_NT = pd.read_csv("Dataframe/dfOperationmetrics_NT.tsv.bz2", sep="\t")#frequence, IC and entroypy of operations unique metric directly assigned
dfTopicmetrics_NT = pd.read_csv("Dataframe/dfTopicmetrics_NT.tsv.bz2", sep="\t")#frequence, IC and entroypy of topics unique metric directly assigned

nbTools = len(dfTool)
nbToolsWithTopic = dfToolTopic["tool"].nunique()
nbToolsWithOperation = dfToolOperation["tool"].nunique()


# Configuration SPARQL end point
endpointURL = "http://localhost:3030/biotoolsEdam/query"
rdfFormat = "turtle"

# Import prefix :
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

# Link tool
biotoolsURI = "https://bio.tools/"
biotoolsOntologyURI = "https://bio.tools/ontology/"
edamURI = "http://edamontology.org/"


def set_file_paths(new_bioschemas_file: str = None, new_edam_file: str = None):
    """
    Optionally update global file paths for bioschemas and edam files.
    If not provided, defaults are kept.
    """
    global bioschemas_file, edam_file

    if new_bioschemas_file:
        bioschemas_file = os.path.abspath(new_bioschemas_file)
    if new_edam_file:
        edam_file = os.path.abspath(new_edam_file)


def displaySparqlResults(results):
    """
    Displays as HTML the result of a SPARQLWrapper query in a Jupyter notebook.

        Parameters:
            results (dictionnary): the result of a call to SPARQLWrapper.query().convert()
    """
    variableNames = results["head"]["vars"]
    # tableCode = '<table><tr><th>{}</th></tr><tr>{}</tr></table>'.format('</th><th>'.join(variableNames), '</tr><tr>'.join('<td>{}</td>'.format('</td><td>'.join([row[vName]['value'] for vName in variableNames]))for row in results["results"]["bindings"]))
    tableCode = "<table><tr><th>{}</th></tr><tr>{}</tr></table>".format(
        "</th><th>".join(variableNames),
        "</tr><tr>".join(
            "<td>{}</td>".format(
                "</td><td>".join(
                    [
                        row[vName]["value"] if vName in row.keys() else "&nbsp;"
                        for vName in variableNames
                    ]
                )
            )
            for row in results["results"]["bindings"]
        ),
    )
    IPython.display.display(IPython.display.HTML(tableCode))


def sparql_results_to_dataframe(results):
    """
    Convert SPARQL JSON results to pandas DataFrame
    """
    # Extract variable names and bindings
    variable_names = results["head"]["vars"]
    bindings = results["results"]["bindings"]

    # Create a list to store rows
    rows = []

    # Process each row
    for row in bindings:
        row_data = {}
        for var_name in variable_names:
            # if variable exists in row keys, get its value, else use None
            if var_name in row.keys():
                row_data[var_name] = row[var_name]["value"]
            else:
                row_data[var_name] = None
        rows.append(row_data)

    # Create DataFrame
    df = pd.DataFrame(rows, columns=variable_names)

    return df


def get_edam_version(endpointURL, prefixes):
    query = """
    SELECT ?ontology ?versionIRI (REPLACE(STR(?versionIRI), 'http://edamontology.org/', '') AS ?versionNumber)
    WHERE {
      ?ontology rdf:type owl:Ontology .
      ?ontology owl:versionIRI ?versionIRI .
    }
    """
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    # Extract the version number from the first result
    version_str = results["results"]["bindings"][0]["versionNumber"]["value"]
    return float(version_str)


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
        # graph = Digraph(graph_attr={'rankdir': 'BT'})
        graph = pgv.AGraph(directed=True, rankdir="BT")

    entityIdent = entityURI.replace(edamURI, "").replace("edam:", "")

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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            startClassIdent = result["subConcept"]["value"].replace(edamURI, "")
            startClassLabel = result["subConceptLabel"]["value"] + (
                "\n(" + startClassIdent + ")" if displayIdentifier else ""
            )
            endClassIdent = result["superConcept"]["value"].replace(edamURI, "")
            endClassLabel = result["superConceptLabel"]["value"] + (
                "\n(" + endClassIdent + ")" if displayIdentifier else ""
            )

            # graph.node(startClassIdent, label=startClassLabel, shape="box")
            # graph.node(endClassIdent, label=endClassLabel, shape="box")
            # graph.edge(startClassIdent, endClassIdent, arrowhead="onormal")
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            startClassIdent = result["subConcept"]["value"].replace(edamURI, "")
            startClassLabel = result["subConceptLabel"]["value"] + (
                "\n(" + startClassIdent + ")" if displayIdentifier else ""
            )
            endClassIdent = result["superConcept"]["value"].replace(edamURI, "")
            endClassLabel = result["superConceptLabel"]["value"] + (
                "\n(" + endClassIdent + ")" if displayIdentifier else ""
            )

            # graph.node(startClassIdent, label=startClassLabel, shape="box")
            # graph.node(endClassIdent, label=endClassLabel, shape="box")
            # graph.edge(startClassIdent, endClassIdent, arrowhead="onormal")
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
    """
    Execute SPARQL query to find EDAM concepts with neighbors and save as DataFrame/TSV.BZ2

    Parameters:
    - endpoint_url: SPARQL endpoint URL
    - prefixes: SPARQL prefixes string
    - output_file: Optional output filename for TSV.BZ2 file

    Returns:
    - pandas DataFrame with query results
    """

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

    # Execute SPARQL query
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    # Convert to pandas DataFrame using the adapted function
    df = sparql_results_to_dataframe(results)
    return df


def get_edam_chained_neighbors_dataframe(endpointURL, prefixes):
    """
    Execute SPARQL query to find EDAM concepts with chained neighbors and save as DataFrame/TSV.BZ2

    Parameters:
    - endpoint_url: SPARQL endpoint URL
    - prefixes: SPARQL prefixes string
    - output_file: Optional output filename for TSV.BZ2 file

    Returns:
    - pandas DataFrame with chained neighbor relationships
    """

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

    # Execute SPARQL query
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    # Convert to pandas DataFrame
    df = sparql_results_to_dataframe(results)

    return df


def getEntityDescriptionGraph(
    entityURI,
    graph=None,
    direction="ancestors",
    displayIdentifier=False,
    highlightEntity=False,
):
    """Return a graph representing the neighbors of an entity.

    Keyword arguments:
    entityURI -- the URI for the entity
    graph -- the graph in which the hierarchy is added. A new graph is created if the value is None. (default: None)
    direction -- should the hierarchy concern the ancestors and/or the descendants of the entity. Possible values: "ancestors", "descendants", "both" (default: "ancestors")
    displayIdentifier -- should the nodes also display their URI (default: False)
    highlightEntity -- should the entity be highlighted (default:False)
    """
    entityIdent = entityURI.replace(edamURI, "").replace("edam:", "")
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
    conceptStyle = {}
    conceptStyle["Class"] = "filled"
    conceptStyle["Tool"] = "filled"
    conceptStyle["Topic"] = "filled"
    conceptStyle["Operation"] = "rounded,filled"
    conceptStyle["Data"] = "filled"
    conceptStyle["Format"] = "rounded,filled"

    conceptShape = {}
    conceptShape["Class"] = "box"
    conceptShape["Tool"] = "oval"
    conceptShape["Topic"] = "box"
    conceptShape["Operation"] = "box"
    conceptShape["Data"] = "hexagon"
    conceptShape["Format"] = "parallelogram"

    conceptFillColor = (
        {}
    )  # values from ColorBrewer paster theme https://en.wikipedia.org/wiki/ColorBrewer
    conceptFillColor["Class"] = "white"
    conceptFillColor["Tool"] = "white"
    conceptFillColor["Topic"] = "#ccebc5"  # green
    conceptFillColor["Operation"] = "#b3cde3"  # blue
    conceptFillColor["Data"] = "#fbb4ae"  # red
    conceptFillColor["Format"] = "#fed9a6"  # orange

    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
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
    """Return a graph representing a tool and its EDAM annotations.

    Keyword arguments:
    toolURI -- the URI for the tool
    graph -- the graph in which the tool and its annotations are added. A new graph is created if the value is None. (default: None)
    showTopics -- should the topics annotating the tool be considered (default: True)
    showOperations -- should the operations annotating the tool be considered (default: True)
    showDeprecatedAnnotations -- should the deprecated topics and operations annotating the tool be considered (default: False)
    highlightDirectAnnotations -- should the topics or operations annotated directly be highliigthed (default:False)
    """
    if graph is None:
        graph = pgv.AGraph(directed=True, rankdir="BT")

    toolIdent = toolURI.replace(biotoolsURI, "")

    if toolURI.startswith("http"):
        toolURI = "<" + toolURI + ">"

    directAnnotationColor = "red" if highlightDirectAnnotations else "black"

    conceptStyle = {}
    conceptStyle["Tool"] = "filled"
    conceptStyle["Topic"] = "filled"
    conceptStyle["TopicDeprecated"] = conceptStyle["Topic"] + ",dashed"
    conceptStyle["TopicAlternative"] = conceptStyle["Topic"] + ",dotted"
    conceptStyle["Operation"] = "rounded,filled"
    conceptStyle["OperationDeprecated"] = conceptStyle["Operation"] + ",dashed"
    conceptStyle["OperationAlternative"] = conceptStyle["Operation"] + ",dotted"

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
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        # print("{}\t{}".format(toolIdent, result["toolLabel"]["value"]))
        if not graph.has_node(toolIdent):
            clusterTools = graph.get_subgraph(name="cluster_tools")
            if clusterTools is None:
                clusterTools = graph.add_subgraph(
                    name="cluster_tools", rankdir="same", style="invis"
                )  # style="invis"
            clusterTools.add_node(
                toolIdent,
                label="{}".format(result["toolLabel"]["value"]),
                shape="ellipse",
                color="blue",
                nodeType="Tool",
                style=conceptStyle["Tool"],
                fillcolor="#ffffff",
            )

    directConcepts = []
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            # print("\t{}\t{}".format(conceptIdent, result["conceptLabel"]["value"]))
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
        query = (
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        directConcepts = []
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            # print("\t{}\t{}".format(conceptIdent, result["conceptLabel"]["value"]))
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
        query = (
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        directConcepts = []
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
        # conceptType = "TopicDeprecated"
        query = (
            """
SELECT DISTINCT ?conceptURI ?conceptLabel ?conceptAlternative ?conceptAlternativeLabel
WHERE {
  VALUES ?tool { """
            + toolURI
            + """ }

  ?tool sc:applicationSubCategory ?conceptURI .
  #?conceptURI rdf:type owl:Class .
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            conceptType = "TopicDeprecated"
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            # print("\t{}\t{}".format(conceptIdent, result["conceptLabel"]["value"]))
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
            if "conceptAlternative" in result.keys():
                conceptType = "TopicAlternative"
                alternativeIdent = result["conceptAlternative"]["value"].replace(
                    edamURI, ""
                )
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
                    sparqlHierarchy = SPARQLWrapper(endpointURL)
                    sparqlHierarchy.setQuery(prefixes + queryHierarchy)
                    sparqlHierarchy.setReturnFormat(JSON)
                    resultsHierarchy = sparqlHierarchy.query().convert()
                    for resultHierarchy in resultsHierarchy["results"]["bindings"]:
                        subConceptIdent = resultHierarchy["subConceptURI"][
                            "value"
                        ].replace(edamURI, "")
                        superConceptIdent = resultHierarchy["superConceptURI"][
                            "value"
                        ].replace(edamURI, "")
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

        # conceptType = "OperationDeprecated"
        query = (
            """
SELECT DISTINCT ?conceptURI ?conceptLabel ?conceptAlternative ?conceptAlternativeLabel
WHERE {
  VALUES ?tool { """
            + toolURI
            + """ }

  ?tool sc:featureList ?conceptURI .
  #?conceptURI rdf:type owl:Class .
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            conceptType = "OperationDeprecated"
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            # print("\t{}\t{}".format(conceptIdent, result["conceptLabel"]["value"]))
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
            if "conceptAlternative" in result.keys():
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
                    sparqlHierarchy = SPARQLWrapper(endpointURL)
                    sparqlHierarchy.setQuery(prefixes + queryHierarchy)
                    sparqlHierarchy.setReturnFormat(JSON)
                    resultsHierarchy = sparqlHierarchy.query().convert()
                    for resultHierarchy in resultsHierarchy["results"]["bindings"]:
                        subConceptIdent = resultHierarchy["subConceptURI"][
                            "value"
                        ].replace(edamURI, "")
                        superConceptIdent = resultHierarchy["superConceptURI"][
                            "value"
                        ].replace(edamURI, "")
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


def getToolLabel(toolURI):
    """Return the label of a tool (or empty string if no label is present).

    Keyword arguments:
    toolURI -- the URI for the tool
    """

    if toolURI.startswith("http"):
        toolURI = "<" + toolURI + ">"
    query = (
        """
SELECT DISTINCT ?tool ?toolLabel
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
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.queryAndConvert()
    return results["results"]["bindings"][0]["toolLabel"]["value"]


def getToolURIByLabel(toolLabel):
    """Return the URI of a tool designated by its label (or None).

    Keyword arguments:
    toolLabel -- the label for the tool
    """

    query = (
        """
SELECT DISTINCT ?tool ?toolLabel
WHERE {
  VALUES ?toolLabel { \""""
        + toolLabel
        + """\" }

  ?tool rdf:type sc:SoftwareApplication .
  ?tool sc:name ?toolLabel .
}
"""
    )
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.queryAndConvert()
    return (
        None
        if len(results["results"]["bindings"]) == 0
        else results["results"]["bindings"][0]["tool"]["value"]
    )


def getToolTopics(toolURI, transitive=False):
    """Return the list of the (URI, label) tuples for the topics associated to a tool.

    Keyword arguments:
    toolURI -- the URI for the tool
    transitive -- also consider the ancestors of the topics directly associated to the tool (default: False)
    """

    if toolURI.startswith("http"):
        toolURI = "<" + toolURI + ">"
    transitiveClause = "/(rdfs:subClassOf*)" if transitive else ""
    query = (
        """
SELECT DISTINCT ?tool ?topic ?topicLabel
WHERE {
  VALUES ?tool { """
        + toolURI
        + """ }

  ?tool sc:applicationSubCategory"""
        + transitiveClause
        + """ ?topic .
  ?topic rdf:type owl:Class .
  FILTER NOT EXISTS { ?topic rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?topic rdfs:label ?tLabel }
  BIND(COALESCE(?tLabel, "") AS ?topicLabel)
}
"""
    )
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.queryAndConvert()
    toolTopics = [
        (result["topic"]["value"], result["topicLabel"]["value"])
        for result in results["results"]["bindings"]
    ]
    return toolTopics


def getToolOperations(toolURI, transitive=False):
    """Return the list of the (URI, label) tuples for the operations associated to a tool.

    Keyword arguments:
    toolURI -- the URI for the tool
    transitive -- also consider the ancestors of the operations directly associated to the tool (default: False)
    """

    if toolURI.startswith("http"):
        toolURI = "<" + toolURI + ">"
    transitiveClause = "/(rdfs:subClassOf*)" if transitive else ""
    query = (
        """
SELECT DISTINCT ?tool ?operation ?operationLabel
WHERE {
  VALUES ?tool { """
        + toolURI
        + """ }

  ?tool sc:featureList"""
        + transitiveClause
        + """ ?operation .
  ?operation rdf:type owl:Class .
  FILTER NOT EXISTS { ?operation rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?operation rdfs:label ?oLabel }
  BIND(COALESCE(?oLabel, "") AS ?operationLabel)
}
"""
    )
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.queryAndConvert()
    toolOperations = [
        (result["operation"]["value"], result["operationLabel"]["value"])
        for result in results["results"]["bindings"]
    ]
    return toolOperations


def getToolsCommonTopics(listToolURI, transitive=False):
    """Return the list of the (URI, label) tuples for the topics associated to all the tools of a list.

    Keyword arguments:
    listToolURI -- list of the URIs for the tools
    transitive -- also consider the ancestors of the topics directly associated to a tool (default: False)
    """
    commonConcepts = []
    if len(listToolURI) > 0:
        commonConcepts = set(getToolTopics(listToolURI[0], transitive=transitive))
    for toolURI in listToolURI[1:]:
        currentConcepts = set(getToolTopics(toolURI, transitive=transitive))
        commonConcepts = commonConcepts.intersection(currentConcepts)
    return list(commonConcepts)


def getToolsCommonOperations(listToolURI, transitive=False):
    """Return the list of the (URI, label) tuples for the operations associated to all the tools of a list.

    Keyword arguments:
    listToolURI -- list of the URIs for the tools
    transitive -- also consider the ancestors of the operations directly associated to a tool (default: False)
    """
    commonConcepts = []
    if len(listToolURI) > 0:
        commonConcepts = set(getToolOperations(listToolURI[0], transitive=transitive))
    for toolURI in listToolURI[1:]:
        currentConcepts = set(getToolOperations(toolURI, transitive=transitive))
        commonConcepts = commonConcepts.intersection(currentConcepts)
    return list(commonConcepts)


def addToolsAndAnnotationsToGraph(
    listToolURI,
    graph=None,
    showTopics=True,
    showOperations=True,
    highlightDirectAnnotations=False,
    highlightIntersection=False,
):
    """Return a graph representing tools and their EDAM annotations.

    Keyword arguments:
      listToolURI -- list of the URIs for the tools
      graph -- the graph in which the tool and its annotations are added. A new graph is created if the value is None.
      showTopics -- should the topics annotating the tool be considered
      showOperations -- should the operations annotating the tool be considered
      highlightDirectAnnotations -- should the topics or operations annotated directly be highlighted
      highlightIntersection -- should the common topics or operations common to all the tools be highlighted
    """
    if graph is None:
        graph = pgv.AGraph(directed=True, rankdir="BT")

    for toolURI in listToolURI:
        addToolAndAnnotationsToGraph(
            toolURI,
            graph=graph,
            showTopics=showTopics,
            showOperations=showOperations,
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
            node.attr["penwidth"] = "3"  # Thicker border
            node.attr["style"] = "filled,bold"  # Bold border and filled style
            node.attr["fillcolor"] = "red"  # Red fill to highlight
            node.attr["highlightIntersection"] = "True"  # Mark this node as highlighted

    return graph


def getScoreColorRGB(scoreValue, scoreMaxValue, color="red"):
    """Return the RGB color (in hex) associated to a score, varying from white (0 score)
    to a target color at the maximum score.

    Keyword arguments:
    scoreValue -- the score to consider (positive number)
    scoreMaxValue -- the maximum value for the score (positive number, >= scoreValue)
    color -- name of the color for the gradient. Possible values are:
             "red", "green", "blue", "orange", "yellow", "pink", "grey".
             (default: "red")
    """
    # Avoid division by 0.
    scoreMaxValue = max(1, scoreMaxValue)
    fraction = scoreValue / scoreMaxValue

    # Define the target color values.
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
        target = (255, 192, 203)  # Light pink (alternative: hot pink (255,105,180))
    elif color == "grey":
        target = (128, 128, 128)
    else:
        return "#ffffff"

    # Interpolate from white (255,255,255) to the target color.
    r = int(255 - fraction * (255 - target[0]))
    g = int(255 - fraction * (255 - target[1]))
    b = int(255 - fraction * (255 - target[2]))

    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def colorGraphNodesAccordingToScore(
    graph, dictTopicScore, dictOperationScore, color="red"
):
    """Modify nodes color according to a score associated to the node. Nodes that are not associated to a score are unaffected.

    Keyword arguments:
    graph -- the graph to be modified
    dictTopicScore -- a dictionary {nodeIdent -> score} for topic nodes
    dictOperationScore -- a dictionary {nodeIdent -> score} for operation nodes
    color -- name of the color that will vary according to scores. Possible values are "red", "green" or "blue" (default: "red")
    """
    topicMaxValue = max(dictTopicScore.values(), default=0)
    operationMaxValue = max(dictOperationScore.values(), default=0)
    for currentNode in graph.nodes_iter():
        if "nodeType" in currentNode.attr.keys():
            currentNodeType = currentNode.attr["nodeType"]
            # currentNode.attr['style'] = conceptStyle[currentNodeType]
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
    """Return a score for a tool computed as the sum of the scores of its annotations (0 if the tool has no annotation).

    Keyword arguments:
    toolURI -- the URI for the tool
    transitive -- also consider the ancestors of the topics or operations directly associated to the tool (default: False)
    dictTopicScore -- a dictionary {nodeIdent -> score} for topic nodes (default: {})
    dictOperationScore -- a dictionary {nodeIdent -> score} for operation nodes (default: {})
    """
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


def getMutualInformation(
    toolURI1: str,
    toolURI2: str,
    dfToolTopicTransitive,
    dfToolOperationTransitive,
    nbToolsWithTopic: int,
    dictAnnotationPairsMutualInformation: dict | None = None,
):
    """
    Compute Mutual Information between two tools based on their
    transitive EDAM Topic + Operation annotations.

    Parameters
    ----------
    toolURI1 : str
    toolURI2 : str
    dfToolTopicTransitive : pd.DataFrame
        DataFrame with columns ["tool", "topic"].
    dfToolOperationTransitive : pd.DataFrame
        DataFrame with columns ["tool", "operation"].
    nbToolsWithTopic : int
        Total number of tools (normalizing denominator).
    dictAnnotationPairsMutualInformation : dict, optional
        Cache for mutual information computation.

    Returns
    -------
    float
        Mutual information score.
    """

    mi = 0.0

    # ----------------------------
    # Get annotation lists
    # ----------------------------
    annotT1 = (
        dfToolTopicTransitive[dfToolTopicTransitive["tool"] == toolURI1][
            "topic"
        ].to_list()
        + dfToolOperationTransitive[dfToolOperationTransitive["tool"] == toolURI1][
            "operation"
        ].to_list()
    )

    annotT2 = (
        dfToolTopicTransitive[dfToolTopicTransitive["tool"] == toolURI2][
            "topic"
        ].to_list()
        + dfToolOperationTransitive[dfToolOperationTransitive["tool"] == toolURI2][
            "operation"
        ].to_list()
    )

    # Pre-cache annotation -> tool sets for speed
    topic2tools = {
        topic: set(
            dfToolTopicTransitive[dfToolTopicTransitive["topic"] == topic]["tool"]
        )
        for topic in set(annotT1 + annotT2)
    }

    # ----------------------------
    # Compute MI
    # ----------------------------
    for a1 in annotT1:

        toolsA1 = topic2tools[a1]
        pa1 = len(toolsA1) / nbToolsWithTopic

        for a2 in annotT2:

            # Check cached value
            if dictAnnotationPairsMutualInformation is not None:
                if (a1, a2) in dictAnnotationPairsMutualInformation:
                    mi += dictAnnotationPairsMutualInformation[(a1, a2)]
                    continue

            toolsA2 = topic2tools[a2]
            pa2 = len(toolsA2) / nbToolsWithTopic
            pa1a2 = len(toolsA1.intersection(toolsA2)) / nbToolsWithTopic

            if pa1 * pa2 * pa1a2 != 0:
                result = pa1a2 * np.log2(pa1a2 / (pa1 * pa2))
            else:
                result = 0.0

            mi += result

            # Save cached values (symmetric)
            if dictAnnotationPairsMutualInformation is not None:
                dictAnnotationPairsMutualInformation[(a1, a2)] = result
                dictAnnotationPairsMutualInformation[(a2, a1)] = result

    return mi


def get_nb_tools() -> int:
    """
    Execute SPARQL query to count distinct SoftwareApplication tools.

    Uses global variables:
      - endpointURL
      - prefixes

    Returns
    -------
    int
        Number of distinct tools.
    """
    query = """
    SELECT (COUNT(DISTINCT ?tool) AS ?nbTools)
    WHERE {
      ?tool rdf:type sc:SoftwareApplication .
    }
    """

    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    nb_tools = int(results["results"]["bindings"][0]["nbTools"]["value"])
    return nb_tools


def get_tools_dataframe() -> pd.DataFrame:
    """
    Execute SPARQL query to get SoftwareApplication tools and their labels.

    Returns
    -------
    pd.DataFrame
        DataFrame with tools and their labels.
    """
    query = """
    SELECT DISTINCT ?tool ?toolLabel
    WHERE {
      ?tool rdf:type sc:SoftwareApplication .
      ?tool sc:name ?toolLabel .
    }
    """

    dfTool = sparqldataframe.query(endpointURL, prefixes + query)
    dfTool.to_csv("Dataframe/dfTool.tsv.bz2", sep="\t", index=False)
    return dfTool


def get_tools_topics_dataframe() -> pd.DataFrame:
    """
    Execute SPARQL query to get SoftwareApplication tools and their topics.

    Returns
    -------
    pd.DataFrame
        DataFrame with tools, topics, and topic labels.
    """
    query = """
    SELECT DISTINCT ?tool ?topic ?topicLabel
    WHERE {
      ?tool rdf:type sc:SoftwareApplication .
      ?tool sc:applicationSubCategory ?topic .
      ?topic rdf:type owl:Class .
      FILTER NOT EXISTS { ?topic rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?topic rdfs:label ?tLabel }
      BIND(COALESCE(?tLabel, "") AS ?topicLabel)
    }
    """

    dfToolTopic = sparqldataframe.query(endpointURL, prefixes + query)
    dfToolTopic.to_csv("Dataframe/dfToolTopic.tsv.bz2", sep="\t", index=False)
    return dfToolTopic


def get_tools_topics_transitive_dataframe() -> pd.DataFrame:
    """
    Execute SPARQL query to get SoftwareApplication tools and their topics (including ancestors).

    Returns
    -------
    pd.DataFrame
        DataFrame with tools, topics, and topic labels (transitive closure).
    """
    query = """
    SELECT DISTINCT ?tool ?topic ?topicLabel
    WHERE {
      ?tool rdf:type sc:SoftwareApplication .
      ?tool sc:applicationSubCategory/(rdfs:subClassOf*) ?topic .
      ?topic rdf:type owl:Class .
      FILTER NOT EXISTS { ?topic rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?topic rdfs:label ?tLabel }
      BIND(COALESCE(?tLabel, "") AS ?topicLabel)
    }
    """

    dfToolTopicTransitive = sparqldataframe.query(endpointURL, prefixes + query)
    dfToolTopicTransitive.to_csv(
        "Dataframe/dfToolTopicTransitive.tsv.bz2", sep="\t", index=False
    )
    return dfToolTopicTransitive


def get_tools_operations_label_dataframe() -> pd.DataFrame:
    """
    Execute SPARQL query to get SoftwareApplication tools and their operations.

    Returns
    -------
    pd.DataFrame
        DataFrame with tools, operations, and operation labels.
    """
    query = """
    SELECT DISTINCT ?tool ?operation ?operationLabel
    WHERE {
      ?tool rdf:type sc:SoftwareApplication .
      ?tool sc:featureList ?operation .
      ?operation rdf:type owl:Class .
      FILTER NOT EXISTS { ?operation rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?operation rdfs:label ?oLabel }
      BIND(COALESCE(?oLabel, "") AS ?operationLabel)
    }
    """

    dfToolOperation = sparqldataframe.query(endpointURL, prefixes + query)
    dfToolOperation.to_csv("Dataframe/dfToolOperation.tsv.bz2", sep="\t", index=False)
    return dfToolOperation


def get_tools_operations_transitive_dataframe() -> pd.DataFrame:
    """
    Execute SPARQL query to get SoftwareApplication tools and their operations (including ancestors).

    Returns
    -------
    pd.DataFrame
        DataFrame with tools, operations, and operation labels (transitive closure).
    """
    query = """
    SELECT DISTINCT ?tool ?operation ?operationLabel
    WHERE {
      ?tool rdf:type sc:SoftwareApplication .
      ?tool sc:featureList/(rdfs:subClassOf*) ?operation .
      ?operation rdf:type owl:Class .
      FILTER NOT EXISTS { ?operation rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?operation rdfs:label ?oLabel }
      BIND(COALESCE(?oLabel, "") AS ?operationLabel)
    }
    """

    dfToolOperationTransitive = sparqldataframe.query(endpointURL, prefixes + query)
    dfToolOperationTransitive.to_csv(
        "Dataframe/dfToolOperationTransitive.tsv.bz2", sep="\t", index=False
    )
    return dfToolOperationTransitive


def get_dftools_with_nbTopics_nbOperations(
    dfTool: pd.DataFrame,
    dfToolTopicTransitive: pd.DataFrame,
    dfToolOperationTransitive: pd.DataFrame,
    output_path: str = "Dataframe/dftools_nbTopics_nbOperations.tsv.bz2",
) -> pd.DataFrame:
    """
    Generate an updated dfTool dataframe including nbTopics and nbOperations counts,
    and save it as a .tsv.bz2 file.

    Parameters
    ----------
    dfTool : pd.DataFrame
        Base dataframe containing tool information (must include a 'tool' column).
    dfToolTopicTransitive : pd.DataFrame
        Dataframe mapping tools to topics (must include a 'tool' column).
    dfToolOperationTransitive : pd.DataFrame
        Dataframe mapping tools to operations (must include a 'tool' column).
    output_path : str, optional
        Path to save the resulting dataframe, default is "dftools.tsv.bz2".

    Returns
    -------
    pd.DataFrame
        The updated dfTool dataframe with 'nbTopics' and 'nbOperations' columns.
    """

    # Compute number of topics per tool
    dfToolNbTopics = (
        dfToolTopicTransitive.groupby(by="tool")
        .size()
        .reset_index(name="nbTopics")
        .sort_values(by="nbTopics", ascending=False)
    )

    # Compute number of operations per tool
    dfToolNbOperations = (
        dfToolOperationTransitive.groupby(by="tool")
        .size()
        .reset_index(name="nbOperations")
        .sort_values(by="nbOperations", ascending=False)
    )

    # Join to main dfTool
    dfTool = dfTool.join(dfToolNbTopics.set_index("tool"), on="tool")
    dfTool = dfTool.join(dfToolNbOperations.set_index("tool"), on="tool")

    # Fill missing values and cast types
    dfTool["nbTopics"] = dfTool["nbTopics"].fillna(0).astype(int)
    dfTool["nbOperations"] = dfTool["nbOperations"].fillna(0).astype(int)

    # Save to compressed TSV
    dfTool.to_csv(output_path, sep="\t", index=False, compression="bz2")

    return dfTool


def generate_df_redundancy_topic(
    output_path: str = "Dataframe/dfToolTopic_redundancy.tsv.bz2",
) -> pd.DataFrame:
    """
    Generate the df_redundancy_topic DataFrame by running a SPARQL query using
    pre-defined global variables (endpointURL, prefixes, edamURI), and save it
    as a compressed .tsv.bz2 file.

    Parameters
    ----------
    output_path : str, optional
        Path to save the resulting DataFrame. Default is 'Dataframe/dfToolTopic_redundancy.tsv.bz2'.

    Returns
    -------
    pd.DataFrame
        The redundancy topic DataFrame.
    """

    # Define SPARQL query
    redundancyQuery = """
    SELECT DISTINCT ?tool ?redundantDirectTopic ?redundantDirectTopicLabel ?directTopic ?directTopicLabel
    WHERE {
      ?tool sc:applicationSubCategory ?redundantDirectTopic .
      ?redundantDirectTopic rdf:type owl:Class .
      FILTER NOT EXISTS { ?redundantDirectTopic rdfs:subClassOf? owl:DeprecatedClass }

      ?tool sc:applicationSubCategory ?directTopic .
      ?directTopic rdf:type owl:Class .
      FILTER NOT EXISTS { ?directTopic rdfs:subClassOf? owl:DeprecatedClass }

      ?directTopic rdfs:subClassOf+ ?redundantDirectTopic .

      OPTIONAL { ?directTopic rdfs:label ?tLabel }
      BIND(COALESCE(?tLabel, "") AS ?directTopicLabel)

      OPTIONAL { ?redundantDirectTopic rdfs:label ?rtLabel }
      BIND(COALESCE(?rtLabel, "") AS ?redundantDirectTopicLabel)
    }
    """

    # Run SPARQL query
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + redundancyQuery)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    # Extract results into a list of dictionaries
    data = [
        {
            "Tool": result["tool"]["value"].replace(edamURI, ""),
            "Direct Topic ID": edamURI
            + result.get("directTopic", {}).get("value", "").replace(edamURI, ""),
            "Direct Topic Label": result.get("directTopicLabel", {}).get("value", ""),
            "Redundant Topic ID": edamURI
            + result.get("redundantDirectTopic", {})
            .get("value", "")
            .replace(edamURI, ""),
            "Redundant Topic Label": result.get("redundantDirectTopicLabel", {}).get(
                "value", ""
            ),
        }
        for result in results["results"]["bindings"]
    ]

    # Create DataFrame
    df_redundancy_topic = pd.DataFrame(data)

    # Save to compressed TSV
    df_redundancy_topic.to_csv(output_path, sep="\t", index=False, compression="bz2")

    return df_redundancy_topic


def generate_df_redundancy_operation(
    output_path: str = "Dataframe/dfToolOperation_redundancy.tsv.bz2",
) -> pd.DataFrame:
    """
    Run a SPARQL query using pre-defined variables (endpointURL, prefixes, edamURI)
    and save the redundancy operation results as a compressed .tsv.bz2 file.
    """

    operationRedundancyQuery = """
    SELECT DISTINCT ?tool ?redundantDirectOperation ?redundantDirectOperationLabel ?directOperation ?directOperationLabel
    WHERE {
      ?tool sc:featureList ?redundantDirectOperation .
      ?redundantDirectOperation rdf:type owl:Class .
      FILTER NOT EXISTS { ?redundantDirectOperation rdfs:subClassOf? owl:DeprecatedClass }

      ?tool sc:featureList ?directOperation .
      ?directOperation rdf:type owl:Class .
      FILTER NOT EXISTS { ?directOperation rdfs:subClassOf? owl:DeprecatedClass }

      ?directOperation rdfs:subClassOf+ ?redundantDirectOperation .

      OPTIONAL { ?directOperation rdfs:label ?dLabel }
      BIND(COALESCE(?dLabel, "") AS ?directOperationLabel)

      OPTIONAL { ?redundantDirectOperation rdfs:label ?rdLabel }
      BIND(COALESCE(?rdLabel, "") AS ?redundantDirectOperationLabel)
    }
    """

    # Run SPARQL query
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes + operationRedundancyQuery)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    # Parse results
    data = [
        {
            "Tool": result["tool"]["value"].replace(edamURI, ""),
            "Direct Operation ID": edamURI
            + result.get("directOperation", {}).get("value", "").replace(edamURI, ""),
            "Direct Operation Label": result.get("directOperationLabel", {}).get(
                "value", ""
            ),
            "Redundant Operation ID": edamURI
            + result.get("redundantDirectOperation", {})
            .get("value", "")
            .replace(edamURI, ""),
            "Redundant Operation Label": result.get(
                "redundantDirectOperationLabel", {}
            ).get("value", ""),
        }
        for result in results["results"]["bindings"]
    ]

    # Create DataFrame
    df_redundancy_operation = pd.DataFrame(data)

    # Save to compressed TSV
    df_redundancy_operation.to_csv(
        output_path, sep="\t", index=False, compression="bz2"
    )

    return df_redundancy_operation


def generate_df_topic_no_redundancy(
    dfToolTopic_path="Dataframe/dfToolTopic.tsv.bz2",
    df_redundancy_topic_path="Dataframe/dfToolTopic_redundancy.tsv.bz2",
    output_path="Dataframe/df_topic_no_redundancy.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute df_topic_no_redundancy using already-generated files.
    """

    # Load both dataframes instead of using globals
    dfToolTopic = pd.read_csv(dfToolTopic_path, sep="\t", compression="bz2")
    df_redundancy_topic = pd.read_csv(
        df_redundancy_topic_path, sep="\t", compression="bz2"
    )

    # Build a set of redundant (Tool, Redundant Topic ID)
    redundant_pairs = set(
        zip(df_redundancy_topic["Tool"], df_redundancy_topic["Redundant Topic ID"])
    )

    # Apply filtering
    df_topic_no_redundancy = dfToolTopic[
        ~dfToolTopic[["tool", "topic"]].apply(tuple, axis=1).isin(redundant_pairs)
    ]

    df_topic_no_redundancy.to_csv(output_path, sep="\t", index=False, compression="bz2")

    return df_topic_no_redundancy


def generate_df_operation_no_redundancy(
    dfToolOperation_path="Dataframe/dfToolOperation.tsv.bz2",
    df_redundancy_operation_path="Dataframe/dfToolOperation_redundancy.tsv.bz2",
    output_path="Dataframe/df_operation_no_redundancy.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute df_operation_no_redundancy using already-generated files.
    """

    dfToolOperation = pd.read_csv(dfToolOperation_path, sep="\t", compression="bz2")
    df_redundancy_operation = pd.read_csv(
        df_redundancy_operation_path, sep="\t", compression="bz2"
    )

    redundant_pairs = set(
        zip(
            df_redundancy_operation["Tool"],
            df_redundancy_operation["Redundant Operation ID"],
        )
    )

    df_operation_no_redundancy = dfToolOperation[
        ~dfToolOperation[["tool", "operation"]]
        .apply(tuple, axis=1)
        .isin(redundant_pairs)
    ]

    df_operation_no_redundancy.to_csv(
        output_path, sep="\t", index=False, compression="bz2"
    )

    return df_operation_no_redundancy


def generate_dfTool_transitive(
    dfTool_path="Dataframe/dfTool.tsv.bz2",
    dfToolTopicTransitive_path="Dataframe/dfToolTopicTransitive.tsv.bz2",
    dfToolOperationTransitive_path="Dataframe/dfToolOperationTransitive.tsv.bz2",
    output_path="Dataframe/dfTool_Transitive.tsv.bz2",
) -> pd.DataFrame:
    """
    Generate dfTool with nbTopics and nbOperations based on transitive closures,
    using already-generated TSV files instead of global variables.
    """

    # Load required dataframes
    dfTool = pd.read_csv(dfTool_path, sep="\t", compression="bz2")
    dfToolTopicTransitive = pd.read_csv(
        dfToolTopicTransitive_path, sep="\t", compression="bz2"
    )
    dfToolOperationTransitive = pd.read_csv(
        dfToolOperationTransitive_path, sep="\t", compression="bz2"
    )

    # Number of topic matches (transitive)
    dfToolNbTopics = (
        dfToolTopicTransitive.groupby("tool").size().reset_index(name="nbTopics")
    )

    # Number of operation matches (transitive)
    dfToolNbOperations = (
        dfToolOperationTransitive.groupby("tool")
        .size()
        .reset_index(name="nbOperations")
    )

    # Join with dfTool
    dfTool_T = dfTool.copy()
    dfTool_T = dfTool_T.join(dfToolNbTopics.set_index("tool"), on="tool")
    dfTool_T = dfTool_T.join(dfToolNbOperations.set_index("tool"), on="tool")

    # Fill missing values
    dfTool_T["nbTopics"] = dfTool_T["nbTopics"].fillna(0).astype(int)
    dfTool_T["nbOperations"] = dfTool_T["nbOperations"].fillna(0).astype(int)

    dfTool_T.to_csv(output_path, sep="\t", index=False, compression="bz2")
    return dfTool_T


def generate_dfTool_no_transitive(
    dfTool_path="Dataframe/dfTool.tsv.bz2",
    dfToolTopic_path="Dataframe/dfToolTopic.tsv.bz2",
    dfToolOperation_path="Dataframe/dfToolOperation.tsv.bz2",
    output_path="Dataframe/dfTool_NoTransitive.tsv.bz2",
) -> pd.DataFrame:
    """
    Generate dfTool with nbTopics and nbOperations without transitive closure,
    using already-generated TSV files instead of global variables.
    """

    # Load required dataframes
    dfTool = pd.read_csv(dfTool_path, sep="\t", compression="bz2")
    dfToolTopic = pd.read_csv(dfToolTopic_path, sep="\t", compression="bz2")
    dfToolOperation = pd.read_csv(dfToolOperation_path, sep="\t", compression="bz2")

    # Number of topic matches (non-transitive)
    dfToolNbTopics = dfToolTopic.groupby("tool").size().reset_index(name="nbTopics")

    # Number of operation matches (non-transitive)
    dfToolNbOperations = (
        dfToolOperation.groupby("tool").size().reset_index(name="nbOperations")
    )

    # Join with dfTool
    dfTool_NT = dfTool.copy()
    dfTool_NT = dfTool_NT.join(dfToolNbTopics.set_index("tool"), on="tool")
    dfTool_NT = dfTool_NT.join(dfToolNbOperations.set_index("tool"), on="tool")

    # Fill missing values
    dfTool_NT["nbTopics"] = dfTool_NT["nbTopics"].fillna(0).astype(int)
    dfTool_NT["nbOperations"] = dfTool_NT["nbOperations"].fillna(0).astype(int)

    dfTool_NT.to_csv(output_path, sep="\t", index=False, compression="bz2")
    return dfTool_NT


def generate_dfTool_no_transitive_no_redundancy(
    dfTool_path="Dataframe/dfTool.tsv.bz2",
    df_topic_no_redundancy_path="Dataframe/df_topic_no_redundancy.tsv.bz2",
    df_operation_no_redundancy_path="Dataframe/df_operation_no_redundancy.tsv.bz2",
    output_path="Dataframe/dfTool_NoTransitive_NoRedundancy.tsv.bz2",
) -> pd.DataFrame:

    dfTool = pd.read_csv(dfTool_path, sep="\t", compression="bz2")
    df_topic_no_redundancy = pd.read_csv(
        df_topic_no_redundancy_path, sep="\t", compression="bz2"
    )
    df_operation_no_redundancy = pd.read_csv(
        df_operation_no_redundancy_path, sep="\t", compression="bz2"
    )

    # nbTopics
    dfToolNbTopics = (
        df_topic_no_redundancy.groupby("tool").size().reset_index(name="nbTopics")
    )

    # nbOperations
    dfToolNbOperations = (
        df_operation_no_redundancy.groupby("tool")
        .size()
        .reset_index(name="nbOperations")
    )

    dfTool_NR = dfTool.copy()
    dfTool_NR = dfTool_NR.join(dfToolNbTopics.set_index("tool"), on="tool")
    dfTool_NR = dfTool_NR.join(dfToolNbOperations.set_index("tool"), on="tool")

    dfTool_NR["nbTopics"] = dfTool_NR["nbTopics"].fillna(0).astype(int)
    dfTool_NR["nbOperations"] = dfTool_NR["nbOperations"].fillna(0).astype(int)

    dfTool_NR.to_csv(output_path, sep="\t", index=False, compression="bz2")

    return dfTool_NR


def get_dfToolTopic_NotOWLClass() -> pd.DataFrame:
    """
    Retrieve tools whose assigned topic is NOT an owl:Class.

    Uses global variables:
      - endpointURL
      - prefixes
      - sparqldataframe

    Returns
    -------
    pd.DataFrame
        DataFrame named dfToolTopic_NotOWLClass containing:
        ['tool', 'topic']
    """

    query = """
    # tools with topic that is not an owl class
    SELECT DISTINCT ?tool ?topic
    WHERE {
      ?tool rdf:type sc:SoftwareApplication .
      ?tool sc:applicationSubCategory/(rdfs:subClassOf*) ?topic .
      FILTER NOT EXISTS {
        ?topic rdf:type owl:Class .
      }
    }
    """

    # Execute SPARQL and produce DataFrame
    dfToolTopic_NotOWLClass = sparqldataframe.query(endpointURL, prefixes + query)

    # Save to compressed TSV file
    dfToolTopic_NotOWLClass.to_csv(
        "Dataframe/dfToolTopic_NotOWLClass.tsv.bz2",
        sep="\t",
        index=False,
        compression="bz2",
    )

    return dfToolTopic_NotOWLClass


def get_dfTool_ObsoleteOperation() -> pd.DataFrame:
    """
    Retrieve tools that have obsolete operation annotations with suggestions covered by other valid annotations.

    Returns
    -------
    pd.DataFrame
        DataFrame named dfTool_ObsoleteOperation containing:
        ['tool']
    """

    query = """
# tools with an obsolete operation annotation that has suggestions covered by another valid annotation
# operation: compatible (rdfs:subClassOf* => 108) or more precise (rdfs:subClassOf+ => 70) than a valid operation

SELECT DISTINCT ?tool #?validItem ?alternativeItem
WHERE {
  VALUES ?annotationType { sc:featureList } # Operation

  ?tool rdf:type sc:SoftwareApplication .
  ?tool ?annotationType ?validItem .
  ?tool ?annotationType ?deprecatedItem .
  
  { ?deprecatedItem rdfs:subClassOf owl:DeprecatedClass }
  UNION
  { ?deprecatedItem owl:deprecated true . }
  UNION
  { ?deprecatedItem owl:deprecated "true" . }
  UNION 
  { ?deprecatedItem owl:deprecated "True" . }
  
  ?deprecatedItem oboInOwl:consider ?alternativeItem .
  ?alternativeItem rdfs:subClassOf+ ?validItem .
}
"""

    # Execute SPARQL and return DataFrame
    dfTool_ObsoleteOperation = sparqldataframe.query(endpointURL, prefixes + query)

    # Save to compressed TSV
    dfTool_ObsoleteOperation.to_csv(
        "Dataframe/dfTool_ObsoleteOperation.tsv.bz2",
        sep="\t",
        index=False,
        compression="bz2",
    )

    return dfTool_ObsoleteOperation


def compute_topic_metrics(
    dfToolTopicTransitive_path="Dataframe/dfToolTopicTransitive.tsv.bz2",
    output_path="Dataframe/dfTopicmetrics.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute topic metrics including Information Content (IC) and entropy
    Parameters:
    -----------
    dfToolTopicTransitive_path : str
        Path to the transitive tool-topic dataframe
    output_path : str
        Path to save the topic metrics dataframe

    Returns:
    --------
    pd.DataFrame
        Dataframe containing topic metrics
    """

    # Read input dataframes
    dfToolTopicTransitive = pd.read_csv(
        dfToolTopicTransitive_path, sep="\t", compression="bz2"
    )

    # Calculate number of tools per topic
    dfTopicNbTools = (
        dfToolTopicTransitive[["tool", "topic"]]
        .groupby(by="topic")
        .size()
        .reset_index(name="nbTools")
        .sort_values(by="nbTools", ascending=False)
    )

    # Create topic dataframe with unique topics and their labels
    dfTopic = (
        dfToolTopicTransitive[["topic", "topicLabel"]]
        .drop_duplicates(subset=["topic", "topicLabel"], keep="first")
        .reset_index(drop=True)
    )

    # Join with tool counts
    dfTopic = dfTopic.join(dfTopicNbTools.set_index("topic"), on="topic")

    # Calculate total number of tools with topics
    nbToolsWithTopic = dfToolTopicTransitive["tool"].nunique()

    # Calculate metrics
    dfTopic["frequence"] = dfTopic["nbTools"] / nbToolsWithTopic
    dfTopic["IC"] = -np.log2(dfTopic["frequence"])
    dfTopic["entropy"] = dfTopic["frequence"] * dfTopic["IC"]

    # Save results
    dfTopic.to_csv(output_path, sep="\t", index=False, compression="bz2")
    dfTopicmetrics = dfTopic

    return dfTopicmetrics


def compute_topic_metrics_NT(
    dfToolTopic_path="Dataframe/dfToolTopic.tsv.bz2",
    output_path="Dataframe/dfTopicmetrics_NT.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute topic metrics without transitive closure including Information Content (IC) and entropy,
    and save the results to global variables.

    Parameters:
    -----------
    dfToolTopic_path : str
        Path to the tool-topic dataframe (without transitive closure)
    output_path : str
        Path to save the topic metrics dataframe

    Returns:
    --------
    pd.DataFrame
        Dataframe containing topic metrics without transitive closure
    """

    # Read input dataframes
    dfToolTopic = pd.read_csv(dfToolTopic_path, sep="\t", compression="bz2")

    # Calculate total number of tools with topics
    nbToolsWithTopic = dfToolTopic["tool"].nunique()

    # Calculate number of tools per topic
    dfTopicNbTools = (
        dfToolTopic[["tool", "topic"]]
        .groupby(by="topic")
        .size()
        .reset_index(name="nbTools")
        .sort_values(by="nbTools", ascending=False)
    )

    # Create topic dataframe with unique topics and their labels
    dfTopic = (
        dfToolTopic[["topic", "topicLabel"]]
        .drop_duplicates(subset=["topic", "topicLabel"], keep="first")
        .reset_index(drop=True)
    )

    # Join with tool counts
    dfTopic = dfTopic.join(dfTopicNbTools.set_index("topic"), on="topic")

    # Calculate metrics
    dfTopic["frequence"] = dfTopic["nbTools"] / nbToolsWithTopic
    dfTopic["IC"] = -np.log2(dfTopic["frequence"])
    dfTopic["entropy"] = dfTopic["frequence"] * dfTopic["IC"]

    # Save results
    dfTopic.to_csv(output_path, sep="\t", index=False, compression="bz2")
    dfTopicmetrics_NT = dfTopic

    return dfTopicmetrics_NT


def compute_operation_metrics(
    dfToolOperationTransitive_path="Dataframe/dfToolOperationTransitive.tsv.bz2",
    output_path="Dataframe/dfOperationmetrics.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute operation metrics including Information Content (IC) and entropy.

    Parameters:
    -----------
    dfToolOperationTransitive_path : str
        Path to the transitive tool-operation dataframe
    dfTool_path : str
        Path to the tool dataframe (to get total number of tools with operations)
    output_path : str
        Path to save the operation metrics dataframe

    Returns:
    --------
    pd.DataFrame
        Dataframe containing operation metrics
    """

    # Read input dataframes
    dfToolOperationTransitive = pd.read_csv(
        dfToolOperationTransitive_path, sep="\t", compression="bz2"
    )

    # Calculate total number of tools with operations
    nbToolsWithOperation = dfToolOperationTransitive["tool"].nunique()

    # Calculate number of tools per operation
    dfOperationNbTools = (
        dfToolOperationTransitive[["tool", "operation"]]
        .groupby(by="operation")
        .size()
        .reset_index(name="nbTools")
        .sort_values(by="nbTools", ascending=False)
    )

    # Create operation dataframe with unique operations and their labels
    dfOperation = (
        dfToolOperationTransitive[["operation", "operationLabel"]]
        .drop_duplicates(subset=["operation", "operationLabel"], keep="first")
        .reset_index(drop=True)
    )

    # Join with tool counts
    dfOperation = dfOperation.join(
        dfOperationNbTools.set_index("operation"), on="operation"
    )

    # Calculate metrics
    dfOperation["frequence"] = dfOperation["nbTools"] / nbToolsWithOperation
    dfOperation["IC"] = -np.log2(dfOperation["frequence"])
    dfOperation["entropy"] = dfOperation["frequence"] * dfOperation["IC"]

    # Save results
    dfOperation.to_csv(output_path, sep="\t", index=False, compression="bz2")
    dfOperationmetrics = dfOperation

    return dfOperationmetrics


def compute_operation_metrics_NT(
    dfToolOperation_path="Dataframe/dfToolOperation.tsv.bz2",
    output_path="Dataframe/dfOperationmetrics_NT.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute operation metrics without transitive closure including Information Content (IC) and entropy.

    Parameters:
    -----------
    dfToolOperation_path : str
        Path to the tool-operation dataframe (without transitive closure)
    output_path : str
        Path to save the operation metrics dataframe

    Returns:
    --------
    pd.DataFrame
        Dataframe containing operation metrics without transitive closure
    """

    # Read input dataframes
    dfToolOperation = pd.read_csv(dfToolOperation_path, sep="\t", compression="bz2")

    # Calculate total number of tools with operations
    nbToolsWithOperation = dfToolOperation["tool"].nunique()

    # Calculate number of tools per operation
    dfOperationNbTools = (
        dfToolOperation[["tool", "operation"]]
        .groupby(by="operation")
        .size()
        .reset_index(name="nbTools")
        .sort_values(by="nbTools", ascending=False)
    )

    # Create operation dataframe with unique operations and their labels
    dfOperation = (
        dfToolOperation[["operation", "operationLabel"]]
        .drop_duplicates(subset=["operation", "operationLabel"], keep="first")
        .reset_index(drop=True)
    )

    # Join with tool counts
    dfOperation = dfOperation.join(
        dfOperationNbTools.set_index("operation"), on="operation"
    )

    # Calculate metrics
    dfOperation["frequence"] = dfOperation["nbTools"] / nbToolsWithOperation
    dfOperation["IC"] = -np.log2(dfOperation["frequence"])
    dfOperation["entropy"] = dfOperation["frequence"] * dfOperation["IC"]

    # Save results
    dfOperation.to_csv(output_path, sep="\t", index=False, compression="bz2")
    dfOperationmetrics_NT = dfOperation

    return dfOperationmetrics_NT


def compute_tool_metrics_with_transitive(
    dfToolTopicTransitive_path="Dataframe/dfToolTopicTransitive.tsv.bz2",
    dfToolOperationTransitive_path="Dataframe/dfToolOperationTransitive.tsv.bz2",
    dfTopic_metrics_path="Dataframe/dfTopicmetrics.tsv.bz2",
    dfOperation_metrics_path="Dataframe/dfOperationmetrics.tsv.bz2",
    dfTool_path="Dataframe/dfTool.tsv.bz2",
    output_path="Dataframe/dfToolallmetrics.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute combined topic and operation metrics for tools using transitive relationships.

    Parameters:
    -----------
    dfToolTopicTransitive_path : str
        Path to the transitive tool-topic dataframe
    dfToolOperationTransitive_path : str
        Path to the transitive tool-operation dataframe
    dfTopic_metrics_path : str
        Path to the topic metrics dataframe
    dfOperation_metrics_path : str
        Path to the operation metrics dataframe
    dfTool_path : str
        Path to the tool dataframe
    output_path : str
        Path to save the combined tool metrics dataframe

    Returns:
    --------
    pd.DataFrame
        Dataframe containing tools with combined topic and operation metrics
    """

    # Read input dataframes
    dfToolTopicTransitive = pd.read_csv(
        dfToolTopicTransitive_path, sep="\t", compression="bz2"
    )
    dfToolOperationTransitive = pd.read_csv(
        dfToolOperationTransitive_path, sep="\t", compression="bz2"
    )
    dfTopic = pd.read_csv(dfTopic_metrics_path, sep="\t", compression="bz2")
    dfOperation = pd.read_csv(dfOperation_metrics_path, sep="\t", compression="bz2")
    dfTool = pd.read_csv(dfTool_path, sep="\t", compression="bz2")

    # Calculate topic scores
    df_topic_scores = (
        dfToolTopicTransitive.join(
            dfTopic[["topic", "IC", "entropy"]].set_index("topic"), on="topic"
        )[["tool", "IC"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"IC": "topicScore"})
        .reset_index()
    )

    dfTool = dfTool.join(df_topic_scores.set_index("tool"), on="tool")
    dfTool["topicScore"] = dfTool["topicScore"].fillna(0)

    # Calculate operation scores
    df_operation_scores = (
        dfToolOperationTransitive.join(
            dfOperation[["operation", "IC", "entropy"]].set_index("operation"),
            on="operation",
        )[["tool", "IC"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"IC": "operationScore"})
        .reset_index()
    )

    dfTool = dfTool.join(df_operation_scores.set_index("tool"), on="tool")
    dfTool["operationScore"] = dfTool["operationScore"].fillna(0)

    # Calculate total score
    dfTool["score"] = dfTool["topicScore"] + dfTool["operationScore"]

    # Calculate topic entropy
    df_topic_entropy = (
        dfToolTopicTransitive.join(
            dfTopic[["topic", "IC", "entropy"]].set_index("topic"), on="topic"
        )[["tool", "entropy"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"entropy": "topicEntropy"})
        .reset_index()
    )

    dfTool = dfTool.join(df_topic_entropy.set_index("tool"), on="tool")
    dfTool["topicEntropy"] = dfTool["topicEntropy"].fillna(0)

    # Calculate operation entropy
    df_operation_entropy = (
        dfToolOperationTransitive.join(
            dfOperation[["operation", "IC", "entropy"]].set_index("operation"),
            on="operation",
        )[["tool", "entropy"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"entropy": "operationEntropy"})
        .reset_index()
    )

    dfTool = dfTool.join(df_operation_entropy.set_index("tool"), on="tool")
    dfTool["operationEntropy"] = dfTool["operationEntropy"].fillna(0)

    # Calculate total entropy
    dfTool["entropy"] = dfTool["topicEntropy"] + dfTool["operationEntropy"]

    # Save results
    dfTool.to_csv(output_path, sep="\t", index=False, compression="bz2")
    dfToolallmetrics = dfTool

    return dfToolallmetrics


def compute_tool_metrics_non_transitive(
    dfToolTopic_path="Dataframe/dfToolTopic.tsv.bz2",
    dfToolOperation_path="Dataframe/dfToolOperation.tsv.bz2",
    dfTopic_metrics_NT_path="Dataframe/dfTopicmetrics_NT.tsv.bz2",
    dfOperation_metrics_NT_path="Dataframe/dfOperationmetrics_NT.tsv.bz2",
    dfTool_path="Dataframe/dfTool.tsv.bz2",
    output_path="Dataframe/dfToolallmetrics_NT.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute combined topic and operation metrics for tools using non-transitive relationships.

    Parameters:
    -----------
    dfToolTopic_path : str
        Path to the non-transitive tool-topic dataframe
    dfToolOperation_path : str
        Path to the non-transitive tool-operation dataframe
    dfTopic_metrics_NT_path : str
        Path to the non-transitive topic metrics dataframe
    dfOperation_metrics_NT_path : str
        Path to the non-transitive operation metrics dataframe
    dfTool_path : str
        Path to the tool dataframe
    output_path : str
        Path to save the combined tool metrics dataframe (non-transitive)

    Returns:
    --------
    pd.DataFrame
        Dataframe containing tools with combined topic and operation metrics (non-transitive)
    """

    # Read input dataframes
    dfToolTopic = pd.read_csv(dfToolTopic_path, sep="\t", compression="bz2")
    dfToolOperation = pd.read_csv(dfToolOperation_path, sep="\t", compression="bz2")
    dfTopicmetrics_NT = pd.read_csv(
        dfTopic_metrics_NT_path, sep="\t", compression="bz2"
    )
    dfOperationmetrics_NT = pd.read_csv(
        dfOperation_metrics_NT_path, sep="\t", compression="bz2"
    )
    dfTool = pd.read_csv(dfTool_path, sep="\t", compression="bz2")

    # Calculate topic scores (non-transitive)
    df_topic_scores = (
        dfToolTopic.join(
            dfTopicmetrics_NT[["topic", "IC", "entropy"]].set_index("topic"), on="topic"
        )[["tool", "IC"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"IC": "topicScore"})
        .reset_index()
    )

    dfTool = dfTool.join(df_topic_scores.set_index("tool"), on="tool")
    dfTool["topicScore"] = dfTool["topicScore"].fillna(0)

    # Calculate operation scores (non-transitive)
    df_operation_scores = (
        dfToolOperation.join(
            dfOperationmetrics_NT[["operation", "IC", "entropy"]].set_index(
                "operation"
            ),
            on="operation",
        )[["tool", "IC"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"IC": "operationScore"})
        .reset_index()
    )

    dfTool = dfTool.join(df_operation_scores.set_index("tool"), on="tool")
    dfTool["operationScore"] = dfTool["operationScore"].fillna(0)

    # Calculate total score
    dfTool["score"] = dfTool["topicScore"] + dfTool["operationScore"]

    # Calculate topic entropy (non-transitive)
    df_topic_entropy = (
        dfToolTopic.join(
            dfTopicmetrics_NT[["topic", "IC", "entropy"]].set_index("topic"), on="topic"
        )[["tool", "entropy"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"entropy": "topicEntropy"})
        .reset_index()
    )

    dfTool = dfTool.join(df_topic_entropy.set_index("tool"), on="tool")
    dfTool["topicEntropy"] = dfTool["topicEntropy"].fillna(0)

    # Calculate operation entropy (non-transitive)
    df_operation_entropy = (
        dfToolOperation.join(
            dfOperationmetrics_NT[["operation", "IC", "entropy"]].set_index(
                "operation"
            ),
            on="operation",
        )[["tool", "entropy"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"entropy": "operationEntropy"})
        .reset_index()
    )

    dfTool = dfTool.join(df_operation_entropy.set_index("tool"), on="tool")
    dfTool["operationEntropy"] = dfTool["operationEntropy"].fillna(0)

    # Calculate total entropy
    dfTool["entropy"] = dfTool["topicEntropy"] + dfTool["operationEntropy"]

    # Save results
    dfTool.to_csv(output_path, sep="\t", index=False, compression="bz2")
    dfToolallmetrics_NT = dfTool

    return dfToolallmetrics_NT


def get_tool_url(tool_name: str) -> str:
    if tool_name.startswith("https://bio.tools/"):
        return tool_name
    return f"https://bio.tools/{tool_name}"


def normalize_tool_input(tools):
    if isinstance(tools, str):
        tools = [tools]
    return [get_tool_url(t) for t in tools]


def _resolve_annotation_type(value: str) -> str:
    """Resolve shorthand to full annotation type."""
    mapping = {
        "T": "Topic",
        "O": "Operation",
        "Topic": "Topic",
        "Operation": "Operation"
    }
    value = value.capitalize() if len(value) > 1 else value.upper()
    if value not in mapping:
        raise ValueError(f"Invalid annotation type: {value}")
    return mapping[value]


# -----------------------------------------
# FETCH ANNOTATION LOGIC (supports multi-type)
# -----------------------------------------

def fetch_annotations(tools, annotation_types=("Topic",), transitive=True, with_label=True):
    tools = normalize_tool_input(tools)

    # Normalize annotation types, remove duplicates
    resolved_types = [
        _resolve_annotation_type(a) for a in annotation_types
    ]
    resolved_types = list(dict.fromkeys(resolved_types))

    result = {tool: {} for tool in tools}

    # Determine which dataframe to use
    for ann_type in resolved_types:

        if ann_type == "Topic":
            df = dfToolTopicTransitive if transitive else dfToolTopic
            col_uri, col_label = "topic", "topicLabel"

        elif ann_type == "Operation":
            df = dfToolOperationTransitive if transitive else dfToolOperation
            col_uri, col_label = "operation", "operationLabel"

        # Collect annotations per tool
        for tool in tools:
            filtered = df[df["tool"] == tool]

            if with_label:
                annotations = [
                    {"URI": row[col_uri], "label": row[col_label]}
                    for _, row in filtered.iterrows()
                ]
            else:
                annotations = [
                    {"URI": row[col_uri]}
                    for _, row in filtered.iterrows()
                ]

            result[tool][ann_type] = annotations

    return result


# -----------------------------------------
# OUTPUT: JSON
# -----------------------------------------

def to_json(annotations: Dict) -> str:
    return json.dumps({"annotation": annotations}, indent=2)


# -----------------------------------------
# OUTPUT: TURTLE
# -----------------------------------------

def _format_as_turtle(results: Dict, include_labels: bool = True) -> str:

    ttl = [
        "@prefix biotools: <https://bio.tools/ontology/> .",
        "@prefix bsc: <http://bioschemas.org/> .",
        "@prefix bsct: <http://bioschemas.org/types/> .",
        "@prefix dcterms: <http://purl.org/dc/terms/> .",
        "@prefix edam: <http://edamontology.org/> .",
        "@prefix sc: <http://schema.org/> .",
        "@prefix schema: <https://schema.org/> .",
        "@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .",
        "@prefix ex: <http://example.org/> .",
        "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .",
        ""
    ]

    for tool_uri, ann_dict in results.items():
        tool_id = tool_uri.replace("https://bio.tools/", "")
        ttl.append(f"ex:Tool_{tool_id} a sc:SoftwareApplication ;")
        ttl.append(f"    sc:url <{tool_uri}> ;")

        triples = []

        # Topics
        for ann in ann_dict.get("Topic", []):
            ann_id = ann["URI"].split(":")[-1]
            triples.append(f"    sc:applicationSubCategory ex:{ann_id} ;")

        # Operations
        for ann in ann_dict.get("Operation", []):
            ann_id = ann["URI"].split(":")[-1]
            triples.append(f"    sc:featureList ex:{ann_id} ;")

        if triples:
            ttl.extend(triples)
            ttl[-1] = ttl[-1].rstrip(" ;") + " ."
        else:
            ttl[-1] = ttl[-1].rstrip(" ;") + " ."

        ttl.append("")

        # Annotation class definitions
        for ann_type, ann_list in ann_dict.items():
            for ann in ann_list:
                ann_id = ann["URI"].split(":")[-1]
                ttl.append(f"ex:{ann_id} a ex:{ann_type} ;")
                if include_labels and "label" in ann:
                    ttl.append(f'    rdfs:label "{ann["label"]}" .')
                else:
                    ttl[-1] = ttl[-1].rstrip(" ;") + " ."
                ttl.append("")
                
    """
    test_turtle = "\n".join(ttl)
    G = rdflib.Graph()
    G.parse(data = test_turtle, format = "turtle")
    """
    return "\n".join(ttl)


# -----------------------------------------
# OUTPUT: SPARQL
# -----------------------------------------

def _format_as_sparql(tool_urls: List[str], annotation_types: List[str], transitive: bool) -> str:

    tools_values = " ".join([f"<{url}>" for url in tool_urls])
    want_topic = any(a.lower() == "topic" for a in annotation_types)
    want_operation = any(a.lower() == "operation" for a in annotation_types)

    # Construct transitive topic hierarchy query
    if transitive and want_topic:
        return f"""PREFIX sc: <http://schema.org/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

CONSTRUCT {{
    ?tool sc:applicationSubCategory ?superTopic .
}}
WHERE {{
    VALUES ?tool {{ {tools_values} }}
    ?tool rdf:type sc:SoftwareApplication .
    ?tool sc:applicationSubCategory ?topic .
    ?topic rdfs:subClassOf* ?superTopic .
    ?superTopic rdf:type owl:Class .
}}"""

    # SELECT variant (non-transitive)
    query = f"""PREFIX sc: <http://schema.org/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

SELECT ?tool ?topic ?topicLabel ?operation ?operationLabel
WHERE {{
    VALUES ?tool {{ {tools_values} }}
    ?tool rdf:type sc:SoftwareApplication .
"""

    if want_topic:
        query += """
    OPTIONAL {
        ?tool sc:applicationSubCategory ?topic .
        OPTIONAL { ?topic rdfs:label ?topicLabel . }
    }
"""

    if want_operation:
        query += """
    OPTIONAL {
        ?tool sc:featureList ?operation .
        OPTIONAL { ?operation rdfs:label ?operationLabel . }
    }
"""

    query += "}"

    return query


DF_TOOL_NO_TRANS = pd.read_csv("Dataframe/dfTool_NoTransitive.tsv.bz2", sep="\t", compression='bz2')
DF_TOOL_TOPICS_OPS = pd.read_csv("Dataframe/dftools_nbTopics_nbOperations.tsv.bz2", sep="\t", compression='bz2')

def get_tool_metrics(tool: str, heritage: bool = True, metric: str = "all") -> dict:
    """
    Fetch metrics for a given tool with selectable type: 'ic', 'entropy', 'count', or 'all'.

    Parameters
    ----------
    tool : str
        Tool name or bio.tools URI.
    heritage : bool
        Whether to use inherited (transitive) metrics (default: True).
    metric : str
        Metric type to return ('ic', 'entropy', 'count', or 'all').

    Returns
    -------
    dict
        A dictionary with the tool’s metrics.
    """
    tool_url = tool if tool.startswith("https://bio.tools/") else f"https://bio.tools/{tool}"
    
    # Choose the appropriate dataframes based on heritage mode
    df_metrics = dfToolallmetrics if heritage else dfToolallmetrics_NT
    df_counts = DF_TOOL_TOPICS_OPS if heritage else DF_TOOL_NO_TRANS

    row_metrics = df_metrics[df_metrics['tool'] == tool_url]
    row_counts = df_counts[df_counts['tool'] == tool_url]

    if row_metrics.empty and row_counts.empty:
        raise ValueError(f"Tool not found: {tool_url}")

    metrics = row_metrics.iloc[0] if not row_metrics.empty else {}
    counts = row_counts.iloc[0] if not row_counts.empty else {}

    result = {"Tool": tool_url}

    # IC (Information Content)
    if metric in ("ic", "all"):
        result["topicScore"] = float(metrics.get("topicScore", 0))
        result["operationScore"] = float(metrics.get("operationScore", 0))
    
    # Entropy
    if metric in ("entropy", "all"):
        result["entropy"] = float(metrics.get("entropy", 0))

    # Count of annotations
    if metric in ("count", "all"):
        result["nbTopics"] = int(counts.get("nbTopics", 0))
        result["nbOperations"] = int(counts.get("nbOperations", 0))

    return result


def format_tool_annotations(metrics: dict) -> dict:
    """
    Format the tool metrics into the EDAM annotation dict structure.
    """
    return {
        "annotation": {
            "Metrics": metrics
        }
    }


def fetch_annotations_with_metrics(
    tools,
    annotation_types=("Topic", "Operation"),
    heritage=True,
    with_label=True,
    metric="all",
    include_annotations=True
):
    """
    Fetch EDAM annotations and/or tool metrics together in one structured dict.

    Parameters
    ----------
    tools : list[str] or str
        Tool name(s) or URLs.
    annotation_types : tuple[str]
        Annotation types to include ('Topic', 'Operation', etc.).
    heritage : bool
        Whether to use inherited (transitive) annotations (default: True).
    with_label : bool
        Whether to include labels in annotations.
    metric : str
        Metric type to include ('ic', 'entropy', 'count', 'all').
    include_annotations : bool
        Whether to include EDAM annotations (default: True).

    Returns
    -------
    dict
        { "https://bio.tools/<tool>" : {
            "annotation": {
              "Topic": [...],
              "Operation": [...],
              "Data": [],
              "Format": [],
              "Metrics": {...}
            }
        }}
    """
    tools = normalize_tool_input(tools)
    combined = {}

    # Fetch annotations only if requested
    annotation_data = {}
    if include_annotations:
        annotation_data = fetch_annotations(
            tools,
            annotation_types=annotation_types,
            transitive=heritage,  # kept for backward compat inside fetch_annotations()
            with_label=with_label
        )

    # Merge metrics and annotations per tool
    for tool in tools:
        try:
            metrics = get_tool_metrics(tool, heritage=heritage, metric=metric)
        except ValueError:
            metrics = {}

        # Safely populate annotation categories
        topic_anns = annotation_data.get(tool, {}).get("Topic", []) if include_annotations else []
        op_anns = annotation_data.get(tool, {}).get("Operation", []) if include_annotations else []
        data_anns = annotation_data.get(tool, {}).get("Data", []) if include_annotations else []
        format_anns = annotation_data.get(tool, {}).get("Format", []) if include_annotations else []

        combined[tool] = {
            "annotation": {
                "Topic": topic_anns,
                "Operation": op_anns,
                "Data": data_anns,
                "Format": format_anns,
                "Metrics": metrics
            }
        }

    return combined