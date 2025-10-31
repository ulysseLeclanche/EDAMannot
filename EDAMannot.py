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

current_dir = os.getcwd()
neighbor_dir = os.path.join(current_dir, "edam") 
bioschemas_file = os.path.join(neighbor_dir, "bioschemas-dump_05_01_2025.ttl")
edam_file = os.path.join(neighbor_dir, "EDAM_1.25.owl")

dfTool = pd.read_csv("Dataframe/dfTool.tsv.bz2", sep="\t")
dfToolTopic = pd.read_csv("Dataframe/dfToolTopic.tsv.bz2", sep="\t")
dfToolTopicTransitive = pd.read_csv("Dataframe/dfToolTopicTransitive.tsv.bz2", sep="\t")
dfToolOperation = pd.read_csv("Dataframe/dfToolOperation.tsv.bz2", sep="\t")
dfToolOperationTransitive = pd.read_csv("Dataframe/dfToolOperationTransitive.tsv.bz2", sep="\t")
nbTools = len(dfTool)
nbToolsWithTopic = dfToolTopic['tool'].nunique()
nbToolsWithOperation = dfToolOperation['tool'].nunique()


#Configuration SPARQL end point 
endpointURL = "http://localhost:3030/biotoolsEdam/query"
rdfFormat = "turtle"

#Import prefix : 
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

#Link tool
biotoolsURI = "https://bio.tools/"
biotoolsOntologyURI = "https://bio.tools/ontology/"
edamURI = "http://edamontology.org/"

def displaySparqlResults(results):
    '''
    Displays as HTML the result of a SPARQLWrapper query in a Jupyter notebook.
    
        Parameters:
            results (dictionnary): the result of a call to SPARQLWrapper.query().convert()
    '''
    variableNames = results['head']['vars']
    #tableCode = '<table><tr><th>{}</th></tr><tr>{}</tr></table>'.format('</th><th>'.join(variableNames), '</tr><tr>'.join('<td>{}</td>'.format('</td><td>'.join([row[vName]['value'] for vName in variableNames]))for row in results["results"]["bindings"]))
    tableCode = '<table><tr><th>{}</th></tr><tr>{}</tr></table>'.format('</th><th>'.join(variableNames), '</tr><tr>'.join('<td>{}</td>'.format('</td><td>'.join([row[vName]['value'] if vName in row.keys() else "&nbsp;" for vName in variableNames]))for row in results["results"]["bindings"]))
    IPython.display.display(IPython.display.HTML(tableCode))
    

def getHierarchyGraph(entityURI, graph=None, direction="ancestors", displayIdentifier=False, highlightEntity=False):
    """Return a graph representing the hierarchy of (in)direct superclasses or subclasses for an entity.

    Keyword arguments:
    entityURI -- the URI for the entity
    graph -- the graph in which the hierarchy is added. A new graph is created if the value is None. (default: None)
    direction -- should the hierarchy concern the ancestors and/or the descendants of the entity. Possible values: "ancestors", "descendants", "both" (default: "ancestors")
    displayIdentifier -- should the nodes also display their URI (default: False)
    highlightEntity -- should the entity be highlighted (default:False)
    """
    if graph is None:
        #graph = Digraph(graph_attr={'rankdir': 'BT'})
        graph = pgv.AGraph(directed=True, rankdir="BT")
    
    entityIdent = entityURI.replace(edamURI, '').replace("edam:", '')
    
    entityType = "Class"
    if entityIdent.startswith("topic_"):
        entityType = "Topic"
    elif entityIdent.startswith("operation_"):
        entityType = "Operation"
    
    conceptStyle = {}
    conceptStyle['Class'] = "filled"
    conceptStyle['Tool'] = "filled"
    conceptStyle['Topic'] = "filled"
    conceptStyle['Operation'] = "rounded,filled"
    
    if entityURI.startswith("http"):
        entityURI = "<" + entityURI + ">"
    
    if (direction == "ancestors") or (direction == "both"):
        query = """
SELECT DISTINCT ?subConcept ?subConceptLabel ?superConcept ?superConceptLabel
WHERE {
  VALUES ?concept { """ + entityURI + """ }
 
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes+query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            startClassIdent = result['subConcept']['value'].replace(edamURI, '')
            startClassLabel = result['subConceptLabel']['value'] + ("\n(" + startClassIdent + ")" if displayIdentifier else "")
            endClassIdent = result['superConcept']['value'].replace(edamURI, '')
            endClassLabel = result['superConceptLabel']['value'] + ("\n(" + endClassIdent + ")" if displayIdentifier else "")
        
            #graph.node(startClassIdent, label=startClassLabel, shape="box")
            #graph.node(endClassIdent, label=endClassLabel, shape="box")
            #graph.edge(startClassIdent, endClassIdent, arrowhead="onormal")
            graph.add_node(startClassIdent, label=startClassLabel, shape="box", color="black", nodeType=entityType, style=conceptStyle[entityType], fillcolor="#ffffff")
            graph.add_node(endClassIdent, label=endClassLabel, shape="box", color="black", nodeType=entityType, style=conceptStyle[entityType], fillcolor="#ffffff")
            graph.add_edge(startClassIdent, endClassIdent, arrowhead="onormal")
    
    if (direction == "descendants") or (direction == "both"):
        query = """
SELECT DISTINCT ?subConcept ?subConceptLabel ?superConcept ?superConceptLabel
WHERE {
  VALUES ?concept { """ + entityURI + """ }
  
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes+query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            startClassIdent = result['subConcept']['value'].replace(edamURI, '')
            startClassLabel = result['subConceptLabel']['value'] + ("\n(" + startClassIdent + ")" if displayIdentifier else "")
            endClassIdent = result['superConcept']['value'].replace(edamURI, '')
            endClassLabel = result['superConceptLabel']['value'] + ("\n(" + endClassIdent + ")" if displayIdentifier else "")
        
            #graph.node(startClassIdent, label=startClassLabel, shape="box")
            #graph.node(endClassIdent, label=endClassLabel, shape="box")
            #graph.edge(startClassIdent, endClassIdent, arrowhead="onormal")
            graph.add_node(startClassIdent, label=startClassLabel, shape="box", color="black", nodeType=entityType, style=conceptStyle[entityType], fillcolor="#ffffff")
            graph.add_node(endClassIdent, label=endClassLabel, shape="box", color="black", nodeType=entityType, style=conceptStyle[entityType], fillcolor="#ffffff")
            graph.add_edge(startClassIdent, endClassIdent, arrowhead="onormal")
    
    if highlightEntity:
        if graph.has_node(entityIdent):
            graph.get_node(entityIdent).attr['color'] = "red"
        else:
            graph.add_node(entityIdent, color="red")
    return graph

def addToolAndAnnotationsToGraph(toolURI, graph=None, showTopics=True, showOperations=True, showDeprecatedAnnotations=False, highlightDirectAnnotations=False):
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
    conceptStyle['Tool'] = "filled"
    conceptStyle['Topic'] = "filled"
    conceptStyle['TopicDeprecated'] = conceptStyle['Topic']+",dashed"
    conceptStyle['TopicAlternative'] = conceptStyle['Topic']+",dotted"
    conceptStyle['Operation'] = "rounded,filled"
    conceptStyle['OperationDeprecated'] = conceptStyle['Operation']+",dashed"
    conceptStyle['OperationAlternative'] = conceptStyle['Operation']+",dotted"
    
    query = """
SELECT DISTINCT ?toolLabel 
WHERE {
  VALUES ?tool { """ + toolURI + """ }

  ?tool rdf:type sc:SoftwareApplication .
  OPTIONAL { ?tool sc:name ?tLabel }
  BIND(COALESCE(?tLabel, "") AS ?toolLabel)
}
"""
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes+query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        #print("{}\t{}".format(toolIdent, result["toolLabel"]["value"]))
        if not graph.has_node(toolIdent):
            clusterTools = graph.get_subgraph(name="cluster_tools")
            if clusterTools is None:
                clusterTools = graph.add_subgraph(name="cluster_tools", rankdir="same", style="invis") # style="invis"
            clusterTools.add_node(toolIdent, label="{}".format(result["toolLabel"]["value"]), shape="ellipse", color="blue", nodeType="Tool", style=conceptStyle["Tool"], fillcolor="#ffffff")
    
    directConcepts = []
    if showTopics:
        directConcepts = []
        conceptType = "Topic"
        query = """
SELECT DISTINCT ?conceptURI ?conceptLabel
WHERE {
  VALUES ?tool { """ + toolURI + """ }

  ?tool sc:applicationSubCategory ?conceptURI .
  ?conceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?conceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?conceptURI rdfs:label ?cLabel }
  BIND(COALESCE(?cLabel, "") AS ?conceptLabel)
}
"""
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes+query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            #print("\t{}\t{}".format(conceptIdent, result["conceptLabel"]["value"]))
            if not graph.has_node(conceptIdent):
                graph.add_node(conceptIdent, label="{}\n({})".format(result["conceptLabel"]["value"], conceptIdent), nodeType=conceptType, shape="box", color=directAnnotationColor, style=conceptStyle[conceptType], fillcolor="#ffffff")
            graph.add_edge(toolIdent, conceptIdent, arrowhead="vee", color="blue", fontcolor="blue", style="dashed")
    
    for conceptURI in directConcepts:
        query = """
SELECT DISTINCT ?subConceptURI ?subConceptLabel ?superConceptURI ?superConceptLabel 
WHERE {
  VALUES ?conceptURI { <""" + conceptURI + """> }
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes+query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        directConcepts = []
        for result in results["results"]["bindings"]:
            subConceptIdent = result["subConceptURI"]["value"].replace(edamURI, "")
            superConceptIdent = result["superConceptURI"]["value"].replace(edamURI, "")
            if not graph.has_node(subConceptIdent):
                graph.add_node(subConceptIdent, label="{}\n{}".format(result["subConceptLabel"]["value"], subConceptIdent), shape="box", color="black", nodeType=conceptType, style=conceptStyle[conceptType], fillcolor="#ffffff")
            if not graph.has_node(superConceptIdent):
                graph.add_node(superConceptIdent, label="{}\n{}".format(result["superConceptLabel"]["value"], superConceptIdent), shape="box", color="black", nodeType=conceptType, style=conceptStyle[conceptType], fillcolor="#ffffff")
            graph.add_edge(subConceptIdent, superConceptIdent, arrowhead="onormal")
            
            
    if showOperations:
        directConcepts = []
        conceptType = "Operation"
        query = """
SELECT DISTINCT ?conceptURI ?conceptLabel
WHERE {
  VALUES ?tool { """ + toolURI + """ }

  ?tool sc:featureList ?conceptURI .
  ?conceptURI rdf:type owl:Class .
  FILTER NOT EXISTS { ?conceptURI rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?conceptURI rdfs:label ?cLabel }
  BIND(COALESCE(?cLabel, "") AS ?conceptLabel)
}
"""
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes+query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            #print("\t{}\t{}".format(conceptIdent, result["conceptLabel"]["value"]))
            if not graph.has_node(conceptIdent):
                graph.add_node(conceptIdent, label="{}\n({})".format(result["conceptLabel"]["value"], conceptIdent), nodeType=conceptType, shape="box", color=directAnnotationColor, style=conceptStyle[conceptType], fillcolor="#ffffff")
            graph.add_edge(toolIdent, conceptIdent, arrowhead="vee", color="blue", fontcolor="blue", style="dashed")
    
    for conceptURI in directConcepts:
        query = """
SELECT DISTINCT ?subConceptURI ?subConceptLabel ?superConceptURI ?superConceptLabel 
WHERE {
  VALUES ?conceptURI { <""" + conceptURI + """> }
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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes+query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        directConcepts = []
        for result in results["results"]["bindings"]:
            subConceptIdent = result["subConceptURI"]["value"].replace(edamURI, "")
            superConceptIdent = result["superConceptURI"]["value"].replace(edamURI, "")
            if not graph.has_node(subConceptIdent):
                graph.add_node(subConceptIdent, label="{}\n{}".format(result["subConceptLabel"]["value"], subConceptIdent), shape="box", color="black", nodeType=conceptType, style=conceptStyle[conceptType], fillcolor="#ffffff")
            if not graph.has_node(superConceptIdent):
                graph.add_node(superConceptIdent, label="{}\n{}".format(result["superConceptLabel"]["value"], superConceptIdent), shape="box", color="black", nodeType=conceptType, style=conceptStyle[conceptType], fillcolor="#ffffff")
            graph.add_edge(subConceptIdent, superConceptIdent, arrowhead="onormal")
    
    if showDeprecatedAnnotations:
        #conceptType = "TopicDeprecated"
        query = """
SELECT DISTINCT ?conceptURI ?conceptLabel ?conceptAlternative ?conceptAlternativeLabel
WHERE {
  VALUES ?tool { """ + toolURI + """ }

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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes+query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            conceptType = "TopicDeprecated"
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            #print("\t{}\t{}".format(conceptIdent, result["conceptLabel"]["value"]))
            if not graph.has_node(conceptIdent):
                graph.add_node(conceptIdent, label="{}\n({})".format(result["conceptLabel"]["value"], conceptIdent), nodeType=conceptType, shape="box", color="grey", style=conceptStyle[conceptType], fillcolor="#ffffff")
            graph.add_edge(toolIdent, conceptIdent, arrowhead="vee", color="grey", fontcolor="grey", style="dashed")
            if "conceptAlternative" in result.keys():
                conceptType = "TopicAlternative"
                alternativeIdent = result["conceptAlternative"]["value"].replace(edamURI, "")
                if not graph.has_node(alternativeIdent):
                    graph.add_node(alternativeIdent, label="{}\n({})".format(result["conceptAlternativeLabel"]["value"], alternativeIdent), nodeType=conceptType, shape="box", color="grey", style=conceptStyle[conceptType], fillcolor="#ffffff")
                    queryHierarchy = """
SELECT DISTINCT ?subConceptURI ?subConceptLabel ?superConceptURI ?superConceptLabel 
WHERE {
  VALUES ?conceptURI { <""" + alternativeURI + """> }
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
                    sparqlHierarchy = SPARQLWrapper(endpointURL)
                    sparqlHierarchy.setQuery(prefixes+queryHierarchy)
                    sparqlHierarchy.setReturnFormat(JSON)
                    resultsHierarchy = sparqlHierarchy.query().convert()
                    for resultHierarchy in resultsHierarchy["results"]["bindings"]:
                        subConceptIdent = resultHierarchy["subConceptURI"]["value"].replace(edamURI, "")
                        superConceptIdent = resultHierarchy["superConceptURI"]["value"].replace(edamURI, "")
                        if not graph.has_node(subConceptIdent):
                            graph.add_node(subConceptIdent, label="{}\n{}".format(resultHierarchy["subConceptLabel"]["value"], subConceptIdent), shape="box", color="grey", nodeType=conceptType, style=conceptStyle[conceptType], fillcolor="#ffffff")
                        if not graph.has_node(superConceptIdent):
                            graph.add_node(superConceptIdent, label="{}\n{}".format(resultHierarchy["superConceptLabel"]["value"], superConceptIdent), shape="box", color="grey", nodeType=conceptType, style=conceptStyle[conceptType], fillcolor="#ffffff")
                        if not graph.has_edge(subConceptIdent, superConceptIdent):
                            graph.add_edge(subConceptIdent, superConceptIdent, arrowhead="onormal", color="grey", style="dotted")
                
                graph.add_edge(conceptIdent, alternativeIdent, arrowhead="vee", color="grey", fontcolor="grey", style="dotted")
    
        
        #conceptType = "OperationDeprecated"
        query = """
SELECT DISTINCT ?conceptURI ?conceptLabel ?conceptAlternative ?conceptAlternativeLabel
WHERE {
  VALUES ?tool { """ + toolURI + """ }

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
        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes+query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            conceptType = "OperationDeprecated"
            directConcepts.append(result["conceptURI"]["value"])
            conceptIdent = result["conceptURI"]["value"].replace(edamURI, "")
            #print("\t{}\t{}".format(conceptIdent, result["conceptLabel"]["value"]))
            if not graph.has_node(conceptIdent):
                graph.add_node(conceptIdent, label="{}\n({})".format(result["conceptLabel"]["value"], conceptIdent), nodeType=conceptType, shape="box", color="grey", style=conceptStyle[conceptType], fillcolor="#ffffff")
            graph.add_edge(toolIdent, conceptIdent, arrowhead="vee", color="grey", fontcolor="grey", style="dashed")
            if "conceptAlternative" in result.keys():
                conceptType = "OperationAlternative"
                alternativeURI = result["conceptAlternative"]["value"]
                alternativeIdent = alternativeURI.replace(edamURI, "")
                if not graph.has_node(alternativeIdent):
                    graph.add_node(alternativeIdent, label="{}\n({})".format(result["conceptAlternativeLabel"]["value"], alternativeIdent), nodeType=conceptType, shape="box", color="grey", style=conceptStyle[conceptType], fillcolor="#ffffff")
                    queryHierarchy = """
SELECT DISTINCT ?subConceptURI ?subConceptLabel ?superConceptURI ?superConceptLabel 
WHERE {
  VALUES ?conceptURI { <""" + alternativeURI + """> }
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
                    sparqlHierarchy = SPARQLWrapper(endpointURL)
                    sparqlHierarchy.setQuery(prefixes+queryHierarchy)
                    sparqlHierarchy.setReturnFormat(JSON)
                    resultsHierarchy = sparqlHierarchy.query().convert()
                    for resultHierarchy in resultsHierarchy["results"]["bindings"]:
                        subConceptIdent = resultHierarchy["subConceptURI"]["value"].replace(edamURI, "")
                        superConceptIdent = resultHierarchy["superConceptURI"]["value"].replace(edamURI, "")
                        if not graph.has_node(subConceptIdent):
                            graph.add_node(subConceptIdent, label="{}\n{}".format(resultHierarchy["subConceptLabel"]["value"], subConceptIdent), shape="box", color="grey", nodeType=conceptType, style=conceptStyle[conceptType], fillcolor="#ffffff")
                        if not graph.has_node(superConceptIdent):
                            graph.add_node(superConceptIdent, label="{}\n{}".format(resultHierarchy["superConceptLabel"]["value"], superConceptIdent), shape="box", color="grey", nodeType=conceptType, style=conceptStyle[conceptType], fillcolor="#ffffff")
                        if not graph.has_edge(subConceptIdent, superConceptIdent):
                            graph.add_edge(subConceptIdent, superConceptIdent, arrowhead="onormal", color="grey", style="dotted")
                graph.add_edge(conceptIdent, alternativeIdent, arrowhead="vee", color="grey", fontcolor="grey", style="dotted")
    
    return graph

def getToolLabel(toolURI):
    """Return the label of a tool (or empty string if no label is present).

    Keyword arguments:
    toolURI -- the URI for the tool
    """
    
    if toolURI.startswith("http"):
        toolURI = "<" + toolURI + ">"
    query = """
SELECT DISTINCT ?tool ?toolLabel
WHERE {
  VALUES ?tool { """ + toolURI + """ }

  ?tool rdf:type sc:SoftwareApplication .
  OPTIONAL { ?tool sc:name ?tLabel }
  BIND(COALESCE(?tLabel, "") AS ?toolLabel)
}
"""
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes+query)
    sparql.setReturnFormat(JSON)
    results = sparql.queryAndConvert()
    return results["results"]["bindings"][0]["toolLabel"]["value"]

def getToolURIByLabel(toolLabel):
    """Return the URI of a tool designated by its label (or None).

    Keyword arguments:
    toolLabel -- the label for the tool
    """
    
    query = """
SELECT DISTINCT ?tool ?toolLabel
WHERE {
  VALUES ?toolLabel { \"""" + toolLabel + """\" }

  ?tool rdf:type sc:SoftwareApplication .
  ?tool sc:name ?toolLabel .
}
"""
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes+query)
    sparql.setReturnFormat(JSON)
    results = sparql.queryAndConvert()
    return None if len(results["results"]["bindings"]) == 0 else results["results"]["bindings"][0]["tool"]["value"]

def getToolTopics(toolURI, transitive=False):
    """Return the list of the (URI, label) tuples for the topics associated to a tool.

    Keyword arguments:
    toolURI -- the URI for the tool
    transitive -- also consider the ancestors of the topics directly associated to the tool (default: False)
    """
    
    if toolURI.startswith("http"):
        toolURI = "<" + toolURI + ">"
    transitiveClause = "/(rdfs:subClassOf*)" if transitive else ""
    query = """
SELECT DISTINCT ?tool ?topic ?topicLabel
WHERE {
  VALUES ?tool { """ + toolURI + """ }

  ?tool sc:applicationSubCategory""" + transitiveClause + """ ?topic .
  ?topic rdf:type owl:Class .
  FILTER NOT EXISTS { ?topic rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?topic rdfs:label ?tLabel }
  BIND(COALESCE(?tLabel, "") AS ?topicLabel)
}
"""
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes+query)
    sparql.setReturnFormat(JSON)
    results = sparql.queryAndConvert()
    toolTopics = [(result["topic"]["value"], result["topicLabel"]["value"]) for result in results["results"]["bindings"]]
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
    query = """
SELECT DISTINCT ?tool ?operation ?operationLabel
WHERE {
  VALUES ?tool { """ + toolURI + """ }

  ?tool sc:featureList""" + transitiveClause + """ ?operation .
  ?operation rdf:type owl:Class .
  FILTER NOT EXISTS { ?operation rdfs:subClassOf? owl:DeprecatedClass }
  OPTIONAL { ?operation rdfs:label ?oLabel }
  BIND(COALESCE(?oLabel, "") AS ?operationLabel)
}
"""
    sparql = SPARQLWrapper(endpointURL)
    sparql.setQuery(prefixes+query)
    sparql.setReturnFormat(JSON)
    results = sparql.queryAndConvert()
    toolOperations = [(result["operation"]["value"], result["operationLabel"]["value"]) for result in results["results"]["bindings"]]
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


def addToolsAndAnnotationsToGraph(listToolURI, graph=None, showTopics=True, showOperations=True, highlightDirectAnnotations=False, highlightIntersection=False):
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
            highlightDirectAnnotations=highlightDirectAnnotations
        )
    
    if highlightIntersection:
        commonConcepts = []
        if showTopics:
            commonConcepts += getToolsCommonTopics(listToolURI, transitive=True)
        if showOperations:
            commonConcepts += getToolsCommonOperations(listToolURI, transitive=True)
    
        for (currentConcept, currentLabel) in commonConcepts:
            currentConceptIdent = currentConcept.replace(edamURI, "")
            node = graph.get_node(currentConceptIdent)
            node.attr['color'] = "red"
            node.attr['penwidth'] = "3"              # Thicker border
            node.attr['style'] = "filled,bold"         # Bold border and filled style
            node.attr['fillcolor'] = "red"             # Red fill to highlight
            node.attr['highlightIntersection'] = "True"  # Mark this node as highlighted
    
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

def colorGraphNodesAccordingToScore(graph, dictTopicScore, dictOperationScore, color="red"):
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
            currentNodeType = currentNode.attr['nodeType']
            #currentNode.attr['style'] = conceptStyle[currentNodeType]
            if currentNodeType == "Topic":
                if currentNode in dictTopicScore.keys():
                    graph.get_node(currentNode).attr['fillcolor'] = getScoreColorRGB(dictTopicScore[currentNode], topicMaxValue, color=color)
            elif currentNodeType == "Operation":
                if currentNode in dictOperationScore.keys():
                    graph.get_node(currentNode).attr['fillcolor'] = getScoreColorRGB(dictOperationScore[currentNode], operationMaxValue, color=color)


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


def getMutualInformation(toolURI1, toolURI2, dictAnnotationPairsMutualInformation=None):
    mi = 0.
    annotT1 = dfToolTopicTransitive[dfToolTopicTransitive['tool'] == toolURI1]['topic'].to_list()
    annotT1 += dfToolOperationTransitive[dfToolOperationTransitive['tool'] == toolURI1]['operation'].to_list()
    annotT2 = dfToolTopicTransitive[dfToolTopicTransitive['tool'] == toolURI2]['topic'].to_list()
    annotT2 += dfToolOperationTransitive[dfToolOperationTransitive['tool'] == toolURI2]['operation'].to_list()
    for a1 in annotT1:
        toolsA1 = set(dfToolTopicTransitive[dfToolTopicTransitive['topic'] == a1]['tool'])
        pa1 = len(toolsA1) / nbToolsWithTopic
        
        for a2 in annotT2:
            if dictAnnotationPairsMutualInformation is not None:
                if (a1, a2) in dictAnnotationPairsMutualInformation.keys():
                    mi += dictAnnotationPairsMutualInformation[(a1, a2)]
                    continue
            toolsA2 = set(dfToolTopicTransitive[dfToolTopicTransitive['topic'] == a2]['tool'])
            pa2 = len(toolsA2) / nbToolsWithTopic
            pa1a2 = len(toolsA1.intersection(toolsA2)) / nbToolsWithTopic

            #if  pa1 * pa2 * pa1a2 != 0:
            #    mi += pa1a2 * np.log2( pa1a2 / (pa1 * pa2) )
            result = pa1a2 * np.log2( pa1a2 / (pa1 * pa2) ) if  pa1 * pa2 * pa1a2 != 0 else 0
            mi += result
            if dictAnnotationPairsMutualInformation is not None:
                dictAnnotationPairsMutualInformation[(a1, a2)] = result
                dictAnnotationPairsMutualInformation[(a2, a1)] = result
    return mi