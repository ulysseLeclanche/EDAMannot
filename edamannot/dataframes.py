"""
This module contains the DataFrame-oriented part of the original EDAMannot.py:
- tool counts and base tool DataFrames;
- direct and transitive topic/operation DataFrames;
- topic/operation redundancy detection;
- non-redundant annotation DataFrames;
- tool count aggregates;
- deprecated-item diagnostic DataFrames.

SPARQL access is delegated to edamannot.sparql.
"""

import pandas as pd

from .config import edamURI
from .sparql import execute_query, query_dataframe


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
      FILTER(STRSTARTS(STR(?tool), "https://bio.tools/"))
    }
    """

    results = execute_query(query)

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
      FILTER(STRSTARTS(STR(?tool), "https://bio.tools/"))
    }
    """

    dfTool = query_dataframe(query)
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
      FILTER(STRSTARTS(STR(?tool), "https://bio.tools/"))
      FILTER NOT EXISTS { ?topic rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?topic rdfs:label ?tLabel }
      BIND(COALESCE(?tLabel, "") AS ?topicLabel)
    }
    """

    dfToolTopic = query_dataframe(query)
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
      FILTER(STRSTARTS(STR(?tool), "https://bio.tools/"))
      FILTER NOT EXISTS { ?topic rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?topic rdfs:label ?tLabel }
      BIND(COALESCE(?tLabel, "") AS ?topicLabel)
    }
    """

    dfToolTopicTransitive = query_dataframe(query)
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
      FILTER(STRSTARTS(STR(?tool), "https://bio.tools/"))
      FILTER NOT EXISTS { ?operation rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?operation rdfs:label ?oLabel }
      BIND(COALESCE(?oLabel, "") AS ?operationLabel)
    }
    """

    dfToolOperation = query_dataframe(query)
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
      FILTER(STRSTARTS(STR(?tool), "https://bio.tools/"))
      FILTER NOT EXISTS { ?operation rdfs:subClassOf? owl:DeprecatedClass }
      OPTIONAL { ?operation rdfs:label ?oLabel }
      BIND(COALESCE(?oLabel, "") AS ?operationLabel)
    }
    """

    dfToolOperationTransitive = query_dataframe(query)
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
    results = execute_query(redundancyQuery)

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
    results = execute_query(operationRedundancyQuery)

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


def get_dfDeprecatedItems(
    output_path: str = "Dataframe/dfDeprecatedItems.tsv.bz2",
) -> pd.DataFrame:
    """
    Retrieve all items marked as owl:deprecated and save as compressed TSV.
    """
    query = """
    SELECT DISTINCT ?deprecatedItem
    WHERE {
      ?deprecatedItem owl:deprecated true .
    }
    ORDER BY ?deprecatedItem
    """

    results = execute_query(query)

    data = [
        {"Deprecated Item": r["deprecatedItem"]["value"]}
        for r in results["results"]["bindings"]
    ]
    dfDeprecatedItems = pd.DataFrame(data)

    dfDeprecatedItems.to_csv(output_path, sep="\t", index=False, compression="bz2")
    return dfDeprecatedItems


def get_dfDeprecatedSuggestedItems(
    output_path: str = "Dataframe/dfDeprecatedSuggestedItems.tsv.bz2",
) -> pd.DataFrame:
    """
    Retrieve deprecated items along with their suggested replacements.
    """
    query = """
    SELECT DISTINCT ?deprecatedItem ?suggestedItem
    WHERE {
      {
        ?deprecatedItem rdfs:subClassOf owl:DeprecatedClass .
      } UNION {
        ?deprecatedItem owl:deprecated true .
      } UNION {
        ?deprecatedItem owl:deprecated "true" .
      } UNION {
        ?deprecatedItem owl:deprecated "True" .
      }
      ?deprecatedItem oboInOwl:consider ?suggestedItem .
    }
    ORDER BY ?deprecatedItem
    """

    results = execute_query(query)

    data = [
        {
            "Deprecated Item": r["deprecatedItem"]["value"],
            "Suggested Item": r["suggestedItem"]["value"],
        }
        for r in results["results"]["bindings"]
    ]
    dfDeprecatedSuggestedItems = pd.DataFrame(data)

    dfDeprecatedSuggestedItems.to_csv(
        output_path, sep="\t", index=False, compression="bz2"
    )
    return dfDeprecatedSuggestedItems


def get_dfToolsWithSomeDeprecatedTopic(
    output_path: str = "Dataframe/dfToolsWithSomeDeprecatedTopic.tsv.bz2",
) -> pd.DataFrame:
    """
    Retrieve the list of tools annotated with at least one deprecated topic.
    """
    query = """
    SELECT DISTINCT ?tool
    WHERE {
      ?tool rdf:type sc:SoftwareApplication .
      ?tool sc:applicationSubCategory/(rdfs:subClassOf*) ?deprecatedTopic .
      ?deprecatedTopic owl:deprecated True .
    }
    ORDER BY ?tool
    """

    results = execute_query(query)

    data = [{"Tool": r["tool"]["value"]} for r in results["results"]["bindings"]]
    dfToolsWithSomeDeprecatedTopic = pd.DataFrame(data)

    dfToolsWithSomeDeprecatedTopic.to_csv(
        output_path, sep="\t", index=False, compression="bz2"
    )
    return dfToolsWithSomeDeprecatedTopic


def get_dfToolsWithSomeDeprecatedOperation(
    output_path: str = "Dataframe/dfToolsWithSomeDeprecatedOperation.tsv.bz2",
) -> pd.DataFrame:
    """
    Retrieve the list of tools annotated with at least one deprecated operation.
    """
    query = """
    SELECT DISTINCT ?tool
    WHERE {
      ?tool rdf:type sc:SoftwareApplication .
      ?tool sc:featureList/(rdfs:subClassOf*) ?deprecatedItem .
      ?deprecatedItem owl:deprecated True .
    }
    ORDER BY ?tool
    """

    results = execute_query(query)

    data = [{"Tool": r["tool"]["value"]} for r in results["results"]["bindings"]]
    dfToolsWithSomeDeprecatedOperation = pd.DataFrame(data)

    dfToolsWithSomeDeprecatedOperation.to_csv(
        output_path, sep="\t", index=False, compression="bz2"
    )
    return dfToolsWithSomeDeprecatedOperation




__all__ = [
    "get_nb_tools",
    "get_tools_dataframe",
    "get_tools_topics_dataframe",
    "get_tools_topics_transitive_dataframe",
    "get_tools_operations_label_dataframe",
    "get_tools_operations_transitive_dataframe",
    "get_dftools_with_nbTopics_nbOperations",
    "generate_df_redundancy_topic",
    "generate_df_redundancy_operation",
    "generate_df_topic_no_redundancy",
    "generate_df_operation_no_redundancy",
    "generate_dfTool_transitive",
    "generate_dfTool_no_transitive",
    "generate_dfTool_no_transitive_no_redundancy",
    "get_dfDeprecatedItems",
    "get_dfDeprecatedSuggestedItems",
    "get_dfToolsWithSomeDeprecatedTopic",
    "get_dfToolsWithSomeDeprecatedOperation",
]
