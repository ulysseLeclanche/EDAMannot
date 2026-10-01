"""
SPARQL utilities for EDAMannot.

This module centralizes the low-level SPARQL access that was previously
implemented repeatedly throughout EDAMannot.py.

Responsibilities
----------------
- Create and configure SPARQLWrapper clients.
- Execute SPARQL queries and return JSON results.
- Execute queries using SPARQLWrapper.queryAndConvert().
- Execute SPARQL queries through sparqldataframe.
- Convert SPARQL JSON results to pandas DataFrames.
- Display SPARQL results in a Jupyter/IPython environment.
- Retrieve the EDAM ontology version.

The historical public helper names are intentionally preserved so the
refactoring can be introduced incrementally without changing the behavior of
the current EDAMannot.py API.
"""

from typing import Any, Optional

import IPython
import pandas as pd
import sparqldataframe
from SPARQLWrapper import JSON, SPARQLWrapper

from .config import endpointURL, prefixes


# ---------------------------------------------------------------------------
# LOW-LEVEL SPARQL CLIENT
# ---------------------------------------------------------------------------

def create_sparql_client(
    endpoint_url: Optional[str] = None,
) -> SPARQLWrapper:
    """
    Create a configured SPARQLWrapper client.

    Parameters
    ----------
    endpoint_url:
        SPARQL endpoint to use. When omitted, the endpoint configured in
        edamannot.config.endpointURL is used.

    Returns
    -------
    SPARQLWrapper
        Configured SPARQLWrapper client.
    """
    return SPARQLWrapper(endpoint_url or endpointURL)


# ---------------------------------------------------------------------------
# QUERY EXECUTION
# ---------------------------------------------------------------------------

def execute_query(
    query: str,
    endpoint_url: Optional[str] = None,
    query_prefixes: Optional[str] = None,
) -> dict[str, Any]:
    """
    Execute a SPARQL SELECT query and return converted JSON results.

    This replaces the repeated pattern found throughout EDAMannot.py:

        sparql = SPARQLWrapper(endpointURL)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()

    Parameters
    ----------
    query:
        SPARQL query body, without prefixes.
    endpoint_url:
        Optional SPARQL endpoint override.
    query_prefixes:
        Optional prefix block. Defaults to config.prefixes.

    Returns
    -------
    dict
        SPARQL JSON results.
    """
    sparql = create_sparql_client(endpoint_url)
    sparql.setQuery((query_prefixes if query_prefixes is not None else prefixes) + query)
    sparql.setReturnFormat(JSON)
    return sparql.query().convert()


def execute_query_and_convert(
    query: str,
    endpoint_url: Optional[str] = None,
    query_prefixes: Optional[str] = None,
) -> dict[str, Any]:
    """
    Execute a SPARQL query using SPARQLWrapper.queryAndConvert().

    This helper preserves the behavior of the original functions that used
    queryAndConvert(), such as getToolLabel(), getToolURIByLabel(),
    getToolTopics(), and getToolOperations().

    Parameters
    ----------
    query:
        SPARQL query body, without prefixes.
    endpoint_url:
        Optional SPARQL endpoint override.
    query_prefixes:
        Optional prefix block. Defaults to config.prefixes.

    Returns
    -------
    dict
        SPARQL JSON results.
    """
    sparql = create_sparql_client(endpoint_url)
    sparql.setQuery((query_prefixes if query_prefixes is not None else prefixes) + query)
    sparql.setReturnFormat(JSON)
    return sparql.queryAndConvert()


def query_dataframe(
    query: str,
    endpoint_url: Optional[str] = None,
    query_prefixes: Optional[str] = None,
) -> pd.DataFrame:
    """
    Execute a SPARQL query through sparqldataframe.

    This preserves the mechanism currently used by the main DataFrame
    extraction functions in EDAMannot.py.

    Parameters
    ----------
    query:
        SPARQL query body, without prefixes.
    endpoint_url:
        Optional SPARQL endpoint override.
    query_prefixes:
        Optional prefix block. Defaults to config.prefixes.

    Returns
    -------
    pandas.DataFrame
        Query results as a DataFrame.
    """
    endpoint = endpoint_url or endpointURL
    prefix_block = query_prefixes if query_prefixes is not None else prefixes

    return sparqldataframe.query(
        endpoint,
        prefix_block + query,
    )


# ---------------------------------------------------------------------------
# RESULT CONVERSION
# ---------------------------------------------------------------------------

def sparql_results_to_dataframe(results: dict[str, Any]) -> pd.DataFrame:
    """
    Convert SPARQL JSON results to a pandas DataFrame.

    This is the original EDAMannot helper extracted without changing its
    output structure.

    Parameters
    ----------
    results:
        Result of SPARQLWrapper.query().convert().

    Returns
    -------
    pandas.DataFrame
        DataFrame whose columns are the SPARQL result variables.
    """
    variable_names = results["head"]["vars"]
    bindings = results["results"]["bindings"]

    rows = []

    for row in bindings:
        row_data = {}

        for var_name in variable_names:
            if var_name in row:
                row_data[var_name] = row[var_name]["value"]
            else:
                row_data[var_name] = None

        rows.append(row_data)

    return pd.DataFrame(rows, columns=variable_names)


# ---------------------------------------------------------------------------
# JUPYTER / IPYTHON DISPLAY
# ---------------------------------------------------------------------------

def displaySparqlResults(results: dict[str, Any]) -> None:
    """
    Display SPARQL JSON results as an HTML table in Jupyter/IPython.

    This keeps the historical function name used by EDAMannot.py.

    Parameters
    ----------
    results:
        Result of SPARQLWrapper.query().convert().
    """
    variable_names = results["head"]["vars"]

    table_code = "<table><tr><th>{}</th></tr><tr>{}</tr></table>".format(
        "</th><th>".join(variable_names),
        "</tr><tr>".join(
            "<td>{}</td>".format(
                "</td><td>".join(
                    [
                        row[var_name]["value"]
                        if var_name in row
                        else "&nbsp;"
                        for var_name in variable_names
                    ]
                )
            )
            for row in results["results"]["bindings"]
        ),
    )

    IPython.display.display(IPython.display.HTML(table_code))


# ---------------------------------------------------------------------------
# EDAM ONTOLOGY INFORMATION
# ---------------------------------------------------------------------------

def get_edam_version(
    endpoint_url: Optional[str] = None,
    query_prefixes: Optional[str] = None,
) -> float:
    """
    Retrieve the EDAM ontology version exposed by the SPARQL endpoint.

    The original EDAMannot.py function accepted endpointURL and prefixes as
    positional arguments. They are now optional so existing explicit calls
    remain possible while new code can simply call get_edam_version().

    Parameters
    ----------
    endpoint_url:
        Optional SPARQL endpoint override.
    query_prefixes:
        Optional prefix block override.

    Returns
    -------
    float
        EDAM ontology version number.

    Raises
    ------
    ValueError
        If the endpoint does not return an EDAM ontology version.
    """
    query = """
    SELECT ?ontology ?versionIRI
           (REPLACE(
               STR(?versionIRI),
               'http://edamontology.org/',
               ''
           ) AS ?versionNumber)
    WHERE {
      ?ontology rdf:type owl:Ontology .
      ?ontology owl:versionIRI ?versionIRI .
    }
    """

    results = execute_query(
        query,
        endpoint_url=endpoint_url,
        query_prefixes=query_prefixes,
    )

    bindings = results["results"]["bindings"]

    if not bindings:
        raise ValueError("No EDAM ontology version was returned by the SPARQL endpoint.")

    version_str = bindings[0]["versionNumber"]["value"]
    return float(version_str)


__all__ = [
    "create_sparql_client",
    "execute_query",
    "execute_query_and_convert",
    "query_dataframe",
    "sparql_results_to_dataframe",
    "displaySparqlResults",
    "get_edam_version",
]
