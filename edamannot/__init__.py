"""
EDAMannot package.

The package-level API is intentionally lazy: importing ``edamannot`` does not
immediately import every submodule (and therefore does not immediately load
optional/heavier dependencies such as Graphviz or SPARQL helpers).
"""

from importlib import import_module
from typing import Final


# Public names exposed at package level.  The mapping is explicit so that the
# dependency direction remains visible and the package does not need to import
# every module eagerly.
_EXPORTS: Final[dict[str, tuple[str, str]]] = {
    # config.py
    "current_dir": ("config", "current_dir"),
    "neighbor_dir": ("config", "neighbor_dir"),
    "notebooks_dir": ("config", "notebooks_dir"),
    "bioschemas_file": ("config", "bioschemas_file"),
    "edam_file": ("config", "edam_file"),
    "edam_neighbors_file": ("config", "edam_neighbors_file"),
    "endpointURL": ("config", "endpointURL"),
    "rdfFormat": ("config", "rdfFormat"),
    "prefixes": ("config", "prefixes"),
    "biotoolsURI": ("config", "biotoolsURI"),
    "biotoolsOntologyURI": ("config", "biotoolsOntologyURI"),
    "edamURI": ("config", "edamURI"),
    "set_file_paths": ("config", "set_file_paths"),

    # sparql.py
    "create_sparql_client": ("sparql", "create_sparql_client"),
    "execute_query": ("sparql", "execute_query"),
    "execute_query_and_convert": ("sparql", "execute_query_and_convert"),
    "query_dataframe": ("sparql", "query_dataframe"),
    "sparql_results_to_dataframe": ("sparql", "sparql_results_to_dataframe"),
    "displaySparqlResults": ("sparql", "displaySparqlResults"),
    "get_edam_version": ("sparql", "get_edam_version"),

    # dataframes.py
    "get_nb_tools": ("dataframes", "get_nb_tools"),
    "get_tools_dataframe": ("dataframes", "get_tools_dataframe"),
    "get_tools_topics_dataframe": ("dataframes", "get_tools_topics_dataframe"),
    "get_tools_topics_transitive_dataframe": (
        "dataframes",
        "get_tools_topics_transitive_dataframe",
    ),
    "get_tools_operations_label_dataframe": (
        "dataframes",
        "get_tools_operations_label_dataframe",
    ),
    "get_tools_operations_transitive_dataframe": (
        "dataframes",
        "get_tools_operations_transitive_dataframe",
    ),
    "get_dftools_with_nbTopics_nbOperations": (
        "dataframes",
        "get_dftools_with_nbTopics_nbOperations",
    ),
    "generate_df_redundancy_topic": (
        "dataframes",
        "generate_df_redundancy_topic",
    ),
    "generate_df_redundancy_operation": (
        "dataframes",
        "generate_df_redundancy_operation",
    ),
    "generate_df_topic_no_redundancy": (
        "dataframes",
        "generate_df_topic_no_redundancy",
    ),
    "generate_df_operation_no_redundancy": (
        "dataframes",
        "generate_df_operation_no_redundancy",
    ),
    "generate_dfTool_transitive": ("dataframes", "generate_dfTool_transitive"),
    "generate_dfTool_no_transitive": (
        "dataframes",
        "generate_dfTool_no_transitive",
    ),
    "generate_dfTool_no_transitive_no_redundancy": (
        "dataframes",
        "generate_dfTool_no_transitive_no_redundancy",
    ),
    "get_dfDeprecatedItems": ("dataframes", "get_dfDeprecatedItems"),
    "get_dfDeprecatedSuggestedItems": (
        "dataframes",
        "get_dfDeprecatedSuggestedItems",
    ),
    "get_dfToolsWithSomeDeprecatedTopic": (
        "dataframes",
        "get_dfToolsWithSomeDeprecatedTopic",
    ),
    "get_dfToolsWithSomeDeprecatedOperation": (
        "dataframes",
        "get_dfToolsWithSomeDeprecatedOperation",
    ),

    # metrics.py
    "getMutualInformation": ("metrics", "getMutualInformation"),
    "compute_topic_metrics": ("metrics", "compute_topic_metrics"),
    "compute_topic_metrics_NT": ("metrics", "compute_topic_metrics_NT"),
    "compute_operation_metrics": ("metrics", "compute_operation_metrics"),
    "compute_operation_metrics_NT": ("metrics", "compute_operation_metrics_NT"),
    "compute_tool_metrics_with_transitive": (
        "metrics",
        "compute_tool_metrics_with_transitive",
    ),
    "compute_tool_metrics_non_transitive": (
        "metrics",
        "compute_tool_metrics_non_transitive",
    ),
    "get_tool_metrics": ("metrics", "get_tool_metrics"),

    # annotations.py
    "get_tool_url": ("annotations", "get_tool_url"),
    "normalize_tool_input": ("annotations", "normalize_tool_input"),
    "getToolLabel": ("annotations", "getToolLabel"),
    "getToolURIByLabel": ("annotations", "getToolURIByLabel"),
    "getToolTopics": ("annotations", "getToolTopics"),
    "getToolOperations": ("annotations", "getToolOperations"),
    "getToolsCommonTopics": ("annotations", "getToolsCommonTopics"),
    "getToolsCommonOperations": ("annotations", "getToolsCommonOperations"),
    "fetch_annotations": ("annotations", "fetch_annotations"),
    "format_tool_annotations": ("annotations", "format_tool_annotations"),
    "to_json": ("annotations", "to_json"),
    "fetch_annotations_with_metrics": (
        "annotations",
        "fetch_annotations_with_metrics",
    ),

    # inference.py
    "infer_edam_neighbors": ("inference", "infer_edam_neighbors"),
    "edam_uri_to_id": ("inference", "edam_uri_to_id"),
    "merge_annotations": ("inference", "merge_annotations"),
    "infer_neighbors_from_annotations": (
        "inference",
        "infer_neighbors_from_annotations",
    ),

    # graphs.py
    "getHierarchyGraph": ("graphs", "getHierarchyGraph"),
    "get_edam_neighbors_dataframe": ("graphs", "get_edam_neighbors_dataframe"),
    "get_edam_chained_neighbors_dataframe": (
        "graphs",
        "get_edam_chained_neighbors_dataframe",
    ),
    "getEntityDescriptionGraph": ("graphs", "getEntityDescriptionGraph"),
    "addToolAndAnnotationsToGraph": (
        "graphs",
        "addToolAndAnnotationsToGraph",
    ),
    "addToolsAndAnnotationsToGraph": (
        "graphs",
        "addToolsAndAnnotationsToGraph",
    ),
    "getScoreColorRGB": ("graphs", "getScoreColorRGB"),
    "dictTopic": ("graphs", "dictTopic"),
    "dictOperation": ("graphs", "dictOperation"),
    "buildTopicOperationDicts": ("graphs", "buildTopicOperationDicts"),
    "colorGraphNodesAccordingToScore": (
        "graphs",
        "colorGraphNodesAccordingToScore",
    ),
    "getToolScore": ("graphs", "getToolScore"),
}


__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    """Lazily resolve a public package attribute from its owning module."""
    try:
        module_name, attribute_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module 'edamannot' has no attribute {name!r}") from exc

    module = import_module(f"{__name__}.{module_name}")
    value = getattr(module, attribute_name)

    # Cache the resolved attribute on the package so subsequent access does not
    # repeat the import/lookup work.
    globals()[name] = value
    return value
