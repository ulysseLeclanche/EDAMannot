"""
Responsibilities
----------------
- Normalize tool identifiers and bio.tools URLs.
- Resolve Topic/Operation annotation aliases.
- Retrieve tool labels and URIs through Fuseki.
- Retrieve direct/transitive Topic and Operation annotations.
- Compute annotations common to several tools.
- Build the structured annotation JSON representation.
- Combine annotations with tool metrics.

"""

import json
from typing import Dict

import pandas as pd

from .config import biotoolsURI, current_dir
from .metrics import get_tool_metrics
from .sparql import execute_query_and_convert


# ---------------------------------------------------------------------------
# TOOL INPUT NORMALIZATION
# ---------------------------------------------------------------------------

def get_tool_url(tool_name: str) -> str:
    """
    Return a canonical bio.tools URL for a tool name or URL.
    """
    if tool_name.startswith("https://bio.tools/"):
        return tool_name

    return f"{biotoolsURI}{tool_name}"


def normalize_tool_input(tools) -> list[str]:
    """
    Normalize a single tool or a sequence of tools to bio.tools URLs.
    """
    if isinstance(tools, str):
        tools = [tools]

    return [get_tool_url(tool) for tool in tools]


# ---------------------------------------------------------------------------
# ANNOTATION TYPE NORMALIZATION
# ---------------------------------------------------------------------------

def _resolve_annotation_type(value: str) -> str:
    """
    Resolve an annotation type alias to its canonical name.

    Supported values:
        T, Topic
        O, Operation
    """
    mapping = {
        "T": "Topic",
        "O": "Operation",
        "Topic": "Topic",
        "Operation": "Operation",
    }

    value = value.capitalize() if len(value) > 1 else value.upper()

    if value not in mapping:
        raise ValueError(f"Invalid annotation type: {value}")

    return mapping[value]


# ---------------------------------------------------------------------------
# TOOL INFORMATION
# ---------------------------------------------------------------------------

def getToolLabel(toolURI):
    """
    Return the label of a tool.

    The behavior is kept compatible with the original implementation.
    """
    if toolURI.startswith("http"):
        toolURI = f"<{toolURI}>"

    query = """
SELECT DISTINCT ?tool ?toolLabel
WHERE {
  VALUES ?tool { """ + toolURI + """ }

  ?tool rdf:type sc:SoftwareApplication .
  OPTIONAL { ?tool sc:name ?tLabel }
  BIND(COALESCE(?tLabel, "") AS ?toolLabel)
}
"""

    results = execute_query_and_convert(query)

    return results["results"]["bindings"][0]["toolLabel"]["value"]


def getToolURIByLabel(toolLabel):
    """
    Return the URI of a tool designated by its label, or None.
    """
    query = """
SELECT DISTINCT ?tool ?toolLabel
WHERE {
  VALUES ?toolLabel { """" + toolLabel + """" }

  ?tool rdf:type sc:SoftwareApplication .
  ?tool sc:name ?toolLabel .
}
"""

    results = execute_query_and_convert(query)
    bindings = results["results"]["bindings"]

    if len(bindings) == 0:
        return None

    return bindings[0]["tool"]["value"]


# ---------------------------------------------------------------------------
# TOOL ANNOTATIONS
# ---------------------------------------------------------------------------

def getToolTopics(toolURI, transitive=False):
    """
    Return (URI, label) tuples for Topic annotations associated with a tool.

    Parameters
    ----------
    toolURI:
        Full bio.tools URI or an HTTP URI accepted by the endpoint.
    transitive:
        Include ancestors through rdfs:subClassOf* when True.
    """
    if toolURI.startswith("http"):
        toolURI = f"<{toolURI}>"

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

    results = execute_query_and_convert(query)

    return [
        (
            result["topic"]["value"],
            result["topicLabel"]["value"],
        )
        for result in results["results"]["bindings"]
    ]


def getToolOperations(toolURI, transitive=False):
    """
    Return (URI, label) tuples for Operation annotations associated with a tool.

    Parameters
    ----------
    toolURI:
        Full bio.tools URI or an HTTP URI accepted by the endpoint.
    transitive:
        Include ancestors through rdfs:subClassOf* when True.
    """
    if toolURI.startswith("http"):
        toolURI = f"<{toolURI}>"

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

    results = execute_query_and_convert(query)

    return [
        (
            result["operation"]["value"],
            result["operationLabel"]["value"],
        )
        for result in results["results"]["bindings"]
    ]


# ---------------------------------------------------------------------------
# COMMON ANNOTATIONS
# ---------------------------------------------------------------------------

def getToolsCommonTopics(listToolURI, transitive=False):
    """
    Return Topic annotations shared by every tool in the list.
    """
    commonConcepts = []

    if len(listToolURI) > 0:
        commonConcepts = set(
            getToolTopics(
                listToolURI[0],
                transitive=transitive,
            )
        )

    for toolURI in listToolURI[1:]:
        currentConcepts = set(
            getToolTopics(
                toolURI,
                transitive=transitive,
            )
        )
        commonConcepts = commonConcepts.intersection(currentConcepts)

    return list(commonConcepts)


def getToolsCommonOperations(listToolURI, transitive=False):
    """
    Return Operation annotations shared by every tool in the list.
    """
    commonConcepts = []

    if len(listToolURI) > 0:
        commonConcepts = set(
            getToolOperations(
                listToolURI[0],
                transitive=transitive,
            )
        )

    for toolURI in listToolURI[1:]:
        currentConcepts = set(
            getToolOperations(
                toolURI,
                transitive=transitive,
            )
        )
        commonConcepts = commonConcepts.intersection(currentConcepts)

    return list(commonConcepts)


# ---------------------------------------------------------------------------
# DATAFRAME-BACKED ANNOTATION API
# ---------------------------------------------------------------------------

_ANNOTATION_DATAFRAME_PATHS = {
    "Topic": {
        False: "Dataframe/dfToolTopic.tsv.bz2",
        True: "Dataframe/dfToolTopicTransitive.tsv.bz2",
    },
    "Operation": {
        False: "Dataframe/dfToolOperation.tsv.bz2",
        True: "Dataframe/dfToolOperationTransitive.tsv.bz2",
    },
}


def _load_annotation_dataframe(annotation_type: str, heritage: bool) -> pd.DataFrame:
    """
    Load the DataFrame required for one annotation type.

    DataFrames are loaded lazily, unlike the original EDAMannot.py globals.
    """
    relative_path = _ANNOTATION_DATAFRAME_PATHS[annotation_type][heritage]
    path = current_dir + "/" + relative_path

    return pd.read_csv(
        path,
        sep="\t",
        compression="bz2",
    )


def fetch_annotations(
    tools,
    annotation_types=("Topic",),
    heritage=True,
    with_label=True,
):
    """
    Fetch EDAM annotations for one or more tools.

    This preserves the public behavior expected by CLI.py.
    """
    tools = normalize_tool_input(tools)

    resolved_types = [
        _resolve_annotation_type(annotation_type)
        for annotation_type in annotation_types
    ]
    resolved_types = list(dict.fromkeys(resolved_types))

    result = {tool: {} for tool in tools}

    for ann_type in resolved_types:
        df = _load_annotation_dataframe(
            ann_type,
            heritage,
        )

        if ann_type == "Topic":
            col_uri = "topic"
            col_label = "topicLabel"
        else:
            col_uri = "operation"
            col_label = "operationLabel"

        for tool in tools:
            filtered = df[df["tool"] == tool]

            if with_label:
                annotations = [
                    {
                        "URI": row[col_uri],
                        "label": row[col_label],
                    }
                    for _, row in filtered.iterrows()
                ]
            else:
                annotations = [
                    {"URI": row[col_uri]}
                    for _, row in filtered.iterrows()
                ]

            result[tool][ann_type] = annotations

    return result


# ---------------------------------------------------------------------------
# OUTPUT / COMBINED API
# ---------------------------------------------------------------------------

def format_tool_annotations(metrics: dict) -> dict:
    """
    Format metrics using the annotation structure expected by EDAMannot.
    """
    return {"annotation": {"Metrics": metrics}}


def to_json(annotations: Dict) -> str:
    """
    Serialize annotations using the historical EDAMannot JSON structure.
    """
    return json.dumps(
        {"annotation": annotations},
        indent=2,
    )


def fetch_annotations_with_metrics(
    tools,
    annotation_types=("Topic", "Operation"),
    heritage=True,
    with_label=True,
    metric="all",
    include_annotations=True,
):
    """
    Fetch EDAM annotations and tool metrics together.

    The returned structure is compatible with the current QC command in
    CLI.py and with the original EDAMannot.py implementation.
    """
    tools = normalize_tool_input(tools)
    combined = {}

    annotation_data = {}

    if include_annotations:
        annotation_data = fetch_annotations(
            tools,
            annotation_types=annotation_types,
            heritage=heritage,
            with_label=with_label,
        )

    for tool in tools:
        try:
            metrics = get_tool_metrics(
                tool,
                heritage=heritage,
                metric=metric,
            )
        except ValueError:
            metrics = {}

        topic_anns = (
            annotation_data.get(tool, {}).get("Topic", [])
            if include_annotations
            else []
        )
        op_anns = (
            annotation_data.get(tool, {}).get("Operation", [])
            if include_annotations
            else []
        )
        data_anns = (
            annotation_data.get(tool, {}).get("Data", [])
            if include_annotations
            else []
        )
        format_anns = (
            annotation_data.get(tool, {}).get("Format", [])
            if include_annotations
            else []
        )

        combined[tool] = {
            "annotation": {
                "Topic": topic_anns,
                "Operation": op_anns,
                "Data": data_anns,
                "Format": format_anns,
                "Metrics": metrics,
            }
        }

    return combined


__all__ = [
    "get_tool_url",
    "normalize_tool_input",
    "getToolLabel",
    "getToolURIByLabel",
    "getToolTopics",
    "getToolOperations",
    "getToolsCommonTopics",
    "getToolsCommonOperations",
    "fetch_annotations",
    "format_tool_annotations",
    "to_json",
    "fetch_annotations_with_metrics",
]
