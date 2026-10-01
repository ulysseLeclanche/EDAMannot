"""
Metrics computation utilities for EDAMannot.

This module contains the metric-related logic extracted from EDAMannot.py:

- Mutual Information between two tools;
- topic metrics (transitive and non-transitive);
- operation metrics (transitive and non-transitive);
- combined tool metrics (transitive and non-transitive);
- retrieval of metrics for a single tool.

DataFrames are loaded only when a metric function needs them.
"""

from typing import Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# INTERNAL DATAFRAME LOADER
# ---------------------------------------------------------------------------

def _read_dataframe(path: str) -> pd.DataFrame:
    """
    Read one of the compressed EDAMannot DataFrames.

    The original code consistently uses tab-separated BZ2-compressed files.
    Keeping this in one helper makes the metric functions independent from
    module-level DataFrame globals.
    """
    return pd.read_csv(
        path,
        sep="\t",
        compression="bz2",
    )


# ---------------------------------------------------------------------------
# MUTUAL INFORMATION
# ---------------------------------------------------------------------------

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



# ---------------------------------------------------------------------------
# TOPIC METRICS
# ---------------------------------------------------------------------------

def compute_topic_metrics(
    dfToolTopicTransitive_path="Dataframe/dfToolTopicTransitive.tsv.bz2",
    output_path="Dataframe/dfTopicmetrics.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute topic metrics including Information Content (IC) and entropy.

    Parameters
    ----------
    dfToolTopicTransitive_path : str
        Path to the transitive tool-topic dataframe.
    output_path : str
        Path to save the topic metrics dataframe.

    Returns
    -------
    pd.DataFrame
        Dataframe containing topic metrics.
    """
    dfToolTopicTransitive = _read_dataframe(dfToolTopicTransitive_path)

    dfTopicNbTools = (
        dfToolTopicTransitive[["tool", "topic"]]
        .groupby(by="topic")
        .size()
        .reset_index(name="nbTools")
        .sort_values(by="nbTools", ascending=False)
    )

    dfTopic = (
        dfToolTopicTransitive[["topic", "topicLabel"]]
        .drop_duplicates(
            subset=["topic", "topicLabel"],
            keep="first",
        )
        .reset_index(drop=True)
    )

    dfTopic = dfTopic.join(
        dfTopicNbTools.set_index("topic"),
        on="topic",
    )

    nbToolsWithTopic = dfToolTopicTransitive["tool"].nunique()

    dfTopic["frequence"] = dfTopic["nbTools"] / nbToolsWithTopic
    dfTopic["IC"] = -np.log2(dfTopic["frequence"])
    dfTopic["entropy"] = dfTopic["frequence"] * dfTopic["IC"]

    dfTopic.to_csv(
        output_path,
        sep="\t",
        index=False,
        compression="bz2",
    )

    return dfTopic


def compute_topic_metrics_NT(
    dfToolTopic_path="Dataframe/dfToolTopic.tsv.bz2",
    output_path="Dataframe/dfTopicmetrics_NT.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute topic metrics without transitive closure.

    Parameters
    ----------
    dfToolTopic_path : str
        Path to the non-transitive tool-topic dataframe.
    output_path : str
        Path to save the topic metrics dataframe.

    Returns
    -------
    pd.DataFrame
        Dataframe containing non-transitive topic metrics.
    """
    dfToolTopic = _read_dataframe(dfToolTopic_path)

    nbToolsWithTopic = dfToolTopic["tool"].nunique()

    dfTopicNbTools = (
        dfToolTopic[["tool", "topic"]]
        .groupby(by="topic")
        .size()
        .reset_index(name="nbTools")
        .sort_values(by="nbTools", ascending=False)
    )

    dfTopic = (
        dfToolTopic[["topic", "topicLabel"]]
        .drop_duplicates(
            subset=["topic", "topicLabel"],
            keep="first",
        )
        .reset_index(drop=True)
    )

    dfTopic = dfTopic.join(
        dfTopicNbTools.set_index("topic"),
        on="topic",
    )

    dfTopic["frequence"] = dfTopic["nbTools"] / nbToolsWithTopic
    dfTopic["IC"] = -np.log2(dfTopic["frequence"])
    dfTopic["entropy"] = dfTopic["frequence"] * dfTopic["IC"]

    dfTopic.to_csv(
        output_path,
        sep="\t",
        index=False,
        compression="bz2",
    )

    return dfTopic


# ---------------------------------------------------------------------------
# OPERATION METRICS
# ---------------------------------------------------------------------------

def compute_operation_metrics(
    dfToolOperationTransitive_path="Dataframe/dfToolOperationTransitive.tsv.bz2",
    output_path="Dataframe/dfOperationmetrics.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute operation metrics including Information Content (IC) and entropy.

    Parameters
    ----------
    dfToolOperationTransitive_path : str
        Path to the transitive tool-operation dataframe.
    output_path : str
        Path to save the operation metrics dataframe.

    Returns
    -------
    pd.DataFrame
        Dataframe containing operation metrics.
    """
    dfToolOperationTransitive = _read_dataframe(
        dfToolOperationTransitive_path
    )

    nbToolsWithOperation = dfToolOperationTransitive["tool"].nunique()

    dfOperationNbTools = (
        dfToolOperationTransitive[["tool", "operation"]]
        .groupby(by="operation")
        .size()
        .reset_index(name="nbTools")
        .sort_values(by="nbTools", ascending=False)
    )

    dfOperation = (
        dfToolOperationTransitive[["operation", "operationLabel"]]
        .drop_duplicates(
            subset=["operation", "operationLabel"],
            keep="first",
        )
        .reset_index(drop=True)
    )

    dfOperation = dfOperation.join(
        dfOperationNbTools.set_index("operation"),
        on="operation",
    )

    dfOperation["frequence"] = (
        dfOperation["nbTools"] / nbToolsWithOperation
    )
    dfOperation["IC"] = -np.log2(dfOperation["frequence"])
    dfOperation["entropy"] = (
        dfOperation["frequence"] * dfOperation["IC"]
    )

    dfOperation.to_csv(
        output_path,
        sep="\t",
        index=False,
        compression="bz2",
    )

    return dfOperation


def compute_operation_metrics_NT(
    dfToolOperation_path="Dataframe/dfToolOperation.tsv.bz2",
    output_path="Dataframe/dfOperationmetrics_NT.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute operation metrics without transitive closure.

    Parameters
    ----------
    dfToolOperation_path : str
        Path to the non-transitive tool-operation dataframe.
    output_path : str
        Path to save the operation metrics dataframe.

    Returns
    -------
    pd.DataFrame
        Dataframe containing non-transitive operation metrics.
    """
    dfToolOperation = _read_dataframe(dfToolOperation_path)

    nbToolsWithOperation = dfToolOperation["tool"].nunique()

    dfOperationNbTools = (
        dfToolOperation[["tool", "operation"]]
        .groupby(by="operation")
        .size()
        .reset_index(name="nbTools")
        .sort_values(by="nbTools", ascending=False)
    )

    dfOperation = (
        dfToolOperation[["operation", "operationLabel"]]
        .drop_duplicates(
            subset=["operation", "operationLabel"],
            keep="first",
        )
        .reset_index(drop=True)
    )

    dfOperation = dfOperation.join(
        dfOperationNbTools.set_index("operation"),
        on="operation",
    )

    dfOperation["frequence"] = (
        dfOperation["nbTools"] / nbToolsWithOperation
    )
    dfOperation["IC"] = -np.log2(dfOperation["frequence"])
    dfOperation["entropy"] = (
        dfOperation["frequence"] * dfOperation["IC"]
    )

    dfOperation.to_csv(
        output_path,
        sep="\t",
        index=False,
        compression="bz2",
    )

    return dfOperation


# ---------------------------------------------------------------------------
# COMBINED TOOL METRICS
# ---------------------------------------------------------------------------

def compute_tool_metrics_with_transitive(
    dfToolTopicTransitive_path="Dataframe/dfToolTopicTransitive.tsv.bz2",
    dfToolOperationTransitive_path="Dataframe/dfToolOperationTransitive.tsv.bz2",
    dfTopic_metrics_path="Dataframe/dfTopicmetrics.tsv.bz2",
    dfOperation_metrics_path="Dataframe/dfOperationmetrics.tsv.bz2",
    dfTool_path="Dataframe/dfTool.tsv.bz2",
    output_path="Dataframe/dfToolallmetrics.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute combined tool metrics using transitive topic and operation
    relationships.
    """
    dfToolTopicTransitive = _read_dataframe(
        dfToolTopicTransitive_path
    )
    dfToolOperationTransitive = _read_dataframe(
        dfToolOperationTransitive_path
    )
    dfTopic = _read_dataframe(dfTopic_metrics_path)
    dfOperation = _read_dataframe(dfOperation_metrics_path)
    dfTool = _read_dataframe(dfTool_path)

    df_topic_scores = (
        dfToolTopicTransitive.join(
            dfTopic[["topic", "IC", "entropy"]].set_index("topic"),
            on="topic",
        )[["tool", "IC"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"IC": "topicScore"})
        .reset_index()
    )

    dfTool = dfTool.join(
        df_topic_scores.set_index("tool"),
        on="tool",
    )
    dfTool["topicScore"] = dfTool["topicScore"].fillna(0)

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

    dfTool = dfTool.join(
        df_operation_scores.set_index("tool"),
        on="tool",
    )
    dfTool["operationScore"] = dfTool["operationScore"].fillna(0)

    dfTool["score"] = (
        dfTool["topicScore"] + dfTool["operationScore"]
    )

    df_topic_entropy = (
        dfToolTopicTransitive.join(
            dfTopic[["topic", "IC", "entropy"]].set_index("topic"),
            on="topic",
        )[["tool", "entropy"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"entropy": "topicEntropy"})
        .reset_index()
    )

    dfTool = dfTool.join(
        df_topic_entropy.set_index("tool"),
        on="tool",
    )
    dfTool["topicEntropy"] = dfTool["topicEntropy"].fillna(0)

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

    dfTool = dfTool.join(
        df_operation_entropy.set_index("tool"),
        on="tool",
    )
    dfTool["operationEntropy"] = dfTool["operationEntropy"].fillna(0)

    dfTool["entropy"] = (
        dfTool["topicEntropy"] + dfTool["operationEntropy"]
    )

    dfTool.to_csv(
        output_path,
        sep="\t",
        index=False,
        compression="bz2",
    )

    return dfTool


def compute_tool_metrics_non_transitive(
    dfToolTopic_path="Dataframe/dfToolTopic.tsv.bz2",
    dfToolOperation_path="Dataframe/dfToolOperation.tsv.bz2",
    dfTopic_metrics_NT_path="Dataframe/dfTopicmetrics_NT.tsv.bz2",
    dfOperation_metrics_NT_path="Dataframe/dfOperationmetrics_NT.tsv.bz2",
    dfTool_path="Dataframe/dfTool.tsv.bz2",
    output_path="Dataframe/dfToolallmetrics_NT.tsv.bz2",
) -> pd.DataFrame:
    """
    Compute combined tool metrics using non-transitive relationships.
    """
    dfToolTopic = _read_dataframe(dfToolTopic_path)
    dfToolOperation = _read_dataframe(dfToolOperation_path)
    dfTopicmetrics_NT = _read_dataframe(dfTopic_metrics_NT_path)
    dfOperationmetrics_NT = _read_dataframe(dfOperation_metrics_NT_path)
    dfTool = _read_dataframe(dfTool_path)

    df_topic_scores = (
        dfToolTopic.join(
            dfTopicmetrics_NT[["topic", "IC", "entropy"]].set_index("topic"),
            on="topic",
        )[["tool", "IC"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"IC": "topicScore"})
        .reset_index()
    )

    dfTool = dfTool.join(
        df_topic_scores.set_index("tool"),
        on="tool",
    )
    dfTool["topicScore"] = dfTool["topicScore"].fillna(0)

    df_operation_scores = (
        dfToolOperation.join(
            dfOperationmetrics_NT[["operation", "IC", "entropy"]]
            .set_index("operation"),
            on="operation",
        )[["tool", "IC"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"IC": "operationScore"})
        .reset_index()
    )

    dfTool = dfTool.join(
        df_operation_scores.set_index("tool"),
        on="tool",
    )
    dfTool["operationScore"] = dfTool["operationScore"].fillna(0)

    dfTool["score"] = (
        dfTool["topicScore"] + dfTool["operationScore"]
    )

    df_topic_entropy = (
        dfToolTopic.join(
            dfTopicmetrics_NT[["topic", "IC", "entropy"]].set_index("topic"),
            on="topic",
        )[["tool", "entropy"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"entropy": "topicEntropy"})
        .reset_index()
    )

    dfTool = dfTool.join(
        df_topic_entropy.set_index("tool"),
        on="tool",
    )
    dfTool["topicEntropy"] = dfTool["topicEntropy"].fillna(0)

    df_operation_entropy = (
        dfToolOperation.join(
            dfOperationmetrics_NT[["operation", "IC", "entropy"]]
            .set_index("operation"),
            on="operation",
        )[["tool", "entropy"]]
        .groupby(by="tool")
        .sum()
        .rename(columns={"entropy": "operationEntropy"})
        .reset_index()
    )

    dfTool = dfTool.join(
        df_operation_entropy.set_index("tool"),
        on="tool",
    )
    dfTool["operationEntropy"] = dfTool["operationEntropy"].fillna(0)

    dfTool["entropy"] = (
        dfTool["topicEntropy"] + dfTool["operationEntropy"]
    )

    dfTool.to_csv(
        output_path,
        sep="\t",
        index=False,
        compression="bz2",
    )

    return dfTool


# ---------------------------------------------------------------------------
# SINGLE-TOOL METRIC ACCESS
# ---------------------------------------------------------------------------

def get_tool_metrics(
    tool: str,
    heritage: bool = True,
    metric: str = "all",
    df_metrics_transitive_path: str = "Dataframe/dfToolallmetrics.tsv.bz2",
    df_metrics_non_transitive_path: str = "Dataframe/dfToolallmetrics_NT.tsv.bz2",
    df_counts_transitive_path: str = "Dataframe/dftools_nbTopics_nbOperations.tsv.bz2",
    df_counts_non_transitive_path: str = "Dataframe/dfTool_NoTransitive.tsv.bz2",
) -> dict:
    """
    Fetch metrics for a given tool.

    Parameters
    ----------
    tool:
        Tool name or bio.tools URI.
    heritage:
        Whether to use inherited (transitive) metrics.
    metric:
        One of 'ic', 'entropy', 'count', or 'all'.
    df_metrics_transitive_path:
        Transitive tool-metrics DataFrame.
    df_metrics_non_transitive_path:
        Non-transitive tool-metrics DataFrame.
    df_counts_transitive_path:
        Tool counts including transitive topic/operation counts.
    df_counts_non_transitive_path:
        Direct/non-transitive tool counts.

    Returns
    -------
    dict
        Tool metrics in the same structure as the original function.
    """
    tool_url = (
        tool
        if tool.startswith("https://bio.tools/")
        else f"https://bio.tools/{tool}"
    )

    if heritage:
        df_metrics = _read_dataframe(df_metrics_transitive_path)
        df_counts = _read_dataframe(df_counts_transitive_path)
    else:
        df_metrics = _read_dataframe(df_metrics_non_transitive_path)
        df_counts = _read_dataframe(df_counts_non_transitive_path)

    row_metrics = df_metrics[df_metrics["tool"] == tool_url]
    row_counts = df_counts[df_counts["tool"] == tool_url]

    if row_metrics.empty and row_counts.empty:
        raise ValueError(f"Tool not found: {tool_url}")

    metrics = row_metrics.iloc[0] if not row_metrics.empty else {}
    counts = row_counts.iloc[0] if not row_counts.empty else {}

    result = {"Tool": tool_url}

    if metric in ("ic", "all"):
        result["topicScore"] = float(metrics.get("topicScore", 0))
        result["operationScore"] = float(metrics.get("operationScore", 0))
        result["score"] = float(metrics.get("score", 0))

    if metric in ("entropy", "all"):
        result["topicEntropy"] = float(metrics.get("topicEntropy", 0))
        result["operationEntropy"] = float(metrics.get("operationEntropy", 0))
        result["entropy"] = float(metrics.get("entropy", 0))

    if metric in ("count", "all"):
        result["nbTopics"] = int(counts.get("nbTopics", 0))
        result["nbOperations"] = int(counts.get("nbOperations", 0))

    return result


__all__ = [
    "getMutualInformation",
    "compute_topic_metrics",
    "compute_topic_metrics_NT",
    "compute_operation_metrics",
    "compute_operation_metrics_NT",
    "compute_tool_metrics_with_transitive",
    "compute_tool_metrics_non_transitive",
    "get_tool_metrics",
]
