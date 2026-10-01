"""
EDAM neighbor inference utilities.

This module contains the EDAM graph-based inference logic extracted from
EDAMannot.py.

Responsibilities
----------------
- Load the precomputed EDAM neighbor graph.
- Retrieve outgoing and incoming neighbor relations.
- Convert EDAM URIs to EDAM identifiers.
- Merge inferred annotations with existing annotations.
"""

from pathlib import Path
from typing import Iterable

from rdflib import Graph, Namespace

from .config import edamURI, edam_neighbors_file


# ---------------------------------------------------------------------------
# EDAM NAMESPACE
# ---------------------------------------------------------------------------

EDAM = Namespace(edamURI)


# ---------------------------------------------------------------------------
# INTERNAL GRAPH LOADING
# ---------------------------------------------------------------------------

def _load_edam_neighbors_graph(
    graph_path: str | Path = edam_neighbors_file,
) -> Graph:
    """
    Load the precomputed EDAM neighbor graph from Turtle.

    Parameters
    ----------
    graph_path:
        Path to the Turtle file containing EDAM neighbor relations.

    Returns
    -------
    rdflib.Graph
        Parsed RDF graph.
    """
    graph = Graph()
    graph.parse(str(graph_path), format="turtle")
    return graph


# ---------------------------------------------------------------------------
# NEIGHBOR INFERENCE
# ---------------------------------------------------------------------------

def infer_edam_neighbors(
    input_entities,
    include_predicates=None,
    include_inverse=True,
):
    """
    Infer EDAM neighbors for a collection of EDAM identifiers.

    Parameters
    ----------
    input_entities:
        EDAM identifiers such as ``topic_3070`` or ``operation_...``.
    include_predicates:
        Optional iterable containing predicate local names to retain.
    include_inverse:
        Also inspect incoming relations when True.

    Returns
    -------
    dict
        Structure compatible with the original EDAMannot implementation:

        {
            "edam:<entity>": {
                "edam:<predicate>": ["target", ...],
                "edam:inverse_<predicate>": ["source", ...]
            }
        }
    """
    graph = _load_edam_neighbors_graph()
    results = {}

    for entity in input_entities:
        entity_uri = EDAM[entity]
        entity_key = f"edam:{entity}"
        results[entity_key] = {}

        # Outgoing relations.
        for predicate, target in graph.predicate_objects(entity_uri):
            predicate_name = str(predicate).rsplit("/", 1)[-1]

            if (
                include_predicates
                and predicate_name not in include_predicates
            ):
                continue

            results[entity_key].setdefault(
                f"edam:{predicate_name}",
                [],
            ).append(
                str(target).rsplit("/", 1)[-1]
            )

        # Incoming relations.
        if include_inverse:
            for source, predicate in graph.subject_predicates(entity_uri):
                predicate_name = str(predicate).rsplit("/", 1)[-1]

                if (
                    include_predicates
                    and predicate_name not in include_predicates
                ):
                    continue

                results[entity_key].setdefault(
                    f"edam:inverse_{predicate_name}",
                    [],
                ).append(
                    str(source).rsplit("/", 1)[-1]
                )

    return results


# ---------------------------------------------------------------------------
# URI / ANNOTATION HELPERS
# ---------------------------------------------------------------------------

def edam_uri_to_id(uri: str) -> str:
    """
    Extract the final identifier from an EDAM URI.
    """
    return uri.rsplit("/", 1)[-1]


def merge_annotations(existing, inferred):
    """
    Merge two annotation lists without duplicate URIs.

    Parameters
    ----------
    existing:
        Existing annotations, each containing a ``URI`` key.
    inferred:
        Inferred annotations, each containing a ``URI`` key.

    Returns
    -------
    list
        Existing annotations followed by previously unseen inferred
        annotations.
    """
    seen = {annotation["URI"] for annotation in existing}
    merged = list(existing)

    for annotation in inferred:
        if annotation["URI"] not in seen:
            merged.append(annotation)
            seen.add(annotation["URI"])

    return merged


def infer_neighbors_from_annotations(
    annotations,
    ann_type,
    with_label=True,
):
    """
    Infer neighboring EDAM annotations from existing annotations.

    Parameters
    ----------
    annotations:
        Annotation list in the structure returned by
        ``annotations.fetch_annotations()``.
    ann_type:
        Annotation type. Kept as a parameter for compatibility with the
        original EDAMannot API.
    with_label:
        Kept as a compatibility parameter. The original implementation did
        not retrieve labels for inferred annotations.

    Returns
    -------
    list[dict]
        Inferred annotations in ``{"URI": ...}`` format.

    Notes
    -----
    ``ann_type`` and ``with_label`` are intentionally retained even though
    the current inference algorithm does not use them to change the query.
    """
    # Preserve the original behavior: the annotation type and label flag are
    # API parameters, but inference itself operates directly on EDAM URIs.
    del ann_type, with_label

    edam_ids = [
        edam_uri_to_id(annotation["URI"])
        for annotation in annotations
    ]

    neighbors = infer_edam_neighbors(
        edam_ids,
        include_inverse=True,
    )

    inferred = []

    for relations in neighbors.values():
        for targets in relations.values():
            for target in targets:
                inferred.append(
                    {
                        "URI": f"{edamURI}{target}",
                    }
                )

    return inferred


__all__ = [
    "EDAM",
    "infer_edam_neighbors",
    "edam_uri_to_id",
    "merge_annotations",
    "infer_neighbors_from_annotations",
]
