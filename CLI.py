import os
import click
import json


# Import your dataframe-generating functions
import EDAMannot as edam

# -----------------------------------------------------
# CLI GROUP
# -----------------------------------------------------


@click.group()
def cli():
    """
    EDAMannot is a command-line toolbox using ShareFAIR-KG.
    The toolbox offersprocessing, metrics and visualisation features leveraging the EDAM hierarchy.
    """
    pass


# -----------------------------------------------------
# INIT COMMAND
# -----------------------------------------------------


@cli.command(name="init")
def initialize():
    """
    Compute all tables containing metrics and annotation information for the tools available in bio.tools.

    These dataframes are necessary for the toolkit to function.

    Command usage : python3 CLI.py init
    """
    click.echo("=== EDAMannot Initialization ===")

    # ------------------------------------------------------------
    # STEP 1 — Launch Fuseki server
    # ------------------------------------------------------------

    click.echo("\n=== Query bioschemas-file, generate and calculate dataframes ===")
    os.makedirs("Dataframe", exist_ok=True)

    generated_files = []

    # ------------------------------------------------------------
    # PRIMARY EXTRACTIONS — base objects
    # ------------------------------------------------------------
    click.echo("→ Getting number of tools")
    nb = edam.get_nb_tools()
    click.echo(f"  Found {nb} tools")

    click.echo("→ Getting dfTool")
    dfTool = edam.get_tools_dataframe()
    generated_files.append("Dataframe/dfTool.tsv.bz2")

    click.echo("→ Getting dfToolTopic (non-transitive)")
    dfToolTopic = edam.get_tools_topics_dataframe()
    generated_files.append("Dataframe/dfToolTopic.tsv.bz2")

    click.echo("→ Getting dfToolOperation (non-transitive)")
    dfToolOperation = edam.get_tools_operations_label_dataframe()
    generated_files.append("Dataframe/dfToolOperation.tsv.bz2")

    click.echo("→ Getting dfToolTopicTransitive")
    dfToolTopicTransitive = edam.get_tools_topics_transitive_dataframe()
    generated_files.append("Dataframe/dfToolTopicTransitive.tsv.bz2")

    click.echo("→ Getting dfToolOperationTransitive")
    dfToolOperationTransitive = edam.get_tools_operations_transitive_dataframe()
    generated_files.append("Dataframe/dfToolOperationTransitive.tsv.bz2")

    # ------------------------------------------------------------
    # AGGREGATES — topic/operation counts
    # ------------------------------------------------------------
    click.echo("→ Generating dfTool with transitive nbTopics & nbOperations")
    dfTool_T = edam.get_dftools_with_nbTopics_nbOperations(
        dfTool,
        dfToolTopicTransitive,
        dfToolOperationTransitive,
        "Dataframe/dftools_nbTopics_nbOperations.tsv.bz2",
    )
    generated_files.append("Dataframe/dftools_nbTopics_nbOperations.tsv.bz2")

    # ------------------------------------------------------------
    # REDUNDANCY QUERIES
    # ------------------------------------------------------------
    click.echo("→ Generating df_redundancy_topic")
    df_redundancy_topic = edam.generate_df_redundancy_topic()
    generated_files.append("Dataframe/dfToolTopic_redundancy.tsv.bz2")

    click.echo("→ Generating df_redundancy_operation")
    df_redundancy_operation = edam.generate_df_redundancy_operation()
    generated_files.append("Dataframe/dfToolOperation_redundancy.tsv.bz2")

    # ------------------------------------------------------------
    # NON-REDUNDANT TABLES
    # ------------------------------------------------------------
    click.echo("→ Generating df_topic_no_redundancy")
    df_topic_no_redundancy = edam.generate_df_topic_no_redundancy()
    generated_files.append("Dataframe/df_topic_no_redundancy.tsv.bz2")

    click.echo("→ Generating df_operation_no_redundancy")
    df_operation_no_redundancy = edam.generate_df_operation_no_redundancy()
    generated_files.append("Dataframe/df_operation_no_redundancy.tsv.bz2")

    # ------------------------------------------------------------
    # dfTool without transitive and without redundancy
    # ------------------------------------------------------------
    click.echo("→ Generating dfTool_NoTransitive")
    dfTool_NoTransitive = edam.generate_dfTool_no_transitive()
    generated_files.append("Dataframe/dfTool_NoTransitive.tsv.bz2")

    click.echo("→ Generating dfTool_NoTransitive_NoRedundancy")
    dfTool_NoTransitive_NoRedundancy = (
        edam.generate_dfTool_no_transitive_no_redundancy()
    )
    generated_files.append("Dataframe/dfTool_NoTransitive_NoRedundancy.tsv.bz2")

    # ------------------------------------------------------------
    # SPECIAL DIAGNOSTIC QUERIES
    # ------------------------------------------------------------

    click.echo("→ Getting dfDeprecatedItems")
    dfDeprecatedItems = edam.get_dfDeprecatedItems()
    generated_files.append("Dataframe/dfDeprecatedItems.tsv.bz2")

    click.echo("→ Getting dfDeprecatedSuggestedItems")
    dfDeprecatedSuggestedItems = edam.get_dfDeprecatedSuggestedItems()
    generated_files.append("Dataframe/dfDeprecatedSuggestedItems.tsv.bz2")

    click.echo("→ Getting dfToolsWithSomeDeprecatedTopic")
    dfToolsWithSomeDeprecatedTopic = edam.get_dfToolsWithSomeDeprecatedTopic()
    generated_files.append("Dataframe/dfToolsWithSomeDeprecatedTopic.tsv.bz2")

    click.echo("→ Getting dfToolsWithSomeDeprecatedOperation")
    dfToolsWithSomeDeprecatedOperation = edam.get_dfToolsWithSomeDeprecatedOperation()
    generated_files.append("Dataframe/dfToolsWithSomeDeprecatedOperation.tsv.bz2")

    # ------------------------------------------------------------------
    # 7) METRICS TABLES
    # ------------------------------------------------------------------
    click.echo("→ Compute dfTopicmetrics")
    dfTopicmetrics = edam.compute_topic_metrics()
    generated_files.append("Dataframe/dfTopicmetrics.tsv.bz2")

    click.echo("→ Compute dfTopicmetrics_NT")
    dfTopicmetrics_NT = edam.compute_topic_metrics_NT()
    generated_files.append("Dataframe/dfTopicmetrics_NT.tsv.bz2")

    click.echo("→ Compute dfOperationmetrics")
    dfOperationmetrics = edam.compute_operation_metrics()
    generated_files.append("Dataframe/dfOperationmetrics.tsv.bz2")

    click.echo("→ Compute dfOperationmetrics_NT")
    dfOperationmetrics_NT = edam.compute_operation_metrics_NT()
    generated_files.append("Dataframe/dfOperationmetrics_NT.tsv.bz2")

    click.echo("→ Compute dfToolallmetrics")
    dfToolallmetrics = edam.compute_tool_metrics_with_transitive()
    generated_files.append("Dataframe/dfToolallmetrics.tsv.bz2")

    click.echo("→ Compute dfToolallmetrics_NT")
    dfToolallmetrics_NT = edam.compute_tool_metrics_non_transitive()
    generated_files.append("Dataframe/dfToolallmetrics_NT.tsv.bz2")

    # ------------------------------------------------------------
    # SUCCESS
    # ------------------------------------------------------------
    click.echo("\n=== Generated files: ===")
    for f in generated_files:
        click.echo(f"  - {f}")


@click.command(name="QC")
@click.argument("tools", nargs=-1)
@click.option(
    "--heritage",
    "-h",
    is_flag=True,
    default=False,
    help="Use heritage (inherited) metrics and annotations (default: False)",
)
@click.option(
    "--metric",
    "-m",
    default="all",
    type=click.Choice(["ic", "entropy", "count", "all"], case_sensitive=False),
    help="Which metric to fetch: 'ic', 'entropy', 'count', or 'all'",
)
@click.option(
    "--annotations",
    is_flag=True,
    default=True,
    help="Include EDAM annotations in the output (default: True)",
)
@click.option(
    "--no-annotations",
    "-noa",
    is_flag=True,
    default=False,
    help="Exclude EDAM annotations from the output (alias: -noa)",
)
@click.option(
    "--output_format",
    "-f",
    default="json",
    type=click.Choice(["json", "dict"], case_sensitive=False),
    help="Output format",
)
def qc(tools, heritage, metric, annotations, no_annotations, output_format):
    """
    Show selected metrics (topicScore, operationScore, score, entropy, count) with optional EDAM annotations.

    Examples of command usage :

    python3 CLI.py QC https://bio.tools/star --heritage --metric all --output_format json

    or using alias options :

    python3 CLI.py QC star -h -m all -f json

    To exclude annotations from the output, see only direct annotations metrics and see selected metrics only :
    python3 CLI.py QC star -m ic -noa

    - topicScore: Score based on the IC (Information Content) of topics associated with a tool.
    - operationScore: Score based on the IC of operations associated with a tool.
    - score: Sum of the two previous scores.
    - topicEntropy: Total entropy of topics associated with a tool.
    - operationEntropy: Total entropy of operations associated with a tool.
    - entropy: Sum of the previous entropies.
    """

    # If user passed --no-annotations or -noa, override
    include_annotations = not no_annotations and annotations

    results = edam.fetch_annotations_with_metrics(
        tools,
        annotation_types=("Topic", "Operation"),
        heritage=heritage,
        with_label=True,
        metric=metric,
        include_annotations=include_annotations,
    )

    if output_format.lower() == "json":
        click.echo(json.dumps(results, indent=2))
    else:
        click.echo(results)


@cli.command(name="describe")
@click.argument("tools", nargs=-1)
@click.option(
    "--annotation",
    "-a",
    "--annotation_type",
    multiple=True,
    default=["Topic"],
    type=click.Choice(["Topic", "Operation", "T", "O"], case_sensitive=False),
    help="Annotation type(s): Topic (T) and/or Operation (O). Can be used multiple times.",
)
@click.option(
    "--heritage",
    "-h",
    is_flag=True,
    default=False,
    help="Use heritage (inherited) annotations (default: False)",
)
@click.option(
    "--no-label",
    "-nL",
    is_flag=True,
    default=False,
    help="Exclude labels from output (default: include labels)",
)
@click.option(
    "--output_format",
    "-f",
    default="json",
    type=click.Choice(["json", "dict"], case_sensitive=False),
    help="Output format",
)
def describe(tools, annotation, heritage, no_label, output_format):
    """
    Describe tools with EDAM annotations can use heritage annotations.

    Example command usage :

    python3 CLI.py describe https://bio.tools/qiime2 --annotation_type Topic --annotation_type Operation --heritage --output_format json

    or using alias options :

    python3 CLI.py describe qiime2 -a T -a O -h -f json
    """

    annotations = edam.fetch_annotations(
        tools,
        annotation_types=annotation,
        heritage=heritage,
        with_label=not no_label,
    )

    if output_format.lower() == "json":
        click.echo(edam.to_json(annotations))

    elif output_format.lower() == "dict":
        click.echo(annotations)


@cli.command(name="describe-viz")
@click.option(
    "-t", "--show-topics", is_flag=True, default=False, help="Include related topics."
)
@click.option(
    "-o",
    "--show-operations",
    is_flag=True,
    default=False,
    help="Include related operations.",
)
@click.option(
    "-h",
    "--highlight",
    is_flag=True,
    default=False,
    help="Highlight intersection between tools.",
)
@click.option(
    "-d",
    "--show-deprecated",
    is_flag=True,
    default=False,
    help="Include deprecated annotations.",
)
@click.option(
    "--title",
    multiple=True,
    required=True,
    help="Tool name(s) (e.g. bwa, qiime2). Can be used multiple times.",
)
@click.option(
    "--output_format",
    "-f",
    type=click.Choice(["SVG", "PNG", "PDF", "CSV"], case_sensitive=False),
    default="SVG",
    help="Graph output format.",
)
@click.option(
    "-O",
    "--output",
    type=click.Path(writable=True),
    help="Output filename without extension.",
)
@click.option(
    "--color-by",
    "-cby",
    type=click.Choice(["none", "count", "ic", "entropy"], case_sensitive=False),
    default="none",
    help="Color topics and operations according to a metric.",
)
@click.option(
    "--color-channel",
    "-cc",
    type=click.Choice(
        ["red", "green", "blue", "orange", "yellow", "pink", "grey"],
        case_sensitive=False,
    ),
    default="red",
    help="Color channel used for score-based node coloring.",
)
def describe_graph(
    show_topics,
    show_operations,
    highlight,
    show_deprecated,
    title,
    output_format,
    output,
    color_by,
    color_channel,
):
    """
    Generate a graph to describe one or more tools using direct and legacy EDAM annotations.
    You can color the annotations according to metrics (Count, IC, entropy).

     Example command usage :

     python3 CLI.py describe-viz --show-topics --show-operations --highlight
     --show-deprecated --title bwa --title qiime2  --color-by count --color-channel red
     --output_format SVG --output bwa_qiime2_common_graph

     or using alias options :

     python3 CLI.py describe-viz -o -h -d --title bwa --title qiime2 -cby count -cc red -f SVG -O bwa_qiime2_common_graph
    """

    BIOTOOLS_URI = "https://bio.tools/"
    tool_uris = [BIOTOOLS_URI + t for t in title]

    # -------- Single tool case --------
    if len(tool_uris) == 1:
        uri = tool_uris[0]

        graph = edam.addToolAndAnnotationsToGraph(
            uri,
            graph=None,
            showTopics=show_topics,
            showOperations=show_operations,
            showDeprecatedAnnotations=show_deprecated,
            highlightDirectAnnotations=True,
        )

        # Coloring according to metric
        if color_by.lower() != "none":
            metric_map = {
                "count": "nbTools",
                "ic": "IC",
                "entropy": "entropy",
            }

            metric_col = metric_map[color_by.lower()]

            # Build the dicts for this metric
            dictTopicScore, dictOperationScore = edam.buildTopicOperationDicts(
                metric_col
            )

            # Apply coloring to graph
            edam.colorGraphNodesAccordingToScore(
                graph,
                dictTopicScore,
                dictOperationScore,
                color=color_channel,
            )

        # CSV mode
        if output_format.lower() == "csv":
            topics = edam.getToolTopics(uri)
            operations = edam.getToolOperations(uri)
            filename = f"{output or 'tool_description'}.csv"
            with open(filename, "w") as f:
                f.write("Type,Annotation\n")
                for t in topics:
                    f.write(f"Topic,{t}\n")
                for o in operations:
                    f.write(f"Operation,{o}\n")

            click.echo(f"Saved {filename}")
            return

        # Graph mode
        data = graph.draw(prog="dot", format=output_format.lower())
        filename = f"{output or 'tool_graph'}.{output_format.lower()}"
        with open(filename, "wb") as f:
            f.write(data)

        click.echo(f"Graph saved as {filename}")
        return

    # -------- Multiple tools case --------
    graph = edam.addToolsAndAnnotationsToGraph(
        tool_uris,
        graph=None,
        showTopics=show_topics,
        showOperations=show_operations,
        highlightIntersection=highlight,
        highlightDirectAnnotations=False,
    )

    # Coloring according to metric
    if color_by.lower() != "none":
        metric_map = {
            "count": "nbTools",
            "ic": "IC",
            "entropy": "entropy",
        }
        metric_col = metric_map[color_by.lower()]

        dictTopicScore, dictOperationScore = edam.buildTopicOperationDicts(metric_col)

        edam.colorGraphNodesAccordingToScore(
            graph,
            dictTopicScore,
            dictOperationScore,
            color=color_channel,
        )

    data = graph.draw(prog="dot", format=output_format.lower())
    filename = f"{output or 'common_graph'}.{output_format.lower()}"
    with open(filename, "wb") as f:
        f.write(data)

    click.echo(f"Graph saved as {filename}")


cli.add_command(qc)

# -----------------------------------------------------
# MAIN
# -----------------------------------------------------

if __name__ == "__main__":
    cli()
