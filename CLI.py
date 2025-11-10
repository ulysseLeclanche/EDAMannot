import os
import time
import click
import subprocess
import time
import signal

# Import your dataframe-generating functions
import EDAMannot as edam

# -----------------------------------------------------
# CLI GROUP
# -----------------------------------------------------


@click.group()
def cli():
    """EDAM Annot CLI."""
    pass


# -----------------------------------------------------
# Fuseki launcher
# -----------------------------------------------------
@cli.group()
def fuseki():
    """Manage Apache Fuseki server."""
    pass

FUSEKI_PID_FILE = ".fuseki.pid"
FUSEKI_LOG_FILE = "fuseki.log"
@fuseki.command("start")
@click.option("--bioschemas", type=click.Path(exists=True), default=None,
              help="Custom Bioschemas TTL file.")
@click.option("--edam", type=click.Path(exists=True), default=None,
              help="Custom EDAM OWL file.")
def fuseki_start(bioschemas, edam):
    """Start the Fuseki server."""
    process = launch_fuseki_server(bioschemas, edam)
    click.echo(f"Fuseki server STARTED (PID={process.pid})")

@fuseki.command("stop")
def fuseki_stop():
    """Stop Fuseki server."""
    if not os.path.exists(FUSEKI_PID_FILE):
        click.echo("Fuseki is not running.")
        return

    with open(FUSEKI_PID_FILE) as f:
        pid = int(f.read().strip())

    try:
        os.kill(pid, signal.SIGTERM)
        time.sleep(1)
        click.echo(f"Fuseki server STOPPED (PID={pid})")
    except Exception:
        click.echo("Failed to terminate Fuseki (maybe already dead).")

    # cleanup
    if os.path.exists(FUSEKI_PID_FILE):
        os.remove(FUSEKI_PID_FILE)

def launch_fuseki_server(bioschemas_file: str = None, edam_file: str = None):
    """
    Launch Fuseki server with default or user-provided files.
    Writes PID to .fuseki.pid
    Redirects output to fuseki.log
    """

    # 1. If already running → reuse
    if os.path.exists(FUSEKI_PID_FILE):
        with open(FUSEKI_PID_FILE) as f:
            pid = int(f.read().strip())
        try:
            os.kill(pid, 0)  # check process exists
            click.echo(f"Fuseki already running (PID={pid})")
            return subprocess.Popen([], pid=pid)  # dummy object
        except Exception:
            pass  # stale pid → start new one

    # 2. Defaults
    default_bioschemas = "edam/bioschemas-dump_05_01_2025.ttl"
    default_edam = "edam/EDAM_1.25.owl"

    bioschemas_file = bioschemas_file or default_bioschemas
    edam_file = edam_file or default_edam

    bioschemas_file = os.path.abspath(bioschemas_file)
    edam_file = os.path.abspath(edam_file)

    click.echo(f"Using Bioschemas : {bioschemas_file}")
    click.echo(f"Using EDAM       : {edam_file}")

    # 3. FUSEKI_HOME check
    fuseki_home = os.environ.get("FUSEKI_HOME")
    if not fuseki_home:
        raise click.ClickException("FUSEKI_HOME is not set.")

    fuseki_cmd = [
        f"{fuseki_home}/fuseki-server",
        f"--file={bioschemas_file}",
        f"--file={edam_file}",
        "/biotoolsEdam",
    ]

    click.echo("\nLaunching Fuseki → fuseki.log")
    click.echo(" ".join(fuseki_cmd))

    log_file = open(FUSEKI_LOG_FILE, "w")

    process = subprocess.Popen(fuseki_cmd, stdout=log_file, stderr=log_file)

    # Write PID
    with open(FUSEKI_PID_FILE, "w") as pidf:
        pidf.write(str(process.pid))

    # Give Fuseki time to boot
    time.sleep(5)

    return process



# -----------------------------------------------------
# INIT COMMAND
# -----------------------------------------------------


@cli.command(name="init")
@click.option("--bioschemas-file", "-b", help="Path to Bioschemas .ttl file")
@click.option("--edam-file", "-e", help="Path to EDAM .owl file")
def initialize(bioschemas_file, edam_file):
    """
    Initialize Fuseki server and generate all EDAMannot dataframes.
    """
    click.echo("=== EDAMannot Initialization ===")

    # ------------------------------------------------------------
    # STEP 1 — Launch Fuseki server
    # ------------------------------------------------------------
    fuseki = launch_fuseki_server(bioschemas_file=bioschemas_file, edam_file=edam_file)

    click.echo("\n=== Query bioschemas-file, generate and calculate dataframes ===")
    os.makedirs("Dataframe", exist_ok=True)

    generated_files = []

    try:
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
        click.echo("→ Getting dfToolTopic_NotOWLClass")
        dfToolTopic_NotOWLClass = edam.get_dfToolTopic_NotOWLClass()
        generated_files.append("Dataframe/dfToolTopic_NotOWLClass.tsv.bz2")

        click.echo("→ Getting dfTool_ObsoleteOperation")
        dfTool_ObsoleteOperation = edam.get_dfTool_ObsoleteOperation()
        generated_files.append("Dataframe/dfTool_ObsoleteOperation.tsv.bz2")

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

    except Exception as e:
        click.echo("\n ERROR while generating dataframes!", err=True)
        click.echo(str(e))
        fuseki.terminate()
        raise click.Abort()

    click.echo("\nShut down Fuseki server")
    fuseki.terminate()


@cli.command(name="describe")
@click.argument("tools", nargs=-1)
@click.option(
    "--annotation", "-a", "--annotation_type",
    multiple=True,
    default=["Topic"],
    type=click.Choice(["Topic", "Operation", "T", "O"], case_sensitive=False),
    help="Annotation type(s): Topic (T) and/or Operation (O). Can be used multiple times."
)
@click.option("--transitive", "-t", is_flag=True, default=False,
              help="Use transitive annotations (default: False)")
@click.option(
    "--no-label", "-nL",
    is_flag=True,
    default=False,
    help="Exclude labels from output (default: include labels)"
)
@click.option(
    "--output_format", "-f", default="json",
    type=click.Choice(["json", "dict", "ttl", "sparql"], case_sensitive=False),
    help="Output format"
)
def describe(tools, annotation, transitive, no_label, output_format):
    """
    Describe tools with EDAM annotations.

    Example:
    python EDAMannot.py describe qiime2 -a T -a O -t -f json
    """
    # Ensure Fuseki is running
    launch_fuseki_server()

    annotations = edam.fetch_annotations(
        tools,
        annotation_types=annotation,
        transitive=transitive,
        with_label=not no_label
    )

    if output_format.lower() == "json":
        click.echo(edam.to_json(annotations))

    elif output_format.lower() == "dict":
        click.echo(annotations)

    elif output_format.lower() == "ttl":
        click.echo(edam._format_as_turtle(annotations, include_labels=not no_label))

    elif output_format.lower() == "sparql":
        tool_urls = edam.normalize_tool_input(tools)
        atypes = [edam._resolve_annotation_type(a) for a in annotation]
        query = edam._format_as_sparql(tool_urls, atypes, transitive)
        click.echo(query)
        
@cli.command(name="describe-graph")
@click.option("-t", "--show-topics", is_flag=True, default=False,
              help="Include related topics.")
@click.option("-o", "--show-operations", is_flag=True, default=False,
              help="Include related operations.")
@click.option("-h", "--highlight", is_flag=True, default=False,
              help="Highlight intersection between tools.")
@click.option("-n", "--show-deprecated", is_flag=True, default=False,
              help="Include deprecated annotations.")
@click.option("--title", multiple=True, required=True,
              help="Tool name(s) (e.g. bwa, qiime2). Can be used multiple times.")
@click.option("--output_format",
              type=click.Choice(["SVG", "PNG", "PDF", "CSV"], case_sensitive=False),
              default="SVG",
              help="Graph output format.")
@click.option("-O", "--output", type=click.Path(writable=True),
              help="Output filename without extension.")
def describe_graph(show_topics, show_operations, highlight,
                   show_deprecated, title, output_format, output):
    """
    Describe one or more tools using EDAM annotations and produce a graph or CSV.

    Example:
        cli.py describe-graph -t -o -h -n --title bwa --title qiime2 --output_format CSV
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

    data = graph.draw(prog="dot", format=output_format.lower())
    filename = f"{output or 'common_graph'}.{output_format.lower()}"
    with open(filename, "wb") as f:
        f.write(data)

    click.echo(f"Graph saved as {filename}")
# -----------------------------------------------------
# MAIN
# -----------------------------------------------------

if __name__ == "__main__":
    cli()
