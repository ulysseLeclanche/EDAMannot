import os
import time
import click
import subprocess

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


def launch_fuseki_server(bioschemas_file: str = None, edam_file: str = None):
    """
    Launch Fuseki server with default or user-provided files.
    Output of Fuseki is redirected to fuseki.log.
    """

    # Default paths
    default_bioschemas = "edam/bioschemas-dump_05_01_2025.ttl"
    default_edam = "edam/EDAM_1.25.owl"

    bioschemas_file = bioschemas_file or default_bioschemas
    edam_file = edam_file or default_edam

    # Normalize paths
    bioschemas_file = os.path.abspath(bioschemas_file)
    edam_file = os.path.abspath(edam_file)

    click.echo(f"Using Bioschemas file: {bioschemas_file}")
    click.echo(f"Using EDAM file:        {edam_file}")

    fuseki_home = os.environ.get("FUSEKI_HOME")
    if not fuseki_home:
        raise click.ClickException("FUSEKI_HOME is not set in environment variables.")

    fuseki_cmd = [
        f"{fuseki_home}/fuseki-server",
        f"--file={bioschemas_file}",
        f"--file={edam_file}",
        "/biotoolsEdam",
    ]

    click.echo("\nLaunching Fuseki server → fuseki.log")
    click.echo(" ".join(fuseki_cmd))

    # --- redirect to log file ---
    log_file = open("fuseki.log", "w")

    process = subprocess.Popen(fuseki_cmd, stdout=log_file, stderr=log_file)

    time.sleep(5)  # let server initialize

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
    # STEP 1 — Launch Fuseki server (with defaults if none given)
    # ------------------------------------------------------------
    fuseki = launch_fuseki_server(bioschemas_file=bioschemas_file, edam_file=edam_file)

    click.echo("\n=== STEP 2 — Generating Dataframes ===")
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

        # ------------------------------------------------------------
        # SUCCESS
        # ------------------------------------------------------------
        click.echo("\n=== STEP 3 — COMPLETED ===")
        click.echo("Generated files:")
        for f in generated_files:
            click.echo(f"  - {f}")

    except Exception as e:
        click.echo("\n ERROR while generating dataframes!", err=True)
        click.echo(str(e))
        fuseki.terminate()
        raise click.Abort()

    click.echo("\nShut down Fuseki server")
    fuseki.terminate()


# -----------------------------------------------------
# MAIN
# -----------------------------------------------------

if __name__ == "__main__":
    cli()
