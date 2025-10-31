import click
from EDAMannot import (
    addToolAndAnnotationsToGraph,
    addToolsAndAnnotationsToGraph,
    getToolURIByLabel,
    getToolTopics,
    getToolOperations,
)

# Default base URI for BioTools
BIOTOOLS_URI = "https://bio.tools/"

@click.group()
def cli():
    """Command-line interface for EDAMannot."""
    pass

# ---- Command 1: Graph of a single tool ----
@click.command()
@click.argument("tool_name")
@click.option("-t", "--show-topics", is_flag=True, default=False, help="Include related topics.")
@click.option("-o", "--show-operations", is_flag=True, default=False, help="Include related operations.")
@click.option("-d", "--show-deprecated", is_flag=True, default=False, help="Include deprecated annotations.")
@click.option("-f", "--format", type=click.Choice(["svg", "png", "pdf"]), default="svg", help="Output format.")
@click.option("-O", "--output", type=click.Path(writable=True), help="Output filename without extension.")
def graph_tool(tool_name, show_topics, show_operations, show_deprecated, output, format):
    """Builds a hierarchical graph for a given tool."""
    uri = BIOTOOLS_URI + tool_name
    graph = addToolAndAnnotationsToGraph(
        uri,
        graph=None,
        showTopics=show_topics,
        showOperations=show_operations,
        showDeprecatedAnnotations=show_deprecated,
        highlightDirectAnnotations=True,
    )

    if output:
        data = graph.draw(prog="dot", format=format)
        with open(f"{output}.{format}", "wb") as f:
            f.write(data)
        click.echo(f"Graph saved as {output}.{format}")
    else:
        click.echo(graph)


# ---- Command 2: Graph for multiple tools ----
@click.command()
@click.argument("tools", nargs=-1)
@click.option("-t", "--show-topics", is_flag=True, default=True, help="Include related topics.")
@click.option("-o", "--show-operations", is_flag=True, default=True, help="Include related operations.")
@click.option("-i", "--highlight-intersection", is_flag=True, default=False, help="Highlight intersection between tools.")
@click.option("-f", "--format", type=click.Choice(["svg", "png", "pdf"]), default="svg", help="Output format.")
@click.option("-O", "--output", type=click.Path(writable=True), help="Output filename without extension.")
def common_graph(tools, show_topics, show_operations, highlight_intersection, output, format):
    """Build a common graph for multiple tools."""
    uris = [BIOTOOLS_URI + t for t in tools]
    graph = addToolsAndAnnotationsToGraph(
        uris,
        graph=None,
        showTopics=show_topics,
        showOperations=show_operations,
        highlightDirectAnnotations=False,
        highlightIntersection=highlight_intersection,
    )

    if output:
        data = graph.draw(prog="dot", format=format)
        with open(f"{output}.{format}", "wb") as f:
            f.write(data)
        click.echo(f"Graph saved as {output}.{format}")
    else:
        click.echo(graph)


# ---- Command 3: Topics ----
@click.command()
@click.argument("tool_uri")
@click.option("-t", "--transitive", is_flag=True, default=False, help="Include transitive topics.")
def topics(tool_uri, transitive):
    """List topics associated with a tool."""
    results = getToolTopics(tool_uri, transitive=transitive)
    for r in results:
        click.echo(r)

# ---- Command 4: Operations ----
@click.command()
@click.argument("tool_uri")
@click.option("-t", "--transitive", is_flag=True, default=False, help="Include transitive operations.")
def operations(tool_uri, transitive):
    """List operations associated with a tool."""
    results = getToolOperations(tool_uri, transitive=transitive)
    for r in results:
        click.echo(r)

# Attach commands to the main CLI group
cli.add_command(graph_tool)
cli.add_command(common_graph)
cli.add_command(topics)
cli.add_command(operations)

if __name__ == "__main__":
    cli(auto_envvar_prefix="BIOTOOLS")  # normal invocation