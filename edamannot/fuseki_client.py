"""
Fuseki client for the EDAMannot enriched-annotation workflow.

This implementation is designed for the existing EDAMannot Fuseki setup:

    http://localhost:3030/sharefair

The dataset already contains the EDAM / ShareFAIR knowledge base in its
default graph and is started with:

    ./fuseki-server --loc=/fuseki/databases/ --update /sharefair

EDAMannot uses two temporary named graphs inside that same dataset:

    urn:edamannot:input
    urn:edamannot:enriched

The RDF is kept in Fuseki. Python does not parse the complete TTL with
rdflib. Input and output are streamed through the HTTP Graph Store endpoint.

First-stage workflow implemented here:

    input.ttl
        -> upload to input named graph
        -> copy input graph to enriched named graph
        -> export enriched graph as Turtle

No EDAM enrichment is performed yet. The next step can add the four
SPARQL UPDATE operations for Topic / Operation / Data / Format.
"""

from __future__ import annotations

import argparse
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


DEFAULT_FUSEKI_BASE = "http://localhost:3030/sharefair"

INPUT_GRAPH = "urn:edamannot:input"
ENRICHED_GRAPH = "urn:edamannot:enriched"

TTL_CONTENT_TYPE = "text/turtle"
TURTLE_ACCEPT = "text/turtle"


class FusekiError(RuntimeError):
    """Raised when a Fuseki HTTP operation fails."""


class FusekiClient:
    """Minimal HTTP client for the existing Fuseki dataset."""

    def __init__(
        self,
        base_url: str = DEFAULT_FUSEKI_BASE,
        timeout: int = 600,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.query_url = f"{self.base_url}/query"
        self.update_url = f"{self.base_url}/update"
        self.data_url = f"{self.base_url}/data"
        self.timeout = timeout

    def _request(
        self,
        url: str,
        *,
        method: str = "GET",
        data: bytes | None = None,
        content_type: str | None = None,
        accept: str | None = None,
    ) -> bytes:
        headers: dict[str, str] = {}

        if content_type:
            headers["Content-Type"] = content_type
        if accept:
            headers["Accept"] = accept

        request = urllib.request.Request(
            url,
            data=data,
            headers=headers,
            method=method,
        )

        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                return response.read()

        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            raise FusekiError(
                f"Fuseki HTTP {exc.code} for {method} {url}\n{body}"
            ) from exc

        except urllib.error.URLError as exc:
            raise FusekiError(
                f"Unable to contact Fuseki at {url}: {exc.reason}"
            ) from exc

    @staticmethod
    def _graph_url(data_url: str, graph: str) -> str:
        params = urllib.parse.urlencode({"graph": graph})
        return f"{data_url}?{params}"

    def ping(self) -> None:
        """
        Check that the existing /sharefair SPARQL query endpoint is reachable.

        The query is deliberately sent with GET. This matches the endpoint
        exposed by the current Fuseki command-line configuration.
        """
        params = urllib.parse.urlencode({"query": "ASK {}"})
        url = f"{self.query_url}?{params}"

        self._request(
            url,
            method="GET",
            accept="application/sparql-results+json",
        )

    def clear_graph(self, graph: str) -> None:
        """
        Remove a named graph if it exists.

        We deliberately use SPARQL UPDATE with DROP SILENT instead of
        DELETE /data?graph=..., because Fuseki returns HTTP 404 when the
        requested named graph does not exist yet.
        """
        update = f"DROP SILENT GRAPH <{graph}>"
        self.run_update(update)

    def upload_ttl(
        self,
        input_path: str | Path,
        graph: str = INPUT_GRAPH,
    ) -> None:
        """
        Replace a named graph with a Turtle file.

        The file is sent directly to Fuseki and is never parsed into an
        rdflib Graph.
        """
        path = Path(input_path)

        if not path.is_file():
            raise FileNotFoundError(f"Turtle file not found: {path}")

        url = self._graph_url(self.data_url, graph)

        with path.open("rb") as source:
            request = urllib.request.Request(
                url,
                data=source,
                headers={"Content-Type": TTL_CONTENT_TYPE},
                method="PUT",
            )

            try:
                with urllib.request.urlopen(
                    request,
                    timeout=self.timeout,
                ):
                    pass

            except urllib.error.HTTPError as exc:
                body = exc.read().decode("utf-8", errors="replace")
                raise FusekiError(
                    f"Fuseki HTTP {exc.code} while uploading {path}\n{body}"
                ) from exc

            except urllib.error.URLError as exc:
                raise FusekiError(
                    f"Unable to upload {path} to Fuseki: {exc.reason}"
                ) from exc

    def run_update(self, update: str) -> None:
        """Execute a SPARQL UPDATE request."""
        self._request(
            self.update_url,
            method="POST",
            data=update.encode("utf-8"),
            content_type="application/sparql-update",
        )

    def copy_graph(
        self,
        source_graph: str = INPUT_GRAPH,
        target_graph: str = ENRICHED_GRAPH,
    ) -> None:
        """
        Copy a named graph inside the existing /sharefair dataset.

        The source and target remain isolated from the default graph.
        """
        source = f"<{source_graph}>"
        target = f"<{target_graph}>"

        update = f"""
        INSERT {{
            GRAPH {target} {{
                ?s ?p ?o
            }}
        }}
        WHERE {{
            GRAPH {source} {{
                ?s ?p ?o
            }}
        }}
        """

        self.run_update(update)

    def export_ttl(
        self,
        output_path: str | Path,
        graph: str = ENRICHED_GRAPH,
    ) -> None:
        """
        Export one named graph as Turtle.

        The HTTP response is streamed directly to disk.
        """
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)

        url = self._graph_url(self.data_url, graph)

        request = urllib.request.Request(
            url,
            headers={"Accept": TURTLE_ACCEPT},
            method="GET",
        )

        try:
            with urllib.request.urlopen(
                request,
                timeout=self.timeout,
            ) as response:
                with path.open("wb") as output:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        output.write(chunk)

        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            raise FusekiError(
                f"Fuseki HTTP {exc.code} while exporting graph {graph}\n{body}"
            ) from exc

        except urllib.error.URLError as exc:
            raise FusekiError(
                f"Unable to export graph {graph} from Fuseki: {exc.reason}"
            ) from exc


def prepare_enrichment_graphs(
    client: FusekiClient,
    input_path: str | Path,
) -> None:
    """
    Prepare the two temporary named graphs.

    Steps:
        1. Check /sharefair/query.
        2. Remove previous EDAMannot temporary graphs.
        3. Upload the source TTL to INPUT_GRAPH.
        4. Copy INPUT_GRAPH to ENRICHED_GRAPH.
    """
    client.ping()

    client.clear_graph(INPUT_GRAPH)
    client.clear_graph(ENRICHED_GRAPH)

    client.upload_ttl(input_path, INPUT_GRAPH)
    client.copy_graph(INPUT_GRAPH, ENRICHED_GRAPH)


def cleanup_enrichment_graphs(
    client: FusekiClient,
) -> None:
    """Remove EDAMannot temporary named graphs."""
    try:
        client.clear_graph(INPUT_GRAPH)
    finally:
        client.clear_graph(ENRICHED_GRAPH)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Test the Fuseki integration used by EDAMannot enrichment."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input Turtle file.",
    )

    parser.add_argument(
        "--output",
        help="Output Turtle file.",
    )

    parser.add_argument(
        "--fuseki",
        default=DEFAULT_FUSEKI_BASE,
        help=f"Fuseki dataset URL (default: {DEFAULT_FUSEKI_BASE}).",
    )

    parser.add_argument(
        "--keep-graphs",
        action="store_true",
        help="Keep the temporary named graphs after the run.",
    )

    args = parser.parse_args()

    client = FusekiClient(args.fuseki)

    try:
        print(f"Fuseki dataset: {client.base_url}")
        print(f"Input graph:     {INPUT_GRAPH}")
        print(f"Enriched graph:  {ENRICHED_GRAPH}")

        prepare_enrichment_graphs(
            client,
            args.input,
        )

        print("Fuseki connection: OK")
        print("Input TTL: loaded")
        print("Enriched graph: prepared")

        if args.output:
            client.export_ttl(
                args.output,
                ENRICHED_GRAPH,
            )
            print(f"Exported TTL: {args.output}")

    finally:
        if not args.keep_graphs:
            cleanup_enrichment_graphs(client)
            print("Temporary graphs removed.")


if __name__ == "__main__":
    main()
