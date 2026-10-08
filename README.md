# EDAMannot

EDAMannot is a command-line toolbox leveraging the **ShareFAIR-KG** knowledge graph from the [ShareFAIR Knowledge Base](https://zenodo.org/records/17737888). It provides processing, metrics, and visualization features based on the [**EDAM ontology**](https://edamontology.org/page).  
The toolbox is particularly useful for EDAM annotations of tools registered in the [bio.tools registry](https://bio.tools/).
This knowledge base and toolbox are supported by the [ShareFAIR](https://projet.liris.cnrs.fr/sharefair/) project WP2.

---

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Examples](#examples)
- [License](#license)
- [Authors and Contact](#authors-and-contact)

---

## Installation

1. **Clone the EDAMannot repository**

```bash
git clone https://github.com/ulysseLeclanche/EDAMannot.git
```

```bash
cd EDAMannot
```

2. **Create a conda environment from the `environment.yml` file**

```bash
conda env create -f environment.yml
```

```bash
conda activate EDAMannot
```

3. **Connect to the ShareFAIR-KG knowledge graph**

Clone the repository and its submodules:

```bash
git clone --recurse-submodules https://gitlab.liris.cnrs.fr/sharefair/knowledge_base_workflow_annotations/ShareFAIR-KG
```

To locally deploy the knowledge base, run `start_knowlegde_base_server.sh` file:
```bash
./start_knowlegde_base_server.sh
```

You can follow the detailed deployment instructions here: [ShareFAIR-KG deployment guide](https://gitlab.liris.cnrs.fr/sharefair/knowledge_base_workflow_annotations/ShareFAIR-KG#deploy-and-query-the-knowlegde-base)

Note: Full functionality of EDAMannot requires a deployed ShareFAIR-KG instance.

## Usage
Once the knowledge graph is deployed, you can use the toolbox from the `EDAMannot` folder. 
You must be in the EDAMannot folder to use it: 

```bash
cd EDAMannot
```

```bash
python3 CLI.py --help
```
It is recommended to initialize the toolbox (especially after updating the knowledge graph):
```bash
python3 CLI.py init
```

## Features

EDAMannot provides the following main commands:

`describe` – Retrieve direct and inherited EDAM annotations for one or more bio.tools tools.

`describe-viz` – Generate a graph representing tool annotations, with colors indicating information content and highlighting intersections.

`QC` – Compute annotation quality metrics, including annotation counts, frequency, informative content (IC), and Shannon entropy.

`enriched_annotation` – Enrich a Bioschemas TTL file with inferred EDAM annotations for Topics, Operations, Data, and Formats using EDAM neighbor relationships.

## Examples
All commands support `--help` for detailed options and examples of use:
```bash
python3 CLI.py command_name --help
```

Here are examples of uses for each of the commands:

### describe

```bash
python3 CLI.py describe https://bio.tools/multiqc --annotation_type Topic --annotation_type Operation --heritage --output_format json
```
or using alias options :
```bash
python3 CLI.py describe qiime2 -a T -a O -h -f json
```

### describe-viz

```bash
python3 CLI.py describe-viz --show-topics --show-operations --highlight
    --show-deprecated --title bwa --title qiime2  --color-by count --color-channel red
    --output_format SVG --output bwa_qiime2_common_graph
```
or using alias options :
```bash
python3 CLI.py describe-viz -o -h -d --title bwa --title qiime2 -cby count -cc red -f SVG -O bwa_qiime2_common_graph
```

### QC

```bash
python3 CLI.py QC https://bio.tools/star --heritage --metric all --output_format json
```
or using alias options :
```bash
python3 CLI.py QC star -h -m all -f json
```

### enriched_annotation

The command automatically builds the EDAM neighbor graph locally, then enriches the Bioschemas dump:

```bash
python3 CLI.py enriched_annotation
```

The command uses two local Fuseki steps on port `3031`:

1. Load `data/EDAM.owl` and generate `data/edam_neighbors.ttl`.
2. Load `data/edam_neighbors.ttl` and enrich `data/bioschemas-dump.ttl`.

The output is written to:

```text
data/bioschemas-dump_enriched.ttl
```

To run the Fuseki steps manually:

```bash
$FUSEKI_HOME/fuseki-server --port=3031 --file="data/EDAM.owl" /edam
python3 edamannot/build_edam_neighbors.py
```

Then stop the first server (Ctrl+C) and run:

```bash
$FUSEKI_HOME/fuseki-server --port=3031 --file="data/edam_neighbors.ttl" /edam
```

## License

This project is licensed under the GNU GPLv3 License. See the [LICENSE](LICENSE) file for details.

## Authors and Contact

- Ulysse LE CLANCHE¹  
- Olivier DAMERON²  
- Alban GAIGNARD³  

### Affiliations

¹ Université Rennes, Inria, CNRS, IRISA—UMR 6074, Rennes 35000, France  
² Université Rennes, Inria, CNRS, IRISA—UMR 6074, Rennes 35000, France  
³ Nantes Université, CNRS, INSERM, l’institut du thorax, F-44000 Nantes, France
