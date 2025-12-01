# EDAMannot

EDAMannot is a command-line toolbox leveraging the **ShareFAIR-KG** knowledge graph. It provides processing, metrics, and visualization features based on the **EDAM ontology**.  
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

1. **Clone the EDAMannot repository**:

```bash
git clone https://github.com/YOUR_USERNAME/EDAMannot.git
cd EDAMannot
```

2. **Create a conda environment from the `environment.yml` file**:

```bash
conda env create -f environment.yml
conda activate EDAMannot
```

3. **Connect to the ShareFAIR-KG knowledge graph**:
Clone the repository and its submodules:

```bash
git clone --recurse-submodules https://gitlab.liris.cnrs.fr/sharefair/knowledge_base_workflow_annotations/knowledge_base.git
```

Follow the deployment instructions here: [ShareFAIR-KG deployment guide](https://gitlab.liris.cnrs.fr/sharefair/knowledge_base_workflow_annotations/knowledge_base#deploy-and-query-the-knowlegde-base)

Note: Full functionality of EDAMannot requires a deployed ShareFAIR-KG instance.

## Usage
Once the knowledge graph is deployed, you can use the toolbox from the `EDAMannot` folder:

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

## Examples
All commands support `--help` for detailed options and examples of use:
```bash
python3 CLI.py command_name --help
```

## License

This project is licensed under the GNU GPLv3 License. See the [LICENSE](LICENSE) file for details.

## Authors and Contact

- Ulysse LE CLANCHE
- Olivier DAMERON
- Alban GAIGNARD

Affiliations: Université de Rennes, IRISA de Rennes, Université de Nantes