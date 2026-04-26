# MultiPrimer Designer

_MultiPrimer Designer_ is a bioinformatics desktop application for designing PCR primers targeting human genetic variants. It provides a graphical interface built on Python/tkinter, covering the full workflow from variant input to a ready-to-use HTML report.

It offers:

* 🧬 **Variant loading & validation**\
  Load genetic variants from CSV or XLSX files using HGVS notation (coding `c.` or genomic `g.`). Variants are automatically validated against MANE Select transcripts via NCBI and Ensembl APIs.

    <details><summary>Supported input format</summary>

    A plain CSV file with three columns (comma, semicolon, or tab-separated):

    ```
    Gene,Transcript,Position
    BRCA1,NM_007294.4,c.5266dup
    TP53,NM_000546.6,c.215C>G
    CFTR,NM_000492.4,c.1521_1523delCTT
    ```

    - **Gene** — gene symbol (e.g. `BRCA1`, `TP53`)
    - **Transcript** — RefSeq accession with version (e.g. `NM_007294.4`)
    - **Position** — HGVS variant (substitutions, deletions, insertions, duplications, complex indels)

    </details>

* 🗺️ **Coordinate mapping**\
  Automatically converts CDS (`c.`) coordinates to genomic (`g.`) coordinates using Ensembl REST API, handling strand orientation, exon boundaries, and both GRCh38 and GRCh37/hg19 reference assemblies.

* 🔬 **Primer design**\
  Wraps the [Primer3](https://primer3.org/) engine to design PCR primers with amplicons 100–500 bp (configurable). Default primer constraints: Tm 57–64 °C, GC% 40–60%, length 18–27 bp.

* 🌍 **Population variant filtering**\
  Primers spanning common SNPs are automatically filtered out. Supports 10 gnomAD/Ensembl populations (global, AFR, EAS, EUR, SAS, AMI, AMR, ASJ, FIN, OTH) with a configurable MAF threshold (default: 0.5%).

* 🔍 **Homology & specificity analysis**\
  Runs BLAST+ (`blastn`) against the human reference genome (GRCh38) to detect off-target binding sites, pseudogenes, and gene duplications. Optional BWA-based specificity checking is also available.

    <details><summary>BLAST+ database</summary>

    The GRCh38 reference genome BLAST database (~900 MB compressed / ~3 GB unpacked) is downloaded automatically the first time homology analysis is run. BLAST+ binaries are installed by `install.py` (see Installation).

    </details>

* 📄 **HTML report export**\
  Generates a professional HTML report containing primer sequences, Tm, GC%, amplicon size, target position, homology hits, and full methodology documentation.

---

## Requirements

- **Python 3.8+** (3.10+ recommended)
- Internet connection (NCBI and Ensembl APIs for sequence retrieval and validation)
- ~3 GB disk space for the BLAST genome database (downloaded on first use)

---

## Installation

The unified `install.py` script handles everything automatically: virtual environment, Python dependencies, MANE transcript database, and BLAST+ binaries.

### macOS

```bash
git clone https://github.com/Adv20202/MultiPrimer_Designer.git
cd MultiPrimer_Designer
python3 install.py
```

BLAST+ is installed via Homebrew (`brew install blast`). If Homebrew is not present, the installer will prompt you.

### Ubuntu / Debian

```bash
git clone https://github.com/Adv20202/MultiPrimer_Designer.git
cd MultiPrimer_Designer
python3 install.py
```

BLAST+ is installed via `apt` (`sudo apt install ncbi-blast+`). The installer will ask for your sudo password.

### Fedora / RHEL

```bash
git clone https://github.com/Adv20202/MultiPrimer_Designer.git
cd MultiPrimer_Designer
python3 install.py
```

BLAST+ is installed via `dnf`/`yum`.

### Windows

```bash
git clone https://github.com/Adv20202/Primer-Designer--GH-.git
cd "Primer-Designer--GH-"
python install.py
```

> Make sure Python is added to PATH (check **"Add Python to PATH"** during installation).

BLAST+ is downloaded from NCBI FTP and unpacked to `tools/blast/bin/` — no administrator rights required.

### What `install.py` does

- Creates a Python virtual environment (`venv/`)
- Installs Python dependencies: `primer3-py`, `openpyxl`, `sv-ttk`
- Downloads the MANE transcript database from NCBI (~10 MB)
- Installs BLAST+ for your platform (see above)
- Generates a platform-specific launch script (`run_primer_designer.sh` / `run_primer_designer.bat`)

---

## Running the application

### macOS / Linux

```bash
./run_primer_designer.sh
```

### Windows

Double-click `run_primer_designer.bat`, or from a terminal:

```bash
run_primer_designer.bat
```

### Directly (all platforms)

```bash
# macOS / Linux
venv/bin/python main.py

# Windows
venv\Scripts\python.exe main.py
```

---

## Quick start

The `data/` directory contains ready-to-use example files:

| File | Contents |
|------|----------|
| `data/example_1.csv` | 5 variants (NBN, CHEK2, BRCA1, APC, ATM) — fast smoke test |
| `data/example_2.csv` | ~60 variants across APC, ATM, BRCA1, BRCA2 |
| `data/example_3.csv` | ~100 variants across CDH1, CDKN2A, CHEK2, MLH1, MSH2, MSH6, MUTYH, NBN, PALB2, PMS2, PTEN, RAD51C, RAD51D, STK11, TP53 |

**Step-by-step workflow:**

1. Launch the application (`./run_primer_designer.sh`)
2. Click **Browse** and select e.g. `data/example_1.csv`
3. Click **Load and Validate** — variants appear in the table
4. Wait for the **"Ready"** status bar message (sequence caching)
5. Click a variant in the table to preview the sequence context
6. Click **Design Primers** — results appear after a moment
7. Click **Export Report** — an HTML report is generated

---

## Evaluation artifacts

For reference and reproducibility (e.g. for use in scientific publications), the repository ships with the example input files together with the actual logs and HTML reports produced by running them through the application. **All runs were performed using the default settings shipped in `config.json` — no parameters were tuned.**

| Input | Generated report | Run log |
|-------|------------------|---------|
| [`data/example_1.csv`](data/example_1.csv) | [`reports/primer_design_report.html`](reports/primer_design_report.html) | [`logs/primer_designer_20260425_204300.log`](logs/primer_designer_20260425_204300.log) |
| [`data/example_2.csv`](data/example_2.csv) | [`reports/primer_design_report_2.html`](reports/primer_design_report_2.html) | [`logs/primer_designer_20260425_204300.log`](logs/primer_designer_20260425_204300.log) |
| [`data/example_3.csv`](data/example_3.csv) | [`reports/primer_design_report_3.html`](reports/primer_design_report_3.html) | [`logs/primer_designer_20260425_204300.log`](logs/primer_designer_20260425_204300.log) |

- **Example variant inputs** are in [`data/`](data/).
- **Generated reports** (HTML, openable in any browser) are in [`reports/`](reports/).
- **Run log** (full pipeline trace: validation, coordinate mapping, Primer3 design, BLAST homology) is in [`logs/`](logs/).
- **Unit tests** are in [`tests/`](tests/) and can be executed with `venv/bin/python -m pytest tests/`.

---

## Configuration

All settings are stored in `config.json` at the project root and can be edited directly.

<details><summary>Primer3 parameters</summary>

```json
"primer3": {
  "primer_min_size": 18,
  "primer_max_size": 27,
  "primer_opt_size": 22,
  "primer_min_tm": 57.0,
  "primer_max_tm": 64.0,
  "primer_opt_tm": 60.0,
  "primer_min_gc": 40.0,
  "primer_max_gc": 60.0,
  "product_min_size": 100,
  "product_max_size": 500
}
```

</details>

<details><summary>Population variant filtering</summary>

```json
"population": {
  "available_populations": ["global", "afr", "ami", "amr", "asj", "eas", "fin", "nfe", "sas", "oth"],
  "default_maf_threshold": 0.005
}
```

</details>

<details><summary>BLAST+ / homology analysis</summary>

```json
"blast": {
  "blastn_path": "blastn",
  "blast_db_path": "",
  "threads": 4,
  "evalue": 1e-10,
  "min_percent_identity": 85.0,
  "min_aligned_length": 50
}
```

`blast_db_path` is filled in automatically after the first database download. To use a pre-existing BLAST database, set the path manually.

</details>

<details><summary>API & general settings</summary>

```json
"api": {
  "ncbi_api_key": "",
  "ncbi_email": "",
  "request_timeout": 30,
  "max_retries": 3,
  "retry_delay": 1.0
},
"default_assembly": "GRCh38",
"min_distance_from_exon_junction": 60,
"min_distance_from_variant": 40,
"log_level": "INFO"
```

An NCBI API key is optional but raises the rate limit from 3 to 10 requests/second. Providing an email (`ncbi_email`) is recommended per NCBI policy.

</details>

---

## Manual BLAST+ installation (if auto-install fails)

| Platform | Command |
|----------|---------|
| macOS | `brew install blast` |
| Ubuntu/Debian | `sudo apt install ncbi-blast+` |
| Windows | Download from [NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and add to PATH, or set `blastn_path` in `config.json` |

---

## Project structure

```
MultiPrimer Designer/
├── main.py                     # Entry point
├── install.py                  # Cross-platform installer
├── requirements.txt
├── config.json                 # Runtime configuration
├── run_primer_designer.sh      # Launch script (macOS/Linux)
├── src/
│   ├── core/                   # Variant models, validation, coordinate mapping, grouping
│   ├── api/                    # NCBI, Ensembl, gnomAD, MANE, PrimerBLAST clients
│   ├── primer/                 # Primer3 designer, BLAST+ homology, BWA specificity
│   ├── gui/                    # tkinter main window, dark theme, sequence visualization
│   └── utils/                  # Config, logging, file parsing, HTML report generator
├── tests/
├── data/                       # Example CSVs, API cache, genome BLAST database
└── logs/
```

---

## License

MIT — see [LICENSE](LICENSE).
