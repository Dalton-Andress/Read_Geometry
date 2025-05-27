# Molecular Coordinate Extractor

A Python script for extracting molecular coordinates from various quantum chemistry file formats. It's designed to be flexible and provide output in human-readable or CSV formats.

## Description

This script parses output and input files from common quantum chemistry software packages to extract atomic coordinates (element symbol, X, Y, Z). It can handle files from Gaussian and MOLPRO, automatically detecting the file type based on its extension.

Key functionalities include:
* Extraction from Gaussian `.com` (input) and `.log` (output) files.
* Extraction from MOLPRO `.in` (input) and `.out` (output) files.
* Automatic detection of normal termination in Gaussian log files to prioritize punch string coordinates, with fallback to standard orientation.
* Generation of a simple molecular formula from extracted coordinates.
* Multiple output formats: a clean table format (default) and CSV.
* Option for compact output (summary line only).
* Support for processing multiple files and glob patterns.
* Debug mode for verbose output during parsing.

## Features

* **Multi-format Support**: Handles Gaussian and MOLPRO input/output files.  
* **Intelligent Gaussian Log Parsing**: Checks for normal termination to decide between punch string or standard orientation coordinates.  
* **Molecular Formula Generation**: Automatically calculates and displays a molecular formula.
* **Flexible Output**:
    * Clean, human-readable table format.  
    * CSV (Comma-Separated Values) format for easy import into spreadsheets or other programs.  
* **Compact Mode**: Option to display only a summary line (file, atom count, formula) without detailed coordinates.  
* **Batch Processing**: Process multiple files at once, including support for glob patterns (e.g., `*.log`).
* **Debug Mode**: Provides detailed insight into the parsing process for troubleshooting.  

## Supported File Types

* **Gaussian**:
    * `.com` (input files)  
    * `.log` (output/log files)  
* **MOLPRO**:
    * `.in` (input files)  
    * `.out` (output files)  

## Prerequisites

* Python 3.7+ (due to usage of dataclasses and modern type hinting, though `pathlib` is standard from 3.4+)

## Installation

1.  Ensure you have Python 3.7 or newer installed.
2.  Download the `ReadGeom.py` script to your local machine.
3.  Make the script executable (optional, but convenient):
    ```bash
    chmod +x ReadGeom.py
    ```

## Usage

The script is run from the command line.

```bash
python ReadGeom.py [OPTIONS] FILE_OR_PATTERN [FILE_OR_PATTERN ...]
```
Or, if made executable:
```bash
./ReadGeom.py [OPTIONS] FILE_OR_PATTERN [FILE_OR_PATTERN ...]
```

### Command-Line Options

```
usage: ReadGeom.py [-h] [-d] [-c] [--format {table,csv}]
                              FILE_OR_PATTERN [FILE_OR_PATTERN ...]

Extract molecular coordinates from quantum chemistry files. Supports Gaussian
(.com/.log) and MOLPRO (.in/.out) files.

positional arguments:
  FILE_OR_PATTERN       Quantum chemistry file(s) or glob patterns to process

options:
  -h, --help            show this help message and exit
  -d, --debug           Enable debug output (prints to stderr)
  -c, --compact         Show only summary line (no detailed coordinates)
  --format {table,csv}  Output format for coordinates (default: table)

Supported file types:
  Gaussian: .com (input), .log (output)
  MOLPRO:   .in (input), .out (output)

Examples:
  python ReadGeom.py molecule.com
  python ReadGeom.py calculation.log --format csv
  python ReadGeom.py molpro_input.in --compact
  python ReadGeom.py input.com output.log molpro.out
  python ReadGeom.py --debug molecule.com
  python ReadGeom.py -d *.com *.log *.in *.out
  python ReadGeom.py --compact --format csv *.com
```

## Examples

1.  **Extract coordinates from a single Gaussian log file (default table format):**
    ```bash
    python ReadGeom.py my_calc.log
    ```

2.  **Extract coordinates from multiple MOLPRO output files into CSV format:**
    ```bash
    python ReadGeom.py job1.out job2.out --format csv > coordinates.csv
    ```

3.  **Process all `.com` files in the current directory and show compact output:**
    ```bash
    python ReadGeom.py *.com -c
    ```

4.  **Enable debug mode for a specific file:**
    ```bash
    python ReadGeom.py problematic_file.in -d
    ```

## Output Formats

### Table Format (Default)

When using the default table format (`--format table` or if no format is specified), the output for each file will be:

```
Coordinates from path/to/your_file.log | Atoms: N | Formula: XnYm...

Element      X                Y                Z              
H            0.0000000000     0.0000000000     0.7689000000   
C            0.0000000000     0.0000000000     -0.3250000000  
...

```
* A blank line precedes the first file's output.
* If `--compact` (`-c`) is used, only the summary line (e.g., `Coordinates from ... | Atoms: ... | Formula: ...`) is printed.
* When processing multiple files, a `---` splitter (followed by a blank line if not in compact mode) separates the output of each file.

### CSV Format

When using `--format csv`, the output will be:

```csv
# File: path/to/your_file.log, Atoms: N, Formula: XnYm...
Element,X,Y,Z
H,0.0000000000,0.0000000000,0.7689000000
C,0.0000000000,0.0000000000,-0.3250000000
...
```
* The summary line is a comment (`#`).
* If `--compact` (`-c`) is used, only the commented summary line is printed.

## Error Handling and Debugging

* The script attempts to handle common file errors (e.g., file not found, unreadable files) gracefully.
* Warnings for unsupported file types or issues during parsing are printed to `stderr`.
* Using the `-d` or `--debug` flag will print detailed parsing step information to `stderr`, which can be helpful for diagnosing issues with specific files.
* The script exits with code `0` on success, `1` if no valid input files are found, and `2` for other unexpected errors.

## License

This project is licensed under the MIT License.
