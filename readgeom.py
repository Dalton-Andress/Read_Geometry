#!/usr/bin/env python3
"""
Molecular File Coordinate Extractor

This script extracts molecular coordinates from various quantum chemistry files:

Gaussian files:
- Input files (.com): Extracts coordinates from the molecular specification section.
- Log files (.log): Extracts coordinates from the punch string (if calculation terminated 
  normally) or from the last "Standard orientation" section found in the file.

MOLPRO files:
- Input files (.in) and output files (.out): Extracts coordinates from geometry 
  blocks, typically defined by `geom={...}` or `geometry={...}`.

For Gaussian log files, coordinates from the punch string appear after the third 
double backslash (\\\\) in the archive entry and are only extracted if "Normal termination" 
is found on (or near) the last line of the file.

"""

# Standard library imports
import argparse  # For parsing command-line arguments
import re        # For regular expression operations
import sys       # For system-specific parameters and functions (like stderr, exit)

# Third-party or specific imports (though dataclasses is standard from 3.7+)
from dataclasses import dataclass # For creating simple classes to store data
from pathlib import Path        # For object-oriented filesystem paths
from typing import List, Optional, Dict # For type hinting

# Define a data class to hold atomic coordinate information.
# This improves readability and type safety when passing coordinate data.
@dataclass
class AtomCoordinate:
    """Represents an atom with its element symbol and 3D coordinates (x, y, z)."""
    element: str  # Element symbol (e.g., 'H', 'C')
    x: float      # X-coordinate
    y: float      # Y-coordinate
    z: float      # Z-coordinate

class MolecularFileParser:
    """
    Parser for quantum chemistry files to extract molecular coordinates.
    Supports:
    - Gaussian input files (.com) and log files (.log)
    - MOLPRO input files (.in) and output files (.out)
    """
    
    # A mapping of atomic numbers to their corresponding element symbols.
    # Used to convert atomic numbers found in some output files (e.g., Gaussian std orient)
    # to standard element symbols.
    ELEMENT_MAP: Dict[int, str] = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
        9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
        16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
        23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
        30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
        37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
        44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
        51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
        58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
        65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
        72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
        79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
        86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
        93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',
        100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db',
        106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg',
        112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts',
        118: 'Og'
    }
    
    def __init__(self, debug: bool = False):
        """
        Initialize the parser.
        Args:
            debug (bool): If True, enables printing of debug messages to stderr.
        """
        self.debug = debug # Flag to control debug output
    
    def _debug_print(self, message: str) -> None:
        """Helper method to print messages only if debug mode is enabled."""
        if self.debug:
            print(message, file=sys.stderr) # Debug messages go to standard error
    
    def _read_file_content(self, file_path: Path) -> Optional[str]:
        """
        Reads the entire content of a file.
        Args:
            file_path (Path): The path to the file to be read.
        Returns:
            Optional[str]: The content of the file as a string, or None if an error occurs.
        """
        try:
            # Attempt to read the file with UTF-8 encoding (common for text files)
            return file_path.read_text(encoding='utf-8')
        except FileNotFoundError:
            # Handle cases where the file does not exist
            self._debug_print(f"Error: File '{file_path}' not found during read attempt.")
            return None
        except Exception as e:
            # Handle other potential errors during file reading (e.g., permissions)
            self._debug_print(f"Error reading file '{file_path}': {e}")
            return None
    
    def _check_normal_termination(self, file_path: Path) -> bool:
        """
        Checks if a Gaussian log file indicates normal termination.
        This is done by looking for "Normal termination" in the last few lines of the file.
        Args:
            file_path (Path): The path to the Gaussian log file.
        Returns:
            bool: True if "Normal termination" is found, False otherwise.
        """
        try:
            # Open the file in read mode
            with file_path.open('r', encoding='utf-8') as file:
                file.seek(0, 2)  # Move the file pointer to the end of the file
                file_size = file.tell() # Get the size of the file
                
                # Define a buffer size to read only a chunk from the end of the file,
                # as the termination message is expected there.
                buffer_size = min(1000, file_size) # Read at most 1000 bytes or file size
                file.seek(max(0, file_size - buffer_size)) # Seek backwards from the end
                content_at_end = file.read() # Read the tail end of the file
                
                # Process the lines from the read chunk
                lines = content_at_end.strip().split('\n')
                last_line = lines[-1] if lines else "" # Get the very last non-empty line
                
                # Check if the specific termination string is present
                return "Normal termination" in last_line
        except Exception as e:
            # Handle errors during file operations for termination check
            self._debug_print(f"Error checking termination status for '{file_path}': {e}")
            return False # Assume not normal termination if an error occurs
    
    def _find_punch_string(self, content: str) -> Optional[str]:
        """
        Finds and extracts the complete punch string (archive entry) from Gaussian log file content.
        The punch string is typically located before "The archive entry for this job was punched."
        Args:
            content (str): The full content of the Gaussian log file.
        Returns:
            Optional[str]: The extracted punch string (concatenated lines), or None if not found.
        """
        # Marker indicating the end of the punch string section
        punch_marker = "The archive entry for this job was punched."
        if punch_marker not in content:
            self._debug_print("Info: Punch marker not found in file content.")
            return None # Punch string section not present
        
        punch_end_index = content.find(punch_marker) # Find the position of the marker
        # Consider only the content before this marker for finding the punch string
        lines_before_marker = content[:punch_end_index].strip().split('\n')
        
        punch_lines: List[str] = [] # To store lines belonging to the punch string
        found_punch_start_delimiter = False # Flag to indicate we are inside the punch string
        
        # Iterate backwards through the lines before the marker
        # The punch string lines are characterized by containing backslashes.
        for line in reversed(lines_before_marker):
            line = line.strip() # Remove leading/trailing whitespace
            if not line: # If an empty line is encountered
                if found_punch_start_delimiter:
                    # If we were already collecting punch lines, an empty line signifies the start (when reading backwards)
                    break 
                continue # Otherwise, skip multiple empty lines before the actual punch string start
                
            if '\\' in line: # Punch string lines typically contain backslashes
                punch_lines.append(line)
                found_punch_start_delimiter = True
            elif found_punch_start_delimiter:
                # If we were collecting punch lines and encounter a line without a backslash,
                # it means we've reached above the punch string block.
                break 
        
        if not punch_lines:
            self._debug_print("Info: No punch lines with backslashes found before marker.")
            return None # No valid punch lines found
        
        punch_lines.reverse() # Reverse to restore the original order of lines
        return ''.join(punch_lines) # Concatenate all punch lines into a single string
    
    def _get_element_symbol(self, atomic_number: int) -> str:
        """
        Converts an atomic number to its element symbol using the ELEMENT_MAP.
        If the atomic number is not in the map, it returns the number as a string.
        Args:
            atomic_number (int): The atomic number to convert.
        Returns:
            str: The element symbol or the atomic number as a string if not found.
        """
        return self.ELEMENT_MAP.get(atomic_number, str(atomic_number))
    
    def _parse_coordinate_line(self, line: str) -> Optional[AtomCoordinate]:
        """
        Parses a single line from a Gaussian "Standard orientation" block to extract atomic coordinates.
        Expected format: CenterNumber AtomicNumber AtomicType X Y Z
        Example:          1          6           0      0.000000    0.000000    0.000000
        Args:
            line (str): The line to parse.
        Returns:
            Optional[AtomCoordinate]: An AtomCoordinate object if parsing is successful, None otherwise.
        """
        parts = line.split() # Split the line by whitespace
        # Check if the line has enough parts (at least 6 for the expected format)
        if len(parts) >= 6: 
            try:
                # Extract atomic number and coordinates based on their typical column positions
                atomic_num = int(parts[1])  # Second part is atomic number
                x_coord = float(parts[3])   # Fourth part is X coordinate
                y_coord = float(parts[4])   # Fifth part is Y coordinate
                z_coord = float(parts[5])   # Sixth part is Z coordinate
                
                element_symbol = self._get_element_symbol(atomic_num) # Convert atomic number to symbol
                return AtomCoordinate(element=element_symbol, x=x_coord, y=y_coord, z=z_coord)
            except (ValueError, IndexError) as e: # Handle errors if parts are not numbers or missing
                self._debug_print(
                    f"Warn: Could not parse std orient line '{line[:40]}...': {e}"
                )
        return None # Return None if parsing fails
    
    def _parse_input_coordinate_line(self, line: str) -> Optional[AtomCoordinate]:
        """
        Parses a coordinate line typically found in Gaussian .com or MOLPRO .in/.out files.
        Expected format: ElementSymbol X Y Z
        Example: C 0.000000 0.000000 0.000000
        Args:
            line (str): The line to parse.
        Returns:
            Optional[AtomCoordinate]: An AtomCoordinate object if successful, None otherwise.
        """
        parts = line.split() # Split the line by whitespace
        # Check if the line has enough parts (at least 4 for Element, X, Y, Z)
        if len(parts) >= 4:
            try:
                element_symbol = parts[0].strip() # First part is the element symbol
                # Validate that the element symbol is purely alphabetic (e.g., C, H, Si)
                # This might be too strict for dummy atoms (X, X1) or custom labels.
                if not re.match(r'^[A-Za-z]+$', element_symbol):
                    self._debug_print(f"Warn: Invalid element format '{element_symbol}' in line '{line}'")
                    return None # Skip if element format is not as expected
                
                # Extract coordinates
                x_coord = float(parts[1])
                y_coord = float(parts[2])
                z_coord = float(parts[3])
                return AtomCoordinate(element=element_symbol, x=x_coord, y=y_coord, z=z_coord)
            except (ValueError, IndexError) as e: # Handle conversion or indexing errors
                self._debug_print(
                    f"Warn: Could not parse input line '{line[:40]}...': {e}"
                )
        return None # Return None if parsing fails
    
    def _extract_from_standard_orientation(
            self, content: str
    ) -> List[AtomCoordinate]:
        """
        Extracts coordinates from the *last* "Standard orientation" block in Gaussian log content.
        Args:
            content (str): The full content of the Gaussian log file.
        Returns:
            List[AtomCoordinate]: A list of AtomCoordinate objects, or an empty list if none found/parsed.
        """
        # Regex to find "Standard orientation:" case-insensitively
        pattern = r'(?i)standard\s+orientation:'
        matches = list(re.finditer(pattern, content)) # Find all occurrences
        
        if not matches:
            self._debug_print("Info: No 'Standard orientation' section found.")
            return [] # Return empty list if no such section
        
        self._debug_print(
            f"Found {len(matches)} 'Std orient' section(s), using last."
        )
        
        last_match = matches[-1] # Get the last match object
        section_start_index = last_match.end() # Content starts after the matched pattern
        
        # Extract the part of the content from the start of the last standard orientation section
        section_content = content[section_start_index:]
        lines_in_section = section_content.split('\n') # Split into lines
        
        coordinates_found: List[AtomCoordinate] = []
        # Counter for '----' delimiter lines which structure the standard orientation block
        dash_line_count = 0 
        
        for line in lines_in_section:
            line = line.strip() # Clean the line
            
            # Detect the dashed lines used as separators
            if '-----' in line and len(line) > 20: # A substantial dashed line
                dash_line_count += 1 
                self._debug_print(f"Found dash line {dash_line_count}: {line[:40]}...")
                continue # Move to the next line
            
            if not line: # Skip any empty lines within the section
                continue
            
            # After the first set of dashes (dash_line_count == 1) are column headers
            # like "Center Number Atomic Type Coordinates (Angstroms)"
            if dash_line_count == 1: 
                # Heuristic check for common header keywords to skip these lines
                if any(header_keyword in line for header_keyword in ['Center', 'Number', 'Atomic', 'Coordinates', 'Type']):
                    self._debug_print(f"Skipping header line: {line}")
                    continue
            
            # The actual coordinate data appears after the second set of dashes (dash_line_count == 2)
            if dash_line_count == 2:
                coord = self._parse_coordinate_line(line) # Try to parse the line
                if coord:
                    coordinates_found.append(coord) # Add successfully parsed coordinate
            
            # The third set of dashes usually signifies the end of the coordinate list
            elif dash_line_count >= 3:
                self._debug_print("Found third dash line, stopping std orient parse.")
                break # Stop parsing this section
        
        self._debug_print(
            f"Extracted {len(coordinates_found)} atoms from 'Std orient'."
        )
        return coordinates_found
    
    def _extract_from_punch(self, content: str) -> List[AtomCoordinate]:
        """
        Extracts coordinates from the Gaussian punch string (archive entry).
        The coordinates are typically in the 4th section of the punch string,
        formatted as Element,X,Y,Z entries separated by backslashes.
        Args:
            content (str): The full content of the Gaussian log file.
        Returns:
            List[AtomCoordinate]: A list of AtomCoordinate objects, or an empty list.
        """
        punch_string = self._find_punch_string(content) # Get the raw punch string
        if not punch_string:
            self._debug_print("Info: No punch string found to extract coordinates from.")
            return [] # Return empty if no punch string
            
        self._debug_print("✓ Punch string located for coordinate extraction.")
        
        # Clean the punch string: remove all whitespace, then split by '\\\\'
        # (double backslash is the delimiter between major sections of the punch string)
        punch_no_space = re.sub(r'\s+', '', punch_string)
        punch_sections = punch_no_space.split('\\\\') 
        
        # Expected sections: Version\NImag\HF\Stoichiometry,Charge,Multiplicity\Coordinates\...\@
        # Coordinates are typically in the 4th main section (index 3).
        if len(punch_sections) < 4:
            self._debug_print(f"Error: Insufficient sections ({len(punch_sections)}) in punch string. Expected at least 4.")
            return []
        
        # The coordinate data itself is a string like "C,0.,0.,0.\H,..."
        coordinate_data_string = punch_sections[3] 
        if not coordinate_data_string:
            self._debug_print("Error: Empty coordinate section in punch string.")
            return []
        
        self._debug_print("✓ Punch coordinate section identified.")
        
        # Individual atom entries are separated by a single backslash within this section
        atom_entries_list = coordinate_data_string.split('\\')
        
        coordinates_found: List[AtomCoordinate] = []
        for entry_str in atom_entries_list:
            if not entry_str: continue # Skip if an entry is empty (e.g., due to trailing backslash)
            
            parts = entry_str.split(',') # Split "Element,X,Y,Z"
            if len(parts) == 4: # Expecting exactly 4 parts
                try:
                    element_symbol = parts[0].strip()
                    # Validate element symbol from punch string
                    if not re.match(r'^[A-Za-z]+(?:-[A-Za-z0-9]+)?$', element_symbol): 
                        # Allows for element symbols like 'C', 'H', or 'C-VTZP' if isotopes/types are encoded.
                        # For simplicity, a stricter '^[A-Za-z]+$' is often used if only plain elements are expected.
                        # Current check: ^[A-Za-z]+$ (from _parse_input_coordinate_line context)
                        if not re.match(r'^[A-Za-z]+$', element_symbol):
                            self._debug_print(f"Warn: Invalid element '{element_symbol}' in punch entry '{entry_str}'")
                            continue
                    
                    x_coord = float(parts[1])
                    y_coord = float(parts[2])
                    z_coord = float(parts[3])
                    coordinates_found.append(AtomCoordinate(element_symbol, x_coord, y_coord, z_coord))
                except ValueError as e: # Handle error if X,Y,Z are not valid floats
                    self._debug_print(f"Warn: Could not parse numeric parts of punch entry '{entry_str}': {e}")
            else:
                 self._debug_print(f"Warn: Malformed punch entry '{entry_str}', expected 4 parts (El,X,Y,Z). Found {len(parts)}.")
        return coordinates_found
    
    def _extract_from_input_file(self, content: str) -> List[AtomCoordinate]:
        """
        Extracts coordinates from a Gaussian input file (.com).
        Relies on finding a charge/multiplicity line, then parsing subsequent lines as XYZ coordinates.
        Args:
            content (str): The full content of the Gaussian .com file.
        Returns:
            List[AtomCoordinate]: A list of AtomCoordinate objects.
        """
        lines = content.split('\n') # Split content into lines
        coordinates_found: List[AtomCoordinate] = []
        coord_section_start_line_idx = -1 # Index of the line *after* chg/mult
        
        # Attempt to find the charge and multiplicity line
        # This line typically contains two integers (e.g., "0 1").
        for i, line_content in enumerate(lines):
            current_line_stripped = line_content.strip()
            # A common pattern for charge/multiplicity line: two integers separated by space(s)
            if re.match(r'^[+-]?\d+\s+[+-]?\d+\s*$', current_line_stripped):
                coord_section_start_line_idx = i + 1 # Coordinates start on the next line
                self._debug_print(f"Found chg/mult line at index {i}: '{current_line_stripped}'")
                break
            # Fallback: if after a few lines, we see something that looks like an older format or just numbers.
            # This part of the original code had a 'pass', indicating a stricter reliance on the above regex.
            # For robustness, if the above fails, this could be expanded.
            # Current strictness: requires the two-integer chg/mult line.

        if coord_section_start_line_idx == -1:
            self._debug_print("Info: No standard chg/mult line found. Cannot reliably parse .com coordinates.")
            return [] # If no chg/mult line, parsing is unreliable for this method

        # Parse lines as coordinates starting from after the chg/mult line
        for i in range(coord_section_start_line_idx, len(lines)):
            current_line_stripped = lines[i].strip()
            
            # An empty line signifies the end of the coordinate block in .com files
            if not current_line_stripped: 
                self._debug_print("Empty line encountered after chg/mult, assuming end of coordinate block.")
                break
            
            # Try to parse the line as an "Element X Y Z" coordinate
            coord = self._parse_input_coordinate_line(current_line_stripped)
            if coord:
                coordinates_found.append(coord)
            else:
                # If a line is not empty and not a valid coordinate,
                # it's assumed to be the start of a new section (e.g., basis sets, variables).
                self._debug_print(f"Stopping .com coordinate parsing at non-coordinate line: '{current_line_stripped[:40]}...'")
                break # Stop parsing coordinates
                
        self._debug_print(f"Extracted {len(coordinates_found)} atoms from Gaussian input file.")
        return coordinates_found
  
    def _extract_from_molpro_file(self, content: str) -> List[AtomCoordinate]:
        """
        Extracts coordinates from a MOLPRO file (input .in or output .out).
        Looks for geometry blocks like `geom={...}` or `geometry={...}`.
        Args:
            content (str): The full content of the MOLPRO file.
        Returns:
            List[AtomCoordinate]: A list of AtomCoordinate objects.
        """
        coordinates_found: List[AtomCoordinate] = []
        # Regex to find geometry blocks. It captures the content inside the curly braces.
        # Allows for optional unit specifications (angstrom, bohr, au) before 'geom' or 'geometry'.
        # Handles variations like 'geom = {', 'geometry{', etc.
        # `([^}]*)` captures everything until the first closing brace. DOTALL makes `.` match newlines.
        geom_block_pattern = re.compile(
            r'(?:angstrom|bohr|au)?\s*;?\s*(?:geometry|geom)\s*=\s*\{([^}]*)\}', 
            re.IGNORECASE | re.DOTALL 
        )
        
        # Iterate through all found geometry blocks in the file
        # The logic takes coordinates from the *first* block that successfully yields any coordinates.
        for match in geom_block_pattern.finditer(content):
            block_content_str = match.group(1) # Get the content within the braces
            self._debug_print("Found MOLPRO geometry block via regex.")
            
            geom_lines = block_content_str.strip().split('\n') # Split block content into lines
            current_block_coords: List[AtomCoordinate] = [] # Coords for this specific block
            
            for line_content in geom_lines:
                line_stripped = line_content.strip()
                # Skip empty lines or comment lines (MOLPRO comments often start with !, *, or #)
                if not line_stripped or line_stripped.startswith(('#', '!', '*')): 
                    continue
                # Skip lines that are likely unit declarations or other keywords within the block
                if line_stripped.lower() in ['angstrom', 'bohr', 'au', 'symmetry'] or '=' in line_stripped:
                    self._debug_print(f"Skipping MOLPRO non-coordinate line: {line_stripped}")
                    continue
                
                # Try to parse the line as an "Element X Y Z" coordinate
                coord = self._parse_input_coordinate_line(line_stripped)
                if coord:
                    current_block_coords.append(coord)
                else:
                    # Line within geometry block that doesn't look like a coordinate
                    self._debug_print(f"Warn: MOLPRO line in geom block not parsed as coord: '{line_stripped}'")
            
            # If this block yielded any coordinates, use them and stop searching further blocks
            if current_block_coords:
                coordinates_found = current_block_coords
                self._debug_print(f"Using coordinates from first successfully parsed MOLPRO block ({len(coordinates_found)} atoms).")
                break # Stop after the first successful block
        
        if not coordinates_found:
            self._debug_print("Info: No coordinates extracted from MOLPRO file using primary regex. Check file structure or for empty geometry blocks.")
        return coordinates_found
    
    def _is_gaussian_input_file(self, file_path: Path) -> bool:
        """Checks if the file extension suggests a Gaussian input file."""
        return file_path.suffix.lower() == '.com'
    
    def _is_gaussian_log_file(self, file_path: Path) -> bool:
        """Checks if the file extension suggests a Gaussian log/output file."""
        return file_path.suffix.lower() == '.log'
    
    def _is_molpro_file(self, file_path: Path) -> bool:
        """Checks if the file extension suggests a MOLPRO file (.in for input, .out for output)."""
        return file_path.suffix.lower() in ['.in', '.out'] # Updated MOLPRO extensions
    
    def _detect_file_type(self, file_path: Path) -> str:
        """
        Detects the type of quantum chemistry file based on its extension.
        Returns:
            str: A string identifier for the file type ('gaussian_input', 'gaussian_log', 'molpro', 'unknown').
        """
        if self._is_gaussian_input_file(file_path):
            return 'gaussian_input'
        elif self._is_gaussian_log_file(file_path):
            return 'gaussian_log'
        elif self._is_molpro_file(file_path):
            return 'molpro'
        else:
            return 'unknown' # Default if no specific type is matched
    
    def extract_coordinates(self, file_path: Path) -> List[AtomCoordinate]:
        """
        Main method to extract molecular coordinates from a given quantum chemistry file.
        It detects the file type and calls the appropriate internal extraction method.
        Args:
            file_path (Path): The path to the quantum chemistry file.
        Returns:
            List[AtomCoordinate]: A list of AtomCoordinate objects, or an empty list on failure.
        """
        self._debug_print(f"Starting extraction for: {file_path}")
        content = self._read_file_content(file_path) # Read the file content
        if content is None:
            # Error reading file (e.g., not found, no permissions)
            # _read_file_content would have printed a debug message.
            return [] 
        
        file_type = self._detect_file_type(file_path) # Determine file type
        self._debug_print(f"Detected file type: {file_type} for {file_path}")
        
        coordinates_extracted: List[AtomCoordinate] = [] # Initialize list for results
        
        # Branch based on detected file type
        if file_type == 'gaussian_input':
            coordinates_extracted = self._extract_from_input_file(content)
        elif file_type == 'gaussian_log':
            # For Gaussian logs, check for normal termination to decide extraction strategy
            normal_termination = self._check_normal_termination(file_path)
            if normal_termination:
                self._debug_print("✓ Gaussian log: Normal termination. Trying punch string.")
                coordinates_extracted = self._extract_from_punch(content)
                # Fallback for punch string: if it exists but parsing it yields no coordinates,
                # or if the punch string itself isn't found, try standard orientation.
                if not coordinates_extracted:
                    if self._find_punch_string(content) is not None:
                         self._debug_print("Info: Punch string found but no coords parsed, trying standard orientation as fallback.")
                    else:
                         self._debug_print("Info: Punch string not found, trying standard orientation.")
                    coordinates_extracted = self._extract_from_standard_orientation(content)
            else: # If not normal termination, directly use standard orientation
                self._debug_print("⚠ Gaussian log: Normal termination not found. Using standard orientation.")
                coordinates_extracted = self._extract_from_standard_orientation(content)
        elif file_type == 'molpro':
            coordinates_extracted = self._extract_from_molpro_file(content)
        elif file_type == 'unknown':
            # If the file type is unknown, print a warning (to stderr) and return empty.
            # This message helps users understand why a file might not be processed.
            print(
                f"Warning: Unsupported or unrecognized file type for '{file_path}'. "
                f"Supported types are Gaussian (.com, .log) and MOLPRO (.in, .out).",
                file=sys.stderr
            )
            return [] # Explicitly return empty list for unknown types
            
        # Final debug message about extraction result for this file
        if coordinates_extracted:
            self._debug_print(f"✓ Successfully parsed {len(coordinates_extracted)} atoms from {file_path}.")
        else:
            self._debug_print(f"✗ No coordinates could be extracted from {file_path} interpreted as {file_type}.")
        return coordinates_extracted

class MolecularFormulaGenerator:
    """
    Generates a standard molecular formula string (e.g., CH4, C2H6O) 
    from a list of atomic coordinates.
    """
    
    @staticmethod
    def generate_formula(coordinates: List[AtomCoordinate]) -> str:
        """
        Generates a compact molecular formula.
        Elements are typically ordered C, then H, then alphabetically.
        Counts of 1 are omitted (e.g., CH4 instead of C1H4).
        Args:
            coordinates (List[AtomCoordinate]): A list of AtomCoordinate objects.
        Returns:
            str: The generated molecular formula string.
        """
        if not coordinates: # Handle empty coordinate list
            return ""
            
        atom_count: Dict[str, int] = {} # Dictionary to store counts of each element
        # Count occurrences of each element
        for coord in coordinates:
            element_symbol = coord.element
            # Normalize element symbol (e.g., 'c' -> 'C', 'cl' -> 'Cl')
            if len(element_symbol) > 1: # For multi-character symbols like 'Cl'
                element_symbol = element_symbol[0].upper() + element_symbol[1:].lower()
            else: # For single-character symbols like 'C'
                element_symbol = element_symbol.upper()
            atom_count[element_symbol] = atom_count.get(element_symbol, 0) + 1 # Increment count
        
        # Define a preferred order for common organic elements (C, then H, then others alphabetically)
        # This is a common convention (Hill system for C and H).
        preferred_order = {'C': 0, 'H': 1} 
        # Sort elements: first by preferred order, then alphabetically for the rest
        sorted_elements = sorted(
            atom_count.keys(), 
            key=lambda el_symbol: (preferred_order.get(el_symbol, float('inf')), el_symbol)
        )

        formula_parts: List[str] = [] # List to build the formula string parts
        for element in sorted_elements:
            count = atom_count[element]
            # Append element symbol and count (if count > 1)
            formula_parts.append(f"{element}{count if count > 1 else ''}")
        return ''.join(formula_parts) # Join parts into final formula string

class CoordinateFormatter:
    """
    Formats and prints molecular coordinates and related information to the console.
    Supports different output formats (e.g., table, CSV).
    """
    
    @staticmethod
    def print_coordinates(
        coordinates: List[AtomCoordinate], 
        file_path: Optional[Path] = None, # Original file path for context in output
        show_details: bool = True,        # Flag to control printing of detailed coordinates
        output_format: str = 'table'      # Desired output format ('table' or 'csv')
    ) -> None:
        """
        Prints the molecular coordinates in the specified format.
        Args:
            coordinates: List of AtomCoordinate objects.
            file_path: Optional path of the source file for context.
            show_details: If False (compact mode), only prints summary information.
            output_format: 'table' for human-readable table, 'csv' for CSV format.
        """
        # Handle the case where no coordinates were extracted
        if not coordinates:
            if output_format == 'csv':
                # For CSV, print a commented line indicating no coordinates
                if file_path:
                    print(f"# File: {file_path}, No coordinates found.")
                else: # Should ideally always have file_path if called from main
                    print("# No coordinates found.")
            else: # For 'table' format
                file_name_display = f"Coordinates from {file_path}" if file_path else "Coordinates"
                print(f"{file_name_display} | No coordinates found.")
            return # Exit if no coordinates

        # Generate molecular formula and get atom count for the summary line
        formula = MolecularFormulaGenerator.generate_formula(coordinates)
        num_atoms = len(coordinates)

        # CSV Output Format
        if output_format == 'csv':
            # Print a commented summary line for CSV context
            if file_path:
                print(f"# File: {file_path}, Atoms: {num_atoms}, Formula: {formula}")
            else: # Fallback if file_path is somehow None
                print(f"# Atoms: {num_atoms}, Formula: {formula}")
            
            if show_details: # Only print CSV data if not in compact mode (-c)
                print("Element,X,Y,Z") # CSV Header
                for coord in coordinates:
                    # Format each coordinate line as CSV
                    print(f"{coord.element},{coord.x:.10f},{coord.y:.10f},{coord.z:.10f}")
        
        # Table Output Format (default)
        else: 
            file_name_display = f"Coordinates from {file_path}" if file_path else "Coordinates"
            # Print the summary information line
            print(f"{file_name_display} | Atoms: {num_atoms} | Formula: {formula}")
            
            if show_details: # Only print detailed coordinates if not in compact mode
                print() # Blank line before the coordinate listing
                for coord in coordinates:
                    # Print each coordinate, formatted for table alignment
                    # Element left-aligned in 8 chars, coordinates in 15 chars with 10 decimal places
                    print(
                        f"{coord.element:<8} {coord.x:<15.10f} "
                        f"{coord.y:<15.10f} {coord.z:<15.10f}"
                    )
                print() # Blank line after the coordinate listing

def create_argument_parser() -> argparse.ArgumentParser:
    """
    Creates and configures the command-line argument parser for the script.
    Returns:
        argparse.ArgumentParser: The configured argument parser object.
    """
    parser = argparse.ArgumentParser(
        # Description shown in the help message
        description=(
            "Extract molecular coordinates from quantum chemistry files. "
            "Supports Gaussian (.com/.log) and MOLPRO (.in/.out) files." # Updated MOLPRO extensions
        ),
        # Use RawDescriptionHelpFormatter to preserve formatting of epilog
        formatter_class=argparse.RawDescriptionHelpFormatter,
        # Detailed help text displayed after argument descriptions
        epilog="""
Supported file types:
  Gaussian: .com (input), .log (output)
  MOLPRO:   .in (input), .out (output)

Examples:
  python %(prog)s molecule.com
  python %(prog)s calculation.log --format csv
  python %(prog)s molpro_input.in --compact
  python %(prog)s input.com output.log molpro.out
  python %(prog)s --debug molecule.com
  python %(prog)s -d *.com *.log *.in *.out
  python %(prog)s --compact --format csv *.com
"""
    )
    
    # Define command-line arguments
    # Positional argument: 'files' - one or more file paths or glob patterns
    parser.add_argument(
        'files', 
        nargs='+', # Requires at least one file/pattern
        help='Quantum chemistry file(s) or glob patterns to process'
    )
    # Optional argument: '-d' or '--debug' - enables debug output
    parser.add_argument(
        '-d', '--debug', 
        action='store_true', # Stores True if flag is present
        help='Enable debug output (prints to stderr)'
    )
    # Optional argument: '-c' or '--compact' - enables compact output (summary only)
    parser.add_argument(
        '-c', '--compact', 
        action='store_true',
        help='Show only summary line (no detailed coordinates)'
    )
    # Optional argument: '--format' - specifies output format
    parser.add_argument(
        '--format',
        type=str, # Expected type of the argument value
        choices=['table', 'csv'], # Allowed values for the format
        default='table', # Default value if not specified
        help='Output format for coordinates (default: table)'
    )
    return parser

def main() -> int:
    """
    Main function to drive the script: parses arguments, processes files, and prints results.
    Returns:
        int: Exit code (0 for success, non-zero for errors).
    """
    # Create parser and parse command-line arguments
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Initialize the file parser, passing the debug flag
    file_parser = MolecularFileParser(debug=args.debug)

    # Process input file arguments, which can be direct paths or glob patterns
    file_paths_collected: List[Path] = []
    for pattern_or_file in args.files:
        path_obj = Path(pattern_or_file)
        # First, check if the argument is a direct path to an existing file
        if path_obj.is_file():
            file_paths_collected.append(path_obj)
        else:
            # If not a direct file, treat it as a glob pattern
            # Path().glob(pattern) searches relative to the current working directory by default
            try:
                found_by_glob = [p for p in Path().glob(pattern_or_file) if p.is_file()]
                if found_by_glob:
                    file_paths_collected.extend(found_by_glob)
                # If glob found nothing, and it wasn't a pattern (no wildcards),
                # and the specific path doesn't exist, then it's a missing specific file.
                elif not path_obj.exists() and not any(c in pattern_or_file for c in '*?[]'):
                    # Print a warning for specific files that are not found
                    print(f"Warning: Specified file not found: {pattern_or_file}", file=sys.stderr)
            except re.error as e: # Catch regex errors if glob pattern is malformed
                 print(f"Warning: Invalid glob pattern '{pattern_or_file}': {e}", file=sys.stderr)


    # If no valid files were found after processing all arguments, print error and exit
    if not file_paths_collected:
        print(f"Error: No valid files found matching input: {', '.join(args.files)}", file=sys.stderr)
        return 1 # Exit with error code 1

    # Remove duplicate file paths (e.g., if specified directly and also matched by a glob) and sort
    file_paths = sorted(list(set(file_paths_collected)))

    # Determine if detailed output should be shown based on --compact flag
    show_details = not args.compact
    multiple_files = len(file_paths) > 1 # Flag if multiple files are being processed
    
    # Debug message if processing multiple files
    if multiple_files and args.debug:
        file_parser._debug_print(f"Processing {len(file_paths)} quantum chemistry files:")
        file_parser._debug_print("-" * 50)
    
    # Iterate through each identified file path for processing
    for i, file_path in enumerate(file_paths):
        # Print a blank line before the output of the first file
        if i == 0: 
            print()

        # Debug messages for processing individual files when multiple are present
        if multiple_files and args.debug:
            # resolve() gives the absolute path, useful for clarity if relative paths were used
            file_parser._debug_print(f"\nProcessing: {file_path.resolve()}")
            file_parser._debug_print("-" * 30)
        
        # Extract coordinates from the current file
        coordinates = file_parser.extract_coordinates(file_path)
        
        # Print the extracted coordinates using the specified format
        CoordinateFormatter.print_coordinates(
            coordinates, 
            file_path=file_path, # Pass file_path for context in the output
            show_details=show_details,
            output_format=args.format
        )
        
        # If processing multiple files and using table format, print a splitter
        if multiple_files and args.format == 'table':
            if i < len(file_paths) - 1: # If it's not the last file
                 print("---") # Print the splitter line
                 if show_details: # Only print a blank line after splitter if not in compact mode
                     print() 
    return 0 # Successful execution

# Standard Python idiom to run the main function when the script is executed
if __name__ == "__main__":
    try:
        return_code = main() # Call the main execution function
        sys.exit(return_code) # Exit with the code returned by main
    except KeyboardInterrupt:
        # Handle Ctrl+C gracefully
        print("\nExtraction interrupted by user.", file=sys.stderr)
        sys.exit(130) # Standard exit code for process terminated by Ctrl+C
    except Exception as e:
        # Catch any other unexpected exceptions during script execution
        print(f"An unexpected error occurred in script execution: {e}", file=sys.stderr)
        # If debug mode was likely intended (simple check of args), print full traceback
        if '--debug' in sys.argv or '-d' in sys.argv:
            import traceback
            traceback.print_exc(file=sys.stderr)
        sys.exit(2) # General error exit code for unexpected issues
        