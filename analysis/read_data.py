#!/usr/bin/env python3

from pathlib import Path
import numpy as np


class MassInflationData:
    """Minimal reader for mass_inf ASCII output files."""
    def __init__(self, folder):
        """
        Read mass_inf data output.

        Parameters
        ----------
        folder : str or Path
            Directory containing .dat files to list and read.
        """
        self.folder = Path(folder)
        self.files = self._list_dat_files()

    def _list_dat_files(self):
        """Return sorted list of all .dat files in the given folder (as Path objects)."""
        return sorted(Path(self.folder).glob('*.dat'))

    def get_v(self, filename=None):
        """
        Return the list of v-slice values from a given file, or from the "fields.dat" file if none is specified.

        Parameters
        ----------
        filename : str or None, optional
            Name of the .dat file to read v-slice values from. If None, use the 'fields.dat' file.

        Returns
        -------
        np.array of float
            Sorted array of v-slice values present in the file.
        """
        if filename is None:
            filename = "fields.dat"  # default to 'fields.dat' file

        filepath = self.folder / filename
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")

        with open(filepath) as f:
            lines = f.readlines()

        v_values = []
        for line in lines:
            if line.startswith('# v ='):
                v_val = float(line.split('=')[1].strip())
                v_values.append(v_val)

        return np.array(v_values)

    def get_metadata(self, filename):
        """
        Return the metadata lines (header comments) from the specified .dat file.

        Parameters
        ----------
        filename : str
            Name of the .dat file to extract metadata from.

        Returns
        -------
        list of str
            List with the labels for each column.
        """
        filepath = self.folder / filename
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")

        with open(filepath) as f:
            lines = f.readlines()

        columns = []
        for line in lines[0:10]:  # look in first 10 lines
            if line.startswith('# Columns: '):
                columns = line.split('Columns: ')[1].strip().split(', ')

        return columns

    def read_file(self, filename, col=1):
        """
        Read an output file and extract data

        Parameters
        ----------
        filename : str
            Filename to read (e.g., 'fields.dat')
        col      : int, optional
            Column index to extract (0-based). Default is 1.

        Returns
        -------
        data : ndarray
            Numeric data (skipping comment lines)
        """
        filepath = self.folder / filename

        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")

        data = np.loadtxt(filepath, usecols=col, unpack=True, comments='#')

        return data

    def get_u(self, filename=None):
        """
        Return the array of u values (first column) from the first v-slice in the file.

        Assumes that all v-slices have the same u values (i.e., the u grid is identical for each v block).

        Parameters
        ----------
        filename : str or None, optional
            Name of the .dat file to read from. If None, uses the 'fields.dat' file.

        Returns
        -------
        np.ndarray
            Array of u values from the first v-slice.
        """
        if filename is None:
            filename = "fields.dat"  # default to fields.dat file

        filepath = self.folder / filename
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")

        u_all = self.read_file(filename, col=0)
        v  = self.get_v(filename)
        Nv = len(v)
        Nu = len(u_all) // Nv
        u  = u_all[0:Nu] # take only the first v-slice

        return u

    def get_data(self, filename, col=1):
        """
        Load a data column from a .dat file and reshape it onto the (u, v) grid.

        Parameters
        ----------
        filename : str
            Name of the .dat file to read (e.g., 'fields.dat').
        col : int, optional
            Column index to extract (0-based, default is 1).

        Returns
        -------
        U : ndarray
            2D array of u coordinates (meshgrid, shape (Nv, Nu)).
        V : ndarray
            2D array of v coordinates (meshgrid, shape (Nv, Nu)).
        data_reshaped : ndarray
            2D array of the selected data column, reshaped to (Nv, Nu)

        Notes
        -----
        Note that the returned arrays have shape (Nv, Nu), ie the first index
        corresponds to the v coordinate and the second to u.
        """
        data = self.read_file(filename, col=col)

        u  = self.get_u(filename)
        v  = self.get_v(filename)
        Nu = len(u)
        Nv = len(v)

        data_reshaped = data.reshape(Nv, Nu)

        U, V = np.meshgrid(u, v, indexing='xy')

        return U, V, data_reshaped
