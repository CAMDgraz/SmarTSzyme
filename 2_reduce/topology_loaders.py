#!/bin/python
# -*- coding: utf-8 -*-
"""
@author: Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
         Pedro A. Sánchez-Murcia [pedro.murcia@medunigraz.at]
"""

from collections import OrderedDict
import re
import os

def load_top(topology: str) -> OrderedDict:
    """
    Load topology file

    Parameters
    ----------
    topology : str
        Path to topology file

    Returns
    -------
        : Topology
        Topology objec
            : Topology
            Topology objectt
    """
    # Test if the file exist
    try:
        os.path.isfile(topology)
    except:
        raise OSError(f'Error!!! File doest not exist or is not a file')

    return get_flags(topology)

# ==============================================================================
# Topologies
# ==============================================================================

def get_flags(top: str) -> OrderedDict:
    """
    Retrieve the information per flag in the topology file.

    Returns
    -------
    found_flags : OrderedDict
        Dictionary containing {(FLAG, (first_line, nrows, format))}
    """
    flags_name = []
    flags_data = OrderedDict()
    flags_format = []
    ignore = False

    with open(top, 'r') as parmfile:
        for line in parmfile:
            if line[0] == '%':
                if line.startswith('%FLAG'):
                    tag, flag = line.split()
                    ignore = flag in {'TITLE'}
                    if not ignore:
                        flags_name.append(flag)
                elif line.startswith('%FORMAT'):
                    if not ignore:
                        fformat = line.strip()[8:-1]
                        flags_format.append(fformat)
                elif line.startswith('%COMMENT'):
                    continue
        
            elif not ignore:
                pyformat = fortran_format(flags_format[-1])
                value_width = pyformat['width']
                values = []
                for i in range(0, pyformat['count']):
                    values_ = line[i*value_width:(i + 1)*value_width].strip()
                    values.append(values_)

                # Filter out blank spaces and covert to data type
                values = list(filter(None, values))
                values = [pyformat['type'](x) for x in values]
                try:
                    flags_data[flags_name[-1]].extend(values)
                except KeyError:
                    flags_data[flags_name[-1]] = values

    return flags_data 

def fortran_format(fformat: str) -> OrderedDict:
    """
    Convert fortran format to python dtype

    Parameters
    ----------
    format : str
        Frotran format

    Returns
    -------
    pytype : str
        Python dtype
    """
    separators = re.findall('\\D', fformat)

    if separators[0] == 'A' or separators[0] == 'a':
        pytype = str
    elif separators[0] == 'F' or separators[0] == 'E':
        pytype = float
    elif separators[0] == 'I':
        pytype = int
    else:
        raise NameError(f'Fortran format type "{separators[0]}"'
                        ' is not recognized"')
    if len(separators) == 1:
        splitted = re.split('\\D', fformat)
        count, width = splitted[0], splitted[1]
        pyformat = OrderedDict([('count', int(count)),
                                ('width', int(width)),
                                ('type', pytype)])
    elif len(separators) == 2:
        splitted = re.split('\\D', fformat)
        count, width, decimal = splitted[0], splitted[1], splitted[2]
        pyformat = OrderedDict([('count', int(count)),
                                ('width', int(width)),
                                ('decimal', int(decimal)),
                                ('type', pytype)])
    else:
        raise IndexError('Number of parameters in the fortran format'
                         f' "{fformat}" coul not be handled')
    return pyformat
