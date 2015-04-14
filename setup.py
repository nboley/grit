"""
Copyright (c) 2011-2015 Nathan Boley

This file is part of GRIT.

GRIT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GRIT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GRIT.  If not, see <http://www.gnu.org/licenses/>.
from setuptools import setup, Extension, find_packages
"""

from setuptools import setup, Extension

try:
    from Cython.Setup import cythonize
    extensions = cythonize([
        Extension("grit.sparsify_support_fns", 
                  ["grit/sparsify_support_fns.pyx", ]),
        Extension("grit.call_peaks_support_fns", 
                  ["grit/call_peaks_support_fns.pyx", ])
    ])
except ImportError:
    extensions = [
        Extension("grit.sparsify_support_fns", 
                  ["grit/sparsify_support_fns.c", ]),
        Extension("grit.call_peaks_support_fns", 
                  ["grit/call_peaks_support_fns.c", ])
    ]

config = {
    'include_package_data': True,
    'ext_modules': extensions, 
    'description': 'GRIT',
    'author': 'Nathan Boley',
    'url': 'http://grit-bio.org/',
    'download_url': 'http://grit-bio.org/git/',
    'author_email': 'npboley@gmail.com',
    'version': '2.0.2',
    'packages': ['grit', 
                 'grit.analyze', 
                 'grit.files', 
                 'grit.lib', 
                 'grit.proteomics'],
    'setup_requires': [],
    'install_requires': [ 'scipy', 'numpy', 'networkx', 'pysam' ],
    'scripts': ['./bin/run_grit', "./bin/call_peaks"],
    'name': 'GRIT'
}

if __name__== '__main__':
    setup(**config)
