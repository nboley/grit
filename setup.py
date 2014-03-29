#from distutils.core import setup
#from distutils.extension import Extension

from setuptools import setup, Extension, find_packages

try: import Cython
except ImportError: USE_CYTHON = False
else: USE_CYTHON = True

# if we have cython, then build the pyx file
if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize("grit/sparsify_support_fns.pyx")
else:
    extensions = [
        Extension("sparsify_support_fns", ["grit/sparsify_support_fns.c", ])
    ]


config = {
    'include_package_data': True,
    'ext_modules': extensions,
    'description': 'GRIT',
    'author': 'Nathan Boley',
    'url': 'http://grit-bio.org/',
    'download_url': 'http://grit-bio.org/git/',
    'author_email': 'npboley@gmail.com',
    'version': '1.1.2b',
    'packages': ['grit', 
                 'grit.analyze', 
                 'grit.files', 
                 'grit.lib', 
                 'grit.proteomics'],
    'setup_requires': [ 'Cython', ],
    'install_requires': [ 'scipy', 'numpy', 'Cython', 'networkx', 'pysam' ],
    'scripts': ['./bin/run_grit.py',],
    'name': 'GRIT'
}

if __name__== '__main__':
    setup(**config)
