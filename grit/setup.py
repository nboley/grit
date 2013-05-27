from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

config = {
    'cmdclass': {'build_ext': build_ext},
    'ext_modules': [
        #Extension("sparsify_support_fns", ["grit/sparsify_support_fns.pyx"]),
        Extension("sparsify_support_fns", ["grit/sparsify_support_fns.c"]),
    ],    
    'description': 'GRIT',
    'author': 'Nathan Boley',
    'url': 'http://grit-bio.org/',
    'download_url': 'http://grit-bio.org/git/',
    'author_email': 'npboley@gmail.com',
    'version': '1.0',
    'packages': ['grit'],
    'requires': ['Cython', 'numpy', 'scipy'],
    'scripts': [],
    'name': 'GRIT'
}

for_ubuntu = """
sudo aptitude install cython  python-dev python-numpy python-scipy python-igraph python-pip python-networkx
easy_install pysam
"""


setup(**config)

