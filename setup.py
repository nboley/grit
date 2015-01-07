from setuptools import setup, Extension, find_packages

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
    'version': '2.0.1',
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
