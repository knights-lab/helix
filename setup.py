"""
Copyright 2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from setuptools import setup, find_packages

__author__ = "Knights Lab"
__copyright__ = "Copyright (c) 2020 --, %s" % __author__
__credits__ = ["Benjamin Hillmann"]
__email__ = "hillm096@cs.umn.edu"
__license__ = "AGPL"
__maintainer__ = "Benjamin Hillmann"

long_description = ''

setup(
    name='helix',
    # version=versioneer.get_version(),
    # cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(exclude=[]),
    url='',
    license=__license__,
    author=__author__,
    author_email=__email__,
    description='',
    long_description=long_description,
    keywords='',
    install_requires=[],
    entry_points={
        'console_scripts': []
    },
    extras_require={
        'develop': ['ipython', 'ipdb'],
        'demo': ['jupyter', 'jupyter_client', 'ipython'],
        'doc': ['sphinx'],
    },
)
