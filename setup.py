#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE

from setuptools import setup, find_packages

setup(
    version="0.4.1-dev",
    packages=find_packages(),
    install_requires=[
        'h5py',
        'ctapipe~=0.17.0',
        'ctaplot~=0.6.2',
        'gammapy~=1.0',
        'pyirf~=0.8.0',
        'jupyterlab~=3.5.2',
        'scikit-learn~=1.2.2',
        'astroquery~=0.4.6',
        'seaborn~=0.13.0',
        'tables~=3.8.0',
        'pydantic>=1.4,<2.0',
    ],
    package_data={
        'sst1mpipe': [
            'data/*',
        ],
    },
    extras_require={
        "docs": [
            "sphinx",
            "sphinx-automodapi",
            "sphinx_argparse",
            "sphinx_rtd_theme",
            "numpydoc",
            "nbsphinx",
            "sphinxcontrib-mermaid",
            "sphinx-togglebutton",
            "CTSRepo @ git+https://github.com/cta-sst-1m/CTS.git"
            ],
    },
    tests_require=[
        'pytest',
    ],
    entry_points={
        'console_scripts': [
            'sst1mpipe_dl1_dl2 = sst1mpipe.scripts.sst1mpipe_dl1_dl2:main',
            'sst1mpipe_data_dl1_dl1_stereo = sst1mpipe.scripts.sst1mpipe_data_dl1_dl1_stereo:main',
            'sst1mpipe_data_get_photon_list = sst1mpipe.scripts.sst1mpipe_data_get_photon_list:main',
            'sst1mpipe_data_dl2_dl3 = sst1mpipe.scripts.sst1mpipe_data_dl2_dl3:main',
            'sst1mpipe_get_tunning_params = sst1mpipe.scripts.sst1mpipe_get_tunning_params:main',
            'sst1mpipe_mc_performance = sst1mpipe.scripts.sst1mpipe_mc_performance:main',
            'sst1mpipe_mc_make_irfs = sst1mpipe.scripts.sst1mpipe_mc_make_irfs:main',
            'sst1mpipe_r0_dl1 = sst1mpipe.scripts.sst1mpipe_r0_dl1:main',
            'sst1mpipe_mc_train_rfs = sst1mpipe.scripts.sst1mpipe_mc_train_rfs:main',
            'sst1mpipe_merge_hdf5 = sst1mpipe.scripts.sst1mpipe_merge_hdf5:main',
            'sst1mpipe_pointing_check = sst1mpipe.scripts.sst1mpipe_pointing_check:main',
            'rate_scans_extraction = sst1mpipe.scripts.rate_scans_extraction:main',
            'create_hdu_indexes = sst1mpipe.scripts.create_hdu_indexes:main',
            'sst1mpipe_extract_dl1_distributions = sst1mpipe.scripts.sst1mpipe_extract_dl1_distributions:main'
        ]
    }
)
