#setup.py

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gaiadr3_zeropoint", # Replace with your own username
    version="0.0.3",
    author="Pau Ramos",
    author_email="p.ramos@unistra.fr",
    license='LGPLv3+',
    description="eDR3 zero point functions from Lindegren et al. 2020 implemented in Python.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/icc-ub/gaiadr3_zeropoint",
    packages=setuptools.find_packages(),
    install_requires=[ 'numpy' ],
    package_data={'': ['LICENSE', 'MANIFEST.in'],
            'zero_point': ['coefficients/*.txt']},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
