from setuptools import setup, find_packages

setup(
    name="pyCAST",  # 
    version="1.0.0",
    author="Giorgio Luciano, Ulderico Fugacci, Silvia Biasotti",
    author_email="giorgio.luciano@cnr.it"  #
    description="Pocket detection and validation pipeline based on alpha shapes and ligand contact analysis.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/giorgioluciano/pyCAST",  # 
    packages=find_packages(exclude=["tests*", "docs*"]),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "biopython",
        "trimesh",
        "matplotlib",
        "scikit-learn",
        "tabulate",
        "plyfile",
        "openpyxl",  # Per leggere file Excel
        "PyMOL-open-source",  # 
    ],
    python_requires='>=3.8',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Cambia se usi un'altra licenza
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "pycast-run=run_analysis:main",  # Consente di eseguire con `pycast-run`
        ],
    },
)
