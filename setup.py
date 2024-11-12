import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ScanpyAutoAnalyzer",
    version="0.0.2",
    author="Stefan Lang",
    author_email="Stefan.Lang@med.lu.se",
    description="A simple tool to prepare and run prepared jupyter notebooks.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/stela2502/ScanpyAutoAnalyzer",
    project_urls={
        "Bug Tracker": "https://github.com/stela2502/ScanpyAutoAnalyzer/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts= [ 
        'scripts/ScanpyAnalysis.py',
        'scripts/VelocytoAnalalysis.py'
     ],
    package_dir={"ScanpyAutoAnalyzer": "src/ScanpyAutoAnalyzer"},
    package_data={"ScanpyAutoAnalyzer": [
        "data/ExampleAnalysis.md",
        "data/VelocytoAnalalysis.md"
        ] 
    },
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
