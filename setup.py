import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="topicpy",
    version="0.3.0",
    author="Filippo Valle",
    author_email="filippo.valle@unito.it",
    description="Package to extract information from topic models.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fvalle1/topicpy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ["tensorflow","gseapy","seaborn", "scikit-learn","tensorflow_probability", "numpy", "matplotlib", "pandas"],
    python_requires='>=3.6',
)
