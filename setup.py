import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="topicpy", # Replace with your own username
    version="0.0.1",
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
    install_requires = ["tensorflow","gseapy","seaborn", "sklearn", "pandarallel","tensorflow_probability", "numpy", "matplotlib", "pandas"],
    python_requires='>=3.6',
)
