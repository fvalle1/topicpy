import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="topicpy", # Replace with your own username
    version="1.1.1",
    author="Filippo Valle",
    author_email="filippo.valle@unito.it",
    description="Package to extract information from topic models.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fvalle1/topics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
