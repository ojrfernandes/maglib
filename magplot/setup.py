from setuptools import setup, find_packages

setup(
    name="magplot",
    version="1.0",
    packages=find_packages(),
    install_requires=["matplotlib", "numpy"],
    author="José Fernandes Jr.",
    description="A plotting library for maglib outputs",
)