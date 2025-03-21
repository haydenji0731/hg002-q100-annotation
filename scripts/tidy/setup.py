from setuptools import setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
	name="tidy",
	version="0.0.1",
	author="HJ Ji",
	author_email="hji20@jh.edu",
	description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
	url="https://github.com/haydenji0731/hg002-q100-annotation",
	install_requires=[
        'orfipy',
        'pandas',
        'biopython',
        'pyfastx',
        'setuptools'
    ],
	python_requires='>=3.8',
	packages=['tidy'],
	entry_points={'console_scripts': ['tidy = tidy.run_tidy:main'],},
)