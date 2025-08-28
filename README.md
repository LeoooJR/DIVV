# DIVV

<center><img src="logo.png" width="500" height="500"></center><br>

![Static Badge](https://img.shields.io/badge/python-3.10%2B-yellow?label=python&labelColor=grey&color=yellow)
![Static Badge](https://img.shields.io/badge/htslib-1.21-grey?color=brightgreen)
![GitHub last commit](https://img.shields.io/github/last-commit/LeoooJR/DIVV)

## :bookmark_tabs: Table of Contents <a name="table">
- [About](#about)
- [Getting Started](#getting-started)
- [Contacts](#contacts)

## :pencil: About The Project <a name="about">

<p align="justify">
DIVV is a software tool designed for comparing two VCF files that have been generated from panel sequencing.<br>
It leverages the CyVCF2 library, which is actively maintained, providing superior VCF file parsing compared to PyVCF. Additionally, tasks are parallelized to enhance efficiency. Moreover, a dynamic report is generated, offering users an optimized and user-friendly visualization of the results.
</p>

## :rocket: Getting Started <a name="getting-started">

```bash
# Clone the repository
git clone https://github.com/LeoooJR/DIVV.git

# Changed directory to the cloned repository
cd DIVV
```

### :computer: Local installation

```bash
# Compile htslib
./install.sh ./src

# Get help about DIVV
python src/main.py --help

# Compare VCF files
python src/main.py --vcfs <VCF1> <VCF2> --report --output <PATH>
```

### :whale2: Using docker

```bash
# Build the image
docker build -t divv .

# Run the container
docker run --rm divv --vcfs <VCF1> <VCF2> --report --output <PATH>

# Get help about DIVV
docker run --rm divv --help
```

## 🗨️ Contacts <a name="contacts"></a>

For questions, feel free to reach out to the maintainers through GitHub or connect on LinkedIn.

[Back to top](#top)
