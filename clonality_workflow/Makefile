.ONESHELL:

SHELL = /bin/bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

build-clonalityParsers-env:
	source ~/.bashrc
	conda env create -f workflow/envs/clonalityParsers_env.yml
	$(CONDA_ACTIVATE) clonalityParsers
	Rscript -e 'remotes::install_github("https://github.com/pawel125/clonalityParsers")'

update-clonalityParsers-env:
	$(CONDA_ACTIVATE) clonalityParsers
	Rscript -e 'remotes::install_github("https://github.com/pawel125/clonalityParsers")'

verify-clonalityParsers-env:
	$(CONDA_ACTIVATE) clonalityParsers
	Rscript -e 'library(clonalityParsers)'
