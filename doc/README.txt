## Run requirements ##
# Galaxy
sudo apt-get install mercurial subversion
sudo apt-get install postgresql-client
# plotting requirements
sudo apt-get install python-gnuplot python-rpy python-rpy2;

# Python packages
sudo apt-get install python-biopython  # For sequence wrangling
sudo apt-get install python-lxml  # For download from MRS
sudo apt-get install python-mysqldb  # For OrthoMCL
sudo apt-get install python-poster  # For Life Science Grid Portal
sudo apt-get install python-networkx  # For drawing Phylo trees

# OrthoMCL
sudo apt-get install libdbd-mysql-perl

# TranslatorX
sudo apt-get install muscle

# Other software
MCL
NCBIBlast+
PAML
PHYLIP


## Build requirements ##
sudo apt-get install \
	python-dev \
	python-setuptools \
	python-nose \
	sloccount;

sudo easy_install coverage;
sudo easy_install pylint;
