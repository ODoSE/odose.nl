 
1. Requirements:
	sudo apt-get install \
		python-dev \
		python-setuptools \
		python-mysqldb \
		python-networkx \
		python-nose \
		python-lxml \
		sloccount \
		muscle;
	sudo easy_install poster;  # For Life Science Grid Portal 
	sudo easy_install coverage;
	sudo easy_install httplib2;
	sudo easy_install pylint;
		
	#BioPython
	sudo apt-get build-dep python-biopython;
	#Download & install 1.54 <= biopython

	#Other software
	MCL
	OrthoMCL
	NCBIBlast+
	TranslatorX
	PAML
	PHYLIP
	Galaxy
	
	#Galaxy plotting requirements
	sudo apt-get install \
	   python-gnuplot \
	   python-rpy;
	   