 
1. Requirements:
	sudo apt-get install \
		python-dev \
		python-setuptools \
		python-httplib2 \
		python-mysqldb \
		python-nose \
		pylint \
		sloccount \
		muscle;
	sudo easy_install coverage;
	
	#BioPython
	sudo apt-get build-dep python-biopython;
	#Download & install 1.54 <= biopython

	#Other software
	MCL
	OrthoMCL
	NCBIBlast+
	TranslatorX
	PAML
	DoFE
	Galaxy
	
	#Galaxy plotting requirements
	sudo apt-get install \
	   python-gnuplot \
	   python-rpy;
	   