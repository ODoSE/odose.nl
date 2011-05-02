all: sloccount nose pylint

sloccount:
	sloccount --duplicates --wide --details src/divergence/*.py src/divergence/test/*.py > sloccount.sc

nose:
	nosetests --with-xunit --with-coverage --verbose
	coverage xml

pylint:
	export PYTHONPATH=src/; \
	pylint --max-line-length=120 --disable="E0602,W0511" -f parseable --include-ids=y src/divergence/ > pylint.txt || exit 0
