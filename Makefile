.PHONY: clean doc doc-clean tests check test install build bdist gh-pages
WHEELHOUSE?=./wheelhouse

build:
	python setup.py build

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg

pylint:
	pylint --errors-only pbcoretools/

install:
	@which pip > /dev/null
	@pip freeze|grep 'pbcoretools=='>/dev/null \
      && pip uninstall -y pbcoretools \
      || echo -n ''
	@pip install ./


clean: doc-clean
	rm -rf build/;\
	find . -name "*.egg-info" | xargs rm -rf;\
	rm -rf dist/;\
	find . -name "*.pyc" | xargs rm -f;
	rm -f nosetests.xml

doc:
	sphinx-apidoc -o doc/ pbcoretools/ && cd doc/ && make html
doc-clean:
	cd doc && rm -rf _templates _static _build searchindex.js objects.inv

doctest:
	cd doc && make doctest

DOCTEST_MODULES= pbcoretools/utils.py #pbcoretools/file_utils.py #etc.

quick-test:
	pytest -v --doctest-modules ${DOCTEST_MODULES} #--ignore=.

unit-test:
	pytest -n auto -v --durations=10 --junit-xml=nosetests.xml --cov=./pbcoretools --cov-report=xml:coverage.xml tests/unit/test_*.py
	sed -i -e 's@filename="@filename="./@g' coverage.xml

test: doctest unit-test

tests: test
check: test

GH_PAGES_SOURCES = pbcoretools doc

gh-pages:
	git checkout gh-pages
	rm -rf _static _sources *.js *.html *.inv
	git checkout master $(GH_PAGES_SOURCES)
	cd doc && make html
	mv -fv doc/_build/html/* .
	rm -rf $(GH_PAGES_SOURCES)
	git add --all && git commit -m "Automatic update of gh-pages branch" && git checkout master

pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'pbcoretools=='>/dev/null \
      && pip uninstall -y pbcoretools \
      || echo -n ''
	@pip install --no-index ./

wheel:
	which pip
	pip wheel -v --wheel-dir=${WHEELHOUSE} --no-deps .
	ls -larth ${WHEELHOUSE}
