.PHONY: clean virtualenv upgrade test package dev dist

PYENV = . env/bin/activate;
PYTHON = $(PYENV) python3
EXTRAS_REQS := $(wildcard requirements-*.txt)

package: env
	$(PYTHON) setup.py sdist

build: env
	$(PYTHON) setup.py build
	$(PYENV) pip install -e . -r requirements.txt

test: dev
	$(PYTHON) `which nosetests` $(NOSEARGS)
	$(PYENV) py.test README.rst

extras: env/make.extras
env/make.extras: $(EXTRAS_REQS) | env
	rm -rf env/build
	$(PYENV) for req in $?; do pip install -r $$req; done
	touch $@

clean:
	python3 setup.py clean
	rm -rf dist build *.so
	find . -type f -name "*.pyc" -exec rm {} \;

nuke: clean
	rm -rf *.egg *.egg-info env bin cover coverage.xml nosetests.xml

env virtualenv: env/bin/activate
env/bin/activate: requirements.txt setup.py
	test -f $@ || virtualenv -p python3 --no-site-packages env
	$(PYENV) pip install -e . -r $<
	touch $@
