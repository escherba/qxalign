import re
import os
import itertools
from setuptools import setup, Extension
# from distutils.core import setup, Extension
from glob import glob


# dependency links
SKIP_RE = re.compile(r'^\s*--find-links\s+(.*)$')

# Regex groups: 0: URL part, 1: package name, 2: package version
EGG_RE = re.compile(r'^(.+)#egg=([a-z0-9_.]+)-([a-z0-9_.-]+)$')

# Regex groups: 0: URL part, 1: package name, 2: branch name
URL_RE = re.compile(r'^\s*(https?://[\w\.]+.*/([^\/]+)/archive/)([^\/]+).zip$')

# our custom way of specifying extra requirements in separate text files
EXTRAS_RE = re.compile(r'.*\bextras\.(\w+)\.txt$')


def parse_reqs(reqs):
    """Parse requirements.txt files into lists of requirements and dependencies
    """
    pkg_reqs = []
    dep_links = []
    for req in reqs:
        # find things like
        # --find-links http://packages.livefyre.com/buildout/packages/
        dep_link_info = re.search(SKIP_RE, req)
        if dep_link_info is not None:
            url = dep_link_info.group(1)
            dep_links.append(url)
            continue
        # add packages of form:
        # git+https://github.com/Livefyre/pymaptools#egg=pymaptools-0.0.3
        egg_info = re.search(EGG_RE, req)
        if egg_info is not None:
            url, egg, version = egg_info.group(0, 2, 3)
            pkg_reqs.append(egg + '==' + version)
            dep_links.append(url)
            continue
        # add packages of form:
        # https://github.com/escherba/matplotlib/archive/qs_fix_build.zip
        zip_info = re.search(URL_RE, req)
        if zip_info is not None:
            url, pkg = zip_info.group(0, 2)
            pkg_reqs.append(pkg)
            dep_links.append(url)
            continue
        pkg_reqs.append(req)
    return pkg_reqs, dep_links


def build_extras(glob_pattern):
    """Generate extras_require mapping
    """
    fnames = glob(glob_pattern)
    result = dict()
    dep_links = []
    for fname in fnames:
        extras_match = re.search(EXTRAS_RE, fname)
        if extras_match is not None:
            extras_file = extras_match.group(0)
            extras_name = extras_match.group(1)
            with open(extras_file, 'r') as fhandle:
                result[extras_name], deps = parse_reqs(fhandle.readlines())
                dep_links.extend(deps)
    return result, dep_links


def resource_string(path):
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), path))
    with open(path, 'r') as fh:
        return fh.read()


INSTALL_REQUIRES, INSTALL_DEPS = parse_reqs(
    resource_string('requirements.txt').splitlines())
TESTS_REQUIRE, TESTS_DEPS = parse_reqs(
    resource_string('requirements-tests.txt').splitlines())
EXTRAS_REQUIRE, EXTRAS_DEPS = build_extras('requirements-extras.*.txt')
DEPENDENCY_LINKS = list(set(itertools.chain(
    INSTALL_DEPS,
    TESTS_DEPS,
    EXTRAS_DEPS
)))


setup(
    ext_modules=[
        Extension("qxalign", sources=["qxalign.c", "align454.c"])
    ],
    name="qxalign",
    author="Eugene Scherba",
    author_email="escherba@bu.edu",
    url="https://github.com/escherba/qxalign3",
    license="BSD",
    version="0.0.1",
    description="Module for quality-aware realignment",
    platforms=["Linux", "Windows", "Mac OS X"],
    install_requires=INSTALL_REQUIRES,
    tests_require=TESTS_REQUIRE,
    dependency_links=DEPENDENCY_LINKS,
    zip_safe=True,
    test_suite='nose.collector',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python",
        "Programming Language :: C",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps."
    ],
    long_description=resource_string('README.rst'),
)
