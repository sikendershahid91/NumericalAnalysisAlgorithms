#!/usr/bin/python

from paver.easy import *
import paver.doctools
import os
import glob
import shutil

@task
def run():
    sh('python3 src/kessel_run.py')
    pass

@task
def test():
    sh('nosetests --with-coverage --cover-erase --cover-package=src --cover-html test')
    pass


@task
def clean():
    for pycfile in glob.glob("*/*/*.pyc"): os.remove(pycfile)
    for pycache in glob.glob("*/__pycache__"): os.removedirs(pycache)
    for pycache in glob.glob("./__pycache__"): shutil.rmtree(pycache)
    for cover in glob.glob("./cover"): shutil.rmtree(cover)
    pass


@task
@needs(['clean', 'test', 'run'])
def default():
    pass

