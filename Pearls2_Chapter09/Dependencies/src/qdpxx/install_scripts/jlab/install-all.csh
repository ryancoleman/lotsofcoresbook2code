#!/bin/tcsh

ssh -n qcdi01    "cd qcd/qdp++; ./install.csh"
ssh -n qcdi02    "cd qcd/qdp++; ./install.csh"
