#!/bin/bash

echo "Run handin template"

echo "Creating the plotting directory if it does not exist"
if [ ! -d "plots" ]; then
  echo "Directory does not exist create it!"
  mkdir plots
fi

echo "Downloading Dataset"
if [ ! -e Vandermonde.txt ]; then
  wget home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt
fi

if [ ! -e 1a.txt ]; then
  echo "Run the script for 1a"
  python3 Poisson.py
fi

if [ ! -e Vandermonde.txt ]; then
  echo "Run the script for 2"
  python3 vandermonde.py
fi

echo "Generating the pdf"

#pdflatex template.tex -quiet=true
pdflatex -interaction=nonstopmode template.tex > tex_output.txt
#bibtex template.aux
#pdflatex template.tex
#pdflatex template.tex


