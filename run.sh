#!/bin/bash

echo "Run handin-1 vangaalen"

echo "Creating the plotting directory if it does not exist"
if [ ! -d "plots" ]; then
  echo "Directory does not exist create it!"
  mkdir plots
fi

if [ ! -e Vandermonde.txt ]; then
  echo "Downloading Dataset"
  wget home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt
fi

if [ ! -e 1a.txt ] || [ ! -e poisson_timing.txt ]; then
  echo "Run the script for 1a"
  python3 Poisson.py
fi

if [ ! -e vandermonde_coefficients.txt ] || [ ! -e vandermonde_timing.txt ] || [ ! -e plots/vandermonde.png ] || [ ! -e plots/vandermonde_itt.png ]; then
  echo "Run the script for 2"
  python3 vandermonde.py
fi

echo "Generating the pdf"

#pdflatex template.tex -quiet=true
pdflatex -interaction=nonstopmode NUR1_vangaalen.tex > tex_output.txt
#bibtex NUR1_vangaalen.aux
pdflatex -interaction=nonstopmode NUR1_vangaalen.tex > tex_output.txt
pdflatex -interaction=nonstopmode NUR1_vangaalen.tex > tex_output.txt