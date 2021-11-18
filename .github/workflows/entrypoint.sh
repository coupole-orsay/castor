#!/bin/sh
export TEXMFHOME=/data/texmfhome
tlmgr repository add ftp://tug.org/historic/systems/texlive/2019/tlnet-final
tlmgr init-usertree
tlmgr --usermode install mdframed zref needspace libertine titling
pandoc --pdf-engine=xelatex doc/doc_TP_coupole.md -o doc/doc_TP_coupole.pdf
