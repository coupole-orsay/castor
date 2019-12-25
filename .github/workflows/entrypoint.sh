#!/bin/sh
export TEXMFHOME=/data/texmfhome
tlmgr init-usertree
tlmgr --usermode install mdframed zref needspace libertine titling
pandoc --pdf-engine=xelatex doc/doc_TP_coupole.md -o doc/doc_TP_coupole.pdf
