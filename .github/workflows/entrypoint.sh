#!/bin/sh
export TEXMFHOME=/data/texmfhome
tug_repository="ftp://tug.org/historic/systems/texlive/2019/tlnet-final"
tlmgr repository add "$tug_repository"
tlmgr option repository "$tug_repository"
tlmgr init-usertree
tlmgr --usermode install mdframed zref needspace libertine titling
pandoc --pdf-engine=xelatex doc/doc_TP_coupole.md -o doc/doc_TP_coupole.pdf
