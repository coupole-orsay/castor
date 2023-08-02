# CASTOR: Codes pour l’ASTronomie à ORsay

Outils pour l’analyse des données de la coupole d’Orsay.

## Documentation

[Documentation pour les TP à la coupole (PDF)][doc-tp-pdf]

[doc-tp-pdf]: https://github.com/coupole-orsay/castor/releases/latest/download/doc_TP_coupole.pdf

## Installation

### Version stable (depuis PyPI; recommandé)

```bash
mkdir castor && cd castor
python3 -m venv venv
source venv/bin/activate
pip3 install castor-orsay
```


### Version de développement (depuis GitHub)

```bash
git clone https://github.com/coupole-orsay/castor
cd castor
python3 -m venv venv
source venv/bin/activate
pip3 install -e .
```


## Utilisation

```bash
cd castor
source venv/bin/activate
castor_<tab>
```


## Licence

Ces outils sont mis à disposition sous licence MIT. Voir `LICENSE.txt`
