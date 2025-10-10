pip3 install -r requirements.txt
git clone https://github.com/BioinfoUNIBA/REDItools2.git
mkdir tools
mv REDItools2 tools
pip install --use-pep517 -r tools/REDItools2/requirements.txt