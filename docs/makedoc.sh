export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

python3 ../code/depthkopy.py i=-1 dochtml newlog

python3 ../code/depthkopy.py --description | sed 's/^/# /' > ../DepthKopy.md
echo >> ../DepthKopy.md

echo '```' >> ../DepthKopy.md
python3 ../code/depthkopy.py --details >> ../DepthKopy.md
echo '```' >> ../DepthKopy.md
echo >> ../DepthKopy.md
echo 'For a better rendering and navigation of this document, please download and open [`./docs/depthkopy.docs.html`](./docs/depthkopy.docs.html), or visit <https://slimsuite.github.io/depthkopy/>.' >> ../DepthKopy.md
echo 'Documentation can also be generated by running DepthKopy with the `dochtml=T` option. (R and pandoc must be installed - see below.)' >> ../DepthKopy.md
echo >> ../DepthKopy.md
echo '## Introduction' >> ../DepthKopy.md
echo >> ../DepthKopy.md
grep -A 10000 Diploidocus depthkopy.docs.Rmd >> ../DepthKopy.md

cp depthkopy.docs.html ../index.html
cp ../DepthKopy.md ../README.md
