cd ../random-generated-models
python modelGenerator.py -size 1000 -numMinActions 2 -type random -numModels 10 -smallestProb 0.01
echo "Created Models"
rm ../case-studies/random-models/RANDOM*
mv generatedModels/*.prism ../case-studies/random-models
echo "Moved"
cd ../scripts
rm -rf currentRun
python secondaryScript.py analyse currentRun