cd ../random-generated-models
python tacasConnectedRandomGenerator.py
python tacasTreeRandomGenerator.py
cd ../sasha-tacas-scripts
python tacasRandomTree.py run && python tacasRandomTree.py read && python tacasRandomTree.py analyse
python tacasRandomTreeConnected.py run && python tacasRandomTreeConnected.py read && python tacasRandomTreeConnected.py analyse