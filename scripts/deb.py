import json
d = json.load(open('c:/github/pycasta/results/coach420_full420/1a26A_pockets.json'))
print(list(d.keys()))
print('volumes:', d.get('pocket_volumes', 'MISSING'))
print('depths:', d.get('pocket_depths', 'MISSING'))