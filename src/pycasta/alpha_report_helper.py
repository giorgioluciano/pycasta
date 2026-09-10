import json
from pathlib import Path
import pandas as pd
from IPython.display import display, Markdown, Image

def show_alpha_report(pdb_name, results_dir, display_pngs=True):
    results = Path(results_dir)
    csvs = sorted(results.glob('*alphasweepmetrics.csv'))
    jsons = sorted(results.glob('*alphadiagnostics.json'))
    if not csvs or not jsons:
        raise FileNotFoundError('Alpha diagnostics artifacts not found')
    csv_path, json_path = csvs[0], jsons[0]
    df = pd.read_csv(csv_path)
    with open(json_path, 'r', encoding='utf-8') as f:
        summary = json.load(f)
    display(Markdown(f"## Alpha report — {pdb_name}"))
    display(Markdown(f"**Base alpha:** {summary.get('basealpha')}  
**Tested alphas:** {summary.get('testedalphas')}"))
    cols = [c for c in ['alpha','numalphatetrahedra','numpreliminarypockets','nummergedpockets','numfinalpockets','meanflowsteps','maxvolume'] if c in df.columns]
    display(df[cols].sort_values('alpha'))
    if display_pngs:
        for key in ['histpng','sweeppng','top3png']:
            p = summary.get('artifacts', {}).get(key)
            if p:
                display(Image(filename=p))
    return df, summary
