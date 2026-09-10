import json
from pathlib import Path

import pandas as pd
from IPython.display import Image, Markdown, display


def load_alpha_artifacts(results_dir):
    results = Path(results_dir)
    csvs = sorted(results.glob('*alphasweepmetrics.csv'))
    jsons = sorted(results.glob('*alphadiagnostics.json'))
    if not csvs or not jsons:
        raise FileNotFoundError('Alpha diagnostics artifacts not found in results_dir')
    csv_path = csvs[0]
    json_path = jsons[0]
    df = pd.read_csv(csv_path)
    with open(json_path, 'r', encoding='utf-8') as f:
        summary = json.load(f)
    return df, summary


def recommend_alpha_from_sweep(df):
    work = df.copy()
    required = ['alpha', 'numfinalpockets', 'maxvolume', 'meanflowsteps']
    missing = [c for c in required if c not in work.columns]
    if missing:
        raise ValueError(f'Missing columns for alpha recommendation: {missing}')
    for col in required:
        work[col] = pd.to_numeric(work[col], errors='coerce').fillna(0)
    work = work.sort_values('alpha').reset_index(drop=True)
    work['vol_norm'] = work['maxvolume'] / max(work['maxvolume'].max(), 1.0)
    work['pocket_norm'] = work['numfinalpockets'] / max(work['numfinalpockets'].max(), 1.0)
    work['flow_norm'] = work['meanflowsteps'] / max(work['meanflowsteps'].max(), 1.0)
    work['delta_volume'] = work['maxvolume'].diff().abs().fillna(0) + work['maxvolume'].diff(-1).abs().fillna(0)
    work['stability'] = 1.0 / (1.0 + work['delta_volume'])
    work['score'] = 0.45 * work['vol_norm'] + 0.20 * work['flow_norm'] + 0.25 * work['stability'] - 0.10 * work['pocket_norm']
    best_idx = work['score'].idxmax()
    best_row = work.loc[best_idx].to_dict()
    return {
        'recommended_alpha': float(best_row['alpha']),
        'score': float(best_row['score']),
        'numfinalpockets': int(best_row['numfinalpockets']),
        'maxvolume': float(best_row['maxvolume']),
        'meanflowsteps': float(best_row['meanflowsteps']),
    }

def recommend_alpha_from_connectivity_sweep(df):
    work = df.copy()

    required = [
        "alpha",
        "alpha_fraction",
        "largest_component_ratio",
        "small_component_fraction",
        "num_alpha_components",
    ]
    missing = [c for c in required if c not in work.columns]
    if missing:
        raise ValueError(f"Missing columns for connectivity alpha recommendation: {missing}")

    for col in required:
        work[col] = pd.to_numeric(work[col], errors="coerce").fillna(0)

    work = work.sort_values("alpha").reset_index(drop=True)

    work["connectivity_score"] = (
        0.60 * work["largest_component_ratio"]
        - 0.20 * work["small_component_fraction"]
        - 0.10 * (work["alpha_fraction"] / max(work["alpha_fraction"].max(), 1e-12))
        - 0.10 * (work["num_alpha_components"] / max(work["num_alpha_components"].max(), 1))
    )

    best_idx = work["connectivity_score"].idxmax()
    best_row = work.loc[best_idx].to_dict()

    return {
        "recommended_alpha": float(best_row["alpha"]),
        "connectivity_score": float(best_row["connectivity_score"]),
        "largest_component_ratio": float(best_row["largest_component_ratio"]),
        "small_component_fraction": float(best_row["small_component_fraction"]),
        "alpha_fraction": float(best_row["alpha_fraction"]),
        "num_alpha_components": int(best_row["num_alpha_components"]),
    }



def show_alpha_report(pdb_name, results_dir, display_pngs=True):
    df, summary = load_alpha_artifacts(results_dir)
    recommendation = recommend_alpha_from_sweep(df)
    display(Markdown(f"## Alpha report — {pdb_name}"))
    display(Markdown(
        f"**Base alpha:** {summary.get('basealpha')}  
"
        f"**Tested alphas:** {summary.get('testedalphas')}  
"
        f"**Recommended alpha:** {recommendation['recommended_alpha']}"
    ))
    cols = [c for c in ['alpha', 'numalphatetrahedra', 'numpreliminarypockets', 'nummergedpockets', 'numfinalpockets', 'meanflowsteps', 'maxvolume'] if c in df.columns]
    display(df[cols].sort_values('alpha'))
    if display_pngs:
        for key in ['histpng', 'sweeppng', 'top3png']:
            path = summary.get('artifacts', {}).get(key)
            if path:
                display(Image(filename=path))
    return {
        'dataframe': df,
        'summary': summary,
        'recommendation': recommendation,
    }
