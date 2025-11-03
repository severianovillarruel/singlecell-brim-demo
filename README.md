# Single-cell Analysis Demo (BRIM-Inspired)

Synthetic scRNA-seq example showing:
- Embedding + clustering (PCA + k-means)
- Simple differential expression (CR_PR vs SD_PD)
- Clear, reproducible structure

**Contents**
- `data/` (created by script)
- `figures/` (created by script)
- `src/generate_data.py` (makes synthetic counts, metadata, DE, figures)
- `notebooks/` (placeholder for a future notebook)

**Quickstart**
```bash
pip install -r requirements.txt
python src/generate_data.py
