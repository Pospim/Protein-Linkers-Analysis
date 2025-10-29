# ELM Analysis - Quick Start

## 🚀 Start Analysis in tmux

```bash
cd ~/Desktop/work/protein_linkers
tmux new-session -s elm_analysis
python run_elm_analysis.py
```

**Detach (close laptop safely):** Press `Ctrl+b` then `d`

## 📊 Check Progress Later

```bash
tmux a-session -s elm_analysis
```

## ⚙️ Common Options

```bash
# Resume from checkpoint
python run_elm_analysis.py --resume

# More frequent checkpoints (every 5 regions)
python run_elm_analysis.py --checkpoint-freq 5

# More workers (doesn't speed up much due to rate limiting)
python run_elm_analysis.py --max-workers 8

# High-confidence motifs only
python run_elm_analysis.py --filter-confidence

# Combine options
python run_elm_analysis.py --resume --checkpoint-freq 5
```

## 📁 Output Files

- `input_2/elm_results.json` - Final results (JSON)
- `input_2/elm_results.csv` - Final results (CSV)
- `input_2/elm_checkpoint.json` - Checkpoint (for resume)
- `input_2/elm_analysis.log` - Detailed log

## 🔍 Monitor Progress

```bash
# In another terminal
tail -f ~/Desktop/work/protein_linkers/input_2/elm_analysis.log
```

## ⏱️ Expected Runtime

- **850 regions** × **60 seconds** = **~14 hours**

## 🆘 Emergency

```bash
# Stop everything
tmux kill-session -t elm_analysis

# Or attach and Ctrl+C
tmux a -t elm_analysis
# Press Ctrl+C, then exit
```

## 📖 Full Documentation

See `TMUX_GUIDE.md` for detailed instructions and troubleshooting.
