# Running ELM Analysis with tmux

This guide shows you how to run the long-running ELM analysis (~14 hours) in a tmux session so you can close your laptop without interrupting the process.

## What is tmux?

tmux is a "terminal multiplexer" that lets you create persistent terminal sessions. Even if you close your laptop, disconnect from SSH, or lose your connection, the session keeps running in the background.

## Quick Start

### 1. Stop the Currently Running Analysis

First, stop your current stuck Jupyter run:
- Go to the Jupyter notebook in VS Code
- Click the "Interrupt" button (square icon) in the notebook toolbar
- Or go to Kernel menu → Interrupt

### 2. Re-run Cell 4 with Domain Merging Fix

Before starting the tmux session, you need to apply the domain overlap fix:
- In the notebook `get_linker_info.ipynb`
- Re-run **Cell 4** (the one with `_merge_overlapping_domains()` and `compute_linker_regions()`)
- This loads the corrected domain merging logic

### 3. Start a tmux Session

Open a terminal and run:

```bash
cd ~/Desktop/work/protein_linkers
tmux new-session -s elm_analysis
```

This creates a new tmux session named "elm_analysis".

### 4. Run the ELM Analysis Script

Inside the tmux session:

```bash
# Basic run with default settings
python run_elm_analysis.py
```

**Available options:**

```bash
# Resume from checkpoint
python run_elm_analysis.py --resume

# Custom number of workers (default: 4)
python run_elm_analysis.py --max-workers 8

# Save checkpoints more frequently (default: every 10 regions)
python run_elm_analysis.py --checkpoint-freq 5

# Use different working directory
python run_elm_analysis.py --workdir ~/my_proteins

# Filter for high-confidence motifs only (LIG, DOC, MOD, TRG)
python run_elm_analysis.py --filter-confidence

# Adjust rate limiting (default: 60 seconds, don't go below ELM's limit!)
python run_elm_analysis.py --rate-limit 60

# Combine multiple options
python run_elm_analysis.py --resume --max-workers 6 --checkpoint-freq 20
```

**Get help:**
```bash
python run_elm_analysis.py --help
```

### 5. Detach from tmux (Close Laptop Safely)

Once the script is running, you can detach from the tmux session:

**Press:** `Ctrl+b` then `d` (release Ctrl+b, then press d)

You'll see: `[detached (from session elm_analysis)]`

**You can now close your laptop!** The analysis continues running in the background.

### 6. Check Progress Later

When you come back, reattach to the session:

```bash
tmux attach-session -s elm_analysis
```

Or use the shortcut:

```bash
tmux a -t elm_analysis
```

### 7. Monitor Progress

While attached, you can see:
- Real-time progress updates: `[123/850] ✓ protein1_domain_1_Name: 5 motifs found`
- Checkpoint saves every 10 regions
- Log file: `~/Desktop/work/protein_linkers/input_2/elm_analysis.log`

To view the log file from another terminal (while tmux is running):

```bash
tail -f ~/Desktop/work/protein_linkers/input_2/elm_analysis.log
```

## Key Features

### Checkpointing
- Progress is saved every 10 completed regions
- Checkpoint file: `~/Desktop/work/protein_linkers/input_2/elm_checkpoint.json`
- If the script crashes or you stop it, resume with: `python run_elm_analysis.py --resume`

### Logging
- All output is logged to: `~/Desktop/work/protein_linkers/input_2/elm_analysis.log`
- Includes timestamps for all events
- Safe to monitor in real-time with `tail -f`

### Output Files
- Final results (JSON): `~/Desktop/work/protein_linkers/input_2/elm_results.json`
- Final results (CSV): `~/Desktop/work/protein_linkers/input_2/elm_results.csv`
- Checkpoint: `~/Desktop/work/protein_linkers/input_2/elm_checkpoint.json`
- Log: `~/Desktop/work/protein_linkers/input_2/elm_analysis.log`

## Useful tmux Commands

| Command | Description |
|---------|-------------|
| `tmux new -s name` | Create new session named "name" |
| `tmux ls` | List all tmux sessions |
| `tmux a -t name` | Attach to session "name" |
| `Ctrl+b d` | Detach from current session |
| `Ctrl+b [` | Enter scroll mode (use arrow keys, q to quit) |
| `tmux kill-session -t name` | Kill session "name" |
| `Ctrl+d` or `exit` | Exit session (only when done!) |

## Common Scenarios

### Scenario 1: Check if still running
```bash
tmux ls
# If you see "elm_analysis" listed, it's running
tmux a -t elm_analysis
```

### Scenario 2: Script finished, want to see results
```bash
tmux a -t elm_analysis
# Review the final output, then exit the session:
exit
```

### Scenario 3: Script crashed, want to resume
```bash
tmux a -t elm_analysis
python run_elm_analysis.py --resume
# Detach again: Ctrl+b d
```

### Scenario 5: Want to speed up checkpointing for safety
```bash
# Save progress every 5 regions instead of 10
python run_elm_analysis.py --checkpoint-freq 5
```

### Scenario 6: Filter for high-confidence motifs only
```bash
# Only analyze LIG, DOC, MOD, TRG motif classes
python run_elm_analysis.py --filter-confidence
```

### Scenario 4: Want to stop everything
```bash
tmux a -t elm_analysis
# Press Ctrl+C to stop the script
exit
# Or just kill the session:
tmux kill-session -t elm_analysis
```

## Estimated Timeline

- **Total regions:** 850
- **Rate limit:** 1 request per 60 seconds
- **Minimum time:** ~850 minutes (14.2 hours)
- **With retries/errors:** ~15-16 hours

**Recommendation:** Start the analysis in the evening before bed. Check progress in the morning!

## Troubleshooting

### "No session found"
The session doesn't exist or was killed. Start a new one with `tmux new -s elm_analysis`.

### Can't see scroll history
Press `Ctrl+b [` to enter scroll mode. Use arrow keys or Page Up/Down. Press `q` to exit.

### Script stopped but tmux is running
Reattach with `tmux a -t elm_analysis`, check the error, and restart with `--resume` if needed.

### Want to run multiple things in tmux
Create multiple sessions with different names:
```bash
tmux new -s analysis1
tmux new -s analysis2
tmux ls  # See all sessions
tmux a -t analysis1  # Switch between them
```

## Alternative: Running in Jupyter via tmux

If you prefer to run the Jupyter notebook in tmux instead:

1. Start tmux: `tmux new -s jupyter`
2. Start Jupyter: `jupyter notebook --no-browser`
3. Copy the URL with token
4. Detach: `Ctrl+b d`
5. In VS Code, connect to the remote Jupyter server using the URL
6. Run your analysis cells
7. The kernel stays running even if you close your laptop!

**Note:** The standalone Python script with checkpointing is safer for long runs.

## Summary

✅ **Before you go:**
1. Stop current Jupyter run
2. Re-run Cell 4 in notebook (domain merging fix)
3. Start tmux: `tmux new -s elm_analysis`
4. Run script: `python run_elm_analysis.py`
5. Detach: `Ctrl+b` then `d`
6. Close laptop safely!

✅ **When you return:**
1. Reattach: `tmux a -t elm_analysis`
2. Check progress
3. When done, review results and `exit`
